"""
Plot subcommand — generate diagnostic and summary plots.

Reads chromosome-level and (optionally) bin-level statistics plus training
loss, then produces histograms, sex-assignment scatter, per-contig boxplots,
per-bin chromosome plots, per-sample CN plots, and model-diagnostic plots.
"""

from __future__ import annotations

import argparse
import logging
import os
from typing import Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages

from gatk_sv_ploidy._util import (
    add_chromosome_labels,
    format_column_name,
    get_chromosome_type,
    save_and_close_plot,
)
from gatk_sv_ploidy._plot_detail import (
    plot_cn_per_bin_chromosome,
    plot_cn_per_contig_boxplot,
    plot_sample_with_variance,
    plot_sex_assignments,
)
from gatk_sv_ploidy._plot_ppd import run_ppd_plots
from gatk_sv_ploidy.data import load_site_data

logger = logging.getLogger(__name__)
_BINQ_FIELD_OPTIONS = ["auto", "BINQ15", "BINQ20", "CALLQ15", "CALLQ20"]


def _resolve_binq_field(
    bin_quality_df: pd.DataFrame,
    requested_field: str,
) -> str:
    if requested_field != "auto":
        return requested_field
    if "BINQ20" in bin_quality_df.columns:
        return "BINQ20"
    return "CALLQ20"


# ── histogram helpers ───────────────────────────────────────────────────────

_CHR_PALETTE = {"Autosomal": "#1f77b4", "chrX": "#ff7f0e", "chrY": "#2ca02c"}


def _apply_plot_depth_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Replace raw-depth summary columns with plot-normalized variants."""
    plot_map = {
        "mean_depth": "plot_mean_depth",
        "std_depth": "plot_std_depth",
        "median_depth": "plot_median_depth",
        "mad_depth": "plot_mad_depth",
    }
    out = df.copy()
    used = False
    for base_col, plot_col in plot_map.items():
        if plot_col in out.columns:
            out[base_col] = out[plot_col]
            used = True
    if used:
        logger.info("Using plot-normalized chromosome depth summaries.")
    return out


def _apply_plot_depth_bin_columns(bin_df: pd.DataFrame | None) -> pd.DataFrame | None:
    """Replace raw per-bin observed depth with plot-normalized depth."""
    if bin_df is None:
        return None
    out = bin_df.copy()
    if "plot_depth" in out.columns:
        out["observed_depth"] = out["plot_depth"]
        logger.info("Using plot-normalized per-bin depth profiles.")
    return out


def _apply_plot_depth_sex_columns(
    sex_df: pd.DataFrame,
    chr_df: pd.DataFrame,
) -> pd.DataFrame:
    """Populate sex-scatter depths from plot-normalized chromosome summaries."""
    out = sex_df.copy()
    chr_depth_df = chr_df[chr_df["chromosome"].isin(["chrX", "chrY"])]
    if chr_depth_df.empty:
        return out

    depth_pivot = chr_depth_df.pivot(
        index="sample",
        columns="chromosome",
        values="median_depth",
    )
    if "chrX" in depth_pivot.columns:
        out["chrX_depth"] = out["sample"].map(depth_pivot["chrX"])
    if "chrY" in depth_pivot.columns:
        out["chrY_depth"] = out["sample"].map(depth_pivot["chrY"])
    return out


def _annotate_ignored_bins(
    bin_df: pd.DataFrame | None,
    ignored_bins_df: pd.DataFrame | None,
) -> pd.DataFrame | None:
    """Merge call-time ignored-bin annotations into per-bin plotting data."""
    if bin_df is None or ignored_bins_df is None:
        return bin_df

    required_cols = {"sample", "chr", "start", "end"}
    missing_cols = required_cols - set(ignored_bins_df.columns)
    if missing_cols:
        missing = ", ".join(sorted(missing_cols))
        raise ValueError(f"Ignored bins file is missing required columns: {missing}")

    keep_cols = [
        col
        for col in [
            "sample",
            "chr",
            "start",
            "end",
            "ignored_in_call",
            "binq_field",
            "binq_value",
            "min_binq",
        ]
        if col in ignored_bins_df.columns
    ]
    annot_df = ignored_bins_df[keep_cols].copy()
    if "ignored_in_call" not in annot_df.columns:
        annot_df["ignored_in_call"] = True

    out = bin_df.merge(
        annot_df,
        on=["sample", "chr", "start", "end"],
        how="left",
        suffixes=("", "__ignored"),
    )
    for col in ("binq_field", "binq_value", "min_binq"):
        ignored_col = f"{col}__ignored"
        if ignored_col not in out.columns:
            continue
        if col in out.columns:
            out[col] = out[col].where(out[col].notna(), out[ignored_col])
        else:
            out[col] = out[ignored_col]
        out = out.drop(columns=[ignored_col])
    ignored_series = out["ignored_in_call"].fillna(False)
    if ignored_series.dtype == object:
        ignored_series = ignored_series.astype(str).str.lower().isin(["true", "1", "yes"])
    out["ignored_in_call"] = ignored_series.astype(bool)

    ignored_fraction = (
        out.groupby(["chr", "start", "end"], sort=False)["ignored_in_call"]
        .mean()
        .rename("ignored_fraction_in_call")
        .reset_index()
    )
    out = out.merge(
        ignored_fraction,
        on=["chr", "start", "end"],
        how="left",
    )
    out["ignored_fraction_in_call"] = out["ignored_fraction_in_call"].fillna(0.0)

    logger.info(
        "Annotated %d ignored sample-bin calls across %d unique bins.",
        int(out["ignored_in_call"].sum()),
        int((ignored_fraction["ignored_fraction_in_call"] > 0).sum()),
    )
    return out


def _annotate_binq_values(
    bin_df: pd.DataFrame | None,
    bin_quality_df: pd.DataFrame | None,
    binq_field: str,
) -> pd.DataFrame | None:
    """Merge per-bin PPD quality values into plotting bin statistics."""
    if bin_df is None or bin_quality_df is None:
        return bin_df

    binq_field = _resolve_binq_field(bin_quality_df, binq_field)
    required_cols = {"chr", "start", "end", binq_field}
    missing_cols = required_cols - set(bin_quality_df.columns)
    if missing_cols:
        missing = ", ".join(sorted(missing_cols))
        raise ValueError(
            f"PPD bin quality file is missing required columns: {missing}"
        )

    annot_df = bin_quality_df[["chr", "start", "end", binq_field]].copy()
    annot_df["binq_field"] = binq_field
    annot_df["binq_value"] = annot_df[binq_field].astype(float)
    annot_df = annot_df.drop(columns=[binq_field])

    out = bin_df.merge(
        annot_df,
        on=["chr", "start", "end"],
        how="left",
        suffixes=("", "__quality"),
    )
    for col in ("binq_field", "binq_value"):
        quality_col = f"{col}__quality"
        if quality_col not in out.columns:
            continue
        if col in out.columns:
            out[col] = out[col].where(out[col].notna(), out[quality_col])
        else:
            out[col] = out[quality_col]
        out = out.drop(columns=[quality_col])

    logger.info(
        "Annotated BINQ values for %d sample-bin rows using %s.",
        int(out["binq_value"].notna().sum()) if "binq_value" in out.columns else 0,
        binq_field,
    )
    return out


def _hist_by_hue(
    data: pd.DataFrame,
    x_col: str,
    hue_col: str,
    title: str,
    xlabel: str,
    output_dir: str,
    filename: str,
    palette: dict,
    bins: int = 50,
    highlight_sample: str = "",
    sample_values: Optional[list] = None,
) -> None:
    """Histogram coloured by *hue_col* with optional sample highlight."""
    fig, ax = plt.subplots(figsize=(10, 6))
    sns.histplot(data=data, x=x_col, hue=hue_col, bins=bins, kde=False,
                 palette=palette, alpha=0.6, edgecolor="black", ax=ax,
                 multiple="stack")

    if highlight_sample and sample_values:
        for v in sample_values:
            ax.axvline(v, color="magenta", linestyle="--", linewidth=2, alpha=0.8)
        ax.plot([], [], color="magenta", linestyle="--", linewidth=2,
                label=" " + highlight_sample)
        ax.legend(loc="best", framealpha=0.9)

    ax.set_xlabel(xlabel)
    ax.set_ylabel("Frequency")
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    save_and_close_plot(output_dir, filename)


# ── chromosome-type histograms ─────────────────────────────────────────────


def plot_histograms_by_chr_type(
    df: pd.DataFrame, output_dir: str, highlight_sample: str = ""
) -> None:
    """Create histograms of every numeric column, grouped by chromosome type.

    Args:
        df: ``chromosome_stats.tsv`` DataFrame.
        output_dir: Base output directory.
        highlight_sample: Sample ID to mark with vertical lines.
    """
    df = df.copy()
    df["chr_type"] = df["chromosome"].apply(get_chromosome_type)

    exclude = {"sample", "chromosome", "chr_type"}
    numeric = [c for c in df.select_dtypes(include=[np.number]).columns if c not in exclude]

    # Prepare highlight values
    hl_vals: dict = {}
    if highlight_sample and highlight_sample in df["sample"].values:
        hdf = df[df["sample"] == highlight_sample]
        for c in numeric:
            hl_vals[c] = hdf[c].tolist()

    # sample_var_map → sample_stdev
    if "sample_var_map" in numeric:
        numeric.remove("sample_var_map")
        df["sample_stdev"] = np.sqrt(df["sample_var_map"])
        sv = np.sqrt(hdf["sample_var_map"]).tolist() if hl_vals else None
        _hist_by_hue(df, "sample_stdev", "chr_type",
                     "Sample Stdev by Chromosome Type", "Sample Stdev",
                     output_dir, "hist_sample_stdev.png", _CHR_PALETTE,
                     highlight_sample=highlight_sample, sample_values=sv)

    for col in numeric:
        _hist_by_hue(
            df, col, "chr_type",
            f"{format_column_name(col)} by Chromosome Type",
            format_column_name(col),
            output_dir, f"hist_{col}.png", _CHR_PALETTE,
            highlight_sample=highlight_sample,
            sample_values=hl_vals.get(col),
        )


# ── aneuploid-only histograms ──────────────────────────────────────────────


def plot_aneuploid_histograms(
    df: pd.DataFrame, output_dir: str, highlight_sample: str = ""
) -> None:
    """Histograms restricted to aneuploid rows.

    Args:
        df: ``chromosome_stats.tsv`` DataFrame.
        output_dir: Base output directory.
        highlight_sample: Sample ID to highlight.
    """
    adf = df[df["is_aneuploid"]].copy()
    if adf.empty:
        logger.warning("No aneuploid samples found — skipping aneuploid histograms")
        return

    adf["chr_type"] = adf["chromosome"].apply(get_chromosome_type)

    # Chromosome count bar chart
    chr_order = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    chr_order = [c for c in chr_order if c in adf["chromosome"].unique()]
    fig, ax = plt.subplots(figsize=(12, 6))
    sns.countplot(data=adf, x="chromosome", order=chr_order, ax=ax,
                  edgecolor="black", alpha=0.7)
    ax.set_xlabel("Chromosome")
    ax.set_ylabel("Count")
    ax.set_title("Aneuploid Samples by Chromosome")
    ax.grid(True, alpha=0.3, axis="y")
    plt.xticks(rotation=45, ha="right")
    save_and_close_plot(output_dir, "hist_aneuploid_chromosome.png")

    hl_vals: dict = {}
    if highlight_sample and highlight_sample in adf["sample"].values:
        hdf = adf[adf["sample"] == highlight_sample]
        for c in (
            "copy_number",
            "mean_cn_probability",
            "plq",
            "mean_depth",
            "median_depth",
            "std_depth",
        ):
            if c in hdf.columns:
                hl_vals[c] = hdf[c].tolist()

    for col in (
        "copy_number",
        "mean_cn_probability",
        "plq",
        "mean_depth",
        "median_depth",
        "std_depth",
    ):
        if col not in adf.columns:
            continue
        _hist_by_hue(
            adf, col, "chr_type",
            f"{format_column_name(col)} (Aneuploid Only)",
            format_column_name(col),
            output_dir, f"hist_aneuploid_{col}.png", _CHR_PALETTE,
            bins=30, highlight_sample=highlight_sample,
            sample_values=hl_vals.get(col),
        )


# ── median-depth distribution ──────────────────────────────────────────────


def plot_median_depth_distributions(
    df: pd.DataFrame, output_dir: str, highlight_sample: str = ""
) -> None:
    """Error-bar plots of median depth ± MAD for every chromosome.

    Args:
        df: ``chromosome_stats.tsv`` DataFrame.
        output_dir: Base output directory.
        highlight_sample: Sample ID to highlight.
    """
    if "median_depth" not in df.columns or "mad_depth" not in df.columns:
        logger.warning("median_depth / mad_depth missing — skipping")
        return

    stats_df = df[
        [
            "sample",
            "chromosome",
            "copy_number",
            "mean_depth",
            "median_depth",
            "mad_depth",
        ]
    ].copy()
    stats_path = os.path.join(output_dir, "median_depth_distribution_stats.tsv")
    stats_df.to_csv(stats_path, sep="\t", index=False)

    chroms = sorted(
        stats_df["chromosome"].unique(),
        key=lambda x: (x in ("chrX", "chrY"),
                       int(x.replace("chr", "")) if x.replace("chr", "").isdigit() else 0, x),
    )

    pdf_path = os.path.join(output_dir, "median_depth_distributions.pdf")
    with PdfPages(pdf_path) as pdf:
        for chrom in chroms:
            cdf = stats_df[stats_df["chromosome"] == chrom].sort_values("median_depth")
            fig, ax = plt.subplots(figsize=(14, 6))

            hl_pos = None
            if highlight_sample and highlight_sample in cdf["sample"].values:
                hl_pos = cdf.index.get_loc(
                    cdf[cdf["sample"] == highlight_sample].index[0]
                )

            ax.errorbar(range(len(cdf)), cdf["median_depth"], yerr=cdf["mad_depth"],
                        fmt="o", markerfacecolor="blue", markeredgecolor="blue",
                        ecolor="black", capsize=2, linestyle="none", alpha=0.7)

            if hl_pos is not None:
                row = cdf.iloc[hl_pos]
                ax.errorbar([hl_pos], [row["median_depth"]], yerr=[row["mad_depth"]],
                            fmt="^", markerfacecolor="magenta", markeredgecolor="black",
                            ecolor="magenta", capsize=4, markersize=10, linewidth=2,
                            linestyle="none", label=" " + highlight_sample)
                ax.legend(loc="best", framealpha=0.9)

            ax.set_xlabel("Sample Index")
            ax.set_ylabel("Median Depth")
            ax.set_title(f"Median Depth for {chrom} with MAD")
            ax.set_xticks(range(len(cdf)))
            ax.set_xticklabels(cdf["sample"], rotation=90, fontsize=6)
            plt.tight_layout()
            pdf.savefig()
            plt.close()

    logger.info("Saved median depth distributions to %s", pdf_path)


# ── training-loss plot ──────────────────────────────────────────────────────


def plot_training_loss(loss_df: pd.DataFrame, output_dir: str) -> None:
    """Line plot of ELBO loss over training epochs.

    Args:
        loss_df: DataFrame with ``epoch`` and ``elbo`` columns.
        output_dir: Base output directory.
    """
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(loss_df["epoch"], loss_df["elbo"])
    ax.set_xlabel("Epoch")
    ax.set_ylabel("ELBO")
    ax.set_title("Training Loss")
    ax.grid(True, alpha=0.3)
    save_and_close_plot(output_dir, "training_loss.png")


# ── bin variance / bias ────────────────────────────────────────────────────


def plot_bin_variance_bias(bin_df: pd.DataFrame, output_dir: str) -> None:
    """Two-panel plot of per-bin bias and stdev.

    Args:
        bin_df: ``bin_stats.tsv.gz`` DataFrame.
        output_dir: Base output directory.
    """
    unique = bin_df.groupby(["chr", "start", "end"]).first().reset_index()
    unique = unique.sort_values(["chr", "start"])

    fig, axes = plt.subplots(2, 1, figsize=(16, 8))
    x = np.arange(len(unique))
    chrs = unique["chr"].values

    axes[0].plot(x, unique["bin_bias"].values, "o", alpha=0.5, markersize=3,
                 color="black", label="Bin bias")
    axes[0].axhline(1.0, color="red", linestyle="--", alpha=0.5, label="No bias")
    axes[0].set_ylabel("Bin mean bias")
    axes[0].set_title("Bin Posteriors")
    axes[0].set_xlim([x.min(), x.max()])

    axes[1].plot(x, np.sqrt(unique["bin_var"].values), "o", alpha=0.5,
                 markersize=3, color="purple", label="Bin stdev")
    axes[1].set_ylabel("Bin stdev")
    axes[1].set_xlim([x.min(), x.max()])

    if "ignored_fraction_in_call" in unique.columns:
        ignored_mask = unique["ignored_fraction_in_call"].to_numpy(dtype=float) > 0
        if np.any(ignored_mask):
            ignored_frac = unique.loc[ignored_mask, "ignored_fraction_in_call"].to_numpy(dtype=float)
            ignored_size = 14.0 + 80.0 * ignored_frac
            axes[0].scatter(
                x[ignored_mask],
                unique.loc[ignored_mask, "bin_bias"].to_numpy(dtype=float),
                s=ignored_size,
                facecolors="none",
                edgecolors="#FFB300",
                linewidths=1.2,
                label="Ignored in call",
                zorder=4,
            )
            axes[1].scatter(
                x[ignored_mask],
                np.sqrt(unique.loc[ignored_mask, "bin_var"].to_numpy(dtype=float)),
                s=ignored_size,
                facecolors="none",
                edgecolors="#FFB300",
                linewidths=1.2,
                label="Ignored in call",
                zorder=4,
            )

    for a in axes:
        add_chromosome_labels(a, chrs)
        if a.get_legend_handles_labels()[0]:
            a.legend(loc="best", framealpha=0.9)

    plt.tight_layout()
    save_and_close_plot(output_dir, "bin_posteriors.png")


# ── CN posterior entropy ────────────────────────────────────────────────────


def plot_cn_posterior_entropy(bin_df: pd.DataFrame, output_dir: str) -> None:
    """Genome-wide plot of per-bin CN posterior Shannon entropy.

    Low-confidence bins (multiple plausible CN states) have high entropy;
    confident bins approach zero.  Systematically elevated entropy in a
    region signals poor model fit.

    Args:
        bin_df: ``bin_stats.tsv.gz`` DataFrame.
        output_dir: Base output directory.
    """
    prob_cols = [f"cn_prob_{i}" for i in range(6)]
    group_cols = prob_cols.copy()
    if "ignored_fraction_in_call" in bin_df.columns:
        group_cols.append("ignored_fraction_in_call")
    probs = bin_df.groupby(["chr", "start", "end"])[group_cols].mean().reset_index()
    probs = probs.sort_values(["chr", "start"])

    p = probs[prob_cols].values
    p = np.clip(p, 1e-10, 1.0)
    entropy = -np.sum(p * np.log2(p), axis=1)

    chrs = probs["chr"].values
    x = np.arange(len(probs))

    fig, ax = plt.subplots(figsize=(16, 4))
    ax.scatter(x, entropy, s=3, alpha=0.5, color="darkred", rasterized=True)
    if "ignored_fraction_in_call" in probs.columns:
        ignored_mask = probs["ignored_fraction_in_call"].to_numpy(dtype=float) > 0
        if np.any(ignored_mask):
            ignored_frac = probs.loc[ignored_mask, "ignored_fraction_in_call"].to_numpy(dtype=float)
            ignored_size = 14.0 + 80.0 * ignored_frac
            ax.scatter(
                x[ignored_mask],
                entropy[ignored_mask],
                s=ignored_size,
                facecolors="none",
                edgecolors="#FFB300",
                linewidths=1.2,
                label="Ignored in call",
                zorder=4,
            )
    ax.set_ylabel("Shannon Entropy (bits)")
    ax.set_title("CN Posterior Entropy per Bin (averaged across samples)")
    ax.set_xlim([x.min(), x.max()])
    ax.grid(True, alpha=0.3, axis="y")
    add_chromosome_labels(ax, chrs)
    if ax.get_legend_handles_labels()[0]:
        ax.legend(loc="best", framealpha=0.9)
    save_and_close_plot(output_dir, "cn_posterior_entropy.png")


# ── chromosome CN heatmap ──────────────────────────────────────────────────


def plot_chromosome_cn_heatmap(
    df: pd.DataFrame, output_dir: str
) -> None:
    """Heatmap of dominant CN per sample × chromosome.

    Colour intensity encodes the CN probability. Useful for spotting
    sample-wide or chromosome-wide patterns.

    Args:
        df: ``chromosome_stats.tsv`` DataFrame.
        output_dir: Base output directory.
    """
    chr_order = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    chr_order = [c for c in chr_order if c in df["chromosome"].unique()]

    cn_pivot = df.pivot(index="sample", columns="chromosome", values="copy_number")
    cn_pivot = cn_pivot.reindex(columns=chr_order)
    prob_pivot = df.pivot(
        index="sample", columns="chromosome", values="mean_cn_probability",
    )
    prob_pivot = prob_pivot.reindex(columns=chr_order)

    # Sort samples: first by chrX CN, then by overall mean CN
    cn_pivot["_sort"] = cn_pivot.mean(axis=1)
    cn_pivot = cn_pivot.sort_values("_sort")
    cn_pivot = cn_pivot.drop("_sort", axis=1)
    prob_pivot = prob_pivot.reindex(cn_pivot.index)

    n_samples = len(cn_pivot)
    show_labels = n_samples <= 60

    fig, ax = plt.subplots(figsize=(14, max(4, 0.25 * n_samples)))
    im = ax.imshow(cn_pivot.values, aspect="auto", cmap="YlOrRd",
                   vmin=0, vmax=5, interpolation="nearest")
    ax.set_xticks(np.arange(len(chr_order)))
    ax.set_xticklabels([c.replace("chr", "") for c in chr_order], rotation=45)
    if show_labels:
        ax.set_yticks(np.arange(n_samples))
        ax.set_yticklabels(cn_pivot.index, fontsize=5)
    else:
        ax.set_yticks([])
    ax.set_xlabel("Chromosome")
    ax.set_ylabel("Sample")
    ax.set_title("Dominant Copy Number per Sample × Chromosome")
    plt.colorbar(im, ax=ax, label="Copy Number")
    plt.tight_layout()
    save_and_close_plot(output_dir, "chromosome_cn_heatmap.png")


def plot_chromosome_plq_heatmap(
    df: pd.DataFrame, output_dir: str
) -> None:
    """Heatmap of chromosome-level PLQ per sample × chromosome."""
    if "plq" not in df.columns:
        logger.warning("plq missing from chromosome stats — skipping PLQ heatmap")
        return

    chr_order = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    chr_order = [c for c in chr_order if c in df["chromosome"].unique()]

    cn_pivot = df.pivot(index="sample", columns="chromosome", values="copy_number")
    cn_pivot = cn_pivot.reindex(columns=chr_order)
    plq_pivot = df.pivot(index="sample", columns="chromosome", values="plq")
    plq_pivot = plq_pivot.reindex(columns=chr_order)

    cn_pivot["_sort"] = cn_pivot.mean(axis=1)
    cn_pivot = cn_pivot.sort_values("_sort")
    cn_pivot = cn_pivot.drop("_sort", axis=1)
    plq_pivot = plq_pivot.reindex(cn_pivot.index)

    n_samples = len(plq_pivot)
    show_labels = n_samples <= 60

    values = np.ma.masked_invalid(plq_pivot.to_numpy(dtype=float))
    cmap = plt.cm.viridis.copy()
    cmap.set_bad(color="#EEEEEE")

    fig, ax = plt.subplots(figsize=(14, max(4, 0.25 * n_samples)))
    im = ax.imshow(values, aspect="auto", cmap=cmap,
                   vmin=0, vmax=99, interpolation="nearest")
    ax.set_xticks(np.arange(len(chr_order)))
    ax.set_xticklabels([c.replace("chr", "") for c in chr_order], rotation=45)
    if show_labels:
        ax.set_yticks(np.arange(n_samples))
        ax.set_yticklabels(plq_pivot.index, fontsize=5)
    else:
        ax.set_yticks([])
    ax.set_xlabel("Chromosome")
    ax.set_ylabel("Sample")
    ax.set_title("Chromosome PLQ per Sample × Chromosome")
    plt.colorbar(im, ax=ax, label="PLQ")
    plt.tight_layout()
    save_and_close_plot(output_dir, "chromosome_plq_heatmap.png")


def plot_binq_genome_profile(
    bin_df: pd.DataFrame,
    output_dir: str,
) -> None:
    """Genome-wide profile of per-bin BINQ values."""
    if "binq_value" not in bin_df.columns:
        logger.warning("binq_value missing from bin stats — skipping BINQ profile")
        return

    group_cols = ["chr", "start", "end", "binq_value"]
    if "binq_field" in bin_df.columns:
        group_cols.append("binq_field")
    if "ignored_fraction_in_call" in bin_df.columns:
        group_cols.append("ignored_fraction_in_call")
    unique = bin_df.groupby(["chr", "start", "end"], sort=False)[group_cols[3:]].first().reset_index()
    unique = unique.sort_values(["chr", "start"])

    values = unique["binq_value"].to_numpy(dtype=float)
    finite_mask = np.isfinite(values)
    if not np.any(finite_mask):
        logger.warning("No finite BINQ values available — skipping BINQ profile")
        return

    field_label = "BINQ"
    if "binq_field" in unique.columns:
        field_series = unique["binq_field"].dropna().astype(str)
        if not field_series.empty:
            field_label = field_series.iloc[0]

    chrs = unique["chr"].to_numpy()
    x = np.arange(len(unique))

    fig, ax = plt.subplots(figsize=(16, 4))
    ax.scatter(
        x[finite_mask],
        values[finite_mask],
        s=5,
        alpha=0.6,
        color="#7E57C2",
        rasterized=True,
    )
    if "ignored_fraction_in_call" in unique.columns:
        ignored_mask = unique["ignored_fraction_in_call"].to_numpy(dtype=float) > 0
        if np.any(ignored_mask):
            ignored_frac = unique.loc[ignored_mask, "ignored_fraction_in_call"].to_numpy(dtype=float)
            ignored_size = 14.0 + 80.0 * ignored_frac
            ax.scatter(
                x[ignored_mask],
                values[ignored_mask],
                s=ignored_size,
                facecolors="none",
                edgecolors="#FFB300",
                linewidths=1.2,
                label="Ignored in call",
                zorder=4,
            )
    ax.set_ylabel(field_label)
    ax.set_title(f"{field_label} per Bin")
    ax.set_ylim([0, 99])
    ax.set_xlim([x.min(), x.max()])
    ax.grid(True, alpha=0.3, axis="y")
    add_chromosome_labels(ax, chrs)
    if ax.get_legend_handles_labels()[0]:
        ax.legend(loc="best", framealpha=0.9)
    save_and_close_plot(output_dir, "binq_genome_profile.png")


# ── ELBO convergence gradient ──────────────────────────────────────────────


def plot_training_loss_with_gradient(
    loss_df: pd.DataFrame, output_dir: str, window: int = 50,
) -> None:
    """Two-panel training diagnostic: ELBO loss + rolling gradient.

    The gradient panel helps assess whether training converged (gradient
    near zero) or is still improving.

    Args:
        loss_df: DataFrame with ``epoch`` and ``elbo`` columns.
        output_dir: Base output directory.
        window: Rolling window size for the gradient estimate.
    """
    epochs = loss_df["epoch"].values
    elbo = loss_df["elbo"].values

    fig, axes = plt.subplots(2, 1, figsize=(12, 7), sharex=True,
                             gridspec_kw={"height_ratios": [2, 1]})

    # Panel 1: ELBO loss
    ax = axes[0]
    ax.plot(epochs, elbo, color="steelblue", linewidth=1)
    ax.set_ylabel("ELBO")
    ax.set_title("Training Loss & Convergence Diagnostic")
    ax.grid(True, alpha=0.3)

    # Panel 2: Rolling gradient (finite differences)
    ax = axes[1]
    if len(elbo) > window:
        gradient = np.convolve(np.diff(elbo), np.ones(window) / window, mode="valid")
        grad_x = epochs[1:len(gradient) + 1]
        ax.plot(grad_x, gradient, color="darkred", linewidth=1)
        ax.axhline(0, color="gray", linestyle="--", alpha=0.5)
    else:
        grad = np.diff(elbo)
        ax.plot(epochs[1:], grad, color="darkred", linewidth=1)
        ax.axhline(0, color="gray", linestyle="--", alpha=0.5)
    ax.set_xlabel("Epoch")
    ax.set_ylabel(f"∂ELBO/∂epoch (rolling {window})")
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    save_and_close_plot(output_dir, "training_loss_gradient.png")


# ── parameter correlation / distribution ───────────────────────────────────


def plot_parameter_diagnostics(bin_df: pd.DataFrame, output_dir: str) -> None:
    """Parameter distribution and correlation diagnostics.

    Three panels:
    1. Bin bias vs bin stdev scatter — correlated values may indicate
       model degeneracy.
    2. Sample variance distribution histogram.
    3. Bin bias distribution histogram.

    Args:
        bin_df: ``bin_stats.tsv.gz`` DataFrame.
        output_dir: Base output directory.
    """
    unique_bins = bin_df.groupby(["chr", "start", "end"]).first().reset_index()
    unique_bins = unique_bins.sort_values(["chr", "start"])
    all_svars = bin_df.groupby("sample")["sample_var"].first().values

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    # Panel 1: bin_bias vs bin_stdev scatter
    ax = axes[0]
    bb = unique_bins["bin_bias"].values
    bs = np.sqrt(unique_bins["bin_var"].values)
    ax.scatter(bb, bs, s=5, alpha=0.4, color="steelblue", rasterized=True)
    ax.set_xlabel("Bin Bias")
    ax.set_ylabel("Bin Stdev")
    ax.set_title("Bin Bias vs Stdev")
    ax.grid(True, alpha=0.3)

    # Panel 2: sample variance histogram
    ax = axes[1]
    ax.hist(np.sqrt(all_svars), bins=30, alpha=0.7, edgecolor="black",
            linewidth=0.5, color="teal")
    ax.set_xlabel("Sample Stdev (√variance)")
    ax.set_ylabel("Count")
    ax.set_title("Sample Variance Distribution")
    ax.grid(True, alpha=0.3)

    # Panel 3: bin bias histogram
    ax = axes[2]
    ax.hist(bb, bins=50, alpha=0.7, edgecolor="black", linewidth=0.5,
            color="coral")
    ax.axvline(1.0, color="red", linestyle="--", alpha=0.6, label="No bias (1.0)")
    ax.set_xlabel("Bin Bias")
    ax.set_ylabel("Count")
    ax.set_title("Bin Bias Distribution")
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    save_and_close_plot(output_dir, "parameter_diagnostics.png")


def plot_site_af_estimates(site_af_df: pd.DataFrame, output_dir: str) -> None:
    """Generate site-level AF estimate diagnostics."""
    required_cols = {
        "chr",
        "input_site_pop_af",
        "naive_bayes_site_pop_af",
        "effective_site_pop_af",
        "pooled_observed_af",
    }
    missing_cols = required_cols - set(site_af_df.columns)
    if missing_cols:
        missing = ", ".join(sorted(missing_cols))
        raise ValueError(
            f"Site AF estimates file is missing required columns: {missing}"
        )

    plot_df = site_af_df.copy()
    plot_df["chr_type"] = plot_df["chr"].apply(get_chromosome_type)

    fig, ax = plt.subplots(figsize=(10, 6))
    sns.histplot(
        data=plot_df,
        x="effective_site_pop_af",
        hue="chr_type",
        bins=50,
        palette=_CHR_PALETTE,
        multiple="stack",
        edgecolor="black",
        alpha=0.6,
        ax=ax,
    )
    ax.set_xlabel("Effective Site Population AF")
    ax.set_ylabel("Count")
    ax.set_title("Effective Site AF Distribution")
    ax.grid(True, alpha=0.3)
    save_and_close_plot(output_dir, "site_af_effective_distribution.png")

    scatter_df = plot_df[np.isfinite(plot_df["pooled_observed_af"])].copy()
    if not scatter_df.empty:
        fig, ax = plt.subplots(figsize=(7, 7))
        sns.scatterplot(
            data=scatter_df,
            x="pooled_observed_af",
            y="naive_bayes_site_pop_af",
            hue="chr_type",
            palette=_CHR_PALETTE,
            s=16,
            alpha=0.5,
            linewidth=0,
            ax=ax,
        )
        ax.plot([0.0, 1.0], [0.0, 1.0], linestyle="--", color="black", alpha=0.6)
        ax.set_xlim(0.0, 1.0)
        ax.set_ylim(0.0, 1.0)
        ax.set_xlabel("Pooled Observed AF")
        ax.set_ylabel("Naive-Bayes Site AF")
        ax.set_title("Pooled Observed AF vs Naive-Bayes Site AF")
        ax.grid(True, alpha=0.3)
        save_and_close_plot(output_dir, "site_af_pooled_vs_naive_bayes.png")

    fig, ax = plt.subplots(figsize=(7, 7))
    sns.scatterplot(
        data=plot_df,
        x="input_site_pop_af",
        y="effective_site_pop_af",
        hue="chr_type",
        palette=_CHR_PALETTE,
        s=16,
        alpha=0.5,
        linewidth=0,
        ax=ax,
    )
    ax.plot([0.0, 1.0], [0.0, 1.0], linestyle="--", color="black", alpha=0.6)
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.set_xlabel("Input Site Population AF")
    ax.set_ylabel("Effective Site Population AF")
    ax.set_title("Input vs Effective Site AF")
    ax.grid(True, alpha=0.3)
    save_and_close_plot(output_dir, "site_af_input_vs_effective.png")


# ── orchestrators ──────────────────────────────────────────────────────────


def _run_ploidy_plots(
    df: pd.DataFrame,
    bin_df: Optional[pd.DataFrame],
    sex_df: pd.DataFrame,
    output_dir: str,
    highlight_sample: str = "",
) -> None:
    """Generate all estimatePloidy-style plots."""
    logger.info("Generating ploidy summary plots ...")
    plot_sex_assignments(sex_df, output_dir, highlight_sample)

    for sex in ("male", "female", "all"):
        for conn in (True, False):
            plot_cn_per_contig_boxplot(
                df,
                output_dir,
                sex_subset=sex,
                connect_samples=conn,
                highlight_sample=highlight_sample,
            )

    if bin_df is not None and not bin_df.empty:
        autosomes = [f"chr{i}" for i in range(1, 23)]
        for chrom in autosomes + ["chrX", "chrY"]:
            if chrom in bin_df["chr"].unique():
                plot_cn_per_bin_chromosome(
                    bin_df,
                    output_dir,
                    chrom,
                    highlight_sample=highlight_sample,
                )
    logger.info("Completed ploidy summary plots.")


def _run_aneuploidy_plots(
    df: pd.DataFrame,
    bin_df: pd.DataFrame,
    loss_df: pd.DataFrame,
    output_dir: str,
    skip_per_sample: bool = False,
    site_data: Optional[dict] = None,
    min_het_alt: int = 3,
) -> None:
    """Generate aneuploidy-detection diagnostic plots."""
    plot_training_loss(loss_df, output_dir)
    plot_training_loss_with_gradient(loss_df, output_dir)
    plot_bin_variance_bias(bin_df, output_dir)
    plot_cn_posterior_entropy(bin_df, output_dir)
    plot_chromosome_cn_heatmap(df, output_dir)
    plot_chromosome_plq_heatmap(df, output_dir)
    plot_binq_genome_profile(bin_df, output_dir)
    plot_parameter_diagnostics(bin_df, output_dir)

    all_vars = bin_df.groupby("sample")["sample_var"].first().values
    sample_groups = {sample: sdf for sample, sdf in bin_df.groupby("sample", sort=False)}
    aneuploid_samples = df[df["is_aneuploid"]]["sample"].unique()
    normal_samples = df[~df["sample"].isin(aneuploid_samples)]["sample"].unique()

    if skip_per_sample:
        logger.info(
            "Skipping per-sample plots (%d aneuploid, %d normal)",
            len(aneuploid_samples),
            len(normal_samples),
        )
        return

    logger.info(
        "Generating per-sample plots (%d aneuploid, %d normal) …",
        len(aneuploid_samples),
        len(normal_samples),
    )

    # Build a sample-name → index mapping for site data lookups
    sample_idx_map: Optional[dict] = None
    if site_data is not None and "sample_ids" in site_data:
        sample_idx_map = {
            str(s): i for i, s in enumerate(site_data["sample_ids"])
        }

    for idx, sid in enumerate(aneuploid_samples, start=1):
        sdata = sample_groups.get(sid)
        if sdata is None:
            continue
        sdf = df[df["sample"] == sid]
        aneu_chrs = [
            (r["chromosome"], r["copy_number"], r["mean_cn_probability"])
            for _, r in sdf[sdf["is_aneuploid"]].iterrows()
        ]
        if idx == 1 or idx % 10 == 0 or idx == len(aneuploid_samples):
            logger.info("Per-sample plots: aneuploid %d/%d (%s)", idx, len(aneuploid_samples), sid)
        plot_sample_with_variance(
            sdata, all_vars, output_dir,
            aneuploid_chrs=aneu_chrs,
            site_data=site_data,
            sample_idx_map=sample_idx_map,
            min_het_alt=min_het_alt,
        )

    for idx, sid in enumerate(normal_samples, start=1):
        sdata = sample_groups.get(sid)
        if sdata is None:
            continue
        if idx == 1 or idx % 10 == 0 or idx == len(normal_samples):
            logger.info("Per-sample plots: normal %d/%d (%s)", idx, len(normal_samples), sid)
        plot_sample_with_variance(
            sdata, all_vars, output_dir,
            site_data=site_data,
            sample_idx_map=sample_idx_map,
            min_het_alt=min_het_alt,
        )


# ── CLI ─────────────────────────────────────────────────────────────────────


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the plot subcommand."""
    p = argparse.ArgumentParser(
        description="Generate diagnostic and summary plots",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("-c", "--chrom-stats", required=True,
                   help="chromosome_stats.tsv (from 'infer')")
    p.add_argument("-o", "--output-dir", required=True)
    p.add_argument("-b", "--bin-stats", default=None,
                   help="bin_stats.tsv.gz (from 'infer')")
    p.add_argument("-t", "--training-loss", default=None,
                   help="training_loss.tsv (from 'infer')")
    p.add_argument("-s", "--sex-assignments", default=None,
                   help="aneuploidy_type_predictions.tsv (from 'call')")
    p.add_argument("--site-data", default=None,
                   help="site_data.npz (from 'preprocess') for per-site AF scatter")
    p.add_argument(
        "--site-af-estimates",
        default=None,
        help="site_af_estimates.tsv.gz (from 'infer') for site-level AF diagnostics",
    )
    p.add_argument("--min-het-alt", type=int, default=3,
                   help="Minimum alt-allele read count to show a site in the AF scatter")
    p.add_argument("--highlight-sample", default="",
                   help="Sample ID to highlight in plots")
    p.add_argument("--skip-per-sample-plots", action="store_true",
                   help="Skip individual sample CN plots")
    p.add_argument("--ppd-bin-summary", default=None,
                   help="ppd_bin_summary.tsv.gz (from 'ppd')")
    p.add_argument("--ppd-chr-summary", default=None,
                   help="ppd_chromosome_summary.tsv (from 'ppd')")
    p.add_argument("--ppd-bin-quality", default=None,
                   help="ppd_bin_quality.tsv (from 'ppd') to add BINQ overlays and diagnostics")
    p.add_argument("--binq-field", choices=_BINQ_FIELD_OPTIONS, default="BINQ20",
                   help="Per-bin quality field to use when --ppd-bin-quality is provided. "
                        "'auto' prefers CALLQ20 when available, otherwise BINQ20")
    p.add_argument("--ignored-bins", default=None,
                   help="ignored_bins.tsv.gz (from 'call') to overlay bins removed during filtered calling")
    return p.parse_args()


def main() -> None:
    """Entry point for ``gatk-sv-ploidy plot``."""
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    df = _apply_plot_depth_columns(pd.read_csv(args.chrom_stats, sep="\t"))
    logger.info("Loaded %d rows for %d samples", len(df), df["sample"].nunique())

    bin_df = None
    if args.bin_stats:
        bin_df = _apply_plot_depth_bin_columns(
            pd.read_csv(args.bin_stats, sep="\t", compression="gzip")
        )
        if args.ppd_bin_quality:
            bin_quality_df = pd.read_csv(args.ppd_bin_quality, sep="\t")
            bin_df = _annotate_binq_values(bin_df, bin_quality_df, args.binq_field)
        if args.ignored_bins:
            ignored_bins_df = pd.read_csv(args.ignored_bins, sep="\t")
            bin_df = _annotate_ignored_bins(bin_df, ignored_bins_df)
    elif args.ignored_bins or args.ppd_bin_quality:
        logger.warning("--ignored-bins and --ppd-bin-quality require --bin-stats for overlay plots — ignoring")

    # ── histograms ──────────────────────────────────────────────────────
    plot_histograms_by_chr_type(df, args.output_dir, args.highlight_sample)
    plot_aneuploid_histograms(df, args.output_dir, args.highlight_sample)
    plot_median_depth_distributions(df, args.output_dir, args.highlight_sample)
    if args.site_af_estimates:
        site_af_df = pd.read_csv(
            args.site_af_estimates,
            sep="\t",
            compression="infer",
        )
        plot_site_af_estimates(site_af_df, args.output_dir)

    # ── ploidy-style plots ──────────────────────────────────────────────
    sex_df = None
    if args.sex_assignments:
        sex_df = pd.read_csv(args.sex_assignments, sep="\t")
    else:
        # Derive minimal sex_df from chromosome stats
        from gatk_sv_ploidy.call import assign_sex_and_aneuploidy_types
        sex_df = assign_sex_and_aneuploidy_types(df)
    sex_df = _apply_plot_depth_sex_columns(sex_df, df)

    _run_ploidy_plots(df, bin_df, sex_df, args.output_dir, args.highlight_sample)

    # ── aneuploidy detection plots ──────────────────────────────────────
    site_data = None
    if args.site_data:
        site_data = load_site_data(args.site_data)

    if args.bin_stats and args.training_loss:
        loss_df = pd.read_csv(args.training_loss, sep="\t")
        _run_aneuploidy_plots(df, bin_df, loss_df, args.output_dir,
                              args.skip_per_sample_plots,
                              site_data=site_data,
                              min_het_alt=args.min_het_alt)
    elif args.bin_stats or args.training_loss:
        logger.warning("Both --bin-stats and --training-loss required for "
                       "aneuploidy detection plots — skipping")

    # ── posterior predictive check plots ─────────────────────────────────
    if args.ppd_bin_summary and args.ppd_chr_summary:
        ppd_bin_df = pd.read_csv(args.ppd_bin_summary, sep="\t",
                                 compression="gzip")
        ppd_chr_df = pd.read_csv(args.ppd_chr_summary, sep="\t")
        run_ppd_plots(ppd_bin_df, ppd_chr_df, args.output_dir,
                      args.highlight_sample)
    elif args.ppd_bin_summary or args.ppd_chr_summary:
        logger.warning("Both --ppd-bin-summary and --ppd-chr-summary "
                       "required for PPD plots — skipping")

    logger.info("Done.")


if __name__ == "__main__":
    main()
