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
from matplotlib.colors import BoundaryNorm, ListedColormap

from gatk_sv_ploidy._util import (
    BINQ_FIELD_OPTIONS,
    add_chromosome_labels,
    baseline_ploidy_label,
    get_chromosome_type,
    resolve_binq_field,
    save_and_close_plot,
)
from gatk_sv_ploidy._plot_detail import (
    plot_cn_per_bin_chromosome,
    plot_cn_per_contig_boxplot,
    plot_sample_with_variance,
    plot_sex_assignments,
)
from gatk_sv_ploidy._plot_ppd import run_ppd_plots
from gatk_sv_ploidy._plot_report import write_plot_report
from gatk_sv_ploidy._plot_style import (
    CHR_TYPE_PALETTE,
    CN_STATE_PALETTE,
    HIGHLIGHT_COLOR,
    TICK_LABEL_SIZE_PT,
    apply_theme,
    double_column_size,
    plot_output_format,
    single_column_size,
)
from gatk_sv_ploidy.data import load_site_data
from gatk_sv_ploidy.infer import load_inference_artifacts
from gatk_sv_ploidy.models import _standardize_multiplicative_bin_factors_numpy

logger = logging.getLogger(__name__)


# ── histogram helpers ───────────────────────────────────────────────────────

_CHR_PALETTE = CHR_TYPE_PALETTE
_DIAGNOSTIC_HISTOGRAM_SPECS: tuple[dict[str, object], ...] = (
    {
        "source_columns": ("sample_overdispersion_map",),
        "title": "Sample Overdispersion by Chromosome Type",
        "xlabel": "Sample Overdispersion",
        "filename": "hist_sample_overdispersion.png",
    },
    {
        "source_columns": ("sample_depth_map",),
        "title": "Sample Depth Map by Chromosome Type",
        "xlabel": "Sample Depth Map",
        "filename": "hist_sample_depth_map.png",
    },
    {
        "source_columns": ("copy_number",),
        "title": "Copy Number by Chromosome Type",
        "xlabel": "Copy Number",
        "filename": "hist_copy_number.png",
    },
    {
        "source_columns": ("plq",),
        "title": "PLQ by Chromosome Type",
        "xlabel": "PLQ",
        "filename": "hist_plq.png",
    },
    {
        "source_columns": ("n_bins_retained",),
        "title": "Bins Retained by Chromosome Type",
        "xlabel": "Bins Retained",
        "filename": "hist_n_bins_retained.png",
    },
    {
        "source_columns": ("frac_bins_retained",),
        "title": "Fraction of Bins Retained by Chromosome Type",
        "xlabel": "Fraction of Bins Retained",
        "filename": "hist_frac_bins_retained.png",
    },
    {
        "source_columns": ("autosomal_baseline_cn",),
        "title": "Autosomal Baseline CN by Chromosome Type",
        "xlabel": "Autosomal Baseline CN",
        "filename": "hist_autosomal_baseline_cn.png",
    },
    {
        "source_columns": ("sample_depth_ratio",),
        "title": "Sample Depth Ratio by Chromosome Type",
        "xlabel": "Sample Depth Ratio",
        "filename": "hist_sample_depth_ratio.png",
    },
    {
        "source_columns": ("global_cn_scale_factor",),
        "title": "Global CN Scale Factor by Chromosome Type",
        "xlabel": "Global CN Scale Factor",
        "filename": "hist_global_cn_scale_factor.png",
    },
)


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
        logger.debug("Using plot-normalized chromosome depth summaries.")
    return out


def _apply_plot_depth_bin_columns(bin_df: pd.DataFrame | None) -> pd.DataFrame | None:
    """Replace raw per-bin observed depth with plot-normalized depth."""
    if bin_df is None:
        return None
    out = bin_df.copy()
    if "plot_depth" in out.columns:
        if "raw_observed_depth" not in out.columns and "observed_depth" in out.columns:
            out["raw_observed_depth"] = out["observed_depth"]
        out["observed_depth"] = out["plot_depth"]
        logger.debug("Using plot-normalized per-bin depth profiles.")
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

    logger.debug(
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

    binq_field = resolve_binq_field(bin_quality_df, binq_field)
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

    logger.debug(
        "Annotated BINQ values for %d sample-bin rows using %s.",
        int(out["binq_value"].notna().sum()) if "binq_value" in out.columns else 0,
        binq_field,
    )
    return out


def _sample_baseline_ploidy_metadata(
    chrom_df: pd.DataFrame,
    sex_df: Optional[pd.DataFrame] = None,
) -> dict[str, dict[str, object]]:
    """Return sample-level baseline ploidy metadata for plotting."""
    metadata: dict[str, dict[str, object]] = {}
    if "sample" not in chrom_df.columns:
        return metadata

    for sample, sdf in chrom_df.groupby("sample", sort=False):
        baseline_cn = 2
        if "autosomal_baseline_cn" in sdf.columns:
            values = pd.to_numeric(
                sdf["autosomal_baseline_cn"],
                errors="coerce",
            ).dropna()
            if not values.empty:
                baseline_cn = int(values.iloc[0])
        metadata[str(sample)] = {
            "baseline_ploidy_type": baseline_ploidy_label(baseline_cn),
            "autosomal_baseline_cn": baseline_cn,
        }

    if sex_df is not None and "sample" in sex_df.columns:
        for _, row in sex_df.iterrows():
            sample = str(row["sample"])
            entry = metadata.setdefault(
                sample,
                {
                    "baseline_ploidy_type": "DIPLOID",
                    "autosomal_baseline_cn": 2,
                },
            )
            if "autosomal_baseline_cn" in row.index and pd.notna(row["autosomal_baseline_cn"]):
                baseline_cn = int(row["autosomal_baseline_cn"])
                entry["autosomal_baseline_cn"] = baseline_cn
                entry["baseline_ploidy_type"] = baseline_ploidy_label(
                    baseline_cn,
                )
            if "baseline_ploidy_type" in row.index and pd.notna(row["baseline_ploidy_type"]):
                entry["baseline_ploidy_type"] = str(row["baseline_ploidy_type"])

    return metadata


def _hist_by_hue(
    data: pd.DataFrame,
    x_col: str,
    hue_col: str,
    title: str,
    xlabel: str,
    output_dir: str,
    filename: str,
    palette: dict,
    bins: int = 30,
    highlight_sample: str = "",
    sample_values: Optional[list] = None,
) -> None:
    """Histogram coloured by *hue_col* with optional sample highlight."""
    fig, ax = plt.subplots(figsize=single_column_size(68))
    sns.histplot(data=data, x=x_col, hue=hue_col, bins=bins, kde=False,
                 palette=palette, alpha=0.6, edgecolor="black", ax=ax,
                 multiple="stack")

    if highlight_sample and sample_values:
        for v in sample_values:
            ax.axvline(v, color=HIGHLIGHT_COLOR, linestyle="--", linewidth=2, alpha=0.8)
        ax.plot([], [], color=HIGHLIGHT_COLOR, linestyle="--", linewidth=2,
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
    """Create curated diagnostic histograms, grouped by chromosome type.

    Args:
        df: ``chromosome_stats.tsv`` DataFrame.
        output_dir: Base output directory.
        highlight_sample: Sample ID to mark with vertical lines.
    """
    df = df.copy()
    df["chr_type"] = df["chromosome"].apply(get_chromosome_type)
    highlight_df = None
    if highlight_sample and highlight_sample in df["sample"].values:
        highlight_df = df[df["sample"] == highlight_sample]

    for spec in _DIAGNOSTIC_HISTOGRAM_SPECS:
        source_col = next(
            (
                col for col in spec["source_columns"]
                if col in df.columns and pd.api.types.is_numeric_dtype(df[col])
            ),
            None,
        )
        if source_col is None:
            continue

        sample_values = None
        if highlight_df is not None:
            sample_values = highlight_df[source_col].tolist()

        _hist_by_hue(
            df,
            source_col,
            "chr_type",
            str(spec["title"]),
            str(spec["xlabel"]),
            output_dir,
            str(spec["filename"]),
            _CHR_PALETTE,
            bins=30,
            highlight_sample=highlight_sample,
            sample_values=sample_values,
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
            fig, ax = plt.subplots(figsize=double_column_size(78))

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
            ax.set_xticklabels(cdf["sample"], rotation=90, fontsize=TICK_LABEL_SIZE_PT)
            plt.tight_layout()
            pdf.savefig(fig, bbox_inches="tight")
            plt.close()

    logger.debug("Saved median depth distributions to %s", pdf_path)


# ── bin variance / bias ────────────────────────────────────────────────────


def plot_bin_variance_bias(bin_df: pd.DataFrame, output_dir: str) -> None:
    """Two-panel plot of per-bin bias and stdev.

    Args:
        bin_df: ``bin_stats.tsv.gz`` DataFrame.
        output_dir: Base output directory.
    """
    unique = bin_df.groupby(["chr", "start", "end"]).first().reset_index()
    unique = unique.sort_values(["chr", "start"])

    fig, axes = plt.subplots(2, 1, figsize=double_column_size(96))
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
                label="Filtered",
                zorder=4,
            )
            axes[1].scatter(
                x[ignored_mask],
                np.sqrt(unique.loc[ignored_mask, "bin_var"].to_numpy(dtype=float)),
                s=ignored_size,
                facecolors="none",
                edgecolors="#FFB300",
                linewidths=1.2,
                label="Filtered",
                zorder=4,
            )

    for a in axes:
        add_chromosome_labels(a, chrs)
        if a.get_legend_handles_labels()[0]:
            a.legend(loc="best", framealpha=0.9)

    plt.tight_layout()
    save_and_close_plot(output_dir, "bin_posteriors.png")


def _coerce_background_factor_matrices(
    map_estimates: dict,
) -> tuple[np.ndarray | None, np.ndarray | None]:
    """Return normalized bin loadings and sample amplitudes when present."""
    if "background_bin_factors" not in map_estimates or "background_sample_factors" not in map_estimates:
        return None, None

    bin_factors = np.asarray(
        map_estimates["background_bin_factors"],
        dtype=np.float64,
    ).squeeze()
    if bin_factors.ndim == 1:
        bin_factors = bin_factors[:, np.newaxis]

    sample_factors = np.asarray(
        map_estimates["background_sample_factors"],
        dtype=np.float64,
    ).squeeze()
    if sample_factors.ndim == 1:
        sample_factors = sample_factors[np.newaxis, :]

    if bin_factors.ndim != 2 or sample_factors.ndim != 2:
        return None, None
    if sample_factors.shape[0] != bin_factors.shape[1]:
        if sample_factors.shape[1] == bin_factors.shape[1]:
            sample_factors = sample_factors.T
        else:
            return None, None

    normalized_bin_factors = bin_factors / np.maximum(
        bin_factors.mean(axis=0, keepdims=True),
        1e-8,
    )
    return normalized_bin_factors, sample_factors


def _coerce_multiplicative_factor_matrices(
    map_estimates: dict,
) -> tuple[np.ndarray | None, np.ndarray | None]:
    """Return standardized multiplicative loadings and sample weights."""
    if (
        "multiplicative_bin_factors" not in map_estimates or
        "multiplicative_sample_factors" not in map_estimates
    ):
        return None, None

    bin_factors = np.asarray(
        map_estimates["multiplicative_bin_factors"],
        dtype=np.float64,
    ).squeeze()
    if bin_factors.ndim == 1:
        bin_factors = bin_factors[:, np.newaxis]

    sample_factors = np.asarray(
        map_estimates["multiplicative_sample_factors"],
        dtype=np.float64,
    ).squeeze()
    if sample_factors.ndim == 1:
        sample_factors = sample_factors[np.newaxis, :]

    if bin_factors.ndim != 2 or sample_factors.ndim != 2:
        return None, None
    if sample_factors.shape[0] != bin_factors.shape[1]:
        if sample_factors.shape[1] == bin_factors.shape[1]:
            sample_factors = sample_factors.T
        else:
            return None, None

    standardized_bin_factors = _standardize_multiplicative_bin_factors_numpy(
        bin_factors,
    )
    return standardized_bin_factors, sample_factors


def _build_factor_plot_specs(map_estimates: dict) -> list[dict[str, object]]:
    """Collect factor matrices and display metadata for available factor types."""
    specs: list[dict[str, object]] = []

    background_bin_loadings, background_sample_factors = _coerce_background_factor_matrices(
        map_estimates,
    )
    if background_bin_loadings is not None and background_sample_factors is not None:
        specs.append(
            {
                "name": "Additive",
                "bin_loadings": background_bin_loadings,
                "sample_weights": background_sample_factors,
                "loading_cmap": "viridis",
                "weight_cmap": "magma",
                "signed": False,
                "loading_label": "Normalized loading",
                "weight_label": "Sample weight",
            }
        )

    multiplicative_bin_loadings, multiplicative_sample_factors = _coerce_multiplicative_factor_matrices(
        map_estimates,
    )
    if (
        multiplicative_bin_loadings is not None and
        multiplicative_sample_factors is not None
    ):
        specs.append(
            {
                "name": "Multiplicative",
                "bin_loadings": multiplicative_bin_loadings,
                "sample_weights": multiplicative_sample_factors,
                "loading_cmap": "coolwarm",
                "weight_cmap": "coolwarm",
                "signed": True,
                "loading_label": "Standardized loading",
                "weight_label": "Sample weight",
            }
        )

    return specs


def _symmetric_color_limits(values: np.ndarray) -> tuple[float, float]:
    """Return symmetric color limits around zero for signed heatmaps."""
    max_abs = float(np.nanmax(np.abs(values)))
    max_abs = max(max_abs, 1e-8)
    return -max_abs, max_abs


def plot_factor_mode_weight_diagnostics(
    map_estimates: dict,
    output_dir: str,
) -> None:
    """Visualize per-mode sample weights for additive and multiplicative factors."""
    factor_specs = _build_factor_plot_specs(map_estimates)
    if not factor_specs:
        logger.debug("Factor weights absent from inference artifacts — skipping mode-weight diagnostics plot.")
        return

    fig, axes = plt.subplots(
        len(factor_specs),
        1,
        figsize=single_column_size(max(55.0, 40.0 * len(factor_specs))),
        squeeze=False,
    )

    for spec_idx, spec in enumerate(factor_specs):
        ax = axes[spec_idx, 0]
        sample_weights = np.asarray(spec["sample_weights"], dtype=np.float64)
        n_factors = sample_weights.shape[0]
        factor_labels = [f"F{idx + 1}" for idx in range(n_factors)]
        box_data = [
            sample_weights[idx, np.isfinite(sample_weights[idx])]
            for idx in range(n_factors)
        ]
        boxplot = ax.boxplot(
            box_data,
            patch_artist=True,
            tick_labels=factor_labels,
            widths=0.65,
        )
        color_map = plt.get_cmap(str(spec["weight_cmap"]))
        for patch_idx, patch in enumerate(boxplot["boxes"]):
            color = color_map((patch_idx + 0.5) / max(n_factors, 1))
            patch.set_facecolor(color)
            patch.set_alpha(0.75)
        for median in boxplot["medians"]:
            median.set_color("black")
            median.set_linewidth(1.5)

        if bool(spec["signed"]):
            ax.axhline(0.0, color="black", linestyle="--", linewidth=1.0, alpha=0.7)
        else:
            ax.set_ylim(bottom=min(0.0, float(np.nanmin(sample_weights)) - 0.01))

        ax.set_title(f"{spec['name']} Mode Weights Across Samples")
        ax.set_xlabel("Mode")
        ax.set_ylabel(str(spec["weight_label"]))
        ax.grid(True, axis="y", alpha=0.25)

    plt.tight_layout()
    save_and_close_plot(output_dir, "factor_mode_weight_diagnostics.png")


def plot_background_factor_diagnostics(
    bin_df: pd.DataFrame,
    map_estimates: dict,
    output_dir: str,
) -> None:
    """Visualize additive and multiplicative factor loadings and sample weights."""
    factor_specs = _build_factor_plot_specs(map_estimates)
    if not factor_specs:
        logger.debug("Factor latents absent from inference artifacts — skipping factor diagnostics plot.")
        return

    unique = (
        bin_df.groupby(["chr", "start", "end"], sort=False)
        .first()
        .reset_index()
    )

    valid_specs: list[dict[str, object]] = []
    for spec in factor_specs:
        bin_loadings = np.asarray(spec["bin_loadings"], dtype=np.float64)
        if len(unique) != bin_loadings.shape[0]:
            logger.warning(
                "Skipping %s factor diagnostics: %d unique bins in bin_stats but %d rows in inference artifacts.",
                spec["name"],
                len(unique),
                bin_loadings.shape[0],
            )
            continue
        valid_specs.append(spec)
    if not valid_specs:
        return

    chromosomes = unique["chr"].astype(str).to_numpy()
    x = np.arange(len(unique))
    chromosome_order = unique["chr"].drop_duplicates().tolist()
    fig_height_mm = max(85.0, 64.0 * len(valid_specs))
    fig, axes = plt.subplots(
        2 * len(valid_specs),
        2,
        figsize=double_column_size(fig_height_mm),
        squeeze=False,
        gridspec_kw={"width_ratios": [3.2, 1.5]},
    )

    for spec_idx, spec in enumerate(valid_specs):
        bin_loadings = np.asarray(spec["bin_loadings"], dtype=np.float64)
        sample_weights = np.asarray(spec["sample_weights"], dtype=np.float64)
        n_factors = bin_loadings.shape[1]
        factor_labels = [f"F{idx + 1}" for idx in range(n_factors)]
        chromosome_means = np.vstack(
            [
                bin_loadings[chromosomes == chrom].mean(axis=0)
                for chrom in chromosome_order
            ]
        )
        sample_strength = np.nanmean(
            np.abs(sample_weights) if bool(spec["signed"]) else sample_weights,
            axis=0,
        )
        sample_order = np.argsort(sample_strength)[::-1]
        row0 = 2 * spec_idx
        row1 = row0 + 1

        loading_kwargs = {}
        weight_kwargs = {}
        if bool(spec["signed"]):
            loading_kwargs["vmin"], loading_kwargs["vmax"] = _symmetric_color_limits(bin_loadings)
            weight_kwargs["vmin"], weight_kwargs["vmax"] = _symmetric_color_limits(sample_weights)

        loadings_im = axes[row0, 0].imshow(
            bin_loadings.T,
            aspect="auto",
            interpolation="nearest",
            cmap=str(spec["loading_cmap"]),
            **loading_kwargs,
        )
        axes[row0, 0].set_title(f"{spec['name']} Bin Loadings")
        axes[row0, 0].set_ylabel("Factor")
        axes[row0, 0].set_yticks(np.arange(n_factors))
        axes[row0, 0].set_yticklabels(factor_labels)
        add_chromosome_labels(axes[row0, 0], chromosomes, x)
        fig.colorbar(
            loadings_im,
            ax=axes[row0, 0],
            fraction=0.03,
            pad=0.02,
            label=str(spec["loading_label"]),
        )

        sample_im = axes[row0, 1].imshow(
            sample_weights[:, sample_order],
            aspect="auto",
            interpolation="nearest",
            cmap=str(spec["weight_cmap"]),
            **weight_kwargs,
        )
        axes[row0, 1].set_title(f"{spec['name']} Sample Weights")
        axes[row0, 1].set_xlabel("Samples sorted by mean weight magnitude")
        axes[row0, 1].set_ylabel("Factor")
        axes[row0, 1].set_yticks(np.arange(n_factors))
        axes[row0, 1].set_yticklabels(factor_labels)
        fig.colorbar(
            sample_im,
            ax=axes[row0, 1],
            fraction=0.05,
            pad=0.03,
            label=str(spec["weight_label"]),
        )

        chromosome_im = axes[row1, 0].imshow(
            chromosome_means.T,
            aspect="auto",
            interpolation="nearest",
            cmap=str(spec["loading_cmap"]),
            **loading_kwargs,
        )
        axes[row1, 0].set_title(f"Mean {spec['name']} Loading by Chromosome")
        axes[row1, 0].set_xlabel("Chromosome")
        axes[row1, 0].set_ylabel("Factor")
        axes[row1, 0].set_xticks(np.arange(len(chromosome_order)))
        axes[row1, 0].set_xticklabels(chromosome_order, rotation=45, ha="right")
        axes[row1, 0].set_yticks(np.arange(n_factors))
        axes[row1, 0].set_yticklabels(factor_labels)
        fig.colorbar(
            chromosome_im,
            ax=axes[row1, 0],
            fraction=0.03,
            pad=0.02,
            label=f"Mean {spec['loading_label'].lower()}",
        )

        box_data = [
            sample_weights[idx, np.isfinite(sample_weights[idx])]
            for idx in range(n_factors)
        ]
        boxplot = axes[row1, 1].boxplot(
            box_data,
            patch_artist=True,
            tick_labels=factor_labels,
            widths=0.65,
        )
        color_map = plt.get_cmap(str(spec["weight_cmap"]))
        for patch_idx, patch in enumerate(boxplot["boxes"]):
            color = color_map((patch_idx + 0.5) / max(n_factors, 1))
            patch.set_facecolor(color)
            patch.set_alpha(0.75)
        for median in boxplot["medians"]:
            median.set_color("black")
            median.set_linewidth(1.5)
        if bool(spec["signed"]):
            axes[row1, 1].axhline(0.0, color="black", linestyle="--", linewidth=1.0, alpha=0.7)
        axes[row1, 1].set_title(f"{spec['name']} Mode Weights")
        axes[row1, 1].set_xlabel("Mode")
        axes[row1, 1].set_ylabel(str(spec["weight_label"]))
        axes[row1, 1].grid(True, axis="y", alpha=0.25)

    plt.tight_layout()
    save_and_close_plot(output_dir, "background_factor_diagnostics.png")
    plot_factor_mode_weight_diagnostics(map_estimates, output_dir)


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

    # Sort samples by sex-chromosome state, then broad CN burden.
    sort_df = pd.DataFrame(index=cn_pivot.index)
    for chrom in ("chrX", "chrY"):
        if chrom in cn_pivot.columns:
            sort_df[chrom] = cn_pivot[chrom]
    sort_df["mean_cn"] = cn_pivot.mean(axis=1)
    sort_df["sample_id"] = cn_pivot.index.astype(str)
    sort_order = sort_df.sort_values(list(sort_df.columns)).index
    cn_pivot = cn_pivot.reindex(sort_order)
    prob_pivot = prob_pivot.reindex(cn_pivot.index)

    n_samples = len(cn_pivot)
    show_labels = n_samples <= 60

    values = np.ma.masked_invalid(cn_pivot.to_numpy(dtype=float))
    cmap = ListedColormap([CN_STATE_PALETTE[i] for i in range(6)])
    cmap.set_bad(color="#F0F3F7")
    norm = BoundaryNorm(np.arange(-0.5, 6.5, 1.0), cmap.N)

    fig, ax = plt.subplots(
        figsize=double_column_size(max(60.0, 22.0 + (2.1 * n_samples)))
    )
    im = ax.imshow(values, aspect="auto", cmap=cmap, norm=norm,
                   interpolation="nearest")
    ax.set_xticks(np.arange(len(chr_order)))
    ax.set_xticklabels([c.replace("chr", "") for c in chr_order], rotation=45)
    if show_labels:
        ax.set_yticks(np.arange(n_samples))
        ax.set_yticklabels(cn_pivot.index, fontsize=TICK_LABEL_SIZE_PT)
    else:
        ax.set_yticks([])
    ax.set_xlabel("Chromosome")
    ax.set_ylabel("Sample")
    ax.set_title("Dominant Copy Number per Sample × Chromosome")
    cbar = plt.colorbar(im, ax=ax, label="Copy Number", ticks=np.arange(6))
    cbar.ax.set_yticklabels([str(i) for i in range(6)])
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

    fig, ax = plt.subplots(
        figsize=double_column_size(max(60.0, 22.0 + (2.1 * n_samples)))
    )
    im = ax.imshow(values, aspect="auto", cmap=cmap,
                   vmin=0, vmax=99, interpolation="nearest")
    ax.set_xticks(np.arange(len(chr_order)))
    ax.set_xticklabels([c.replace("chr", "") for c in chr_order], rotation=45)
    if show_labels:
        ax.set_yticks(np.arange(n_samples))
        ax.set_yticklabels(plq_pivot.index, fontsize=TICK_LABEL_SIZE_PT)
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

    fig, ax = plt.subplots(figsize=double_column_size(58))
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
                label="Filtered",
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

    fig, axes = plt.subplots(2, 1, figsize=single_column_size(96), sharex=True,
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

    fig, ax = plt.subplots(figsize=single_column_size(68))
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
        fig, ax = plt.subplots(figsize=single_column_size(89))
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

    fig, ax = plt.subplots(figsize=single_column_size(89))
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
    map_estimates: Optional[dict] = None,
    sex_df: Optional[pd.DataFrame] = None,
) -> None:
    """Generate aneuploidy-detection diagnostic plots."""
    plot_training_loss_with_gradient(loss_df, output_dir)
    plot_bin_variance_bias(bin_df, output_dir)
    plot_chromosome_cn_heatmap(df, output_dir)
    plot_chromosome_plq_heatmap(df, output_dir)
    plot_binq_genome_profile(bin_df, output_dir)
    if map_estimates is not None:
        plot_background_factor_diagnostics(bin_df, map_estimates, output_dir)

    all_vars = bin_df.groupby("sample")["sample_var"].first().values
    sample_groups = {sample: sdf for sample, sdf in bin_df.groupby("sample", sort=False)}
    baseline_metadata = _sample_baseline_ploidy_metadata(df, sex_df)
    aneuploid_samples = [str(sample) for sample in df[df["is_aneuploid"]]["sample"].unique()]
    non_diploid_samples = [
        sample for sample, metadata in baseline_metadata.items()
        if (
            str(metadata.get("baseline_ploidy_type", "DIPLOID")) != "DIPLOID" or
            int(metadata.get("autosomal_baseline_cn", 2)) != 2
        )
    ]
    highlighted_samples = list(dict.fromkeys(aneuploid_samples + non_diploid_samples))
    normal_samples = [
        str(sample) for sample in sample_groups
        if str(sample) not in set(highlighted_samples)
    ]

    if skip_per_sample:
        logger.info(
            "Skipping per-sample plots (%d highlighted, %d normal)",
            len(highlighted_samples),
            len(normal_samples),
        )
        return

    logger.info(
        "Generating per-sample plots (%d highlighted, %d normal) …",
        len(highlighted_samples),
        len(normal_samples),
    )

    # Build a sample-name → index mapping for site data lookups
    sample_idx_map: Optional[dict] = None
    has_aggregate_af_columns = (
        {"mean_het_af", "n_het_sites"}.issubset(bin_df.columns) and
        bin_df["n_het_sites"].sum() > 0
    )
    if site_data is not None and "sample_ids" in site_data:
        sample_idx_map = {
            str(s): i for i, s in enumerate(site_data["sample_ids"])
        }
    elif site_data is not None:
        logger.warning(
            "site_data.npz lacks sample_ids; per-sample AF panels require "
            "raw counts with unambiguous sample IDs and will be skipped."
        )
    elif has_aggregate_af_columns:
        logger.warning(
            "Per-sample AF panels require --site-data raw counts; aggregate "
            "mean_het_af/n_het_sites columns are ignored."
        )

    for idx, sid in enumerate(highlighted_samples, start=1):
        sdata = sample_groups.get(sid)
        if sdata is None:
            continue
        sdf = df[df["sample"] == sid]
        aneu_chrs = [
            (r["chromosome"], r["copy_number"], r["mean_cn_probability"])
            for _, r in sdf[sdf["is_aneuploid"]].iterrows()
        ]
        metadata = baseline_metadata.get(
            str(sid),
            {"baseline_ploidy_type": "DIPLOID", "autosomal_baseline_cn": 2},
        )
        if idx == 1 or idx % 10 == 0 or idx == len(highlighted_samples):
            logger.debug(
                "Per-sample plots: highlighted %d/%d",
                idx,
                len(highlighted_samples),
            )
        plot_sample_with_variance(
            sdata, all_vars, output_dir,
            aneuploid_chrs=aneu_chrs,
            baseline_ploidy_type=str(metadata.get("baseline_ploidy_type", "DIPLOID")),
            autosomal_baseline_cn=int(metadata.get("autosomal_baseline_cn", 2)),
            site_data=site_data,
            sample_idx_map=sample_idx_map,
            min_het_alt=min_het_alt,
        )

    for idx, sid in enumerate(normal_samples, start=1):
        sdata = sample_groups.get(sid)
        if sdata is None:
            continue
        metadata = baseline_metadata.get(
            str(sid),
            {"baseline_ploidy_type": "DIPLOID", "autosomal_baseline_cn": 2},
        )
        if idx == 1 or idx % 10 == 0 or idx == len(normal_samples):
            logger.debug(
                "Per-sample plots: normal %d/%d",
                idx,
                len(normal_samples),
            )
        plot_sample_with_variance(
            sdata, all_vars, output_dir,
            baseline_ploidy_type=str(metadata.get("baseline_ploidy_type", "DIPLOID")),
            autosomal_baseline_cn=int(metadata.get("autosomal_baseline_cn", 2)),
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
    p.add_argument("--artifacts", default=None,
                   help="inference_artifacts.npz (from 'infer'); auto-detected next to --bin-stats or --chrom-stats when present")
    p.add_argument("-t", "--training-loss", default=None,
                   help="training_loss.tsv (from 'infer')")
    p.add_argument("-s", "--sex-assignments", default=None,
                   help="aneuploidy_type_predictions.tsv (from 'call')")
    p.add_argument("--site-data", default=None,
                   help="site_data.npz (from 'preprocess'); required for raw per-site AF panels")
    p.add_argument(
        "--site-af-estimates",
        default=None,
        help="site_af_estimates.tsv.gz (from 'infer') for site-level AF diagnostics",
    )
    p.add_argument("--min-het-alt", type=int, default=0,
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
                   help="ppd_bin_quality.tsv (from 'ppd') to add per-bin quality overlays and diagnostics")
    p.add_argument("--binq-field", choices=BINQ_FIELD_OPTIONS, default="BINQ20",
                   help="Per-bin quality field to use when --ppd-bin-quality is provided. "
                        "Defaults to BINQ20. 'auto' prefers BINQ20 and falls back to CALLQ20")
    p.add_argument("--ignored-bins", default=None,
                   help="ignored_bins.tsv.gz (from 'call') to overlay bins removed during filtered calling")
    p.add_argument(
        "--pdf",
        action="store_true",
        help="Write figure artifacts as PDF instead of the default PNG. "
             "The multi-page median depth summary remains PDF in either mode.",
    )
    return p.parse_args()


def main() -> None:
    """Entry point for ``gatk-sv-ploidy plot``."""
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    with plot_output_format("pdf" if args.pdf else "png"):
        apply_theme()

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

        artifacts_path = None
        for candidate in (
            args.artifacts,
            os.path.join(os.path.dirname(args.bin_stats), "inference_artifacts.npz") if args.bin_stats else None,
            os.path.join(os.path.dirname(args.chrom_stats), "inference_artifacts.npz"),
        ):
            if candidate and os.path.exists(candidate):
                artifacts_path = candidate
                break
        if args.artifacts and artifacts_path is None:
            logger.warning("Inference artifacts not found at %s — skipping background diagnostics plot.", args.artifacts)
        map_estimates = None
        if artifacts_path is not None:
            map_estimates, _ = load_inference_artifacts(artifacts_path)

        # ── histograms ──────────────────────────────────────────────────────
        plot_histograms_by_chr_type(df, args.output_dir, args.highlight_sample)
        plot_median_depth_distributions(df, args.output_dir, args.highlight_sample)
        if args.site_af_estimates:
            site_af_df = pd.read_csv(
                args.site_af_estimates,
                sep="\t",
                compression="infer",
            )
            plot_site_af_estimates(site_af_df, args.output_dir)

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
                                  min_het_alt=args.min_het_alt,
                                  map_estimates=map_estimates,
                                  sex_df=sex_df)
        elif args.bin_stats or args.training_loss:
            logger.warning("Both --bin-stats and --training-loss required for "
                           "aneuploidy detection plots — skipping")

        write_plot_report(
            args.output_dir,
            chrom_df=df,
            bin_df=bin_df,
            sex_df=sex_df,
            highlight_sample=args.highlight_sample,
        )

        logger.info("Done.")


if __name__ == "__main__":
    main()
