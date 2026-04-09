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

# ── histogram helpers ───────────────────────────────────────────────────────

_CHR_PALETTE = {"Autosomal": "#1f77b4", "chrX": "#ff7f0e", "chrY": "#2ca02c"}


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
        for c in ("copy_number", "mean_cn_probability", "mean_depth",
                   "median_depth", "std_depth"):
            if c in hdf.columns:
                hl_vals[c] = hdf[c].tolist()

    for col in ("copy_number", "mean_cn_probability", "mean_depth",
                "median_depth", "std_depth"):
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

    stats_df = df[["sample", "chromosome", "copy_number", "mean_depth",
                    "median_depth", "mad_depth"]].copy()
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

    for a in axes:
        add_chromosome_labels(a, chrs)

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
    probs = bin_df.groupby(["chr", "start", "end"])[prob_cols].mean().reset_index()
    probs = probs.sort_values(["chr", "start"])

    p = probs[prob_cols].values
    p = np.clip(p, 1e-10, 1.0)
    entropy = -np.sum(p * np.log2(p), axis=1)

    chrs = probs["chr"].values
    x = np.arange(len(probs))

    fig, ax = plt.subplots(figsize=(16, 4))
    ax.scatter(x, entropy, s=3, alpha=0.5, color="darkred", rasterized=True)
    ax.set_ylabel("Shannon Entropy (bits)")
    ax.set_title("CN Posterior Entropy per Bin (averaged across samples)")
    ax.set_xlim([x.min(), x.max()])
    ax.grid(True, alpha=0.3, axis="y")
    add_chromosome_labels(ax, chrs)
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


# ── orchestrators ──────────────────────────────────────────────────────────


def _run_ploidy_plots(
    df: pd.DataFrame,
    bin_df: Optional[pd.DataFrame],
    sex_df: pd.DataFrame,
    output_dir: str,
    highlight_sample: str = "",
) -> None:
    """Generate all estimatePloidy-style plots."""
    plot_sex_assignments(sex_df, output_dir, highlight_sample)

    for sex in ("male", "female", "all"):
        for conn in (True, False):
            plot_cn_per_contig_boxplot(df, output_dir, sex_subset=sex,
                                      connect_samples=conn,
                                      highlight_sample=highlight_sample)

    if bin_df is not None and not bin_df.empty:
        autosomes = [f"chr{i}" for i in range(1, 23)]
        for chrom in autosomes + ["chrX", "chrY"]:
            if chrom in bin_df["chr"].unique():
                plot_cn_per_bin_chromosome(bin_df, output_dir, chrom,
                                          highlight_sample=highlight_sample)


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
    plot_parameter_diagnostics(bin_df, output_dir)

    all_vars = bin_df.groupby("sample")["sample_var"].first().values
    aneuploid_samples = df[df["is_aneuploid"]]["sample"].unique()
    normal_samples = df[~df["sample"].isin(aneuploid_samples)]["sample"].unique()

    if skip_per_sample:
        logger.info("Skipping per-sample plots (%d aneuploid, %d normal)",
                     len(aneuploid_samples), len(normal_samples))
        return

    logger.info("Generating per-sample plots (%d aneuploid, %d normal) …",
                len(aneuploid_samples), len(normal_samples))

    # Build a sample-name → index mapping for site data lookups
    sample_idx_map: Optional[dict] = None
    if site_data is not None and "sample_ids" in site_data:
        sample_idx_map = {
            str(s): i for i, s in enumerate(site_data["sample_ids"])
        }

    for sid in aneuploid_samples:
        sdata = bin_df[bin_df["sample"] == sid].copy()
        sdf = df[df["sample"] == sid]
        aneu_chrs = [
            (r["chromosome"], r["copy_number"], r["mean_cn_probability"])
            for _, r in sdf[sdf["is_aneuploid"]].iterrows()
        ]
        plot_sample_with_variance(
            sdata, all_vars, output_dir,
            aneuploid_chrs=aneu_chrs,
            site_data=site_data,
            sample_idx_map=sample_idx_map,
            min_het_alt=min_het_alt,
        )

    for sid in normal_samples:
        sdata = bin_df[bin_df["sample"] == sid].copy()
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
    return p.parse_args()


def main() -> None:
    """Entry point for ``gatk-sv-ploidy plot``."""
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    df = pd.read_csv(args.chrom_stats, sep="\t")
    logger.info("Loaded %d rows for %d samples", len(df), df["sample"].nunique())

    # ── histograms ──────────────────────────────────────────────────────
    plot_histograms_by_chr_type(df, args.output_dir, args.highlight_sample)
    plot_aneuploid_histograms(df, args.output_dir, args.highlight_sample)
    plot_median_depth_distributions(df, args.output_dir, args.highlight_sample)

    # ── ploidy-style plots ──────────────────────────────────────────────
    bin_df = None
    if args.bin_stats:
        bin_df = pd.read_csv(args.bin_stats, sep="\t", compression="gzip")

    sex_df = None
    if args.sex_assignments:
        sex_df = pd.read_csv(args.sex_assignments, sep="\t")
    else:
        # Derive minimal sex_df from chromosome stats
        from gatk_sv_ploidy.call import assign_sex_and_aneuploidy_types
        sex_df = assign_sex_and_aneuploidy_types(df)

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
