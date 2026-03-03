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
) -> None:
    """Generate aneuploidy-detection diagnostic plots."""
    plot_training_loss(loss_df, output_dir)
    plot_bin_variance_bias(bin_df, output_dir)

    all_vars = bin_df.groupby("sample")["sample_var"].first().values
    aneuploid_samples = df[df["is_aneuploid"]]["sample"].unique()
    normal_samples = df[~df["sample"].isin(aneuploid_samples)]["sample"].unique()

    if skip_per_sample:
        logger.info("Skipping per-sample plots (%d aneuploid, %d normal)",
                     len(aneuploid_samples), len(normal_samples))
        return

    logger.info("Generating per-sample plots (%d aneuploid, %d normal) …",
                len(aneuploid_samples), len(normal_samples))

    for sid in aneuploid_samples:
        sdata = bin_df[bin_df["sample"] == sid].copy()
        sdf = df[df["sample"] == sid]
        aneu_chrs = [
            (r["chromosome"], r["copy_number"], r["mean_cn_probability"])
            for _, r in sdf[sdf["is_aneuploid"]].iterrows()
        ]
        plot_sample_with_variance(sdata, all_vars, output_dir, aneuploid_chrs=aneu_chrs)

    for sid in normal_samples:
        sdata = bin_df[bin_df["sample"] == sid].copy()
        plot_sample_with_variance(sdata, all_vars, output_dir)


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
    p.add_argument("--highlight-sample", default="",
                   help="Sample ID to highlight in plots")
    p.add_argument("--skip-per-sample-plots", action="store_true",
                   help="Skip individual sample CN plots")
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
    if args.bin_stats and args.training_loss:
        loss_df = pd.read_csv(args.training_loss, sep="\t")
        _run_aneuploidy_plots(df, bin_df, loss_df, args.output_dir,
                              args.skip_per_sample_plots)
    elif args.bin_stats or args.training_loss:
        logger.warning("Both --bin-stats and --training-loss required for "
                       "aneuploidy detection plots — skipping")

    logger.info("Done.")


if __name__ == "__main__":
    main()
