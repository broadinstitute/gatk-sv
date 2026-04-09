"""
Posterior predictive check plotting helpers.

Visualisations that compare the posterior predictive distribution to
observed data.  Called from :mod:`gatk_sv_ploidy.plot` when PPD output
files are available.
"""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats as sp_stats

from gatk_sv_ploidy._util import (
    add_chromosome_labels,
    get_chromosome_type,
    save_and_close_plot,
)

logger = logging.getLogger(__name__)

_CHR_PALETTE = {"Autosomal": "#1f77b4", "chrX": "#ff7f0e", "chrY": "#2ca02c"}


# ── observed vs predicted scatter ───────────────────────────────────────────


def plot_obs_vs_predicted(
    ppd_bin_df: pd.DataFrame,
    output_dir: str,
    highlight_sample: str = "",
) -> None:
    """Scatter of observed depth vs PPD mean, coloured by CN state.

    A well-fitting model produces points clustered along the identity line.

    Args:
        ppd_bin_df: ``ppd_bin_summary.tsv.gz`` DataFrame.
        output_dir: Base output directory.
        highlight_sample: Sample to emphasise.
    """
    fig, ax = plt.subplots(figsize=(8, 8))

    cn_colors = {0: "#004D40", 1: "#FFC107", 2: "#1E88E5",
                 3: "#D81B60", 4: "#38006B", 5: "#FF6D00"}

    for cn in sorted(ppd_bin_df["cn_map"].unique()):
        sub = ppd_bin_df[ppd_bin_df["cn_map"] == cn]
        if highlight_sample:
            sub = sub[sub["sample"] != highlight_sample]
        ax.scatter(
            sub["ppd_mean"], sub["observed_depth"],
            s=2, alpha=0.15, color=cn_colors.get(cn, "gray"),
            label=f"CN={cn}", rasterized=True,
        )

    if highlight_sample and highlight_sample in ppd_bin_df["sample"].values:
        hl = ppd_bin_df[ppd_bin_df["sample"] == highlight_sample]
        ax.scatter(
            hl["ppd_mean"], hl["observed_depth"],
            s=8, alpha=0.7, color="magenta", marker="^",
            label=highlight_sample, zorder=5,
        )

    lims = [
        min(ax.get_xlim()[0], ax.get_ylim()[0]),
        max(ax.get_xlim()[1], ax.get_ylim()[1]),
    ]
    ax.plot(lims, lims, "k--", alpha=0.4, linewidth=1, label="Identity")
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.set_xlabel("Predicted Depth (PPD Mean)")
    ax.set_ylabel("Observed Depth")
    ax.set_title("Observed vs Posterior Predictive Depth")
    ax.legend(loc="upper left", fontsize=8, markerscale=3)
    ax.grid(True, alpha=0.3)
    save_and_close_plot(output_dir, "ppd_obs_vs_predicted.png", subdir="ppd")


# ── residual genome-wide profile ───────────────────────────────────────────


def plot_residual_genome_profile(
    ppd_bin_df: pd.DataFrame,
    output_dir: str,
    highlight_sample: str = "",
) -> None:
    """Genome-wide residual profile averaged across samples.

    Bins with large mean residuals indicate systematic model misfit.

    Args:
        ppd_bin_df: ``ppd_bin_summary.tsv.gz`` DataFrame.
        output_dir: Base output directory.
        highlight_sample: Sample to plot individually.
    """
    # Average residual per bin (across samples)
    bin_avg = ppd_bin_df.groupby(["chr", "start", "end"]).agg(
        mean_residual=("residual", "mean"),
        mean_z_score=("z_score", "mean"),
    ).reset_index().sort_values(["chr", "start"])

    chrs = bin_avg["chr"].values
    x = np.arange(len(bin_avg))

    fig, axes = plt.subplots(2, 1, figsize=(16, 8), sharex=True)

    # Panel 1: Mean residual
    ax = axes[0]
    ax.scatter(x, bin_avg["mean_residual"], s=3, alpha=0.5, color="steelblue",
               rasterized=True)
    ax.axhline(0, color="red", linestyle="--", alpha=0.5)
    ax.set_ylabel("Mean Residual (obs − pred)")
    ax.set_title("Genome-wide Residual Profile (averaged across samples)")
    ax.grid(True, alpha=0.3)
    add_chromosome_labels(ax, chrs)

    # Panel 2: Mean z-score
    ax = axes[1]
    ax.scatter(x, bin_avg["mean_z_score"], s=3, alpha=0.5, color="purple",
               rasterized=True)
    ax.axhline(0, color="red", linestyle="--", alpha=0.5)
    ax.axhline(2, color="orange", linestyle=":", alpha=0.4)
    ax.axhline(-2, color="orange", linestyle=":", alpha=0.4)
    ax.set_ylabel("Mean Z-Score")
    ax.grid(True, alpha=0.3)
    add_chromosome_labels(ax, chrs)

    plt.tight_layout()
    save_and_close_plot(output_dir, "ppd_residual_genome.png", subdir="ppd")

    # Highlighted sample profile
    if highlight_sample and highlight_sample in ppd_bin_df["sample"].values:
        hl = ppd_bin_df[ppd_bin_df["sample"] == highlight_sample].sort_values(
            ["chr", "start"]
        )
        chrs_hl = hl["chr"].values
        x_hl = np.arange(len(hl))

        fig, axes = plt.subplots(2, 1, figsize=(16, 8), sharex=True)
        axes[0].scatter(x_hl, hl["residual"], s=3, alpha=0.6, color="steelblue",
                        rasterized=True)
        axes[0].axhline(0, color="red", linestyle="--", alpha=0.5)
        axes[0].set_ylabel("Residual")
        axes[0].set_title(f"Residual Profile: {highlight_sample}")
        axes[0].grid(True, alpha=0.3)
        add_chromosome_labels(axes[0], chrs_hl)

        axes[1].scatter(x_hl, hl["z_score"], s=3, alpha=0.6, color="purple",
                        rasterized=True)
        axes[1].axhline(0, color="red", linestyle="--", alpha=0.5)
        axes[1].axhline(2, color="orange", linestyle=":", alpha=0.4)
        axes[1].axhline(-2, color="orange", linestyle=":", alpha=0.4)
        axes[1].set_ylabel("Z-Score")
        axes[1].grid(True, alpha=0.3)
        add_chromosome_labels(axes[1], chrs_hl)

        plt.tight_layout()
        safe = highlight_sample.replace("/", "_").replace(" ", "_")
        save_and_close_plot(
            output_dir, f"ppd_residual_{safe}.png", subdir="ppd",
        )


# ── calibration histogram ──────────────────────────────────────────────────


def plot_calibration_histogram(
    ppd_bin_df: pd.DataFrame,
    output_dir: str,
) -> None:
    """Histogram of PPD tail probabilities — should be uniform if calibrated.

    Args:
        ppd_bin_df: ``ppd_bin_summary.tsv.gz`` DataFrame.
        output_dir: Base output directory.
    """
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # One-sided tail probability
    ax = axes[0]
    ax.hist(ppd_bin_df["tail_prob"], bins=50, density=True, alpha=0.7,
            edgecolor="black", linewidth=0.3, color="steelblue")
    ax.axhline(1.0, color="red", linestyle="--", alpha=0.6, label="Uniform(0,1)")
    ax.set_xlabel("Tail Probability P(draw ≥ obs)")
    ax.set_ylabel("Density")
    ax.set_title("PPD Calibration: Tail Probabilities")
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Two-sided tail probability
    ax = axes[1]
    ax.hist(ppd_bin_df["two_tail_prob"], bins=50, density=True, alpha=0.7,
            edgecolor="black", linewidth=0.3, color="teal")
    ax.axhline(1.0, color="red", linestyle="--", alpha=0.6, label="Uniform(0,1)")
    ax.set_xlabel("Two-Tail Probability")
    ax.set_ylabel("Density")
    ax.set_title("PPD Calibration: Two-Tailed")
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    save_and_close_plot(output_dir, "ppd_calibration_histogram.png", subdir="ppd")


# ── QQ plot ─────────────────────────────────────────────────────────────────


def plot_qq(
    ppd_bin_df: pd.DataFrame,
    output_dir: str,
) -> None:
    """QQ plot of z-scores against a standard normal, by chromosome type.

    Well-calibrated z-scores follow the identity line on a QQ plot.

    Args:
        ppd_bin_df: ``ppd_bin_summary.tsv.gz`` DataFrame.
        output_dir: Base output directory.
    """
    fig, ax = plt.subplots(figsize=(8, 8))

    df = ppd_bin_df.copy()
    df["chr_type"] = df["chr"].apply(get_chromosome_type)

    for chr_type, color in _CHR_PALETTE.items():
        sub = df[df["chr_type"] == chr_type]
        if sub.empty:
            continue
        z = np.sort(sub["z_score"].values)
        n = len(z)
        theoretical = sp_stats.norm.ppf(np.linspace(0.5 / n, 1 - 0.5 / n, n))
        ax.scatter(theoretical, z, s=2, alpha=0.3, color=color,
                   label=chr_type, rasterized=True)

    lims = [-4, 4]
    ax.plot(lims, lims, "k--", alpha=0.5, linewidth=1, label="Identity")
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.set_xlabel("Theoretical Quantiles (Std Normal)")
    ax.set_ylabel("Observed Z-Score Quantiles")
    ax.set_title("QQ Plot of PPD Z-Scores")
    ax.legend(loc="upper left")
    ax.grid(True, alpha=0.3)
    ax.set_aspect("equal")
    save_and_close_plot(output_dir, "ppd_qq_plot.png", subdir="ppd")


# ── PPD interval coverage ──────────────────────────────────────────────────


def plot_interval_coverage(
    ppd_bin_df: pd.DataFrame,
    output_dir: str,
) -> None:
    """Bar chart comparing nominal vs empirical predictive interval coverage.

    For each nominal level (50%, 80%, 90%, 95%), shows the fraction of
    observations actually falling within the corresponding PPD quantile
    interval.

    Args:
        ppd_bin_df: ``ppd_bin_summary.tsv.gz`` DataFrame.
        output_dir: Base output directory.
    """
    obs = ppd_bin_df["observed_depth"].values
    levels = {
        "50%": ("ppd_q25", "ppd_q75"),
        "90%": ("ppd_q05", "ppd_q95"),
    }

    nominal = []
    empirical = []
    labels = []
    for label, (lo_col, hi_col) in levels.items():
        lo = ppd_bin_df[lo_col].values
        hi = ppd_bin_df[hi_col].values
        inside = np.mean((obs >= lo) & (obs <= hi))
        nom = int(label.rstrip("%")) / 100.0
        nominal.append(nom)
        empirical.append(inside)
        labels.append(label)

    fig, ax = plt.subplots(figsize=(6, 5))
    x = np.arange(len(labels))
    width = 0.35
    ax.bar(x - width / 2, nominal, width, label="Nominal", color="lightgray",
           edgecolor="black")
    ax.bar(x + width / 2, empirical, width, label="Empirical", color="steelblue",
           edgecolor="black")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_xlabel("Predictive Interval")
    ax.set_ylabel("Coverage Fraction")
    ax.set_title("Predictive Interval Coverage")
    ax.legend()
    ax.set_ylim([0, 1.05])
    ax.grid(True, alpha=0.3, axis="y")
    save_and_close_plot(output_dir, "ppd_interval_coverage.png", subdir="ppd")


# ── Bayesian p-value by chromosome ─────────────────────────────────────────


def plot_bayesian_pvalue_heatmap(
    ppd_chr_df: pd.DataFrame,
    output_dir: str,
) -> None:
    """Heatmap of Bayesian p-values (samples × chromosomes).

    Values near 0 or 1 indicate poor model fit for that chromosome/sample.

    Args:
        ppd_chr_df: ``ppd_chromosome_summary.tsv`` DataFrame.
        output_dir: Base output directory.
    """
    pivot = ppd_chr_df.pivot(
        index="sample", columns="chromosome", values="bayesian_pvalue",
    )
    # Sort chromosomes
    chr_order = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    chr_order = [c for c in chr_order if c in pivot.columns]
    pivot = pivot.reindex(columns=chr_order)

    # Sort samples by mean p-value (most extreme first)
    pivot["_mean_dev"] = (pivot.values - 0.5).mean(axis=1)
    pivot = pivot.sort_values("_mean_dev")
    pivot = pivot.drop("_mean_dev", axis=1)

    # Truncate sample labels if too many
    n_samples = len(pivot)
    show_labels = n_samples <= 50

    fig, ax = plt.subplots(figsize=(14, max(4, 0.3 * n_samples)))
    im = ax.imshow(pivot.values, aspect="auto", cmap="RdYlGn",
                   vmin=0, vmax=1, interpolation="nearest")
    ax.set_xticks(np.arange(len(chr_order)))
    ax.set_xticklabels([c.replace("chr", "") for c in chr_order], rotation=45)
    if show_labels:
        ax.set_yticks(np.arange(n_samples))
        ax.set_yticklabels(pivot.index, fontsize=6)
    else:
        ax.set_yticks([])
    ax.set_xlabel("Chromosome")
    ax.set_ylabel("Sample")
    ax.set_title("Bayesian P-Value per Sample × Chromosome")
    plt.colorbar(im, ax=ax, label="Bayesian p-value")
    plt.tight_layout()
    save_and_close_plot(output_dir, "ppd_bayesian_pvalue_heatmap.png", subdir="ppd")


# ── per-chromosome RMSE distribution ───────────────────────────────────────


def plot_chromosome_rmse(
    ppd_chr_df: pd.DataFrame,
    output_dir: str,
) -> None:
    """Boxplot of per-chromosome RMSE across samples.

    Args:
        ppd_chr_df: ``ppd_chromosome_summary.tsv`` DataFrame.
        output_dir: Base output directory.
    """
    chr_order = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    chr_order = [c for c in chr_order if c in ppd_chr_df["chromosome"].unique()]

    fig, ax = plt.subplots(figsize=(14, 6))
    bp_data = []
    for chrom in chr_order:
        vals = ppd_chr_df[ppd_chr_df["chromosome"] == chrom]["rmse"].values
        bp_data.append(vals)

    ax.boxplot(bp_data, positions=np.arange(len(chr_order)), widths=0.6,
               patch_artist=True,
               boxprops=dict(facecolor="steelblue", alpha=0.5),
               medianprops=dict(color="red", linewidth=1.5),
               showfliers=True, flierprops=dict(markersize=3, alpha=0.4))
    ax.set_xticks(np.arange(len(chr_order)))
    ax.set_xticklabels([c.replace("chr", "") for c in chr_order], rotation=45)
    ax.set_xlabel("Chromosome")
    ax.set_ylabel("RMSE")
    ax.set_title("Per-Chromosome RMSE Distribution Across Samples")
    ax.grid(True, alpha=0.3, axis="y")
    plt.tight_layout()
    save_and_close_plot(output_dir, "ppd_chromosome_rmse.png", subdir="ppd")


# ── orchestrator ────────────────────────────────────────────────────────────


def run_ppd_plots(
    ppd_bin_df: pd.DataFrame,
    ppd_chr_df: pd.DataFrame,
    output_dir: str,
    highlight_sample: str = "",
) -> None:
    """Run all posterior predictive check plots.

    Args:
        ppd_bin_df: Per-bin PPD summary DataFrame.
        ppd_chr_df: Per-chromosome PPD summary DataFrame.
        output_dir: Base output directory.
        highlight_sample: Sample ID to highlight.
    """
    logger.info("Generating PPD diagnostic plots …")

    plot_obs_vs_predicted(ppd_bin_df, output_dir, highlight_sample)
    plot_residual_genome_profile(ppd_bin_df, output_dir, highlight_sample)
    plot_calibration_histogram(ppd_bin_df, output_dir)
    plot_qq(ppd_bin_df, output_dir)
    plot_interval_coverage(ppd_bin_df, output_dir)
    plot_bayesian_pvalue_heatmap(ppd_chr_df, output_dir)
    plot_chromosome_rmse(ppd_chr_df, output_dir)

    logger.info("PPD plots complete.")
