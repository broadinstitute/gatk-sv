"""
Internal detailed plotting helpers for the plot subcommand.

Contains the more complex visualisations: per-sample CN/depth plots,
per-contig boxplots, per-bin chromosome plots, and the sex-assignment
scatter.  These are called from :mod:`gatk_sv_ploidy.plot`.
"""

from __future__ import annotations

import logging
import os
from typing import List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from gatk_sv_ploidy._util import add_chromosome_labels

logger = logging.getLogger(__name__)

# CN state palette (matches original R / Python colours)
_CN_COLORS = ["#004D40", "#FFC107", "#1E88E5", "#D81B60", "#38006B", "#FF6D00"]


# ── per-sample plots ───────────────────────────────────────────────────────


def _equal_width_x(chr_array: np.ndarray, n: int) -> np.ndarray:
    """Map bin indices to an x-axis where every chromosome has equal width."""
    changes = np.where(chr_array[:-1] != chr_array[1:])[0] + 1
    bounds = np.concatenate([[0], changes, [n]])
    x = np.zeros(n)
    for i in range(len(bounds) - 1):
        s, e = bounds[i], bounds[i + 1]
        x[s:e] = i + np.linspace(0, 1, e - s, endpoint=False)
    return x


def plot_sample_with_variance(
    sample_data: pd.DataFrame,
    all_sample_vars: np.ndarray,
    output_dir: str,
    aneuploid_chrs: Optional[List[Tuple[str, int, float]]] = None,
) -> None:
    """Three-panel plot: depth + CN, CN posterior stack, sample-var histogram.

    Args:
        sample_data: Bin-level rows for a single sample.
        all_sample_vars: Array of every sample's variance (for histogram).
        output_dir: Base output directory.
        aneuploid_chrs: Optional list of ``(chr, cn, prob)`` tuples for the
            aneuploid chromosomes to highlight.
    """
    name = sample_data["sample"].iloc[0]
    is_aneu = aneuploid_chrs is not None and len(aneuploid_chrs) > 0

    obs = sample_data["observed_depth"].values
    cn_map = sample_data["cn_map"].values
    chrs = sample_data["chr"].values
    probs = sample_data[
        ["cn_prob_0", "cn_prob_1", "cn_prob_2", "cn_prob_3", "cn_prob_4", "cn_prob_5"]
    ].values
    svar = sample_data["sample_var"].iloc[0]

    x = _equal_width_x(chrs, len(obs))
    fig, axes = plt.subplots(3, 1, figsize=(8, 8))
    fig.text(0.1 / 8, 0.98, name, fontsize=12, fontweight="bold", va="top", ha="left")

    # Panel 1 — depth + MAP CN
    ax = axes[0]
    ax.plot(x, cn_map, "o", alpha=0.5, markersize=4, color="black",
            label="Predicted bin copy number")
    ax.plot(x, obs, "r-", linewidth=1, alpha=0.7, label="Normalised read depth")
    ax.legend(loc="lower right", bbox_to_anchor=(1.0, 1.02), ncol=2, borderaxespad=0)
    ax.grid(True, axis="y", alpha=1, linestyle="-", linewidth=1)
    ax.set_ylim([-0.5, 5.5])
    ax.set_xlim([x.min(), x.max()])

    # Panel 2 — CN posterior stackplot
    ax = axes[1]
    ax.stackplot(x, probs.T, labels=[f"CN={i}" for i in range(6)],
                 alpha=0.7, colors=_CN_COLORS)
    ax.set_ylabel("Copy Number Probability")
    ax.legend(loc="lower right", bbox_to_anchor=(1.0, 1.02), ncol=6, borderaxespad=0)
    ax.set_ylim([0, 1])
    ax.set_xlim([x.min(), x.max()])

    # Panel 3 — sample variance histogram
    ax = axes[2]
    ax.hist(np.sqrt(all_sample_vars), bins=30, alpha=0.7, edgecolor="black",
            linewidth=0.5, color="gray")
    ax.axvline(np.sqrt(svar), color="red", linestyle="--", linewidth=2,
               label=f"This sample: {np.sqrt(svar):.3f}")
    ax.set_xlabel("Sample stdev")
    ax.set_ylabel("Count")
    ax.set_yscale("log")
    ax.legend(loc="lower right", bbox_to_anchor=(1.0, 1.02), borderaxespad=0)
    ax.grid(True, alpha=0.3)

    # Highlight aneuploid chromosomes
    if is_aneu:
        aneu_names = {c for c, _, _ in aneuploid_chrs}
        for cname in aneu_names:
            idx = np.where(chrs == cname)[0]
            if len(idx) > 0:
                for a in axes[:2]:
                    a.axvspan(x[idx[0]], x[idx[-1]], alpha=0.15, color="red", zorder=0)

    for a in axes[:2]:
        add_chromosome_labels(a, chrs, x_transformed=x)

    plt.tight_layout()
    dest = os.path.join(output_dir, "sample_plots")
    os.makedirs(dest, exist_ok=True)
    safe = name.replace("/", "_").replace(" ", "_")
    fig.savefig(os.path.join(dest, f"{safe}.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)


# ── per-contig boxplot ──────────────────────────────────────────────────────


def plot_cn_per_contig_boxplot(
    df: pd.DataFrame,
    output_dir: str,
    sex_subset: str = "all",
    connect_samples: bool = True,
    highlight_sample: str = "",
) -> None:
    """Boxplot with jitter showing normalised depth per contig.

    Args:
        df: ``chromosome_stats.tsv`` DataFrame.
        output_dir: Base output directory.
        sex_subset: ``'male'``, ``'female'``, or ``'all'``.
        connect_samples: Draw lines connecting the same sample across contigs.
        highlight_sample: Sample ID to overlay in magenta.
    """
    plot_df = df.copy()

    # Filter by inferred sex
    if sex_subset == "male":
        males = _samples_by_chrX(df, lambda cn: cn < 2)
        plot_df = plot_df[plot_df["sample"].isin(males)]
        suffix = "chrX_lessThan_2copies"
    elif sex_subset == "female":
        females = _samples_by_chrX(df, lambda cn: cn >= 2)
        plot_df = plot_df[plot_df["sample"].isin(females)]
        suffix = "chrX_atLeast_2copies"
    else:
        suffix = "all_samples"

    chr_order = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    chr_order = [c for c in chr_order if c in plot_df["chromosome"].unique()]

    wide = plot_df.pivot(index="sample", columns="chromosome", values="median_depth")
    wide = wide.reindex(columns=chr_order)
    cn_wide = plot_df.pivot(index="sample", columns="chromosome", values="copy_number")
    cn_wide = cn_wide.reindex(columns=chr_order)

    fig, ax = plt.subplots(figsize=(12, 6))
    y_max = max(4, plot_df["median_depth"].max())
    ax.set_ylim(0, y_max)

    # Alternating chromosome shading
    for i in range(0, len(chr_order), 2):
        ax.axvspan(i - 0.5, i + 0.5, color="lightgray", alpha=0.3, zorder=0)

    for y in np.arange(0, y_max + 0.5, 0.5):
        style = ("-", 0.5, 0.8) if y == int(y) else (":", 0.3, 0.5)
        ax.axhline(y, color="gray", linestyle=style[0], alpha=style[1],
                   linewidth=style[2], zorder=1)
    ax.axhline(2, color="black", linestyle="-", linewidth=1.5, zorder=2)

    # Lines connecting samples
    if connect_samples:
        for sample in wide.index:
            if sample == highlight_sample:
                continue
            ax.plot(np.arange(len(chr_order)), wide.loc[sample].values,
                    color="gray", alpha=0.2, linewidth=0.4, zorder=3)

    # Jittered scatter with gain/loss colouring
    rng = np.random.RandomState(0)
    for i, chrom in enumerate(chr_order):
        cdf = plot_df[plot_df["chromosome"] == chrom]
        vals = cdf["median_depth"].values
        cns = cdf["copy_number"].values
        jx = rng.normal(i, 0.1, size=len(vals))

        exp_lo, exp_hi = _expected_cn_range(chrom)
        colours = [
            "blue" if cn > exp_hi else ("red" if cn < exp_lo else "#838393")
            for cn in cns
        ]
        ax.scatter(jx, vals, s=10, alpha=0.5, c=colours, zorder=4)

    # Boxplots
    bp_data = [plot_df[plot_df["chromosome"] == c]["median_depth"].values
               for c in chr_order]
    ax.boxplot(
        bp_data, positions=np.arange(len(chr_order)), widths=0.6,
        patch_artist=False, showcaps=True, showfliers=False,
        boxprops=dict(linewidth=1.5, color="black"),
        whiskerprops=dict(linewidth=1.5, color="black"),
        medianprops=dict(linewidth=1.5, color="black"),
        capprops=dict(linewidth=0), zorder=5,
    )

    if highlight_sample and highlight_sample in wide.index:
        ax.plot(np.arange(len(chr_order)), wide.loc[highlight_sample].values,
                color="magenta", linewidth=2, zorder=6, label=" " + highlight_sample)
        ax.legend(loc="upper right", framealpha=0.9)

    ax.set_xticks(np.arange(len(chr_order)))
    ax.set_xticklabels([c.replace("chr", "") for c in chr_order], rotation=90)
    ax.set_xlabel("Chromosome", fontsize=12)
    ax.set_ylabel("Normalised Depth", fontsize=12)

    n = wide.shape[0]
    ni = int(wide.isnull().any(axis=1).sum())
    ax.text(0.02, 0.98, f"N={n} samples ({ni} incomplete)",
            transform=ax.transAxes, va="top", ha="left", fontsize=10,
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.8))

    plt.tight_layout()
    contour = "with_contours" if connect_samples else "no_contours"
    dest = os.path.join(output_dir, "raw_depth_contig")
    os.makedirs(dest, exist_ok=True)
    fig.savefig(
        os.path.join(dest, f"estimated_CN_per_contig.{suffix}.{contour}.png"),
        dpi=300, bbox_inches="tight",
    )
    plt.close(fig)
    logger.info("Created contig boxplot (%s, %s)", suffix, contour)


# ── per-bin-per-chromosome plot ─────────────────────────────────────────────


def plot_cn_per_bin_chromosome(
    bin_df: pd.DataFrame,
    output_dir: str,
    chromosome: str,
    sex_subset: Optional[str] = None,
    highlight_sample: str = "",
) -> None:
    """Boxplot with jitter showing depth per bin for a single chromosome.

    Args:
        bin_df: ``bin_stats.tsv.gz`` DataFrame.
        output_dir: Base output directory.
        chromosome: e.g. ``'chr1'``, ``'chrX'``.
        sex_subset: ``'male'``, ``'female'``, or ``None``.
        highlight_sample: Sample ID to overlay.
    """
    cdf = bin_df[bin_df["chr"] == chromosome].copy()
    if len(cdf) == 0:
        return

    plot_col = "observed_depth" if chromosome in ("chrX", "chrY") or sex_subset is None else "cn_map"

    if sex_subset == "male":
        males = _samples_by_chrX_bin(bin_df, lambda cn: cn < 2)
        cdf = cdf[cdf["sample"].isin(males)]
        suffix = "chrX_lessThan_2copies"
    elif sex_subset == "female":
        females = _samples_by_chrX_bin(bin_df, lambda cn: cn >= 2)
        cdf = cdf[cdf["sample"].isin(females)]
        suffix = "chrX_atLeast_2copies"
    else:
        suffix = "all_samples"

    cdf = cdf.sort_values("start")
    bins = cdf[["start", "end"]].drop_duplicates().sort_values("start").reset_index(drop=True)
    bins["bin_idx"] = bins.index
    cdf = cdf.merge(bins[["start", "end", "bin_idx"]], on=["start", "end"], how="left")

    fig, ax = plt.subplots(figsize=(12, 6))
    ax.set_ylim(0, 5)

    for y in np.arange(0, 5.5, 0.5):
        style = ("-", 0.5, 0.8) if y == int(y) else (":", 0.3, 0.5)
        ax.axhline(y, color="gray", linestyle=style[0], alpha=style[1], linewidth=style[2])
    ax.axhline(2, color="black", linestyle="-", linewidth=1.5)

    # Sample lines
    rng = np.random.RandomState(0)
    for sample in cdf["sample"].unique():
        if sample == highlight_sample:
            continue
        sdf = cdf[cdf["sample"] == sample].sort_values("bin_idx")
        ax.plot(sdf["bin_idx"], sdf[plot_col], color="gray", alpha=0.2, linewidth=2)

    # Jittered scatter
    bin_indices = sorted(cdf["bin_idx"].unique())
    for bi in bin_indices:
        vals = cdf[cdf["bin_idx"] == bi][plot_col].values
        jx = rng.normal(bi, 0.1, size=len(vals))
        ax.scatter(jx, vals, s=10, alpha=0.5, color="#838393")

    # Boxplots
    bp_data, bp_pos = [], []
    for bi in bin_indices:
        vals = cdf[cdf["bin_idx"] == bi][plot_col].values
        if len(vals) > 0:
            bp_data.append(vals)
            bp_pos.append(bi)
    if bp_data:
        ax.boxplot(
            bp_data, positions=bp_pos, widths=0.6,
            patch_artist=False, showcaps=True, showfliers=False,
            boxprops=dict(linewidth=1.0, color="black"),
            whiskerprops=dict(linewidth=1.0, color="black"),
            medianprops=dict(linewidth=1.0, color="black"),
            capprops=dict(linewidth=0),
        )

    if highlight_sample:
        sdf = cdf[cdf["sample"] == highlight_sample].sort_values("bin_idx")
        if len(sdf) > 0:
            ax.plot(sdf["bin_idx"], sdf[plot_col], color="magenta", linewidth=2,
                    label=" " + highlight_sample)
            ax.legend(loc="upper right", framealpha=0.9)

    ylabel = "Normalised Depth" if plot_col == "observed_depth" else "Estimated Copy Number"
    ax.set_xlabel(f"{chromosome} Position (Binned)", fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_xticks([])

    n = cdf["sample"].nunique()
    ax.text(0.02, 0.98, f"N={n} samples", transform=ax.transAxes, va="top",
            ha="left", fontsize=10,
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.8))

    plt.tight_layout()
    dest = os.path.join(output_dir, "raw_depth_binned")
    os.makedirs(dest, exist_ok=True)
    fig.savefig(
        os.path.join(dest, f"estimated_CN_per_bin.{suffix}.{chromosome}.png"),
        dpi=300, bbox_inches="tight",
    )
    plt.close(fig)
    logger.info("Created per-bin plot: %s (%s)", chromosome, suffix)


# ── sex-assignment scatter ──────────────────────────────────────────────────

_SEX_COLORS = {
    "MALE": "#00BFF4",
    "FEMALE": "#fd8eff",
    "TURNER": "#e02006",
    "TRIPLE_X": "#7B2AB3",
    "KLINEFELTER": "#FF6A09",
    "JACOBS": "#29840f",
    "OTHER": "#8F1336",
}


def plot_sex_assignments(
    sex_df: pd.DataFrame,
    output_dir: str,
    highlight_sample: str = "",
) -> None:
    """Scatter of chrX vs chrY normalised depth, coloured by sex.

    Args:
        sex_df: DataFrame with ``sample``, ``sex``, ``chrX_depth``,
            ``chrY_depth``.
        output_dir: Base output directory.
        highlight_sample: Sample ID to mark with a triangle.
    """
    logger.info("Generating sex-assignment scatter …")
    fig, ax = plt.subplots(figsize=(8, 8))
    lim = 3
    ax.set_xlim(-0.1, lim + 0.1)
    ax.set_ylim(-0.1, lim + 0.1)

    for i in range(lim + 1):
        ax.axhline(i, color="gray", linestyle=":", alpha=0.5)
        ax.axvline(i, color="gray", linestyle=":", alpha=0.5)

    for xc in range(lim + 1):
        for yc in range(lim + 1):
            ax.text(xc, yc, "X" * xc + "Y" * yc, fontsize=16, fontweight="bold",
                    color="lightgray", ha="center", va="center", alpha=0.8)

    for assignment in sex_df["sex"].unique():
        sub = sex_df[sex_df["sex"] == assignment]
        non_hl = sub[sub["sample"] != highlight_sample]
        if len(non_hl):
            ax.scatter(non_hl["chrX_depth"], non_hl["chrY_depth"],
                       c=_SEX_COLORS.get(assignment, "#8F1336"),
                       s=50, alpha=0.7, label=assignment, edgecolors="none")
        hl = sub[sub["sample"] == highlight_sample]
        if len(hl):
            ax.scatter(hl["chrX_depth"], hl["chrY_depth"],
                       c=_SEX_COLORS.get(assignment, "#8F1336"),
                       s=100, alpha=1.0, marker="^", edgecolors="black", linewidths=2)

    if highlight_sample and highlight_sample in sex_df["sample"].values:
        row = sex_df[sex_df["sample"] == highlight_sample].iloc[0]
        ax.axhline(row["chrY_depth"], color="gray", linestyle="-", linewidth=2,
                   alpha=0.7, label=" " + highlight_sample)
        ax.axvline(row["chrX_depth"], color="gray", linestyle="-", linewidth=2,
                   alpha=0.7)

    ax.legend(loc="upper right", framealpha=0.9)
    ax.set_xlabel("chrX Normalised Depth", fontsize=12)
    ax.set_ylabel("chrY Normalised Depth", fontsize=12)

    plt.tight_layout()
    fig.savefig(os.path.join(output_dir, "sex_assignments.png"), dpi=300,
                bbox_inches="tight")
    plt.close(fig)
    logger.info("Saved sex_assignments.png")


# ── private helpers ─────────────────────────────────────────────────────────


def _expected_cn_range(chrom: str) -> Tuple[int, int]:
    """Return (expected_min_cn, expected_max_cn) for *chrom*."""
    if chrom == "chrX":
        return 1, 2
    if chrom == "chrY":
        return 0, 1
    return 2, 2


def _samples_by_chrX(df: pd.DataFrame, predicate) -> list:
    """Return sample IDs whose chrX CN satisfies *predicate*."""
    out = []
    for sample, sdf in df.groupby("sample"):
        cn = sdf[sdf["chromosome"] == "chrX"]["copy_number"].values
        if len(cn) > 0 and predicate(cn[0]):
            out.append(sample)
    return out


def _samples_by_chrX_bin(bin_df: pd.DataFrame, predicate) -> list:
    """Like :func:`_samples_by_chrX` but from bin-level data."""
    wide = bin_df.pivot_table(index="sample", columns="chr",
                              values="cn_map", aggfunc="mean")
    if "chrX" not in wide.columns:
        return list(wide.index)
    return [s for s in wide.index if predicate(wide.loc[s, "chrX"])]
