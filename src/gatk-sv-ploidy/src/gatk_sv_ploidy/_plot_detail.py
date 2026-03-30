"""
Internal detailed plotting helpers for the plot subcommand.

Contains the more complex visualisations: per-sample CN/depth plots,
per-contig boxplots, per-bin chromosome plots, and the sex-assignment
scatter.  These are called from :mod:`gatk_sv_ploidy.plot`.
"""

from __future__ import annotations

import logging
import os
from typing import Dict, List, Optional, Tuple

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


def _draw_site_af_scatter(
    ax: plt.Axes,
    sample_data: pd.DataFrame,
    x: np.ndarray,
    site_data: Dict[str, np.ndarray],
    sample_idx: int,
    min_het_alt: int = 3,
) -> None:
    """Draw individual site AFs as a scatter on *ax*.

    Each valid site in each bin is plotted as a translucent dot at the bin's
    x-coordinate (with small horizontal jitter).  This shows the raw per-site
    data that the model uses, revealing the bimodal AF distribution
    (homozygous-reference near 0 and heterozygous near 0.5).
    """
    # Build lookup: (chr, start, end) -> npz bin index
    bin_chr = site_data["bin_chr"]
    bin_start = site_data["bin_start"]
    bin_end = site_data["bin_end"]
    bin_lookup: Dict[tuple, int] = {}
    for bi in range(len(bin_chr)):
        bin_lookup[(str(bin_chr[bi]), int(bin_start[bi]), int(bin_end[bi]))] = bi

    sa = site_data["site_alt"]     # (n_bins, max_sites, n_samples)
    st = site_data["site_total"]   # (n_bins, max_sites, n_samples)
    sm = site_data["site_mask"]    # (n_bins, max_sites, n_samples)

    # Determine per-bin x-spacing for jitter scale
    if len(x) > 1:
        # median gap between consecutive bins
        dx = np.median(np.diff(x[x > 0])) if np.any(x > 0) else 0.01
    else:
        dx = 0.01
    jitter_half = dx * 0.3

    rng = np.random.RandomState(42)
    scatter_x: list = []
    scatter_af: list = []

    for row_i, (_, row) in enumerate(sample_data.iterrows()):
        key = (str(row["chr"]), int(row["start"]), int(row["end"]))
        bi = bin_lookup.get(key)
        if bi is None:
            continue

        mask = sm[bi, :, sample_idx]
        if not np.any(mask):
            continue

        alt = sa[bi, mask, sample_idx]
        total = st[bi, mask, sample_idx]
        valid = (total > 0) & (alt >= min_het_alt)
        if not np.any(valid):
            continue

        af = alt[valid] / total[valid]
        n_pts = len(af)
        jitter = rng.uniform(-jitter_half, jitter_half, size=n_pts)
        scatter_x.extend(x[row_i] + jitter)
        scatter_af.extend(af)

    if scatter_x:
        ax.scatter(
            scatter_x, scatter_af,
            s=1, alpha=0.15, color="#00897B",
            rasterized=True, zorder=2, linewidths=0,
        )
        logger.info("  AF scatter: %d site points for %s",
                     len(scatter_x), sample_data["sample"].iloc[0])


def plot_sample_with_variance(
    sample_data: pd.DataFrame,
    all_sample_vars: np.ndarray,
    output_dir: str,
    aneuploid_chrs: Optional[List[Tuple[str, int, float]]] = None,
    site_data: Optional[Dict[str, np.ndarray]] = None,
    sample_idx_map: Optional[Dict[str, int]] = None,
    min_het_alt: int = 3,
) -> None:
    """Multi-panel plot: depth + CN, optional AF scatter, CN posterior, sample-var histogram.

    When *site_data* is provided (the ``.npz`` from preprocessing), each
    individual SNP allele fraction is drawn as a translucent dot, giving a
    faithful view of the per-site data the model actually uses.  Falls back
    to a stem plot of per-bin mean het AF when only bin-level columns are
    available.

    Args:
        sample_data: Bin-level rows for a single sample.
        all_sample_vars: Array of every sample's variance (for histogram).
        output_dir: Base output directory.
        aneuploid_chrs: Optional list of ``(chr, cn, prob)`` tuples for the
            aneuploid chromosomes to highlight.
        site_data: Optional dict from ``site_data.npz`` with keys
            ``site_alt``, ``site_total``, ``site_mask``, ``bin_chr``,
            ``bin_start``, ``bin_end``, ``sample_ids``.
        sample_idx_map: Optional mapping of sample name to column index in
            the site-data arrays.
    """
    name = sample_data["sample"].iloc[0]
    is_aneu = aneuploid_chrs is not None and len(aneuploid_chrs) > 0

    # Determine what AF data is available
    has_site_scatter = False
    _si: Optional[int] = None
    if site_data is not None and sample_idx_map is not None:
        _si = sample_idx_map.get(name)
        if _si is not None:
            has_site_scatter = True

    has_af = has_site_scatter or (
        "mean_het_af" in sample_data.columns
        and "n_het_sites" in sample_data.columns
        and sample_data["n_het_sites"].sum() > 0
    )

    obs = sample_data["observed_depth"].values
    cn_map = sample_data["cn_map"].values
    chrs = sample_data["chr"].values
    probs = sample_data[
        ["cn_prob_0", "cn_prob_1", "cn_prob_2", "cn_prob_3", "cn_prob_4", "cn_prob_5"]
    ].values
    svar = sample_data["sample_var"].iloc[0]

    x = _equal_width_x(chrs, len(obs))

    # Layout: tall depth panel, optional tall AF panel, half-height CN-prob
    # and stdev panels.
    if has_af:
        fig, axes = plt.subplots(
            4, 1, figsize=(8, 9),
            gridspec_kw={"height_ratios": [2, 2, 1, 1]},
        )
        ax_depth, ax_af, ax_cn_prob, ax_hist = axes
    else:
        fig, axes = plt.subplots(
            3, 1, figsize=(8, 6),
            gridspec_kw={"height_ratios": [2, 1, 1]},
        )
        ax_depth, ax_cn_prob, ax_hist = axes
        ax_af = None

    fig.text(0.1 / 8, 0.98, name, fontsize=12, fontweight="bold", va="top", ha="left")

    # Panel 1 — depth + MAP CN
    ax = ax_depth
    ax.plot(x, cn_map, "o", alpha=0.5, markersize=4, color="black",
            label="Predicted bin copy number")
    ax.plot(x, obs, "r-", linewidth=1, alpha=0.7, label="Normalised read depth")
    ax.legend(loc="lower right", bbox_to_anchor=(1.0, 1.02), ncol=2, borderaxespad=0)
    ax.grid(True, axis="y", alpha=1, linestyle="-", linewidth=1)
    ax.set_ylim([-0.5, 5.5])
    ax.set_xlim([x.min(), x.max()])

    # Panel 2 — AF (per-site scatter when available, else bin-mean stems)
    if ax_af is not None:
        ax = ax_af
        if has_site_scatter:
            _draw_site_af_scatter(
                ax, sample_data, x, site_data, _si,
                min_het_alt=min_het_alt,
            )
        else:
            # Fallback: bin-level mean het AF stems
            af_vals = sample_data["mean_het_af"].values
            n_sites_vals = sample_data["n_het_sites"].values
            af_valid = (n_sites_vals > 0) & ~np.isnan(af_vals)
            if af_valid.any():
                markers, stems, base = ax.stem(
                    x[af_valid], af_vals[af_valid],
                )
                plt.setp(stems, linewidth=0.8, alpha=0.2, color="#00897B")
                plt.setp(markers, markersize=2, alpha=0.1, color="#00897B")
                plt.setp(base, linewidth=0.5)
        ax.axhline(0.5, color="red", linestyle="--", alpha=0.1,
                   linewidth=1, label="Expected het AF (0.5)")
        ax.set_ylabel("Allele Fraction")
        ax.set_ylim([-0.05, 1.05])
        ax.set_xlim([x.min(), x.max()])
        ax.legend(loc="lower right", bbox_to_anchor=(1.0, 1.02), borderaxespad=0)
        ax.grid(True, axis="y", alpha=0.3)

    # CN posterior stackplot (half height)
    ax = ax_cn_prob
    ax.stackplot(x, probs.T, labels=[f"CN={i}" for i in range(6)],
                 alpha=0.7, colors=_CN_COLORS)
    ax.set_ylabel("Copy Number Probability")
    ax.legend(loc="lower right", bbox_to_anchor=(1.0, 1.02), ncol=6, borderaxespad=0)
    ax.set_ylim([0, 1])
    ax.set_xlim([x.min(), x.max()])

    # Sample variance histogram (half height)
    ax = ax_hist
    ax.hist(np.sqrt(all_sample_vars), bins=30, alpha=0.7, edgecolor="black",
            linewidth=0.5, color="gray")
    ax.axvline(np.sqrt(svar), color="red", linestyle="--", linewidth=2,
               label=f"This sample: {np.sqrt(svar):.3f}")
    ax.set_xlabel("Sample stdev")
    ax.set_ylabel("Count")
    ax.set_yscale("log")
    ax.legend(loc="lower right", bbox_to_anchor=(1.0, 1.02), borderaxespad=0)
    ax.grid(True, alpha=0.3)

    # Highlight aneuploid chromosomes on spatial panels
    spatial_axes = [ax_depth]
    if ax_af is not None:
        spatial_axes.append(ax_af)
    spatial_axes.append(ax_cn_prob)

    if is_aneu:
        aneu_names = {c for c, _, _ in aneuploid_chrs}
        for cname in aneu_names:
            idx = np.where(chrs == cname)[0]
            if len(idx) > 0:
                for a in spatial_axes:
                    a.axvspan(x[idx[0]], x[idx[-1]], alpha=0.15, color="red", zorder=0)

    for a in spatial_axes:
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

    if plot_df.empty:
        logger.info("Skipping contig boxplot for %s: no samples after filtering", suffix)
        return

    chr_order = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    chr_order = [c for c in chr_order if c in plot_df["chromosome"].unique()]

    if not chr_order:
        logger.info("Skipping contig boxplot for %s: no chromosomes available", suffix)
        return

    wide = plot_df.pivot(index="sample", columns="chromosome", values="median_depth")
    wide = wide.reindex(columns=chr_order)
    cn_wide = plot_df.pivot(index="sample", columns="chromosome", values="copy_number")
    cn_wide = cn_wide.reindex(columns=chr_order)

    fig, ax = plt.subplots(figsize=(12, 6))
    finite_depths = plot_df["median_depth"].to_numpy(dtype=float)
    finite_depths = finite_depths[np.isfinite(finite_depths)]
    y_max = max(4, float(finite_depths.max())) if len(finite_depths) > 0 else 4
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
    bp_data = []
    bp_positions = []
    for i, chrom in enumerate(chr_order):
        cdf = plot_df[plot_df["chromosome"] == chrom]
        vals = cdf["median_depth"].to_numpy(dtype=float)
        cns = cdf["copy_number"].to_numpy()
        valid = np.isfinite(vals)
        vals = vals[valid]
        cns = cns[valid]

        if len(vals) == 0:
            continue

        jx = rng.normal(i, 0.1, size=len(vals))

        exp_lo, exp_hi = _expected_cn_range(chrom)
        colours = [
            "blue" if cn > exp_hi else ("red" if cn < exp_lo else "#838393")
            for cn in cns
        ]
        ax.scatter(jx, vals, s=10, alpha=0.5, c=colours, zorder=4)
        bp_data.append(vals)
        bp_positions.append(i)

    if bp_data:
        ax.boxplot(
            bp_data, positions=bp_positions, widths=0.6,
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
