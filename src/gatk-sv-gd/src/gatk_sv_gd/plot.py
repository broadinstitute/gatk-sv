"""
GD CNV Plotting

Generate visualisation plots for GD copy-number variant calls produced by
gd_cnv_call.py (or the combined plot_gd_cnv_output.py).  Shows depth
profiles at GD loci with annotations for segmental duplications, gene
transcripts, breakpoint intervals, and copy-number calls.

Flanking regions are compressed in display space by ``--flank-scale``
(default 0.10) so the locus body dominates the figure width.

Usage:
    python gd_cnv_plot.py \
        --calls gd_cnv_calls.tsv.gz \
        --cn-posteriors cn_posteriors.tsv.gz \
        --gd-table gd_table.tsv \
        --gtf genes.gtf.gz \
        --segdup-bed segdups.bed.gz \
        --output-dir plots/
"""

import argparse
import os
import sys
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from gatk_sv_gd._util import TeeStream, setup_logging
from gatk_sv_gd.annotations import (
    FlankCompressor,
    GapsAnnotation,
    GTFParser,
    SegDupAnnotation,
    _build_gap_positions,
    _build_heatmap_matrix,
    _build_ploidy_lookup,
    _depths_with_gaps,
    _draw_overview_column,
    _sort_samples_by_ploidy,
    draw_annotations_panel,
    get_sample_columns,
)
from gatk_sv_gd.models import GDLocus, GDTable
from gatk_sv_gd.highres import normalize_highres_bins, query_highres_bins


# =============================================================================
# Raw-depth helpers
# =============================================================================


def estimate_lowres_bin_size(df: pd.DataFrame) -> float:
    """Return the median bin size (bp) across all rows of *df*."""
    return float(np.median((df["End"] - df["Start"]).values))


def _build_raw_region_df(
    locus: GDLocus,
    region_start: int,
    region_end: int,
    raw_counts_df: Optional[pd.DataFrame],
    sample_cols: List[str],
    raw_sample_medians: Dict[str, float],
    lowres_median_bin_size: float,
    highres_path: Optional[str] = None,
    size_threshold_factor: float = 10.0,
) -> Optional[pd.DataFrame]:
    """
    Build a normalised raw-depth DataFrame for a locus's visible region,
    optionally substituting high-resolution bins for small intervals/flanks.

    For each GD body interval and each flank, the region's genomic span is
    compared against ``size_threshold_factor × lowres_median_bin_size``.
    When the span is smaller than that threshold **and** *highres_path* is
    provided, bins for that sub-region are queried on demand via pysam
    (reusing :func:`~gatk_sv_gd.highres.query_highres_bins` and
    :func:`~gatk_sv_gd.highres.normalize_highres_bins` from the preprocess
    module), and the corresponding low-resolution bins are spliced out.
    Larger regions – and small regions for which the high-res query returns
    no data – retain their low-resolution bins.

    All returned depths are normalised to the 2.0 = diploid reference scale
    used by the processed ``depth_df``.

    Args:
        locus: GDLocus defining body intervals and chromosome.
        region_start: Leftmost coordinate of the visible region (derived
            from the processed depth_df bin extents).
        region_end: Rightmost coordinate of the visible region.
        raw_counts_df: Full low-resolution raw count matrix held in memory.
            May be ``None`` when no low-res file was loaded.
        sample_cols: Sample column names to include in the output.
        raw_sample_medians: Per-sample autosomal median raw counts from the
            low-resolution file, used to normalise both low-res and
            high-res bins.
        lowres_median_bin_size: Median bin size (bp) of the low-resolution
            file; passed to
            :func:`~gatk_sv_gd.highres.normalize_highres_bins` for
            bin-size-ratio correction.
        highres_path: Path to a bgzipped, tabix-indexed high-resolution
            read-count file (.tsv.gz + .tbi).  ``None`` disables high-res
            substitution entirely.
        size_threshold_factor: Body intervals and flanks whose genomic span
            is less than this multiple of *lowres_median_bin_size* are
            served from the high-res file (default 10).

    Returns:
        Normalised DataFrame with columns ``Chr``, ``Start``, ``End`` and
        one column per sample, sorted by ``Start``.  Returns ``None`` when
        no usable data is available.
    """
    chrom = locus.chrom
    threshold_bp = size_threshold_factor * lowres_median_bin_size

    # Restrict to samples that have a valid per-sample median.
    common_samples = [
        s for s in sample_cols
        if s in raw_sample_medians and raw_sample_medians[s] > 0
    ]
    if not common_samples:
        return None

    # Array form required by normalize_highres_bins.
    column_medians_arr = np.array([raw_sample_medians[s] for s in common_samples])

    # --- Step 1: normalised low-res base covering the full visible region ---
    base_df = pd.DataFrame()
    if raw_counts_df is not None:
        mask = (
            (raw_counts_df["Chr"] == chrom)
            & (raw_counts_df["End"] > region_start)
            & (raw_counts_df["Start"] < region_end)
        )
        raw_sub = raw_counts_df[mask].copy().sort_values("Start")
        if len(raw_sub) > 0:
            usable = [s for s in common_samples if s in raw_sub.columns]
            if usable:
                meds = np.array([raw_sample_medians[s] for s in usable])
                normed_vals = 2.0 * raw_sub[usable].values.astype(float) / meds
                base_df = pd.concat(
                    [
                        raw_sub[["Chr", "Start", "End"]].reset_index(drop=True),
                        pd.DataFrame(normed_vals, columns=usable),
                    ],
                    axis=1,
                )

    if highres_path is None:
        return base_df if len(base_df) > 0 else None

    # --- Step 2: identify small regions to replace with high-res bins ---
    regions_to_check: List[Tuple[int, int, str]] = list(locus.get_intervals())
    if region_start < locus.start:
        regions_to_check.append((region_start, locus.start, "left_flank"))
    if region_end > locus.end:
        regions_to_check.append((locus.end, region_end, "right_flank"))

    highres_parts: List[pd.DataFrame] = []
    excluded_ranges: List[Tuple[int, int]] = []

    for r_start, r_end, r_name in regions_to_check:
        region_size = r_end - r_start
        if region_size >= threshold_bp:
            continue  # large enough — low-res is adequate
        hr_raw = query_highres_bins(
            highres_path, chrom, r_start, r_end, common_samples,
        )
        if len(hr_raw) == 0:
            print(f"    [high-res plot] {r_name}: no bins returned, keeping low-res")
            continue
        normed = normalize_highres_bins(
            hr_raw, common_samples, column_medians_arr, lowres_median_bin_size,
        )
        keep_cols = ["Chr", "Start", "End"] + [
            s for s in common_samples if s in normed.columns
        ]
        highres_parts.append(normed.reset_index()[keep_cols])
        excluded_ranges.append((r_start, r_end))
        print(
            f"    [high-res plot] {r_name} ({region_size:,} bp "
            f"< {threshold_bp:,.0f} bp threshold): {len(hr_raw)} high-res bins"
        )

    if not highres_parts:
        return base_df if len(base_df) > 0 else None

    # --- Step 3: splice out low-res bins covered by high-res regions ---
    if len(base_df) > 0:
        bin_mids = (base_df["Start"] + base_df["End"]) / 2
        keep = pd.Series(True, index=base_df.index)
        for ex_s, ex_e in excluded_ranges:
            keep &= ~((bin_mids >= ex_s) & (bin_mids < ex_e))
        base_df = base_df[keep]

    pieces = ([base_df] if len(base_df) > 0 else []) + highres_parts
    combined = pd.concat(pieces, axis=0, ignore_index=True)
    combined = (
        combined
        .drop_duplicates(subset=["Chr", "Start", "End"])
        .sort_values("Start")
        .reset_index(drop=True)
    )
    return combined if len(combined) > 0 else None


# =============================================================================
# Top-level Plot Functions
# =============================================================================


def plot_locus_overview(
    locus: GDLocus,
    calls_df: pd.DataFrame,
    depth_df: pd.DataFrame,
    gtf: Optional[GTFParser],
    segdup: Optional[SegDupAnnotation],
    output_dir: str,
    padding: int = 50000,
    ploidy_df: Optional[pd.DataFrame] = None,
    ploidy_lookup: Optional[Dict[Tuple[str, str], int]] = None,
    min_gene_label_spacing: float = 0.05,
    raw_counts_df: Optional[pd.DataFrame] = None,
    raw_sample_medians: Optional[Dict[str, float]] = None,
    gaps: Optional[GapsAnnotation] = None,
    flank_scale: float = 0.20,
    lowres_median_bin_size: Optional[float] = None,
    highres_path: Optional[str] = None,
):
    """
    Create overview plot for a GD locus showing all samples.

    When *raw_counts_df* is provided the figure has two columns:
      - **Left**: normalised raw depth (before filtering/rebinning)
      - **Right**: processed rebinned counts from the model
    """
    if not locus.breakpoints:
        print(f"  Warning: No breakpoints defined for locus {locus.cluster}, skipping")
        return

    chrom = locus.chrom
    mask = (
        (depth_df["Cluster"] == locus.cluster) &
        (depth_df["Chr"] == chrom)
    )
    region_df = depth_df[mask].copy().sort_values("Start")

    if len(region_df) == 0:
        print(f"  Warning: No depth data found for locus {locus.cluster} "
              f"({chrom}:{locus.start:,}-{locus.end:,}), skipping plot")
        return

    region_start = int(region_df["Start"].min())
    region_end = int(region_df["End"].max())

    # Build coordinate transform
    xform = FlankCompressor(region_start, region_end,
                            locus.start, locus.end, flank_scale=flank_scale)

    sample_cols = get_sample_columns(region_df)
    if ploidy_lookup is None:
        ploidy_lookup = _build_ploidy_lookup(ploidy_df)

    carriers = set(calls_df[
        (calls_df["cluster"] == locus.cluster) &
        (calls_df["is_carrier"])
    ]["sample"].unique())

    carrier_cols = _sort_samples_by_ploidy(
        [s for s in sample_cols if s in carriers], chrom, ploidy_lookup)
    non_carrier_cols = _sort_samples_by_ploidy(
        [s for s in sample_cols if s not in carriers], chrom, ploidy_lookup)
    n_carriers = len(carrier_cols)
    n_non_carriers = len(non_carrier_cols)

    all_ploidies = sorted(set(
        ploidy_lookup.get((s, chrom), 2) for s in sample_cols
    ))
    n_ploidy_panels = len(all_ploidies)

    # ---- Build raw-normalised DataFrame for the left column (if available) ----
    have_raw = (
        (raw_counts_df is not None or highres_path is not None)
        and raw_sample_medians is not None
        and lowres_median_bin_size is not None
    )
    raw_region_df: Optional[pd.DataFrame] = None
    if have_raw:
        raw_region_df = _build_raw_region_df(
            locus, region_start, region_end,
            raw_counts_df, sample_cols, raw_sample_medians,
            lowres_median_bin_size, highres_path=highres_path,
        )
        if raw_region_df is None or len(raw_region_df) == 0:
            have_raw = False
            raw_region_df = None

    # ---- figure layout ----
    figure_scale = 1.0
    n_cols = 2 if have_raw else 1
    carrier_height = min(6, max(1, n_carriers * 0.15)) if n_carriers > 0 else 0.5
    non_carrier_height = min(6, max(1, n_non_carriers * 0.15)) if n_non_carriers > 0 else 0.5
    depth_panel_height = 2.0
    indiv_depth_panel_height = 2.0
    fig_height = (3 + carrier_height + non_carrier_height
                  + n_ploidy_panels * depth_panel_height
                  + indiv_depth_panel_height)

    height_ratios = (
        [1, carrier_height, non_carrier_height]
        + [depth_panel_height] * n_ploidy_panels
        + [indiv_depth_panel_height]
    )
    n_rows = 3 + n_ploidy_panels + 1
    fig_width = 16 * n_cols

    fig, axes = plt.subplots(
        n_rows, n_cols, figsize=(fig_width * figure_scale, fig_height * figure_scale),
        gridspec_kw={"height_ratios": height_ratios},
        squeeze=False,
    )

    locus_label = f"{locus.cluster} ({chrom}:{locus.start:,}-{locus.end:,})"

    if have_raw:
        _draw_overview_column(
            [axes[r, 0] for r in range(n_rows)],
            raw_region_df, locus, calls_df,
            region_start, region_end,
            [s for s in carrier_cols if s in raw_region_df.columns],
            [s for s in non_carrier_cols if s in raw_region_df.columns],
            all_ploidies, ploidy_lookup,
            [s for s in sample_cols if s in raw_region_df.columns],
            carriers, gtf, segdup,
            col_title=f"Raw normalised — {locus_label}",
            min_gene_label_spacing=min_gene_label_spacing,
            show_colorbar=False,
            gaps=gaps,
            xform=xform,
        )
        _draw_overview_column(
            [axes[r, 1] for r in range(n_rows)],
            region_df, locus, calls_df,
            region_start, region_end,
            carrier_cols, non_carrier_cols,
            all_ploidies, ploidy_lookup,
            sample_cols, carriers, gtf, segdup,
            col_title=f"Processed — {locus_label}",
            min_gene_label_spacing=min_gene_label_spacing,
            show_colorbar=False,
            gaps=gaps,
            xform=xform,
        )
    else:
        _draw_overview_column(
            [axes[r, 0] for r in range(n_rows)],
            region_df, locus, calls_df,
            region_start, region_end,
            carrier_cols, non_carrier_cols,
            all_ploidies, ploidy_lookup,
            sample_cols, carriers, gtf, segdup,
            col_title=locus_label,
            min_gene_label_spacing=min_gene_label_spacing,
            show_colorbar=False,
            gaps=gaps,
            xform=xform,
        )

    plt.tight_layout(h_pad=0.3, w_pad=0.5)
    fig.subplots_adjust(hspace=0.08)

    locus_dir = os.path.join(output_dir, "locus_plots")
    os.makedirs(locus_dir, exist_ok=True)
    filename = f"{locus.cluster.replace('/', '_')}_overview.png"
    plt.savefig(os.path.join(locus_dir, filename), dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Created: locus_plots/{filename}")


def plot_sample_at_locus(
    sample_id: str,
    locus: GDLocus,
    calls_df: pd.DataFrame,
    depth_df: pd.DataFrame,
    gtf: Optional[GTFParser],
    segdup: Optional[SegDupAnnotation],
    output_dir: str,
    padding: int = 50000,
    min_gene_label_spacing: float = 0.05,
    gaps: Optional[GapsAnnotation] = None,
    flank_scale: float = 0.20,
):
    """Create detailed plot for a specific sample at a GD locus."""
    if not locus.breakpoints:
        return

    chrom = locus.chrom
    mask = (
        (depth_df["Cluster"] == locus.cluster) &
        (depth_df["Chr"] == chrom)
    )
    region_df = depth_df[mask].copy().sort_values("Start")

    if len(region_df) == 0 or sample_id not in region_df.columns:
        return

    region_start = int(region_df["Start"].min())
    region_end = int(region_df["End"].max())

    xform = FlankCompressor(region_start, region_end,
                            locus.start, locus.end, flank_scale=flank_scale)

    sample_calls = calls_df[
        (calls_df["cluster"] == locus.cluster) &
        (calls_df["sample"] == sample_id)
    ]

    figure_scale = 1.0
    fig, axes = plt.subplots(3, 1, figsize=(14 * figure_scale, 8 * figure_scale),
                             gridspec_kw={"height_ratios": [1, 2, 1]})

    bin_starts = region_df["Start"].values
    bin_ends = region_df["End"].values
    bin_mids = (bin_starts + bin_ends) / 2
    sample_depth = region_df[sample_id].values

    # Transform bar geometry
    d_bin_mids = xform(bin_mids)
    d_bar_widths = xform(bin_ends) - xform(bin_starts)

    # Panel 1: Annotations
    is_carrier = sample_calls["is_carrier"].any() if len(sample_calls) > 0 else False
    carrier_str = " [CARRIER]" if is_carrier else ""
    title = f"{sample_id} at {locus.cluster}{carrier_str}"
    draw_annotations_panel(axes[0], locus, region_start, region_end, chrom, title, gtf, segdup,
                           gaps=gaps, show_gd_entries=True,
                           min_gene_label_spacing=min_gene_label_spacing, xform=xform)
    axes[0].set_ylabel("Annotations")

    # Panel 2: Depth profile
    ax = axes[1]
    ax.bar(d_bin_mids, sample_depth, width=d_bar_widths * 0.9, alpha=0.7,
           color="steelblue", edgecolor="none")

    intervals = locus.get_intervals()
    for i, (start, end, name) in enumerate(intervals):
        ax.axvspan(xform(start), xform(end), alpha=0.1,
                   color=plt.cm.Set3(i / max(len(intervals), 1)))

    ax.axhline(2.0, color="green", linestyle="-", alpha=0.5, linewidth=1, label="CN=2")
    ax.axhline(1.0, color="orange", linestyle="--", alpha=0.5, linewidth=1, label="CN=1")
    ax.axhline(3.0, color="purple", linestyle="--", alpha=0.5, linewidth=1, label="CN=3")

    for i, (bp_start, bp_end) in enumerate(locus.breakpoints):
        ax.axvspan(xform(bp_start), xform(bp_end), alpha=0.2, color="red", zorder=0)

    ax.set_xlim(0.0, xform.d_end)
    ax.set_ylim(0, 5)
    ax.set_ylabel("Normalized Depth")
    ax.legend(loc="upper right")
    ax.grid(True, alpha=0.3, axis="y")
    xform.format_genomic_ticks(ax, locus.breakpoints)

    # Panel 3: Call summary
    ax = axes[2]
    ax.set_xlim(0.0, xform.d_end)
    ax.set_ylim(0, 1)

    y_pos = 0.5
    for _, call in sample_calls.iterrows():
        if call["is_carrier"]:
            color = "red" if call["svtype"] == "DEL" else "blue"
            alpha = 0.8 if call["is_best_match"] else 0.4
            ds, de = xform(call["start"]), xform(call["end"])
            ax.axvspan(ds, de, alpha=alpha, color=color)
            label = f"{call['svtype']} CN={call['mean_cn']:.2f}"
            if call["is_best_match"]:
                label += " (best)"
            ax.text((ds + de) / 2, y_pos, label,
                    ha="center", va="center", fontsize=9,
                    fontweight="bold" if call["is_best_match"] else "normal")

    ax.set_yticks([])
    ax.set_xlabel(f"Position on {chrom}")
    ax.set_ylabel("Calls")
    xform.format_genomic_ticks(ax, locus.breakpoints)

    plt.tight_layout()

    sample_dir = os.path.join(output_dir, "sample_plots", locus.cluster.replace("/", "_"))
    os.makedirs(sample_dir, exist_ok=True)
    filename = f"{sample_id.replace('/', '_')}.png"
    plt.savefig(os.path.join(sample_dir, filename), dpi=150, bbox_inches="tight")
    plt.close()


def plot_carrier_summary(calls_df: pd.DataFrame, output_dir: str):
    """Create summary plots of carrier calls across all loci."""
    carriers = calls_df[calls_df["is_carrier"]].copy()

    if len(carriers) == 0:
        print("No carriers to plot.")
        return

    figure_scale = 1.0
    fig, axes = plt.subplots(1, 2, figsize=(8 * figure_scale, 4 * figure_scale))

    ax = axes[0]
    locus_counts = carriers.groupby(["cluster", "svtype"])["sample"].nunique().unstack(fill_value=0)
    locus_counts.plot(kind="barh", ax=ax, color={"DEL": "red", "DUP": "blue"}, alpha=0.7)
    ax.set_xlabel("Number of Carriers")
    ax.set_ylabel("Locus")
    ax.set_title("Carriers per GD Locus")
    ax.legend(title="SV Type")

    ax = axes[1]
    all_scores = carriers["log_prob_score"].dropna()
    bins = np.linspace(all_scores.min(), all_scores.max(), 21)
    for svtype, color in [("DEL", "red"), ("DUP", "blue")]:
        sv_data = carriers[carriers["svtype"] == svtype]
        if len(sv_data) > 0:
            ax.hist(sv_data["log_prob_score"], bins=bins, alpha=0.6, color=color,
                    label=f"{svtype} (n={len(sv_data)})", edgecolor="black")
    ax.set_xlabel("Log Probability Score")
    ax.set_ylabel("Count")
    ax.set_title("Log Probability Score Distribution in Carriers")
    ax.legend()

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "carrier_summary.png"), dpi=150, bbox_inches="tight")
    plt.close()
    print("Created: carrier_summary.png")


def plot_confidence_distribution(calls_df: pd.DataFrame, output_dir: str):
    """Plot distribution of confidence scores."""
    figure_scale = 1.0
    fig, axes = plt.subplots(1, 2, figsize=(8 * figure_scale, 4 * figure_scale))

    ax = axes[0]
    carriers = calls_df[calls_df["is_carrier"]]
    non_carriers = calls_df[~calls_df["is_carrier"]]
    all_scores = calls_df["log_prob_score"].dropna()
    bins = np.linspace(all_scores.min(), all_scores.max(), 31)
    ax.hist(non_carriers["log_prob_score"], bins=bins, alpha=0.6, label="Non-carriers",
            color="gray", edgecolor="black")
    ax.hist(carriers["log_prob_score"], bins=bins, alpha=0.6, label="Carriers",
            color="green", edgecolor="black")
    ax.set_xlabel("Log Probability Score")
    ax.set_ylabel("Count")
    ax.set_yscale("log")
    ax.legend()

    ax = axes[1]
    non_carriers_mask = ~calls_df["is_carrier"]
    ax.scatter(calls_df.loc[non_carriers_mask, "mean_depth"],
               calls_df.loc[non_carriers_mask, "log_prob_score"],
               c="gray", alpha=0.1, s=10, label="Non-carriers")
    carriers_mask = calls_df["is_carrier"]
    ax.scatter(calls_df.loc[carriers_mask, "mean_depth"],
               calls_df.loc[carriers_mask, "log_prob_score"],
               c="green", alpha=0.5, s=10, label="Carriers")
    ax.set_xlabel("Mean Depth (normalized)")
    ax.set_ylabel("Confidence")

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "confidence_distribution.png"), dpi=150, bbox_inches="tight")
    plt.close()
    print("Created: confidence_distribution.png")


def create_carrier_pdf(
    calls_df: pd.DataFrame,
    depth_df: pd.DataFrame,
    gd_table: GDTable,
    gtf: Optional[GTFParser],
    segdup: Optional[SegDupAnnotation],
    output_dir: str,
    padding: int = 50000,
    min_gene_label_spacing: float = 0.05,
    raw_counts_df: Optional[pd.DataFrame] = None,
    gaps: Optional[GapsAnnotation] = None,
    flank_scale: float = 0.20,
    lowres_median_bin_size: Optional[float] = None,
    highres_path: Optional[str] = None,
):
    """Create a PDF with plots for all carrier samples."""
    carriers = calls_df[calls_df["is_carrier"]]

    if len(carriers) == 0:
        print("No carriers to include in PDF.")
        return

    # Precompute per-sample genome-wide autosomal median raw counts
    raw_sample_medians: Dict[str, float] = {}
    if raw_counts_df is not None:
        autosomal = [str(c) for c in range(1, 23)] + [f"chr{c}" for c in range(1, 23)]
        raw_auto = raw_counts_df[raw_counts_df["Chr"].isin(autosomal)]
        raw_sample_cols = [c for c in raw_counts_df.columns
                          if c not in ("Chr", "Start", "End", "source_file", "Bin")]
        for s in raw_sample_cols:
            med = raw_auto[s].median()
            if med > 0:
                raw_sample_medians[s] = med
        print(f"  Computed raw autosomal medians for {len(raw_sample_medians)} samples")

    if lowres_median_bin_size is None and raw_counts_df is not None:
        lowres_median_bin_size = estimate_lowres_bin_size(raw_counts_df)

    pdf_path = os.path.join(output_dir, "carrier_plots.pdf")
    print(f"Creating carrier PDF: {pdf_path}")

    with PdfPages(pdf_path) as pdf:
        for cluster in sorted(carriers["cluster"].unique()):
            locus = gd_table.loci.get(cluster)
            if locus is None:
                continue

            if not locus.breakpoints:
                print(f"  Warning: No breakpoints for {cluster}, skipping")
                continue

            cluster_carriers = carriers[carriers["cluster"] == cluster]["sample"].unique()
            print(f"  Adding {len(cluster_carriers)} carriers for {cluster}")

            chrom = locus.chrom
            locus_mask = (
                (depth_df["Cluster"] == cluster) &
                (depth_df["Chr"] == chrom)
            )
            region_df = depth_df[locus_mask].sort_values("Start")
            if len(region_df) == 0:
                continue

            region_start = int(region_df["Start"].min())
            region_end = int(region_df["End"].max())

            # Build coordinate transform for this locus
            xform = FlankCompressor(region_start, region_end,
                                    locus.start, locus.end,
                                    flank_scale=flank_scale)

            cluster_calls_df = calls_df[calls_df["cluster"] == cluster]
            _raw_region = None
            if (
                (raw_counts_df is not None or highres_path is not None)
                and raw_sample_medians
                and lowres_median_bin_size is not None
            ):
                _raw_region = _build_raw_region_df(
                    locus, region_start, region_end,
                    raw_counts_df, list(raw_sample_medians.keys()),
                    raw_sample_medians, lowres_median_bin_size,
                    highres_path=highres_path,
                )

            for sample_id in sorted(cluster_carriers):
                if sample_id not in region_df.columns:
                    continue

                sample_calls = cluster_calls_df[cluster_calls_df["sample"] == sample_id]

                figure_scale = 1.0
                fig, axes = plt.subplots(2, 1, figsize=(12 * figure_scale, 4 * figure_scale),
                                         gridspec_kw={"height_ratios": [1, 2]})

                bin_mids = (region_df["Start"].values + region_df["End"].values) / 2
                sample_depth = region_df[sample_id].values

                # Transform bar geometry
                d_bin_mids = xform(bin_mids)
                d_bar_widths = xform(region_df["End"].values) - xform(region_df["Start"].values)

                # Panel 1: Annotations
                carrier_call = (
                    sample_calls[sample_calls["is_carrier"]].iloc[0]
                    if len(sample_calls[sample_calls["is_carrier"]]) > 0
                    else None
                )
                call_info = (
                    f" - {carrier_call['svtype']} confidence={carrier_call['log_prob_score']:.2f}"
                    if carrier_call is not None else ""
                )
                title = f"{sample_id} at {cluster}{call_info}"
                draw_annotations_panel(
                    axes[0], locus, region_start, region_end, chrom, title,
                    gtf, segdup, gaps=gaps, show_gd_entries=False,
                    min_gene_label_spacing=min_gene_label_spacing, xform=xform,
                )

                # Panel 2: Depth
                ax = axes[1]

                best_call = None
                if len(sample_calls) > 0:
                    carrier_calls_local = sample_calls[sample_calls["is_carrier"]]
                    if len(carrier_calls_local) > 0:
                        best_matches = carrier_calls_local[carrier_calls_local["is_best_match"]]
                        if len(best_matches) > 0:
                            best_call = best_matches.iloc[0]
                        else:
                            best_call = carrier_calls_local.loc[
                                carrier_calls_local["log_prob_score"].idxmax()]

                if best_call is not None:
                    interval_start = best_call["start"]
                    interval_end = best_call["end"]
                    mean_depth = best_call["mean_depth"]
                    svtype = best_call["svtype"]
                    color = "#FF6B6B" if svtype == "DEL" else "#6B9BD1"
                    d_is, d_ie = xform(interval_start), xform(interval_end)
                    ax.axvspan(d_is, d_ie, alpha=0.2, color=color, zorder=1,
                               label=f"{svtype} region")
                    ax.hlines(mean_depth, d_is, d_ie,
                              colors="black", linewidth=2.5, alpha=0.8, zorder=2,
                              label=f"Mean depth={mean_depth:.2f}")

                ax.bar(d_bin_mids, sample_depth, width=d_bar_widths * 0.9, alpha=0.6,
                       color="steelblue", edgecolor="none", zorder=3)

                # Overlay raw depth trace if available
                if _raw_region is not None and sample_id in _raw_region.columns:
                    raw_mids = (_raw_region["Start"].values + _raw_region["End"].values) / 2
                    norm_raw = _raw_region[sample_id].values.astype(float)
                    ax.plot(xform(raw_mids), norm_raw, color="darkorange", linewidth=0.6,
                            alpha=0.7, zorder=4, label="Raw depth")

                ax.axhline(2.0, color="green", linestyle="-", alpha=0.4, linewidth=1.5,
                           label="Reference CN=2", zorder=0)
                ax.axhline(1.0, color="orange", linestyle=":", alpha=0.4, linewidth=1, zorder=0)
                ax.axhline(3.0, color="purple", linestyle=":", alpha=0.4, linewidth=1, zorder=0)

                ax.set_xlim(0.0, xform.d_end)
                ax.set_ylim(0, 5)
                ax.set_xlabel(f"Position on {chrom}")
                ax.set_ylabel("Normalized Depth")
                ax.grid(True, alpha=0.3, axis="y", zorder=0)
                ax.legend(loc="upper right", fontsize=8)
                xform.format_genomic_ticks(ax, locus.breakpoints)

                plt.tight_layout()
                pdf.savefig(fig)
                plt.close()

    print(f"  Saved carrier PDF: {pdf_path}")


# =============================================================================
# Main
# =============================================================================


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Plot GD CNV detection results",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--calls", "-c",
        required=True,
        help="GD CNV calls file (gd_cnv_calls.tsv.gz) produced by gd_cnv_call.py",
    )
    parser.add_argument(
        "--cn-posteriors",
        required=True,
        help="CN posteriors file (cn_posteriors.tsv.gz) with depth values",
    )
    parser.add_argument(
        "--gd-table", "-g",
        required=True,
        help="GD locus definition table (TSV)",
    )
    parser.add_argument(
        "--ploidy-table",
        required=False,
        help="Ploidy estimates TSV (ploidy_estimates.tsv) from gd_cnv_pyro.py. "
             "Columns: sample, contig, median_depth, ploidy. "
             "If not provided, ploidy=2 is assumed for all sample/contig pairs.",
    )
    parser.add_argument(
        "--gtf",
        required=False,
        help="Gene annotation GTF file (can be gzipped)",
    )
    parser.add_argument(
        "--segdup-bed",
        required=False,
        help="Segmental duplication BED file (can be gzipped)",
    )
    parser.add_argument(
        "--gaps-bed",
        required=False,
        help="Reference assembly gaps BED file (can be gzipped).  Regions "
             "are drawn as hatched gray spans on annotation panels.",
    )
    parser.add_argument(
        "--raw-counts",
        required=False,
        help="Raw read-count matrix TSV (can be gzipped).  When provided, "
             "the raw depth profile is superimposed on carrier plots.",
    )
    parser.add_argument(
        "--high-res-counts",
        required=False,
        help="Bgzipped, tabix-indexed high-resolution read-count file "
             "(.tsv.gz + .tbi).  Requires --raw-counts.  For each GD body "
             "interval and flank whose genomic span is less than 10\u00d7 the "
             "estimated low-resolution bin size, bins are queried on demand "
             "via pysam rather than slicing the low-resolution matrix.",
    )
    parser.add_argument(
        "--output-dir", "-o",
        required=True,
        help="Output directory for plots",
    )
    parser.add_argument(
        "--padding",
        type=int,
        default=50000,
        help="Padding around locus boundaries (bp)",
    )
    parser.add_argument(
        "--plot-all-samples",
        action="store_true",
        help="Plot individual sample plots for all samples (not just carriers)",
    )
    parser.add_argument(
        "--sample",
        type=str,
        help="Plot only this specific sample",
    )
    parser.add_argument(
        "--min-gene-label-spacing",
        type=float,
        default=0.05,
        help="Minimum distance between adjacent gene label centres as a fraction "
             "of the plot width (default: 0.05 = 5%%).",
    )
    parser.add_argument(
        "--locus",
        action="append",
        dest="loci",
        metavar="LOCUS",
        help="Restrict plotting to this locus (cluster name).  May be given "
             "multiple times to select several loci.  When omitted all loci "
             "are plotted.",
    )
    parser.add_argument(
        "--flank-scale",
        type=float,
        default=0.20,
        help="Fraction of the total plot width that each flanking region "
             "should occupy.  0.20 means each flank gets 20%% of the "
             "display width, leaving 60%% for the locus body (when both "
             "flanks are present).  Smaller values compress flanks more.  "
             "Use 0 to hide flanks entirely.",
    )
    return parser.parse_args()


def main():
    """Main function."""
    args = parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    log_fh = setup_logging(args.output_dir)

    print(f"Output directory: {args.output_dir}")
    print(f"Flank scale: {args.flank_scale}")

    # Load data
    print("\nLoading data...")

    print(f"  Loading calls: {args.calls}")
    calls_df = pd.read_csv(args.calls, sep="\t", compression="infer")
    print(f"    {len(calls_df)} call records")

    print(f"  Loading CN posteriors: {args.cn_posteriors}")
    cn_posteriors_df = pd.read_csv(args.cn_posteriors, sep="\t", compression="infer")
    print(f"    {len(cn_posteriors_df)} bin-sample records")

    print(f"  Loading GD table: {args.gd_table}")
    gd_table = GDTable(args.gd_table)
    print(f"    {len(gd_table.loci)} loci")

    # Load ploidy table if provided
    ploidy_df = None
    if args.ploidy_table:
        print(f"  Loading ploidy table: {args.ploidy_table}")
        ploidy_df = pd.read_csv(args.ploidy_table, sep="\t")
        print(f"    {len(ploidy_df)} sample/contig ploidy records")

    # Convert cn_posteriors to depth_df format (bins x samples matrix)
    print("\n  Converting posteriors to depth matrix format...")
    depth_df = cn_posteriors_df.pivot(
        index=["cluster", "chr", "start", "end"],
        columns="sample",
        values="depth",
    ).reset_index()
    depth_df = depth_df.rename(columns={
        "cluster": "Cluster", "chr": "Chr", "start": "Start", "end": "End",
    })
    sample_cols_all = [c for c in depth_df.columns
                       if c not in ["Cluster", "Chr", "Start", "End"]]
    print(f"    {len(depth_df)} bin-rows x {len(sample_cols_all)} samples (across all loci)")

    # Optional annotations
    gtf = None
    if args.gtf:
        print(f"\n  Loading GTF: {args.gtf}")
        gtf = GTFParser(args.gtf)

    segdup = None
    if args.segdup_bed:
        print(f"  Loading segdup BED: {args.segdup_bed}")
        segdup = SegDupAnnotation(args.segdup_bed)

    gaps = None
    if args.gaps_bed:
        print(f"  Loading gaps BED: {args.gaps_bed}")
        gaps = GapsAnnotation(args.gaps_bed)

    # Load raw counts matrix if provided
    raw_counts_df = None
    raw_sample_medians: Dict[str, float] = {}
    if args.raw_counts:
        print(f"\n  Loading raw counts: {args.raw_counts}")
        raw_counts_df = pd.read_csv(args.raw_counts, sep="\t", compression="infer")
        if "#Chr" in raw_counts_df.columns:
            raw_counts_df.rename(columns={"#Chr": "Chr"}, inplace=True)
        raw_counts_df["Start"] = raw_counts_df["Start"].astype(int)
        raw_counts_df["End"] = raw_counts_df["End"].astype(int)
        print(f"    {len(raw_counts_df)} bins")
        autosomal = [str(c) for c in range(1, 23)] + [f"chr{c}" for c in range(1, 23)]
        raw_auto = raw_counts_df[raw_counts_df["Chr"].isin(autosomal)]
        raw_sample_cols = [c for c in raw_counts_df.columns
                          if c not in ("Chr", "Start", "End", "source_file", "Bin")]
        for s in raw_sample_cols:
            med = raw_auto[s].median()
            if med > 0:
                raw_sample_medians[s] = med
        print(f"    Computed raw autosomal medians for {len(raw_sample_medians)} samples")

    lowres_median_bin_size: Optional[float] = None
    if raw_counts_df is not None:
        lowres_median_bin_size = estimate_lowres_bin_size(raw_counts_df)
        print(f"    Estimated low-res median bin size: {lowres_median_bin_size:,.0f} bp")

    highres_path: Optional[str] = getattr(args, "high_res_counts", None)
    if highres_path:
        if not args.raw_counts:
            print("ERROR: --high-res-counts requires --raw-counts for per-sample median estimation")
            sys.exit(1)
        print(f"  High-res counts: {highres_path}")

    # ---- Determine which loci to plot ----
    if args.loci:
        loci_to_plot = {}
        for name in args.loci:
            if name in gd_table.loci:
                loci_to_plot[name] = gd_table.loci[name]
            else:
                print(f"  WARNING: --locus '{name}' not found in GD table, skipping")
        if not loci_to_plot:
            print("ERROR: none of the requested loci were found in the GD table")
            sys.exit(1)
        print(f"  Restricting to {len(loci_to_plot)} of {len(gd_table.loci)} loci")
    else:
        loci_to_plot = gd_table.loci

    # Optionally subset calls to selected loci so summary plots reflect the
    # same scope as the per-locus plots.
    if args.loci:
        plot_calls_df = calls_df[calls_df["cluster"].isin(loci_to_plot)].copy()
    else:
        plot_calls_df = calls_df

    # Build a filtered GDTable-like object for create_carrier_pdf
    class _FilteredGDTable:
        def __init__(self, loci_dict):
            self.loci = loci_dict

    filtered_gd_table = _FilteredGDTable(loci_to_plot)

    # Create summary plots
    print("\nCreating summary plots...")
    plot_carrier_summary(plot_calls_df, args.output_dir)
    plot_confidence_distribution(plot_calls_df, args.output_dir)

    # Create locus overview plots
    print("\nCreating locus overview plots...")
    ploidy_lookup = _build_ploidy_lookup(ploidy_df)
    depth_by_cluster = {k: v for k, v in depth_df.groupby("Cluster")}
    for cluster, locus in loci_to_plot.items():
        print(f"  Processing {cluster}...")
        cluster_depth = depth_by_cluster.get(cluster)
        if cluster_depth is None:
            print(f"  Warning: No depth data for {cluster}, skipping")
            continue
        plot_locus_overview(
            locus, calls_df, cluster_depth, gtf, segdup,
            args.output_dir, padding=args.padding,
            ploidy_lookup=ploidy_lookup,
            min_gene_label_spacing=args.min_gene_label_spacing,
            raw_counts_df=raw_counts_df,
            raw_sample_medians=raw_sample_medians if raw_sample_medians else None,
            gaps=gaps,
            flank_scale=args.flank_scale,
            lowres_median_bin_size=lowres_median_bin_size,
            highres_path=highres_path,
        )

    # Create individual sample plots
    if args.plot_all_samples or args.sample:
        print("\nCreating individual sample plots...")
        sample_cols = [c for c in depth_df.columns if c not in ['Chr', 'Start', 'End']]
        carriers = calls_df[calls_df["is_carrier"]]

        for cluster, locus in loci_to_plot.items():
            if args.sample:
                samples_to_plot = [args.sample]
            elif args.plot_all_samples:
                samples_to_plot = sample_cols
            else:
                samples_to_plot = carriers[carriers["cluster"] == cluster]["sample"].unique()

            for sample_id in samples_to_plot:
                plot_sample_at_locus(
                    sample_id, locus, calls_df, depth_df, gtf, segdup,
                    args.output_dir, padding=args.padding,
                    min_gene_label_spacing=args.min_gene_label_spacing,
                    gaps=gaps,
                    flank_scale=args.flank_scale,
                )

    # Create carrier PDF
    print("\nCreating carrier PDF...")
    create_carrier_pdf(
        plot_calls_df, depth_df, filtered_gd_table, gtf, segdup,
        args.output_dir, padding=args.padding,
        min_gene_label_spacing=args.min_gene_label_spacing,
        raw_counts_df=raw_counts_df,
        gaps=gaps,
        flank_scale=args.flank_scale,
        lowres_median_bin_size=lowres_median_bin_size,
        highres_path=highres_path,
    )

    print("\n" + "=" * 80)
    print("Plotting complete!")
    print("=" * 80)


if __name__ == "__main__":
    main()
