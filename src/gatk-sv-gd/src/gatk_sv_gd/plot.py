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
from time import perf_counter
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


EMPTY_CALLS_COLUMNS = [
    "sample",
    "cluster",
    "GD_ID",
    "chrom",
    "start",
    "end",
    "svtype",
    "BP1",
    "BP2",
    "is_terminal",
    "n_bins",
    "mean_depth",
    "sample_ploidy",
    "matched_haplotype",
    "hap_cn_state",
    "matched_seg_start",
    "matched_seg_end",
    "matched_seg_n_bins",
    "matched_interval_bp",
    "interval_coverage",
    "reciprocal_overlap",
    "min_interval_confidence",
    "left_flank_non_event_median",
    "right_flank_non_event_median",
    "min_flank_non_event_confidence",
    "is_carrier",
    "is_best_match",
    "log_prob_score",
    "confidence_score",
    "calling_method",
]


PLOT_TIMING_ENABLED = os.getenv("GATK_SV_GD_PLOT_TIMING", "").strip().lower() in {
    "1", "true", "yes", "on",
}

PDF_MAX_SIGNAL_BINS = 300


def _print_timing(label: str, start_time: float) -> None:
    """Print elapsed wall-clock time for a plotting stage."""
    if PLOT_TIMING_ENABLED:
        print(f"    [timing] {label}: {perf_counter() - start_time:.3f}s")

# =============================================================================
# Viterbi overlay data container
# =============================================================================


def _merge_segments(segments: List[dict]) -> List[dict]:
    """Merge consecutive segments with the same CN state into single spans.

    Input segments must be sorted by ``start``.  Adjacent segments whose
    ``cn_state`` values match are fused: the merged segment inherits the
    earliest ``start`` and latest ``end``.  This guarantees every bin is
    covered by exactly one output segment with no gaps or overlaps between
    same-CN-state runs.
    """
    if not segments:
        return []
    merged: List[dict] = [dict(segments[0])]
    for seg in segments[1:]:
        prev = merged[-1]
        if seg["cn_state"] == prev["cn_state"]:
            prev["end"] = seg["end"]
        else:
            merged.append(dict(seg))
    return merged


class ViterbiOverlayData:
    """Pre-loaded Viterbi paths produced by the calling step.

    Loads the ``viterbi_paths.tsv.gz`` file written by ``call.py`` so that
    the exact same path computed during calling is drawn on carrier plots.

    Supports the new per-haplotype format (``haplotype`` column):
      - haplotype=0: total CN segments
      - haplotype=1: haplotype 1 segments (diploid only)
      - haplotype=2: haplotype 2 segments (diploid only)
    """

    def __init__(self, paths_df: pd.DataFrame):
        # Ensure the cn_state column exists.  Support old files that used
        # a ``mean_cn`` column (float) by rounding to int.
        if "cn_state" not in paths_df.columns:
            if "mean_cn" in paths_df.columns:
                print(
                    "  WARNING: viterbi_paths.tsv.gz uses the old mean_cn "
                    "column.  Re-run gd_cnv_call.py to get the new cn_state "
                    "format.  Falling back to rounding mean_cn to int."
                )
                paths_df = paths_df.copy()
                paths_df["cn_state"] = paths_df["mean_cn"].round().astype(int)
            else:
                raise ValueError(
                    "viterbi_paths.tsv.gz has neither 'cn_state' nor 'mean_cn' "
                    "columns.  Please re-run gd_cnv_call.py to regenerate it."
                )
        if "category" not in paths_df.columns:
            paths_df = paths_df.copy()
            paths_df["category"] = paths_df["cn_state"].apply(
                lambda cn: "DEL" if cn < 2 else ("DUP" if cn > 2 else "REF")
            )

        # Add haplotype column if absent (backward compat with pre-diploid files)
        if "haplotype" not in paths_df.columns:
            paths_df = paths_df.copy()
            paths_df["haplotype"] = 0

        # Index by (sample, cluster, haplotype) for O(1) lookup.
        # Merge consecutive records with the same category into single
        # segments so that each bin is covered by exactly one segment with
        # no gaps or overlaps.
        self._paths: Dict[Tuple[str, str, int], List[dict]] = {}
        for (sample, cluster, hap), grp in paths_df.groupby(
            ["sample", "cluster", "haplotype"]
        ):
            rows = grp.sort_values("start")
            merged = _merge_segments([
                {
                    "start": int(row["start"]),
                    "end": int(row["end"]),
                    "cn_state": int(row["cn_state"]),
                    "category": str(row["category"]),
                }
                for _, row in rows.iterrows()
            ])
            self._paths[(str(sample), str(cluster), int(hap))] = merged

    def get_path(
        self,
        sample_id: str,
        cluster: str,
        haplotype: int = 0,
    ) -> Optional[List[dict]]:
        """Look up the pre-computed Viterbi category segments for one
        sample/cluster/haplotype.

        Args:
            sample_id: Sample identifier.
            cluster: Cluster name.
            haplotype: 0 for total CN (default), 1 or 2 for per-haplotype.

        Returns a list of segment dicts with keys ``start``, ``end``,
        ``cn_state`` (int), and ``category`` (DEL/REF/DUP), or *None*
        if no path was saved for this combination.
        """
        return self._paths.get((str(sample_id), str(cluster), haplotype))

    def has_haplotype_paths(
        self,
        sample_id: str,
        cluster: str,
    ) -> bool:
        """Return True if per-haplotype paths are available for this sample."""
        return any(
            self.get_path(sample_id, cluster, haplotype) is not None
            for haplotype in (1, 2)
        )


def _get_confidence_column(calls_df: pd.DataFrame) -> str:
    """Return the preferred confidence column present in the calls table."""
    if "confidence_score" in calls_df.columns:
        return "confidence_score"
    if "log_prob_score" in calls_df.columns:
        return "log_prob_score"
    raise ValueError(
        "Calls table is missing both 'confidence_score' and 'log_prob_score'."
    )


def _get_confidence_label(confidence_column: str) -> str:
    """Return a human-readable axis label for the confidence column."""
    return "Confidence Score" if confidence_column == "confidence_score" else "Log Probability Score"


def _sanitize_plot_label(label: str) -> str:
    """Return a filesystem-safe plot label for locus and sample outputs."""
    sanitized = str(label)
    for char in ("/", "-", ":"):
        sanitized = sanitized.replace(char, "_")
    return sanitized


def _build_carrier_best_match_mask(calls_df: pd.DataFrame) -> pd.Series:
    """Return the mask used to identify selected carrier calls."""
    mask = calls_df["is_carrier"] == True  # noqa: E712
    if "is_best_match" in calls_df.columns:
        mask = mask & (calls_df["is_best_match"] == True)  # noqa: E712
    return mask


def _compute_raw_sample_medians(
    raw_counts_df: Optional[pd.DataFrame],
) -> Dict[str, float]:
    """Compute per-sample autosomal raw-depth medians for plot overlays."""
    raw_sample_medians: Dict[str, float] = {}
    if raw_counts_df is None:
        return raw_sample_medians

    autosomal = [str(c) for c in range(1, 23)] + [f"chr{c}" for c in range(1, 23)]
    raw_auto = raw_counts_df[raw_counts_df["Chr"].isin(autosomal)]
    raw_sample_cols = [
        c for c in raw_counts_df.columns
        if c not in ("Chr", "Start", "End", "source_file", "Bin")
    ]
    for sample_id in raw_sample_cols:
        med = raw_auto[sample_id].median()
        if med > 0:
            raw_sample_medians[sample_id] = med
    return raw_sample_medians


def _parse_eval_sample_list(value: object) -> List[str]:
    """Parse a comma-delimited eval report sample field into sample IDs."""
    if pd.isna(value):
        return []
    text = str(value).strip()
    if not text:
        return []
    return [sample.strip() for sample in text.split(",") if sample.strip()]


def _build_gd_to_cluster_map(
    loci_by_cluster: Dict[str, GDLocus],
) -> Dict[str, str]:
    """Map each GD_ID to its containing cluster."""
    gd_to_cluster: Dict[str, str] = {}
    for cluster, locus in loci_by_cluster.items():
        for entry in locus.gd_entries:
            gd_to_cluster[str(entry["GD_ID"])] = str(cluster)
    return gd_to_cluster


def _build_eval_pdf_specs(
    eval_report_df: pd.DataFrame,
    calls_df: pd.DataFrame,
    loci_by_cluster: Dict[str, GDLocus],
    confidence_threshold: float,
) -> Dict[str, List[dict]]:
    """Build TP/FP/FN page specs from an eval report and calls table."""
    confidence_column = _get_confidence_column(calls_df)
    carrier_mask = _build_carrier_best_match_mask(calls_df)
    selected_calls = calls_df[
        carrier_mask & (calls_df[confidence_column] >= confidence_threshold)
    ].copy()
    predicted_by_gd = {
        str(gd_id): set(group["sample"].astype(str).unique())
        for gd_id, group in selected_calls.groupby("GD_ID")
    }
    gd_to_cluster = _build_gd_to_cluster_map(loci_by_cluster)

    specs: Dict[str, List[dict]] = {
        "true_positives": [],
        "false_positives": [],
        "false_negatives": [],
    }
    category_labels = {
        "true_positives": "True positive",
        "false_positives": "False positive",
        "false_negatives": "False negative",
    }

    for _, row in eval_report_df.iterrows():
        gd_id = str(row["GD_ID"])
        cluster = gd_to_cluster.get(gd_id)
        if cluster is None:
            continue

        fp_samples = set(_parse_eval_sample_list(row.get("FP_samples", "")))
        fn_samples = set(_parse_eval_sample_list(row.get("FN_samples", "")))
        tp_samples = sorted(predicted_by_gd.get(gd_id, set()) - fp_samples)

        expected_tp = int(row.get("TP", len(tp_samples)))
        if expected_tp != len(tp_samples):
            print(
                f"  WARNING: derived TP count for {gd_id} is {len(tp_samples)} "
                f"but eval report says {expected_tp}. "
                "Make sure the carrier confidence threshold matches the eval run."
            )

        for category_name, sample_ids in (
            ("true_positives", tp_samples),
            ("false_positives", sorted(fp_samples)),
            ("false_negatives", sorted(fn_samples)),
        ):
            for sample_id in sample_ids:
                specs[category_name].append(
                    {
                        "cluster": cluster,
                        "sample": str(sample_id),
                        "gd_id": gd_id,
                        "title_suffix": f"{category_labels[category_name]}: {gd_id}",
                    }
                )

    return specs


def _render_pdf_sample_page(
    pdf: PdfPages,
    sample_id: str,
    cluster: str,
    locus: GDLocus,
    region_df: pd.DataFrame,
    cluster_calls_df: pd.DataFrame,
    confidence_column: str,
    confidence_threshold: float,
    gtf: Optional[GTFParser],
    segdup: Optional[SegDupAnnotation],
    min_gene_label_spacing: float,
    xform: FlankCompressor,
    raw_region_df: Optional[pd.DataFrame],
    minor_baf_region_df: Optional[pd.DataFrame],
    baf_sites_region_df: Optional[pd.DataFrame],
    event_del_region_df: Optional[pd.DataFrame],
    event_dup_region_df: Optional[pd.DataFrame],
    gaps: Optional[GapsAnnotation],
    viterbi_data: Optional[ViterbiOverlayData],
    title_suffix: Optional[str] = None,
) -> bool:
    """Render one carrier-style review page into an output PDF."""
    if sample_id not in region_df.columns:
        return False

    page_start = perf_counter()

    sample_calls = cluster_calls_df[cluster_calls_df["sample"] == sample_id]
    pdf_calls = sample_calls[sample_calls[confidence_column] >= confidence_threshold]

    region_start = int(region_df["Start"].min())
    region_end = int(region_df["End"].max())
    chrom = locus.chrom

    figure_scale = 1.0
    fig, axes = plt.subplots(
        4,
        1,
        figsize=(12 * figure_scale, 7.5 * figure_scale),
        gridspec_kw={"height_ratios": [1, 2, 1, 1]},
    )

    sample_depth = region_df[sample_id].values
    sample_minor_baf, sample_baf_sites = _extract_sample_baf_vectors(
        region_df,
        minor_baf_region_df,
        baf_sites_region_df,
        sample_id,
    )

    carrier_call = None
    if len(pdf_calls) > 0:
        best_matches = pdf_calls[pdf_calls["is_best_match"]]
        if len(best_matches) > 0:
            carrier_call = best_matches.loc[best_matches[confidence_column].idxmax()]
        else:
            carrier_call = pdf_calls.loc[pdf_calls[confidence_column].idxmax()]

    call_info = (
        f" - {carrier_call['svtype']} confidence={carrier_call[confidence_column]:.2f}"
        if carrier_call is not None else ""
    )

    event_probs = None
    event_svtype = None
    if carrier_call is not None:
        event_svtype = str(carrier_call["svtype"])
        event_region_df = event_del_region_df if event_svtype == "DEL" else event_dup_region_df
        event_probs = _extract_sample_event_probabilities(
            region_df,
            event_region_df,
            sample_id,
        )

    plot_region_df = region_df
    plot_sample_depth = np.asarray(sample_depth, dtype=float)
    plot_sample_minor_baf = sample_minor_baf
    plot_sample_baf_sites = sample_baf_sites
    plot_event_probs = event_probs
    if len(region_df) > PDF_MAX_SIGNAL_BINS:
        (
            plot_region_df,
            plot_sample_depth,
            plot_sample_minor_baf,
            plot_sample_baf_sites,
            plot_event_probs,
        ) = _coarsen_pdf_page_signals(
            locus,
            region_df,
            sample_depth,
            sample_minor_baf,
            sample_baf_sites,
            event_probs,
            max_total_bins=PDF_MAX_SIGNAL_BINS,
        )

    plot_bin_mids = (plot_region_df["Start"].values + plot_region_df["End"].values) / 2
    d_bin_mids = xform(plot_bin_mids)
    d_bar_widths = xform(plot_region_df["End"].values) - xform(plot_region_df["Start"].values)

    title = f"{sample_id} at {cluster}{call_info}"
    if title_suffix:
        title = f"{title} [{title_suffix}]"
    stage_start = perf_counter()
    draw_annotations_panel(
        axes[0],
        locus,
        region_start,
        region_end,
        chrom,
        title,
        gtf,
        segdup,
        gaps=gaps,
        show_gd_entries=False,
        min_gene_label_spacing=min_gene_label_spacing,
        xform=xform,
    )
    _print_timing(f"{cluster} {sample_id} pdf annotations", stage_start)

    ax = axes[1]
    stage_start = perf_counter()
    best_call = carrier_call
    if best_call is not None:
        svtype = best_call["svtype"]
        mean_depth = best_call["mean_depth"]
        color = "#FF6B6B" if svtype == "DEL" else "#6B9BD1"

        interval_start = best_call["start"]
        interval_end = best_call["end"]
        d_is, d_ie = xform(interval_start), xform(interval_end)
        ax.axvspan(d_is, d_ie, alpha=0.2, color=color, zorder=1,
                   label=f"{svtype} region")
        ax.hlines(mean_depth, d_is, d_ie, colors="black", linewidth=2.5, alpha=0.8, zorder=2, label=f"Mean depth={mean_depth:.2f}")

    _plot_depth_bars_with_baf(
        ax,
        d_bin_mids,
        d_bar_widths,
        plot_sample_depth,
        minor_baf_values=plot_sample_minor_baf,
        baf_site_counts=plot_sample_baf_sites,
        zorder=3,
    )

    if raw_region_df is not None and sample_id in raw_region_df.columns:
        raw_mids = (raw_region_df["Start"].values + raw_region_df["End"].values) / 2
        norm_raw = raw_region_df[sample_id].values.astype(float)
        ax.plot(xform(raw_mids), norm_raw, color="darkorange", linewidth=0.6,
                alpha=0.7, zorder=4, label="Raw depth")

    if viterbi_data is not None:
        vit_segs = viterbi_data.get_path(sample_id, cluster)
        _plot_viterbi_trace(
            ax,
            vit_segs,
            xform,
            color="magenta",
            label="Viterbi CN",
            linewidth=2.0,
            alpha=0.85,
            zorder=5,
        )

        if viterbi_data.has_haplotype_paths(sample_id, cluster):
            hap_styles = [
                (1, "red", "Hap 1", -0.04),
                (2, "blue", "Hap 2", 0.04),
            ]
            for hap_id, hap_color, hap_label, hap_offset in hap_styles:
                hap_segs = viterbi_data.get_path(sample_id, cluster, hap_id)
                _plot_viterbi_trace(
                    ax,
                    hap_segs,
                    xform,
                    color=hap_color,
                    label=hap_label,
                    linewidth=1.6,
                    alpha=0.8,
                    zorder=6,
                    vertical_offset=hap_offset,
                )

    ax.axhline(2.0, color="green", linestyle="-", alpha=0.4, linewidth=1.5,
               label="Reference CN=2", zorder=0)
    ax.axhline(1.0, color="orange", linestyle=":", alpha=0.4, linewidth=1, zorder=0)
    ax.axhline(3.0, color="purple", linestyle=":", alpha=0.4, linewidth=1, zorder=0)
    ax.set_xlim(0.0, xform.d_end)
    ax.set_ylim(0, 5)
    ax.set_ylabel("Normalized Depth")
    ax.grid(True, alpha=0.3, axis="y", zorder=0)
    ax.legend(loc="upper right", fontsize=8)
    _print_timing(f"{cluster} {sample_id} pdf depth panel", stage_start)

    baf_ax = axes[2]
    stage_start = perf_counter()
    _plot_baf_signal_panel(
        baf_ax,
        d_bin_mids,
        d_bar_widths,
        plot_sample_minor_baf,
        plot_sample_baf_sites,
        xform,
        locus,
        chrom,
        show_xlabel=False,
    )
    _print_timing(f"{cluster} {sample_id} pdf baf panel", stage_start)

    event_ax = axes[3]
    stage_start = perf_counter()
    _plot_event_marginal_panel(
        event_ax,
        d_bin_mids,
        d_bar_widths,
        plot_event_probs,
        xform,
        locus,
        chrom,
        event_svtype,
    )
    _print_timing(f"{cluster} {sample_id} pdf event panel", stage_start)

    stage_start = perf_counter()
    plt.tight_layout()
    _print_timing(f"{cluster} {sample_id} pdf layout", stage_start)
    stage_start = perf_counter()
    pdf.savefig(fig)
    _print_timing(f"{cluster} {sample_id} pdf savefig", stage_start)
    plt.close()
    _print_timing(f"{cluster} {sample_id} pdf total", page_start)
    return True


def _get_event_probability_column(svtype: str) -> str:
    """Return the event-marginal column matching one SV type."""
    if svtype == "DEL":
        return "prob_del_event"
    if svtype == "DUP":
        return "prob_dup_event"
    raise ValueError(f"Unsupported svtype: {svtype}")


def _sorted_viterbi_segments(segments: Optional[List[dict]]) -> List[dict]:
    """Return Viterbi segments sorted left-to-right with stable tie breaks."""
    if not segments:
        return []
    return sorted(
        segments,
        key=lambda seg: (
            int(seg["start"]),
            int(seg["end"]),
            int(seg["cn_state"]),
        ),
    )


def _build_viterbi_trace_coords(
    segments: Optional[List[dict]],
    xform,
    vertical_offset: float = 0.0,
) -> Tuple[List[float], List[float]]:
    """Build a monotonic step trace for Viterbi segments.

    Adjacent segments are drawn as clean vertical steps. Any genomic gaps
    between segments are bridged horizontally using the previous state so the
    trace stays visually continuous across masked bins.
    """
    sorted_segments = _sorted_viterbi_segments(segments)
    if not sorted_segments:
        return [], []

    trace_x: List[float] = []
    trace_y: List[float] = []
    prev_end: Optional[int] = None
    prev_cn: Optional[float] = None

    for seg in sorted_segments:
        seg_start = int(seg["start"])
        seg_end = int(seg["end"])
        seg_cn = float(seg["cn_state"]) + vertical_offset
        start_x = xform(seg_start)
        end_x = xform(seg_end)

        if not trace_x:
            trace_x.extend([start_x, end_x])
            trace_y.extend([seg_cn, seg_cn])
        else:
            if prev_end is not None and seg_start > prev_end:
                trace_x.append(start_x)
                trace_y.append(prev_cn)

            if prev_cn != seg_cn:
                trace_x.append(start_x)
                trace_y.append(seg_cn)

            trace_x.append(end_x)
            trace_y.append(seg_cn)

        prev_end = seg_end
        prev_cn = seg_cn

    return trace_x, trace_y


def _plot_viterbi_trace(
    ax,
    segments: Optional[List[dict]],
    xform,
    color: str,
    label: str,
    linewidth: float,
    alpha: float,
    zorder: int,
    vertical_offset: float = 0.0,
    rasterized: bool = False,
) -> None:
    """Plot one Viterbi step trace using the shared carrier-plot code path."""
    trace_x, trace_y = _build_viterbi_trace_coords(
        segments,
        xform,
        vertical_offset=vertical_offset,
    )
    if not trace_x:
        return
    ax.plot(
        trace_x,
        trace_y,
        color=color,
        linewidth=linewidth,
        alpha=alpha,
        zorder=zorder,
        label=label,
        solid_joinstyle="round",
        solid_capstyle="round",
        rasterized=rasterized,
    )


# =============================================================================
# Raw-depth helpers
# =============================================================================


def estimate_lowres_bin_size(df: pd.DataFrame) -> float:
    """Return the median bin size (bp) across all rows of *df*."""
    return float(np.median((df["End"] - df["Start"]).values))


def _plot_depth_bars_with_baf(
    ax,
    x_positions: np.ndarray,
    bar_widths: np.ndarray,
    depth_values: np.ndarray,
    minor_baf_values: Optional[np.ndarray] = None,
    baf_site_counts: Optional[np.ndarray] = None,
    zorder: int = 3,
    rasterized: bool = False,
) -> None:
    """Draw normalized-depth bars, optionally split by modeled minor BAF.

    When BAF summaries are available for a bin/sample pair, the total depth bar
    is split into minor- and major-allele components using the same minor AF
    summary that the model consumes. Bins lacking usable BAF evidence fall back
    to the legacy single-color depth bar.
    """
    bar_widths = np.asarray(bar_widths, dtype=float) * 0.9
    depth_values = np.asarray(depth_values, dtype=float)

    if minor_baf_values is None or baf_site_counts is None:
        ax.bar(
            x_positions,
            depth_values,
            width=bar_widths,
            alpha=0.6,
            color="steelblue",
            edgecolor="none",
            zorder=zorder,
            label="Normalized depth",
            rasterized=rasterized,
        )
        return

    minor_baf_values = np.asarray(minor_baf_values, dtype=float)
    baf_site_counts = np.asarray(baf_site_counts)
    valid_baf = np.logical_and.reduce([
        np.isfinite(minor_baf_values),
        np.isfinite(depth_values),
        baf_site_counts > 0,
    ])

    if not np.any(valid_baf):
        ax.bar(
            x_positions,
            depth_values,
            width=bar_widths,
            alpha=0.6,
            color="steelblue",
            edgecolor="none",
            zorder=zorder,
            label="Normalized depth",
            rasterized=rasterized,
        )
        return

    clipped_minor_baf = np.clip(minor_baf_values, 0.0, 0.5)
    minor_depth = np.where(valid_baf, depth_values * clipped_minor_baf, 0.0)
    major_depth = np.where(valid_baf, depth_values - minor_depth, 0.0)

    if np.any(~valid_baf):
        ax.bar(
            x_positions[~valid_baf],
            depth_values[~valid_baf],
            width=bar_widths[~valid_baf],
            alpha=0.6,
            color="gray",
            edgecolor="none",
            zorder=zorder,
            label="Depth (no BAF)",
            rasterized=rasterized,
        )

    ax.bar(
        x_positions[valid_baf],
        minor_depth[valid_baf],
        width=bar_widths[valid_baf],
        alpha=0.75,
        color="#F4A261",
        edgecolor="none",
        zorder=zorder,
        label="Minor-allele depth",
        rasterized=rasterized,
    )
    ax.bar(
        x_positions[valid_baf],
        major_depth[valid_baf],
        bottom=minor_depth[valid_baf],
        width=bar_widths[valid_baf],
        alpha=0.65,
        color="#4C78A8",
        edgecolor="none",
        zorder=zorder,
        label="Major-allele depth",
        rasterized=rasterized,
    )


def _aligned_region_sample_vector(
    region_df: pd.DataFrame,
    value_df: Optional[pd.DataFrame],
    sample_id: str,
) -> Optional[np.ndarray]:
    """Return one sample vector aligned to ``region_df`` rows."""
    if value_df is None or sample_id not in value_df.columns:
        return None

    key_cols = ["Cluster", "Chr", "Start", "End"]
    if len(region_df) == len(value_df):
        region_keys = region_df[key_cols].reset_index(drop=True)
        value_keys = value_df[key_cols].reset_index(drop=True)
        if region_keys.equals(value_keys):
            return value_df[sample_id].to_numpy(copy=True)

    aligned = region_df[key_cols].merge(
        value_df[key_cols + [sample_id]],
        on=key_cols,
        how="left",
    )
    return aligned[sample_id].to_numpy(copy=True)


def _extract_sample_baf_vectors(
    region_df: pd.DataFrame,
    minor_baf_df: Optional[pd.DataFrame],
    baf_sites_df: Optional[pd.DataFrame],
    sample_id: str,
) -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
    """Align per-bin BAF summaries to a region/sample depth slice."""
    return (
        _aligned_region_sample_vector(region_df, minor_baf_df, sample_id),
        _aligned_region_sample_vector(region_df, baf_sites_df, sample_id),
    )


def _plot_baf_signal_panel(
    ax,
    x_positions: np.ndarray,
    bar_widths: np.ndarray,
    minor_baf_values: Optional[np.ndarray],
    baf_site_counts: Optional[np.ndarray],
    xform: FlankCompressor,
    locus: GDLocus,
    chrom: str,
    show_xlabel: bool = True,
    rasterized: bool = False,
) -> None:
    """Render the modeled per-bin minor-BAF signal for one sample."""
    ax.set_xlim(0.0, xform.d_end)
    ax.set_ylim(0.0, 1.0)

    for y_value in (1.0 / 3.0, 0.5, 2.0 / 3.0):
        ax.axhline(y_value, color="gray", linestyle=":", linewidth=0.8, alpha=0.25, zorder=0)

    if minor_baf_values is None or baf_site_counts is None:
        ax.text(
            0.5,
            0.5,
            "No BAF signal available",
            transform=ax.transAxes,
            ha="center",
            va="center",
            color="dimgray",
            fontsize=9,
        )
        if show_xlabel:
            ax.set_xlabel(f"Position on {chrom}")
        ax.set_ylabel("Minor BAF")
        xform.format_genomic_ticks(ax, locus.breakpoints)
        return

    minor_baf_values = np.asarray(minor_baf_values, dtype=float)
    baf_site_counts = np.asarray(baf_site_counts)
    bar_widths = np.asarray(bar_widths, dtype=float) * 0.85
    valid_baf = np.logical_and.reduce([
        np.isfinite(minor_baf_values),
        baf_site_counts > 0,
    ])

    if np.any(valid_baf):
        clipped_minor_baf = np.clip(minor_baf_values[valid_baf], 0.0, 1.0)
        ax.bar(
            x_positions[valid_baf],
            clipped_minor_baf,
            width=bar_widths[valid_baf],
            color="#7B61A8",
            alpha=0.25,
            edgecolor="none",
            zorder=1,
            rasterized=rasterized,
        )
        baf_stems = ax.vlines(
            x_positions[valid_baf],
            0.0,
            clipped_minor_baf,
            colors="#7B61A8",
            linewidth=0.7,
            alpha=0.5,
            zorder=1.5,
            rasterized=rasterized,
        )
        if hasattr(baf_stems, "set_capstyle"):
            baf_stems.set_capstyle("round")
        ax.scatter(
            x_positions[valid_baf],
            clipped_minor_baf,
            s=10,
            color="#5B3F8C",
            alpha=0.8,
            linewidths=0,
            zorder=2,
            rasterized=rasterized,
        )
    else:
        ax.text(
            0.5,
            0.5,
            "No BAF-supported bins",
            transform=ax.transAxes,
            ha="center",
            va="center",
            color="dimgray",
            fontsize=9,
        )

    if show_xlabel:
        ax.set_xlabel(f"Position on {chrom}")
    ax.set_ylabel("Minor BAF")
    xform.format_genomic_ticks(ax, locus.breakpoints)


def _extract_sample_event_probabilities(
    region_df: pd.DataFrame,
    event_probability_df: Optional[pd.DataFrame],
    sample_id: str,
) -> Optional[np.ndarray]:
    """Align per-bin event marginals to a region/sample depth slice."""
    return _aligned_region_sample_vector(region_df, event_probability_df, sample_id)


def _plot_event_marginal_panel(
    ax,
    x_positions: np.ndarray,
    bar_widths: np.ndarray,
    event_probabilities: Optional[np.ndarray],
    xform: FlankCompressor,
    locus: GDLocus,
    chrom: str,
    svtype: Optional[str],
    rasterized: bool = False,
) -> None:
    """Render per-bin event marginal probabilities for one sample."""
    ax.set_xlim(0.0, xform.d_end)
    ax.set_ylim(0.0, 1.0)
    ax.axhline(0.5, color="gray", linestyle=":", linewidth=0.8, alpha=0.25, zorder=0)

    if event_probabilities is None or svtype is None:
        ax.text(
            0.5,
            0.5,
            "No event marginals available",
            transform=ax.transAxes,
            ha="center",
            va="center",
            color="dimgray",
            fontsize=9,
        )
        ax.set_xlabel(f"Position on {chrom}")
        ax.set_ylabel("Event Prob.")
        xform.format_genomic_ticks(ax, locus.breakpoints)
        return

    event_probabilities = np.asarray(event_probabilities, dtype=float)
    valid = np.isfinite(event_probabilities)
    color = "#C23B22" if svtype == "DEL" else "#2A6FBB"

    if np.any(valid):
        clipped = np.clip(event_probabilities[valid], 0.0, 1.0)
        ax.bar(
            x_positions[valid],
            clipped,
            width=np.asarray(bar_widths, dtype=float)[valid] * 0.85,
            color=color,
            alpha=0.25,
            edgecolor="none",
            zorder=1,
            rasterized=rasterized,
        )
        ax.plot(
            x_positions[valid],
            clipped,
            color=color,
            linewidth=1.2,
            alpha=0.9,
            zorder=2,
            rasterized=rasterized,
        )
    else:
        ax.text(
            0.5,
            0.5,
            "No finite event marginals",
            transform=ax.transAxes,
            ha="center",
            va="center",
            color="dimgray",
            fontsize=9,
        )

    ax.set_xlabel(f"Position on {chrom}")
    ax.set_ylabel(f"P({svtype} on ≥1 hap)")
    xform.format_genomic_ticks(ax, locus.breakpoints)


def _build_raw_region_df(
    locus: GDLocus,
    region_start: int,
    region_end: int,
    raw_counts_df: Optional[pd.DataFrame],
    sample_cols: List[str],
    raw_sample_medians: Dict[str, float],
    lowres_median_bin_size: float,
    processed_region_df: Optional[pd.DataFrame] = None,
    highres_path: Optional[str] = None,
    size_threshold_factor: float = 10.0,
    highres_max_bins: Optional[int] = None,
) -> Optional[pd.DataFrame]:
    """
    Build a normalised raw-depth DataFrame for a locus's visible region,
    optionally substituting high-resolution bins for small intervals/flanks.

    For each GD body interval and flank, the span is compared against
    ``size_threshold_factor × lowres_median_bin_size``.  When the span is
    smaller than that threshold **and** *highres_path* is provided, bins
    for that sub-region are queried on demand via pysam.

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
        processed_region_df: Processed region DataFrame used in the model /
            plotted depth bars. When provided, this is inspected to detect
            locus-wide high-resolution replacement during preprocess.
        highres_path: Path to a bgzipped, tabix-indexed high-resolution
            read-count file (.tsv.gz + .tbi).  ``None`` disables high-res
            substitution entirely.
        size_threshold_factor: Regions whose genomic span is less than
            this multiple of *lowres_median_bin_size* are served from
            the high-res file (default 10).
        highres_max_bins: Optional cap on the number of queried high-res bins
            retained for plotting. When set, consecutive bins are merged by
            mean before normalisation so plotting can be much faster while
            preserving coarse visual structure.

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

    def _count_interval_bins(df: pd.DataFrame, start: int, end: int) -> int:
        if df is None or len(df) == 0:
            return 0
        mids = (df["Start"].values + df["End"].values) / 2
        return int(((mids >= start) & (mids < end)).sum())

    def _query_and_normalize_highres(query_start: int, query_end: int) -> Optional[pd.DataFrame]:
        hr_raw = query_highres_bins(
            highres_path,
            chrom,
            query_start,
            query_end,
            common_samples,
            max_bins=highres_max_bins,
        )
        if len(hr_raw) == 0:
            return None
        if highres_max_bins is not None:
            print(
                f"    [high-res plot] coarsened query {chrom}:{query_start:,}-{query_end:,} "
                f"to {len(hr_raw)} bins for plotting"
            )
        normed = normalize_highres_bins(
            hr_raw, common_samples, column_medians_arr, lowres_median_bin_size,
        )
        keep_cols = ["Chr", "Start", "End"] + [
            s for s in common_samples if s in normed.columns
        ]
        return (
            normed.reset_index()[keep_cols]
            .sort_values("Start")
            .reset_index(drop=True)
            .copy()
        )

    # --- Step 1: normalised low-res base covering the full visible region ---
    base_df = pd.DataFrame()
    if raw_counts_df is not None:
        mask = (
            (raw_counts_df["Chr"] == chrom) &
            (raw_counts_df["End"] > region_start) &
            (raw_counts_df["Start"] < region_end)
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

    # --- Step 1b: mirror locus-wide high-res replacement from preprocess ---
    # Preprocess switches the *entire locus* to high-res when any body interval
    # is under-covered and the high-res result improves interval coverage. The
    # processed bars in region_df therefore come from high-res across A-C / C-D
    # too, not just the originally under-covered interval. To keep the plotted
    # raw signal comparable, detect that condition and source the full visible
    # region from high-res here as well.
    if processed_region_df is not None and len(base_df) > 0:
        locus_wide_highres = False
        for iv_start, iv_end, iv_name in locus.get_intervals():
            processed_count = _count_interval_bins(processed_region_df, iv_start, iv_end)
            lowres_count = _count_interval_bins(base_df, iv_start, iv_end)
            if processed_count > lowres_count:
                print(
                    f"    [high-res plot] full-locus replacement triggered: "
                    f"processed {iv_name} has {processed_count} bins vs "
                    f"{lowres_count} low-res bins"
                )
                locus_wide_highres = True
                break

        if locus_wide_highres:
            full_region_highres = _query_and_normalize_highres(region_start, region_end)
            if full_region_highres is not None:
                print(
                    f"    [high-res plot] using high-res across full visible region "
                    f"{chrom}:{region_start:,}-{region_end:,} "
                    f"({len(full_region_highres)} bins)"
                )
                return full_region_highres
            print(
                "    [high-res plot] full-region query returned no bins; "
                "falling back to low-res / partial substitution"
            )

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
        normed = _query_and_normalize_highres(r_start, r_end)
        if normed is None:
            print(f"    [high-res plot] {r_name}: no bins returned, keeping low-res")
            continue
        highres_parts.append(normed)
        excluded_ranges.append((r_start, r_end))
        print(
            f"    [high-res plot] {r_name} ({region_size:,} bp "
            f"< {threshold_bp:,.0f} bp threshold): {len(normed)} high-res bins"
        )

    if not highres_parts:
        return base_df if len(base_df) > 0 else None

    # --- Step 3: splice out low-res bins covered by high-res regions ---
    if len(base_df) > 0:
        keep = pd.Series(True, index=base_df.index)
        for ex_s, ex_e in excluded_ranges:
            # True interval overlap: drop low-res bins that overlap any
            # high-res replacement region (no midpoint heuristic).
            keep &= ~((base_df["Start"] < ex_e) & (base_df["End"] > ex_s))
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


def _rebin_region_df(
    df: pd.DataFrame,
    locus: GDLocus,
    max_bins_per_region: int = 100,
) -> pd.DataFrame:
    """Rebin a raw-depth DataFrame so that no single sub-region (body or
    individual flank) exceeds *max_bins_per_region* bins.

    Flanks and body are rebinned independently so that the boundary
    between them is preserved.  Within each sub-region, consecutive
    bins are grouped and their depths averaged, while ``Start`` / ``End``
    take the min / max of each group so the merged bin spans the full
    range.

    Args:
        df: Normalised raw-depth DataFrame with ``Chr``, ``Start``,
            ``End`` and per-sample columns, sorted by ``Start``.
        locus: GDLocus defining body boundaries.
        max_bins_per_region: Target maximum number of bins for each
            sub-region.

    Returns:
        Rebinned DataFrame with the same column schema as *df*.
    """
    if df is None or len(df) == 0:
        return df

    sample_cols = [c for c in df.columns if c not in ("Chr", "Start", "End")]

    mids = (df["Start"].values + df["End"].values) / 2
    left_mask = mids < locus.start
    body_mask = (mids >= locus.start) & (mids < locus.end)
    right_mask = mids >= locus.end

    parts: List[pd.DataFrame] = []
    for mask in (left_mask, body_mask, right_mask):
        sub = df[mask].copy()
        if len(sub) == 0:
            continue
        if len(sub) <= max_bins_per_region:
            parts.append(sub)
            continue
        # Rebin: split into max_bins_per_region groups of roughly equal size
        n = len(sub)
        group_ids = np.arange(n) * max_bins_per_region // n
        grouped = sub.assign(_g=group_ids).groupby("_g", sort=True)
        agg = grouped.agg(
            Chr=("Chr", "first"),
            Start=("Start", "min"),
            End=("End", "max"),
            **{s: (s, "mean") for s in sample_cols},
        ).reset_index(drop=True)
        parts.append(agg)

    if not parts:
        return df
    return pd.concat(parts, ignore_index=True).sort_values("Start").reset_index(drop=True)


def _allocate_segment_bin_targets(segment_counts: List[int], max_total_bins: int) -> List[int]:
    """Allocate a total display-bin budget across left/body/right segments."""
    positive_indices = [index for index, count in enumerate(segment_counts) if count > 0]
    if not positive_indices:
        return [0] * len(segment_counts)

    total_bins = sum(segment_counts)
    if total_bins <= max_total_bins:
        return list(segment_counts)

    raw_targets = [count * max_total_bins / total_bins if count > 0 else 0.0 for count in segment_counts]
    targets = [0] * len(segment_counts)
    for index in positive_indices:
        targets[index] = min(segment_counts[index], max(1, int(np.floor(raw_targets[index]))))

    current_total = sum(targets)
    if current_total < max_total_bins:
        residual_order = sorted(
            positive_indices,
            key=lambda index: (raw_targets[index] - targets[index], segment_counts[index]),
            reverse=True,
        )
        while current_total < max_total_bins:
            advanced = False
            for index in residual_order:
                if targets[index] < segment_counts[index]:
                    targets[index] += 1
                    current_total += 1
                    advanced = True
                    if current_total == max_total_bins:
                        break
            if not advanced:
                break
    elif current_total > max_total_bins:
        shrink_order = sorted(
            positive_indices,
            key=lambda index: (targets[index] - raw_targets[index], targets[index]),
            reverse=True,
        )
        while current_total > max_total_bins:
            reduced = False
            for index in shrink_order:
                if targets[index] > 1:
                    targets[index] -= 1
                    current_total -= 1
                    reduced = True
                    if current_total == max_total_bins:
                        break
            if not reduced:
                break

    return targets


def _rebin_aligned_region_dfs_for_display(
    locus: GDLocus,
    region_df: pd.DataFrame,
    aligned_dfs: List[Optional[pd.DataFrame]],
    max_total_bins: int = PDF_MAX_SIGNAL_BINS,
) -> Tuple[pd.DataFrame, List[Optional[pd.DataFrame]]]:
    """Rebin aligned plot DataFrames together for PDF display only.

    The same grouping is applied to every compatible DataFrame so depth, BAF,
    and event traces remain row-aligned after coarsening.
    """
    if region_df is None or len(region_df) == 0 or len(region_df) <= max_total_bins:
        return region_df, aligned_dfs

    key_cols = [column for column in ("Cluster", "Chr", "Start", "End") if column in region_df.columns]
    mids = (region_df["Start"].values + region_df["End"].values) / 2
    segment_masks = [
        mids < locus.start,
        (mids >= locus.start) & (mids < locus.end),
        mids >= locus.end,
    ]
    segment_counts = [int(mask.sum()) for mask in segment_masks]
    segment_targets = _allocate_segment_bin_targets(segment_counts, max_total_bins)

    def _rebin_one(df: Optional[pd.DataFrame]) -> Optional[pd.DataFrame]:
        if df is None or len(df) == 0:
            return df
        if len(df) != len(region_df):
            return _rebin_region_df(df, locus, max_bins_per_region=max(1, max_total_bins // 3))
        if key_cols:
            region_keys = region_df[key_cols].reset_index(drop=True)
            df_keys = df[key_cols].reset_index(drop=True)
            if not region_keys.equals(df_keys):
                return _rebin_region_df(df, locus, max_bins_per_region=max(1, max_total_bins // 3))

        value_cols = [column for column in df.columns if column not in key_cols]
        rebinned_parts: List[pd.DataFrame] = []
        for mask, target in zip(segment_masks, segment_targets):
            sub = df[mask].copy()
            if len(sub) == 0:
                continue
            if target <= 0 or len(sub) <= target:
                rebinned_parts.append(sub)
                continue
            group_ids = np.arange(len(sub)) * target // len(sub)
            grouped = sub.assign(_g=group_ids).groupby("_g", sort=True)
            agg_spec = {}
            for column in key_cols:
                if column == "Start":
                    agg_spec[column] = (column, "min")
                elif column == "End":
                    agg_spec[column] = (column, "max")
                else:
                    agg_spec[column] = (column, "first")
            agg_spec.update({column: (column, "mean") for column in value_cols})
            rebinned_parts.append(grouped.agg(**agg_spec).reset_index(drop=True))

        if not rebinned_parts:
            return df
        return pd.concat(rebinned_parts, ignore_index=True).sort_values("Start").reset_index(drop=True)

    rebinned_region = _rebin_one(region_df)
    rebinned_aligned = [_rebin_one(df) for df in aligned_dfs]
    return rebinned_region, rebinned_aligned


def _coarsen_pdf_page_signals(
    locus: GDLocus,
    region_df: pd.DataFrame,
    sample_depth: np.ndarray,
    minor_baf_values: Optional[np.ndarray],
    baf_site_counts: Optional[np.ndarray],
    event_probabilities: Optional[np.ndarray],
    max_total_bins: int = PDF_MAX_SIGNAL_BINS,
) -> Tuple[pd.DataFrame, np.ndarray, Optional[np.ndarray], Optional[np.ndarray], Optional[np.ndarray]]:
    """Coarsen one PDF page's sample-specific signal arrays for display only."""
    if len(region_df) <= max_total_bins:
        return region_df, sample_depth, minor_baf_values, baf_site_counts, event_probabilities

    starts = region_df["Start"].to_numpy(dtype=int)
    ends = region_df["End"].to_numpy(dtype=int)
    mids = (starts + ends) / 2
    segment_masks = [
        mids < locus.start,
        (mids >= locus.start) & (mids < locus.end),
        mids >= locus.end,
    ]
    segment_counts = [int(mask.sum()) for mask in segment_masks]
    segment_targets = _allocate_segment_bin_targets(segment_counts, max_total_bins)

    sample_depth = np.asarray(sample_depth, dtype=float)
    minor_baf_arr = None if minor_baf_values is None else np.asarray(minor_baf_values, dtype=float)
    baf_sites_arr = None if baf_site_counts is None else np.asarray(baf_site_counts)
    event_probs_arr = None if event_probabilities is None else np.asarray(event_probabilities, dtype=float)

    rebinned_rows: List[dict] = []
    rebinned_depth: List[float] = []
    rebinned_minor_baf: List[float] = [] if minor_baf_arr is not None else None
    rebinned_baf_sites: List[float] = [] if baf_sites_arr is not None else None
    rebinned_event_probs: List[float] = [] if event_probs_arr is not None else None

    for mask, target in zip(segment_masks, segment_targets):
        indices = np.flatnonzero(mask)
        if len(indices) == 0:
            continue
        if target <= 0 or len(indices) <= target:
            for index in indices:
                rebinned_rows.append({
                    key: region_df.iloc[index][key]
                    for key in region_df.columns
                    if key in ("Cluster", "Chr", "Start", "End")
                })
                rebinned_depth.append(float(sample_depth[index]))
                if rebinned_minor_baf is not None:
                    rebinned_minor_baf.append(float(minor_baf_arr[index]))
                if rebinned_baf_sites is not None:
                    rebinned_baf_sites.append(float(baf_sites_arr[index]))
                if rebinned_event_probs is not None:
                    rebinned_event_probs.append(float(event_probs_arr[index]))
            continue

        group_ids = np.arange(len(indices)) * target // len(indices)
        for group_id in range(target):
            group_indices = indices[group_ids == group_id]
            if len(group_indices) == 0:
                continue
            rebinned_rows.append({
                key: (
                    region_df.iloc[group_indices[0]][key]
                    if key not in ("Start", "End") else
                    int(starts[group_indices].min()) if key == "Start" else int(ends[group_indices].max())
                )
                for key in region_df.columns
                if key in ("Cluster", "Chr", "Start", "End")
            })
            rebinned_depth.append(float(np.nanmean(sample_depth[group_indices])))
            if rebinned_minor_baf is not None:
                group_sites = baf_sites_arr[group_indices] if baf_sites_arr is not None else None
                valid_minor = np.isfinite(minor_baf_arr[group_indices])
                if group_sites is not None:
                    valid_minor &= np.isfinite(group_sites) & (group_sites > 0)
                if np.any(valid_minor):
                    if group_sites is not None:
                        weights = group_sites[group_indices * 0 + valid_minor] if False else group_sites[valid_minor]
                        values = minor_baf_arr[group_indices][valid_minor]
                        weight_sum = float(np.sum(weights))
                        rebinned_minor_baf.append(float(np.sum(values * weights) / weight_sum) if weight_sum > 0 else float(np.nanmean(values)))
                    else:
                        rebinned_minor_baf.append(float(np.nanmean(minor_baf_arr[group_indices][valid_minor])))
                else:
                    rebinned_minor_baf.append(np.nan)
            if rebinned_baf_sites is not None:
                valid_sites = baf_sites_arr[group_indices]
                rebinned_baf_sites.append(float(np.nansum(valid_sites[np.isfinite(valid_sites)])))
            if rebinned_event_probs is not None:
                valid_event = event_probs_arr[group_indices]
                rebinned_event_probs.append(float(np.nanmean(valid_event)))

    plot_region_df = pd.DataFrame(rebinned_rows)
    return (
        plot_region_df,
        np.asarray(rebinned_depth, dtype=float),
        None if rebinned_minor_baf is None else np.asarray(rebinned_minor_baf, dtype=float),
        None if rebinned_baf_sites is None else np.asarray(rebinned_baf_sites, dtype=float),
        None if rebinned_event_probs is None else np.asarray(rebinned_event_probs, dtype=float),
    )


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
    raw_heatmap_viz_bins: int = 400,
    raw_highres_max_bins: int = 800,
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

    overview_start = perf_counter()

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
        (raw_counts_df is not None or highres_path is not None) and
        raw_sample_medians is not None and
        lowres_median_bin_size is not None
    )
    raw_region_df: Optional[pd.DataFrame] = None
    if have_raw:
        stage_start = perf_counter()
        raw_region_df = _build_raw_region_df(
            locus, region_start, region_end,
            raw_counts_df, sample_cols, raw_sample_medians,
            lowres_median_bin_size,
            processed_region_df=region_df,
            highres_path=highres_path,
            highres_max_bins=raw_highres_max_bins,
        )
        if raw_region_df is not None:
            raw_region_df = _rebin_region_df(raw_region_df, locus)
        if raw_region_df is None or len(raw_region_df) == 0:
            have_raw = False
            raw_region_df = None
        _print_timing(f"{locus.cluster} raw-region prep", stage_start)

    # ---- figure layout ----
    figure_scale = 1.0
    n_cols = 2 if have_raw else 1
    carrier_height = min(6, max(1, n_carriers * 0.15)) if n_carriers > 0 else 0.5
    non_carrier_height = min(6, max(1, n_non_carriers * 0.15)) if n_non_carriers > 0 else 0.5
    depth_panel_height = 2.0
    indiv_depth_panel_height = 2.0
    fig_height = (
        3 + carrier_height + non_carrier_height +
        n_ploidy_panels * depth_panel_height +
        indiv_depth_panel_height
    )

    height_ratios = (
        [1, carrier_height, non_carrier_height] +
        [depth_panel_height] * n_ploidy_panels +
        [indiv_depth_panel_height]
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
        stage_start = perf_counter()
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
            heatmap_viz_bins=raw_heatmap_viz_bins,
        )
        _print_timing(f"{locus.cluster} raw column draw", stage_start)

        stage_start = perf_counter()
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
            heatmap_viz_bins=raw_heatmap_viz_bins,
        )
        _print_timing(f"{locus.cluster} processed column draw", stage_start)
    else:
        stage_start = perf_counter()
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
            heatmap_viz_bins=raw_heatmap_viz_bins,
        )
        _print_timing(f"{locus.cluster} single column draw", stage_start)

    stage_start = perf_counter()
    plt.tight_layout(h_pad=0.3, w_pad=0.5)
    fig.subplots_adjust(hspace=0.08)
    _print_timing(f"{locus.cluster} layout", stage_start)

    locus_dir = os.path.join(output_dir, "locus_plots")
    os.makedirs(locus_dir, exist_ok=True)
    filename = f"{_sanitize_plot_label(locus.cluster)}_overview.png"
    stage_start = perf_counter()
    fig.savefig(os.path.join(locus_dir, filename), dpi=100)
    _print_timing(f"{locus.cluster} savefig", stage_start)
    plt.close()
    _print_timing(f"{locus.cluster} overview total", overview_start)
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

            label = f"{call['svtype']}"
            matched_haplotype = call.get("matched_haplotype")
            hap_cn_state = call.get("hap_cn_state")
            if pd.notna(matched_haplotype) and pd.notna(hap_cn_state):
                label += f" H{int(matched_haplotype)} CN={int(hap_cn_state)}"
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

    sample_dir = os.path.join(output_dir, "sample_plots", _sanitize_plot_label(locus.cluster))
    os.makedirs(sample_dir, exist_ok=True)
    filename = f"{_sanitize_plot_label(sample_id)}.png"
    plt.savefig(os.path.join(sample_dir, filename), dpi=150, bbox_inches="tight")
    plt.close()


def plot_carrier_summary(calls_df: pd.DataFrame, output_dir: str):
    """Create summary plots of carrier calls across all loci."""
    carriers = calls_df[calls_df["is_carrier"]].copy()
    confidence_column = _get_confidence_column(calls_df)
    confidence_label = _get_confidence_label(confidence_column)

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
    all_scores = carriers[confidence_column].dropna()
    bins = np.linspace(all_scores.min(), all_scores.max(), 21)
    for svtype, color in [("DEL", "red"), ("DUP", "blue")]:
        sv_data = carriers[carriers["svtype"] == svtype]
        if len(sv_data) > 0:
            ax.hist(sv_data[confidence_column], bins=bins, alpha=0.6, color=color,
                    label=f"{svtype} (n={len(sv_data)})", edgecolor="black")
    ax.set_xlabel(confidence_label)
    ax.set_ylabel("Count")
    ax.set_title(f"{confidence_label} Distribution in Carriers")
    ax.legend()

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "carrier_summary.png"), dpi=150, bbox_inches="tight")
    plt.close()
    print("Created: carrier_summary.png")


def plot_confidence_distribution(calls_df: pd.DataFrame, output_dir: str):
    """Plot distribution of confidence scores."""
    if len(calls_df) == 0:
        print("No calls to plot confidence distribution.")
        return

    confidence_column = _get_confidence_column(calls_df)
    confidence_label = _get_confidence_label(confidence_column)
    figure_scale = 1.0
    fig, axes = plt.subplots(1, 2, figsize=(8 * figure_scale, 4 * figure_scale))

    ax = axes[0]
    carriers = calls_df[calls_df["is_carrier"]]
    non_carriers = calls_df[~calls_df["is_carrier"]]
    all_scores = calls_df[confidence_column].dropna()
    bins = np.linspace(all_scores.min(), all_scores.max(), 31)
    ax.hist(non_carriers[confidence_column], bins=bins, alpha=0.6, label="Non-carriers",
            color="gray", edgecolor="black")
    ax.hist(carriers[confidence_column], bins=bins, alpha=0.6, label="Carriers",
            color="green", edgecolor="black")
    ax.set_xlabel(confidence_label)
    ax.set_ylabel("Count")
    ax.set_yscale("log")
    ax.legend()

    ax = axes[1]
    non_carriers_mask = ~calls_df["is_carrier"]
    ax.scatter(calls_df.loc[non_carriers_mask, "mean_depth"],
               calls_df.loc[non_carriers_mask, confidence_column],
               c="gray", alpha=0.1, s=10, label="Non-carriers")
    carriers_mask = calls_df["is_carrier"]
    ax.scatter(calls_df.loc[carriers_mask, "mean_depth"],
               calls_df.loc[carriers_mask, confidence_column],
               c="green", alpha=0.5, s=10, label="Carriers")
    ax.set_xlabel("Mean Depth (normalized)")
    ax.set_ylabel(confidence_label)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "confidence_distribution.png"), dpi=150, bbox_inches="tight")
    plt.close()
    print("Created: confidence_distribution.png")


def create_carrier_pdf(
    calls_df: pd.DataFrame,
    depth_df: pd.DataFrame,
    loci_by_cluster: Dict[str, GDLocus],
    gtf: Optional[GTFParser],
    segdup: Optional[SegDupAnnotation],
    output_dir: str,
    event_marginals_df: Optional[pd.DataFrame] = None,
    event_del_df: Optional[pd.DataFrame] = None,
    event_dup_df: Optional[pd.DataFrame] = None,
    minor_baf_df: Optional[pd.DataFrame] = None,
    baf_sites_df: Optional[pd.DataFrame] = None,
    padding: int = 50000,
    min_gene_label_spacing: float = 0.05,
    raw_counts_df: Optional[pd.DataFrame] = None,
    gaps: Optional[GapsAnnotation] = None,
    flank_scale: float = 0.20,
    lowres_median_bin_size: Optional[float] = None,
    highres_path: Optional[str] = None,
    viterbi_data: Optional[ViterbiOverlayData] = None,
    confidence_threshold: float = 0.5,
):
    """Create a PDF with plots for calls above the confidence threshold."""
    confidence_column = _get_confidence_column(calls_df)
    carriers_mask = _build_carrier_best_match_mask(calls_df)
    carriers = calls_df[
        carriers_mask & (calls_df[confidence_column] >= confidence_threshold)
    ].copy()

    if len(carriers) == 0:
        print(
            "No calls at or above the carrier PDF confidence threshold "
            f"({confidence_threshold:.2f})."
        )
        return

    raw_sample_medians = _compute_raw_sample_medians(raw_counts_df)
    if raw_sample_medians:
        print(f"  Computed raw autosomal medians for {len(raw_sample_medians)} samples")

    if lowres_median_bin_size is None and raw_counts_df is not None:
        lowres_median_bin_size = estimate_lowres_bin_size(raw_counts_df)

    pdf_path = os.path.join(output_dir, "carrier_plots.pdf")
    print(f"Creating carrier PDF: {pdf_path}")
    pdf_start = perf_counter()

    with PdfPages(pdf_path) as pdf:
        for cluster in sorted(carriers["cluster"].unique()):
            cluster_start = perf_counter()
            locus = loci_by_cluster.get(cluster)
            if locus is None:
                continue

            if not locus.breakpoints:
                print(f"  Warning: No breakpoints for {cluster}, skipping")
                continue

            cluster_carriers = carriers[carriers["cluster"] == cluster]["sample"].unique()
            print(
                f"  Adding {len(cluster_carriers)} high-confidence sample(s) for {cluster} "
                f"(threshold={confidence_threshold:.2f})"
            )

            chrom = locus.chrom
            stage_start = perf_counter()
            locus_mask = (
                (depth_df["Cluster"] == cluster) &
                (depth_df["Chr"] == chrom)
            )
            region_df = depth_df[locus_mask].sort_values("Start")
            baf_region_df = None
            baf_sites_region_df = None
            event_del_region_df = None
            event_dup_region_df = None
            if minor_baf_df is not None:
                baf_region_df = minor_baf_df[locus_mask].sort_values("Start")
            if baf_sites_df is not None:
                baf_sites_region_df = baf_sites_df[locus_mask].sort_values("Start")
            if event_del_df is not None:
                event_del_region_df = event_del_df[locus_mask].sort_values("Start")
            if event_dup_df is not None:
                event_dup_region_df = event_dup_df[locus_mask].sort_values("Start")
            if len(region_df) == 0:
                continue
            if len(region_df) > PDF_MAX_SIGNAL_BINS:
                region_df, rebinned_frames = _rebin_aligned_region_dfs_for_display(
                    locus,
                    region_df,
                    [baf_region_df, baf_sites_region_df, event_del_region_df, event_dup_region_df],
                    max_total_bins=PDF_MAX_SIGNAL_BINS,
                )
                (
                    baf_region_df,
                    baf_sites_region_df,
                    event_del_region_df,
                    event_dup_region_df,
                ) = rebinned_frames
            _print_timing(f"{cluster} carrier region slice", stage_start)

            region_start = int(region_df["Start"].min())
            region_end = int(region_df["End"].max())

            # Build coordinate transform for this locus
            stage_start = perf_counter()
            xform = FlankCompressor(region_start, region_end,
                                    locus.start, locus.end,
                                    flank_scale=flank_scale)
            _print_timing(f"{cluster} carrier xform", stage_start)

            stage_start = perf_counter()
            cluster_calls_df = calls_df[calls_df["cluster"] == cluster]
            _raw_region = None
            if all([
                raw_counts_df is not None or highres_path is not None,
                bool(raw_sample_medians),
                lowres_median_bin_size is not None,
            ]):
                _raw_region = _build_raw_region_df(
                    locus, region_start, region_end,
                    raw_counts_df, list(raw_sample_medians.keys()),
                    raw_sample_medians, lowres_median_bin_size,
                    processed_region_df=region_df,
                    highres_path=highres_path,
                )
                if _raw_region is not None:
                    _raw_region = _rebin_region_df(_raw_region, locus)
            _print_timing(f"{cluster} carrier raw-region prep", stage_start)

            rendered_pages = 0
            for sample_id in sorted(cluster_carriers):
                rendered = _render_pdf_sample_page(
                    pdf,
                    str(sample_id),
                    cluster,
                    locus,
                    region_df,
                    cluster_calls_df,
                    confidence_column,
                    confidence_threshold,
                    gtf,
                    segdup,
                    min_gene_label_spacing,
                    xform,
                    _raw_region,
                    baf_region_df,
                    baf_sites_region_df,
                    event_del_region_df,
                    event_dup_region_df,
                    gaps,
                    viterbi_data,
                )
                if not rendered:
                    continue
                rendered_pages += 1

            _print_timing(
                f"{cluster} carrier cluster total ({rendered_pages} pages)",
                cluster_start,
            )

    print(f"  Saved carrier PDF: {pdf_path}")
    _print_timing("carrier pdf total", pdf_start)


def create_eval_category_pdfs(
    eval_report_df: pd.DataFrame,
    calls_df: pd.DataFrame,
    depth_df: pd.DataFrame,
    loci_by_cluster: Dict[str, GDLocus],
    gtf: Optional[GTFParser],
    segdup: Optional[SegDupAnnotation],
    output_dir: str,
    event_marginals_df: Optional[pd.DataFrame] = None,
    event_del_df: Optional[pd.DataFrame] = None,
    event_dup_df: Optional[pd.DataFrame] = None,
    minor_baf_df: Optional[pd.DataFrame] = None,
    baf_sites_df: Optional[pd.DataFrame] = None,
    padding: int = 50000,
    min_gene_label_spacing: float = 0.05,
    raw_counts_df: Optional[pd.DataFrame] = None,
    gaps: Optional[GapsAnnotation] = None,
    flank_scale: float = 0.20,
    lowres_median_bin_size: Optional[float] = None,
    highres_path: Optional[str] = None,
    viterbi_data: Optional[ViterbiOverlayData] = None,
    confidence_threshold: float = 0.5,
):
    """Create TP/FP/FN review PDFs when an eval report is provided."""
    confidence_column = _get_confidence_column(calls_df)
    raw_sample_medians = _compute_raw_sample_medians(raw_counts_df)
    if raw_sample_medians:
        print(f"  Computed raw autosomal medians for {len(raw_sample_medians)} samples")

    if lowres_median_bin_size is None and raw_counts_df is not None:
        lowres_median_bin_size = estimate_lowres_bin_size(raw_counts_df)

    specs_by_category = _build_eval_pdf_specs(
        eval_report_df,
        calls_df,
        loci_by_cluster,
        confidence_threshold,
    )
    output_paths = {
        "true_positives": os.path.join(output_dir, "true_positives.pdf"),
        "false_positives": os.path.join(output_dir, "false_positives.pdf"),
        "false_negatives": os.path.join(output_dir, "false_negatives.pdf"),
    }

    for category_name, pdf_path in output_paths.items():
        page_specs = specs_by_category.get(category_name, [])
        if not page_specs:
            print(f"No pages for {os.path.basename(pdf_path)}")
            continue

        print(f"Creating PDF: {pdf_path}")
        category_start = perf_counter()
        with PdfPages(pdf_path) as pdf:
            pages_by_cluster: Dict[str, List[dict]] = {}
            for spec in page_specs:
                pages_by_cluster.setdefault(str(spec["cluster"]), []).append(spec)

            for cluster in sorted(pages_by_cluster):
                cluster_start = perf_counter()
                locus = loci_by_cluster.get(cluster)
                if locus is None or not locus.breakpoints:
                    print(f"  Warning: No breakpoints for {cluster}, skipping")
                    continue

                chrom = locus.chrom
                stage_start = perf_counter()
                locus_mask = (
                    (depth_df["Cluster"] == cluster) &
                    (depth_df["Chr"] == chrom)
                )
                region_df = depth_df[locus_mask].sort_values("Start")
                if len(region_df) == 0:
                    continue
                _print_timing(f"{cluster} {category_name} region slice", stage_start)

                baf_region_df = None
                baf_sites_region_df = None
                event_del_region_df = None
                event_dup_region_df = None
                if minor_baf_df is not None:
                    baf_region_df = minor_baf_df[locus_mask].sort_values("Start")
                if baf_sites_df is not None:
                    baf_sites_region_df = baf_sites_df[locus_mask].sort_values("Start")
                if event_del_df is not None:
                    event_del_region_df = event_del_df[locus_mask].sort_values("Start")
                if event_dup_df is not None:
                    event_dup_region_df = event_dup_df[locus_mask].sort_values("Start")
                if len(region_df) > PDF_MAX_SIGNAL_BINS:
                    region_df, rebinned_frames = _rebin_aligned_region_dfs_for_display(
                        locus,
                        region_df,
                        [baf_region_df, baf_sites_region_df, event_del_region_df, event_dup_region_df],
                        max_total_bins=PDF_MAX_SIGNAL_BINS,
                    )
                    (
                        baf_region_df,
                        baf_sites_region_df,
                        event_del_region_df,
                        event_dup_region_df,
                    ) = rebinned_frames

                region_start = int(region_df["Start"].min())
                region_end = int(region_df["End"].max())
                stage_start = perf_counter()
                xform = FlankCompressor(
                    region_start,
                    region_end,
                    locus.start,
                    locus.end,
                    flank_scale=flank_scale,
                )
                _print_timing(f"{cluster} {category_name} xform", stage_start)

                stage_start = perf_counter()
                cluster_calls_df = calls_df[calls_df["cluster"] == cluster]
                raw_region_df = None
                if all([
                    raw_counts_df is not None or highres_path is not None,
                    bool(raw_sample_medians),
                    lowres_median_bin_size is not None,
                ]):
                    raw_region_df = _build_raw_region_df(
                        locus,
                        region_start,
                        region_end,
                        raw_counts_df,
                        list(raw_sample_medians.keys()),
                        raw_sample_medians,
                        lowres_median_bin_size,
                        processed_region_df=region_df,
                        highres_path=highres_path,
                    )
                    if raw_region_df is not None:
                        raw_region_df = _rebin_region_df(raw_region_df, locus)
                _print_timing(f"{cluster} {category_name} raw-region prep", stage_start)

                cluster_specs = sorted(
                    pages_by_cluster[cluster],
                    key=lambda spec: (str(spec["sample"]), str(spec.get("gd_id", ""))),
                )
                print(f"  Adding {len(cluster_specs)} page(s) for {cluster}")
                for spec in cluster_specs:
                    _render_pdf_sample_page(
                        pdf,
                        str(spec["sample"]),
                        cluster,
                        locus,
                        region_df,
                        cluster_calls_df,
                        confidence_column,
                        confidence_threshold,
                        gtf,
                        segdup,
                        min_gene_label_spacing,
                        xform,
                        raw_region_df,
                        baf_region_df,
                        baf_sites_region_df,
                        event_del_region_df,
                        event_dup_region_df,
                        gaps,
                        viterbi_data,
                        title_suffix=spec.get("title_suffix"),
                    )

                _print_timing(
                    f"{cluster} {category_name} cluster total ({len(cluster_specs)} pages)",
                    cluster_start,
                )

        print(f"  Saved PDF: {pdf_path}")
        _print_timing(f"{category_name} pdf total", category_start)


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
    parser.add_argument(
        "--viterbi-paths",
        required=False,
        help="Viterbi paths file (viterbi_paths.tsv.gz) produced by "
             "gd_cnv_call.py.  When provided, the Viterbi CN segmentation "
             "is overlaid on carrier depth plots.",
    )
    parser.add_argument(
        "--event-marginals",
        required=False,
        help="Per-bin event marginal file (event_marginals.tsv.gz) produced by "
             "gd_cnv_call.py. When omitted, plot.py will look for a sibling "
             "file next to --calls.",
    )
    parser.add_argument(
        "--eval-report",
        required=False,
        help="truth_evaluation_report.tsv produced by the eval tool. When provided, "
             "plot.py writes true_positives.pdf, false_positives.pdf, and "
             "false_negatives.pdf instead of carrier_plots.pdf.",
    )
    parser.add_argument(
        "--carrier-confidence-threshold",
        type=float,
        default=0.6,
        help="Minimum confidence required for a call to appear in carrier_plots.pdf.",
    )
    return parser.parse_args()


def main():
    """Main function."""
    args = parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    setup_logging(args.output_dir)

    print(f"Output directory: {args.output_dir}")
    print(f"Flank scale: {args.flank_scale}")

    # Load data
    print("\nLoading data...")

    print(f"  Loading calls: {args.calls}")
    try:
        calls_df = pd.read_csv(args.calls, sep="\t", compression="infer")
    except pd.errors.EmptyDataError:
        print("    Calls file is empty; continuing with a no-calls DataFrame")
        calls_df = pd.DataFrame(columns=EMPTY_CALLS_COLUMNS)
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

    # Guard against duplicate (cluster, chr, start, end, sample) rows that
    # can arise when a boundary bin is shared between adjacent intervals.
    _dup_key = ["cluster", "chr", "start", "end", "sample"]
    n_before = len(cn_posteriors_df)
    cn_posteriors_df = cn_posteriors_df.drop_duplicates(subset=_dup_key)
    n_dropped = n_before - len(cn_posteriors_df)
    if n_dropped:
        print(f"    NOTE: dropped {n_dropped} duplicate bin-sample rows before pivoting")

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

    minor_baf_df = None
    baf_sites_df = None
    if "minor_baf_median" in cn_posteriors_df.columns and "baf_n_sites" in cn_posteriors_df.columns:
        minor_baf_df = cn_posteriors_df.pivot(
            index=["cluster", "chr", "start", "end"],
            columns="sample",
            values="minor_baf_median",
        ).reset_index()
        minor_baf_df = minor_baf_df.rename(columns={
            "cluster": "Cluster", "chr": "Chr", "start": "Start", "end": "End",
        })
        baf_sites_df = cn_posteriors_df.pivot(
            index=["cluster", "chr", "start", "end"],
            columns="sample",
            values="baf_n_sites",
        ).reset_index()
        baf_sites_df = baf_sites_df.rename(columns={
            "cluster": "Cluster", "chr": "Chr", "start": "Start", "end": "End",
        })
        print("    Built aligned minor-BAF matrices for carrier PDF bars")

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
        raw_sample_cols = [
            c for c in raw_counts_df.columns
            if c not in ("Chr", "Start", "End", "source_file", "Bin")
        ]
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

    ploidy_lookup = _build_ploidy_lookup(ploidy_df)

    # Load Viterbi overlay data if paths file provided
    viterbi_data: Optional[ViterbiOverlayData] = None
    if args.viterbi_paths:
        print(f"\nLoading Viterbi paths: {args.viterbi_paths}")
        paths_df = pd.read_csv(
            args.viterbi_paths, sep="\t", compression="infer",
        )
        viterbi_data = ViterbiOverlayData(paths_df)
        print(f"  {len(paths_df)} path records loaded")

    event_marginals_df: Optional[pd.DataFrame] = None
    event_del_df: Optional[pd.DataFrame] = None
    event_dup_df: Optional[pd.DataFrame] = None
    event_marginals_path = args.event_marginals
    if event_marginals_path is None:
        sibling_event_marginals = os.path.join(
            os.path.dirname(args.calls),
            "event_marginals.tsv.gz",
        )
        if os.path.exists(sibling_event_marginals):
            event_marginals_path = sibling_event_marginals

    if event_marginals_path:
        print(f"\nLoading event marginals: {event_marginals_path}")
        event_marginals_df = pd.read_csv(
            event_marginals_path,
            sep="\t",
            compression="infer",
        )
        event_marginals_df = event_marginals_df.rename(columns={
            "cluster": "Cluster",
            "chrom": "Chr",
            "start": "Start",
            "end": "End",
        })
        if "prob_del_event" in event_marginals_df.columns:
            event_del_df = event_marginals_df.pivot(
                index=["Cluster", "Chr", "Start", "End"],
                columns="sample",
                values="prob_del_event",
            ).reset_index()
        if "prob_dup_event" in event_marginals_df.columns:
            event_dup_df = event_marginals_df.pivot(
                index=["Cluster", "Chr", "Start", "End"],
                columns="sample",
                values="prob_dup_event",
            ).reset_index()
        print(f"  {len(event_marginals_df)} event-marginal records loaded")
    else:
        print("\nNo event marginals provided; event panel will show a placeholder.")

    eval_report_df: Optional[pd.DataFrame] = None
    if args.eval_report:
        print(f"\nLoading eval report: {args.eval_report}")
        eval_report_df = pd.read_csv(args.eval_report, sep="\t")
        print(f"  {len(eval_report_df)} eval report rows loaded")

    # Create summary plots
    print("\nCreating summary plots...")
    plot_carrier_summary(plot_calls_df, args.output_dir)
    plot_confidence_distribution(plot_calls_df, args.output_dir)

    # Create locus overview plots
    print("\nCreating locus overview plots...")
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

    if eval_report_df is not None:
        print("\nCreating eval-category PDFs...")
        create_eval_category_pdfs(
            eval_report_df,
            plot_calls_df,
            depth_df,
            loci_to_plot,
            gtf,
            segdup,
            args.output_dir,
            event_marginals_df=event_marginals_df,
            event_del_df=event_del_df,
            event_dup_df=event_dup_df,
            minor_baf_df=minor_baf_df,
            baf_sites_df=baf_sites_df,
            padding=args.padding,
            min_gene_label_spacing=args.min_gene_label_spacing,
            raw_counts_df=raw_counts_df,
            gaps=gaps,
            flank_scale=args.flank_scale,
            lowres_median_bin_size=lowres_median_bin_size,
            highres_path=highres_path,
            viterbi_data=viterbi_data,
            confidence_threshold=args.carrier_confidence_threshold,
        )
    else:
        print("\nCreating carrier PDF...")
        create_carrier_pdf(
            plot_calls_df, depth_df, loci_to_plot, gtf, segdup,
            args.output_dir,
            event_marginals_df=event_marginals_df,
            event_del_df=event_del_df,
            event_dup_df=event_dup_df,
            minor_baf_df=minor_baf_df,
            baf_sites_df=baf_sites_df,
            padding=args.padding,
            min_gene_label_spacing=args.min_gene_label_spacing,
            raw_counts_df=raw_counts_df,
            gaps=gaps,
            flank_scale=args.flank_scale,
            lowres_median_bin_size=lowres_median_bin_size,
            highres_path=highres_path,
            viterbi_data=viterbi_data,
            confidence_threshold=args.carrier_confidence_threshold,
        )

    print("\n" + "=" * 80)
    print("Plotting complete!")
    print("=" * 80)


if __name__ == "__main__":
    main()
