"""
High-resolution bin processing and locus bin collection.

Handles querying tabix-indexed high-resolution count files, normalising
them to match the low-resolution scale, and assembling all bins across
GD loci into a single combined DataFrame.

Also serves as the CLI entry-point for the ``preprocess`` subcommand,
which runs data loading, normalisation, quality filtering, ploidy
estimation, and bin collection *without* model training, writing the
results to disk for downstream consumption by ``infer``.
"""

import argparse
import gzip
import os
from collections import defaultdict
from typing import Dict, Iterator, List, Optional, Tuple

import numpy as np
import pandas as pd
import pysam
from intervaltree import IntervalTree

from gatk_sv_gd import _util
from gatk_sv_gd._util import get_sample_columns, setup_logging
from gatk_sv_gd.models import GDLocus, GDTable, validate_gd_table_for_preprocess
from gatk_sv_gd.depth import ExclusionMask
from gatk_sv_gd.bins import (
    LocusBinMapping,
    assign_bins_to_intervals,
    compute_bin_quality_mask,
    compute_flank_regions_from_bins,
    extract_locus_bins,
    filter_low_quality_bins,
    get_flank_filter_params,
    read_data,
    rebin_locus_intervals,
)
from gatk_sv_gd.highres import normalize_highres_bins, query_highres_bins
from gatk_sv_gd.output import build_ploidy_map, estimate_ploidy, write_locus_metadata


# ── Region parsing helpers ───────────────────────────────────────────


def _parse_region(region_str: str) -> Tuple[str, Optional[int], Optional[int]]:
    """Parse a region string like ``chr1:3000-4000`` or ``chr1``.

    Returns ``(chrom, start, end)`` where *start* and *end* are ``None``
    when only a chromosome is specified.
    """
    if ":" in region_str:
        chrom, coords = region_str.split(":", 1)
        parts = coords.replace(",", "").split("-")
        if len(parts) != 2:
            raise ValueError(
                "Invalid region format: expected chrom:start-end"
            )
        try:
            return chrom, int(parts[0]), int(parts[1])
        except ValueError:
            raise ValueError("Invalid region coordinates") from None
    return region_str, None, None


def _is_chr_x(chrom: str) -> bool:
    """Return True when *chrom* refers to chromosome X."""
    return str(chrom).strip().lower() in {"chrx", "x"}


def _is_chr_y(chrom: str) -> bool:
    """Return True when *chrom* refers to chromosome Y."""
    return str(chrom).strip().lower() in {"chry", "y"}


def _flatten_multi_args(arg_groups: List[List[str]]) -> List[str]:
    """Flatten argparse lists produced by repeated multi-value options."""
    return [value for group in arg_groups for value in group]


def _locus_overlaps_regions(
    locus: GDLocus,
    regions: List[Tuple[str, Optional[int], Optional[int]]],
) -> bool:
    """Return True if *locus* overlaps any of the parsed regions."""
    for chrom, start, end in regions:
        if locus.chrom != chrom:
            continue
        if start is None:
            return True  # whole-chromosome match
        if locus.start < end and locus.end > start:
            return True
    return False


def _filter_and_prepare_locus_bins(
    locus_df: pd.DataFrame,
    locus: GDLocus,
    flank_regions: List[Tuple[int, int, str]],
    left_bound: int,
    right_bound: int,
    max_bins_per_interval: int,
    exclusion_mask: Optional[ExclusionMask] = None,
    flank_exclusion_mask: Optional[ExclusionMask] = None,
    exclusion_threshold: float = 0.5,
    filter_params: Optional[dict] = None,
    exclusion_bypass_regions: Optional[List[Tuple[int, int]]] = None,
    min_rebin_coverage: float = 0.5,
    ploidy_map: Optional[Dict[Tuple[str, str], int]] = None,
    par_mask: Optional[ExclusionMask] = None,
) -> Tuple[pd.DataFrame, Dict[str, List[int]]]:
    """
    Shared helper: apply mask filtering, quality filtering, trimming,
    rebinning, and interval assignment to a locus DataFrame.

    This centralises the per-locus post-extraction processing so the same
    logic is used for both low-resolution and high-resolution bins.

    Args:
        locus_df: DataFrame of bins covering the locus (+ flanks).
        locus: GDLocus definition.
        flank_regions: Pre-computed flank region tuples.
        left_bound: Left edge of active region for trimming.
        right_bound: Right edge of active region for trimming.
        max_bins_per_interval: Maximum bins per interval after rebinning.
        exclusion_mask: Optional ExclusionMask for filtering.
        flank_exclusion_mask: Optional ExclusionMask applied only to
            left/right flank bins, not GD body intervals.
        exclusion_threshold: Overlap fraction threshold for masking.
        filter_params: If provided, a dict with keys ``median_min``,
            ``median_max``, ``mad_max`` to apply
            :func:`filter_low_quality_bins` to *locus_df* before any other
            processing.  Skipped when *None*.
        exclusion_bypass_regions: Genomic ranges where masking is
            skipped (bins with any true interval overlap with a bypass
            range are kept regardless of mask overlap).
        ploidy_map: Optional ``{(sample, chrom): ploidy}`` lookup used
            to adjust per-sample depths to diploid-equivalent before
            computing bin-level median/MAD statistics for quality
            filtering.
        par_mask: Optional pseudoautosomal interval mask. Body bins that
            overlap PAR are temporarily exempted from quality filtering;
            flank bins that overlap PAR are excluded from baseline regions.

    Returns:
        Tuple of (processed DataFrame, interval_bins dict).
    """
    if len(locus_df) == 0:
        return locus_df, {}

    # Trim to active region
    locus_df = locus_df[
        (locus_df["End"] > left_bound) & (locus_df["Start"] < right_bound)
    ].copy()
    print(f"    Bins in active region: {len(locus_df)}")

    if len(locus_df) == 0:
        return locus_df, {}

    # Separate body bins and flank candidates by midpoint
    bin_mids = (locus_df["Start"] + locus_df["End"]) / 2
    is_body = (bin_mids >= locus.start) & (bin_mids < locus.end)
    body_df = locus_df[is_body].copy()
    flank_df = locus_df[~is_body].copy()

    # Body bins: exclusion masking with threshold + bypass
    if exclusion_mask is not None and len(body_df) > 0:
        bypass = exclusion_bypass_regions or []
        _bs = body_df["Start"].values
        _be = body_df["End"].values
        in_bypass = np.zeros(len(body_df), dtype=bool)
        for bs, be in bypass:
            in_bypass |= (_bs < be) & (_be > bs)
        overlaps = np.zeros(len(body_df))
        non_bypass = ~in_bypass
        if non_bypass.any():
            overlaps[non_bypass] = exclusion_mask.get_overlap_fractions_batch(
                locus.chrom, _bs[non_bypass], _be[non_bypass],
            )
        body_keep = in_bypass | (overlaps < exclusion_threshold)
        n_masked = int((~body_keep).sum())
        if n_masked > 0:
            print(f"    Masked {n_masked} body bins due to exclusion overlap")
        body_df = body_df[body_keep].copy()

    # Flank bins: strict zero-overlap exclusion masking
    if exclusion_mask is not None and len(flank_df) > 0:
        flank_ov = exclusion_mask.get_overlap_fractions_batch(
            locus.chrom, flank_df["Start"].values, flank_df["End"].values,
        )
        flank_keep = flank_ov == 0.0
        n_masked = int((~flank_keep).sum())
        if n_masked > 0:
            print(f"    Masked {n_masked} flank bins due to exclusion overlap")
        flank_df = flank_df[flank_keep].copy()

    # Additional flank-only exclusion masking
    if flank_exclusion_mask is not None and len(flank_df) > 0:
        fe_ov = flank_exclusion_mask.get_overlap_fractions_batch(
            locus.chrom, flank_df["Start"].values, flank_df["End"].values,
        )
        fe_keep = fe_ov == 0.0
        n_masked = int((~fe_keep).sum())
        if n_masked > 0:
            print(f"    Masked {n_masked} flank bins due to flank-exclusion overlap")
        flank_df = flank_df[fe_keep].copy()

    if par_mask is not None and len(flank_df) > 0:
        par_ov = par_mask.get_overlap_fractions_batch(
            locus.chrom, flank_df["Start"].values, flank_df["End"].values,
        )
        par_keep = par_ov == 0.0
        n_masked = int((~par_keep).sum())
        if n_masked > 0:
            print(f"    Masked {n_masked} flank bins due to PAR overlap")
        flank_df = flank_df[par_keep].copy()

    # Optional quality filtering
    if filter_params is not None:
        flank_filter_params = get_flank_filter_params(filter_params, locus.chrom)

        def _qf(sub_df, label, quality_params):
            if len(sub_df) == 0:
                return sub_df
            keep, quality_stats = compute_bin_quality_mask(
                sub_df,
                median_min=quality_params["median_min"],
                median_max=quality_params["median_max"],
                mad_max=quality_params["mad_max"],
                ploidy_map=ploidy_map,
                par_mask=par_mask,
            )
            n_filt = quality_stats["filtered"]
            if n_filt > 0:
                print(
                    f"    Quality-filtered {n_filt}/{len(sub_df)} {label} bins "
                    f"(median [{quality_params['median_min']}, {quality_params['median_max']}], "
                    f"MAD <= {quality_params['mad_max']})"
                )
            if quality_stats["par_ignored"] > 0:
                print(f"    Ignored {quality_stats['par_ignored']}/{len(sub_df)} {label} PAR bins")
            return sub_df[keep].copy()
        body_df = _qf(body_df, "body", filter_params)
        flank_df = _qf(flank_df, "flank", flank_filter_params)

    # Combine body + flank bins
    parts = [p for p in [body_df, flank_df] if len(p) > 0]
    if parts:
        locus_df = pd.concat(parts).sort_values("Start").reset_index(drop=True)
    else:
        locus_df = pd.DataFrame(columns=locus_df.columns)
    print(f"    Combined bins (body + flanks): {len(locus_df)}")

    if len(locus_df) == 0:
        return locus_df, {}

    # Rebin
    if max_bins_per_interval > 0 and len(locus_df) > 0:
        n_before = len(locus_df)
        locus_df = rebin_locus_intervals(locus_df, locus, max_bins_per_interval, flank_regions,
                                         min_rebin_coverage=min_rebin_coverage)
        if len(locus_df) < n_before:
            print(f"    Bins after rebinning: {len(locus_df)} (reduced from {n_before})")

    # Assign to intervals
    interval_bins = assign_bins_to_intervals(locus_df, locus, flank_regions) if len(locus_df) > 0 else {}

    # Drop breakpoint ranges
    bp_range_indices = interval_bins.pop("breakpoint_ranges", [])
    if bp_range_indices:
        print(f"    Masking {len(bp_range_indices)} breakpoint-range bin(s)")
        locus_df = locus_df.drop(index=bp_range_indices)

    return locus_df, interval_bins


def _select_highres_interval_replacements(
    undercovered: List[Tuple[str, int]],
    hr_interval_bins: Dict[str, List[int]],
    min_bins_per_interval: Optional[int] = None,
) -> List[Tuple[str, int, int]]:
    """Return body intervals whose high-res replacement improves bin count.

    High-resolution rescue should be decided per under-covered body interval,
    not by replacing the entire locus wholesale. Only intervals whose
    processed high-res result yields more modeled bins than the current
    low-resolution result are selected for replacement. When
    *min_bins_per_interval* is provided, the high-res result must also meet
    that usual hard-failure threshold.
    """
    replacements: List[Tuple[str, int, int]] = []
    for name, orig_count in undercovered:
        hr_count = len(hr_interval_bins.get(name, []))
        meets_minimum = min_bins_per_interval is None or hr_count >= min_bins_per_interval
        if hr_count > orig_count and meets_minimum:
            replacements.append((name, orig_count, hr_count))
    return replacements


def _merge_highres_interval_replacements(
    locus_df: pd.DataFrame,
    hr_locus_df: pd.DataFrame,
    interval_bins: Dict[str, List[int]],
    hr_interval_bins: Dict[str, List[int]],
    replacement_intervals: List[str],
    locus: GDLocus,
    flank_regions: List[Tuple[int, int, str]],
) -> Tuple[pd.DataFrame, Dict[str, List[int]]]:
    """Replace selected body intervals with high-res bins and reassign bins.

    The replacement is intentionally interval-scoped: untouched body intervals
    and flanks remain on the existing low-resolution representation.
    """
    if not replacement_intervals:
        return locus_df, interval_bins

    orig_drop_indices = sorted({
        idx
        for name in replacement_intervals
        for idx in interval_bins.get(name, [])
    })
    hr_add_indices = sorted({
        idx
        for name in replacement_intervals
        for idx in hr_interval_bins.get(name, [])
    })

    remaining_df = locus_df.drop(index=orig_drop_indices) if orig_drop_indices else locus_df
    hr_replacement_df = hr_locus_df.loc[hr_add_indices] if hr_add_indices else hr_locus_df.iloc[0:0].copy()

    merged_df = (
        pd.concat([remaining_df, hr_replacement_df], axis=0)
        .sort_values("Start")
        .reset_index(drop=True)
    )

    merged_interval_bins = assign_bins_to_intervals(merged_df, locus, flank_regions)
    bp_range_indices = merged_interval_bins.pop("breakpoint_ranges", [])
    if bp_range_indices:
        print(f"    [highres] masked {len(bp_range_indices)} breakpoint-range bin(s) after merge")
        merged_df = merged_df.drop(index=bp_range_indices).reset_index(drop=True)
        merged_interval_bins = assign_bins_to_intervals(merged_df, locus, flank_regions)
        merged_interval_bins.pop("breakpoint_ranges", None)

    return merged_df, merged_interval_bins


def collect_all_locus_bins(
    df: pd.DataFrame,
    gd_table: GDTable,
    exclusion_mask: Optional[ExclusionMask],
    flank_exclusion_mask: Optional[ExclusionMask] = None,
    par_mask: Optional[ExclusionMask] = None,
    exclusion_threshold: float = 0.5,
    locus_padding: int = 0,
    min_bins_per_interval: int = 3,
    max_bins_per_interval: int = 10,
    highres_counts_path: Optional[str] = None,
    column_medians: Optional[np.ndarray] = None,
    lowres_median_bin_size: Optional[float] = None,
    filter_params: Optional[dict] = None,
    exclusion_bypass_threshold: float = 0.9,
    min_rebin_coverage: float = 0.5,
    min_flank_bases: int = 50000,
    max_flank_bases: int = 1_000_000,
    min_flank_bins: int = 10,
    min_flank_coverage: float = 0.1,
    regions: Optional[List[Tuple[str, Optional[int], Optional[int]]]] = None,
    ploidy_map: Optional[Dict[Tuple[str, str], int]] = None,
) -> Tuple[pd.DataFrame, List[LocusBinMapping], Dict[str, GDLocus]]:
    """
    Collect all bins across all GD loci into a single DataFrame.

    Non-NAHR loci are automatically skipped (only NAHR loci are
    processed).  Each body interval (region between adjacent breakpoints)
    must contain at least *min_bins_per_interval* bins after all
    processing (hard failure).

    When a body interval is ≥ *exclusion_bypass_threshold* overlapped by
    the exclusion mask, masking is skipped
    for bins in that interval so that coverage is preserved.

    When a high-resolution counts file is provided (``highres_counts_path``),
    loci with any under-covered interval are automatically re-extracted
    from the high-res file.

    Args:
        df: DataFrame with all (low-resolution, normalised) bins.
        gd_table: GDTable with locus definitions.
        exclusion_mask: Optional ExclusionMask for filtering.
        flank_exclusion_mask: Optional ExclusionMask applied only to
            left/right flanks, not GD body intervals.
        par_mask: Optional pseudoautosomal interval mask. PAR bins are
            excluded from flanking regions because they do not provide a
            valid non-PAR chrX/chrY baseline under contig-wide ploidy
            handling. Body bins overlapping PAR remain temporarily exempted
            from quality filtering.
        exclusion_threshold: Minimum overlap fraction with the mask to
            exclude a bin.
        locus_padding: Padding around locus boundaries.
        min_bins_per_interval: Hard-failure minimum: each body interval
            must have at least this many bins after all processing.
        max_bins_per_interval: Maximum bins per interval after rebinning
            (0 = no rebinning).
        highres_counts_path: Optional path to a bgzipped, tabix-indexed
            high-resolution read-count file.  When set, loci with any
            under-covered interval are re-queried at this resolution. If a
            high-res body interval remains under-covered after normal
            rebinning, preprocessing retries that high-res locus with 10x
            more bins per interval and fails if the interval still misses
            *min_bins_per_interval*.
        column_medians: Per-sample autosomal median *raw* counts from the
            low-res file.  Required when *highres_counts_path* is set.
        lowres_median_bin_size: Median bin size (bp) in the low-res file.
            Required when *highres_counts_path* is set.
        filter_params: Quality-filter thresholds dict (``median_min``,
            ``median_max``, ``mad_max``) applied to high-res bins.
        exclusion_bypass_threshold: If a body interval is ≥ this fraction
            overlapped by the mask, masking is skipped for bins in
            that interval (default 0.9 = 90%).
        min_flank_bases: Minimum cumulative base pairs each flank must
            cover regardless of locus size (default 50 000).
        max_flank_bases: Maximum cumulative bp target per flank.  Caps
            locus-size so that very large loci do not force flanks to
            consume all available clean bins (default 1 000 000).
        min_flank_bins: Minimum number of bins each flank must contain
            regardless of base-pair thresholds (default 10).
        min_flank_coverage: Minimum fraction of the effective bp target
            that a flank's accumulated bin coverage must reach.  Flanks
            below this threshold are rejected (default 0.1 = 10%).
        regions: Optional list of parsed region tuples ``(chrom, start,
            end)`` to restrict processing to.  Only loci overlapping at
            least one region are included.  ``None`` means process all.
        ploidy_map: Optional ``{(sample, chrom): ploidy}`` lookup used
            to adjust per-sample depths to diploid-equivalent before
            computing bin-level median/MAD statistics for quality
            filtering on sex chromosomes.

    Returns:
        Tuple of:
        - Combined DataFrame with all locus bins
        - List of LocusBinMapping objects tracking bin assignments
        - Dict of included loci (cluster -> GDLocus)
    """
    all_locus_dfs = []
    all_mappings = []
    included_loci = {}
    current_idx = 0

    sample_cols = get_sample_columns(df)

    print(f"\n{'=' * 80}")
    print("COLLECTING BINS ACROSS ALL GD LOCI")
    print(f"{'=' * 80}")

    # Pre-build a per-chromosome cache of filtered (masked + quality) bins.
    # Flank computation walks outward through these; using the full chromosome
    # guarantees filtered bins are found even when multi-megabase masked deserts
    # surround the locus.
    #
    # Quality filtering (median/MAD) is applied here so that the flank
    # extension only counts high-quality bins toward coverage targets.
    # This prevents the flank from terminating prematurely on a cluster
    # of low-quality bins that will later be removed.
    print("\nBuilding per-chromosome filtered bin cache...")
    # Build quality-filter params dict — always, so flank extension can
    # exclude low-quality bins during accumulation.
    _flank_filter_params: Optional[dict] = None
    if filter_params is not None:
        _flank_filter_params = filter_params
    else:
        # When no explicit filter_params were provided (no high-res path),
        # build one from the same thresholds used by the genome-wide filter
        # so that flank extension is quality-aware even for low-res data.
        # The caller can skip this by setting all thresholds to extremes.
        pass  # will be passed as None; quality filtering already applied to df

    chrom_filtered: Dict[str, pd.DataFrame] = {}
    cache_total_bins = 0
    cache_input_bins = 0
    cache_excluded_bins = 0
    cache_par_excluded_bins = 0
    cache_quality_filtered_bins = 0
    for chrom, chrom_df in df.groupby("Chr"):
        keep = np.ones(len(chrom_df), dtype=bool)
        n_excluded = 0
        n_quality = 0
        n_par_excluded = 0
        par_overlap = np.zeros(len(chrom_df), dtype=bool)
        chrom_flank_filter_params = get_flank_filter_params(_flank_filter_params, chrom)
        if exclusion_mask is not None:
            overlaps = exclusion_mask.get_overlap_fractions_batch(
                chrom, chrom_df["Start"].values, chrom_df["End"].values
            )
            # Strict zero-overlap filtering: the chrom_filtered cache is
            # used exclusively for flank boundary computation.  Flanks
            # establish baseline depth and must be free of any exclusion-
            # region contamination, so only truly clean bins are kept.
            exclusion_keep = overlaps == 0.0
            n_excluded = int((~exclusion_keep).sum())
            keep &= exclusion_keep
        if par_mask is not None:
            par_overlap = par_mask.get_overlap_fractions_batch(
                chrom, chrom_df["Start"].values, chrom_df["End"].values,
            ) > 0.0
            n_par_excluded = int(par_overlap.sum())
            keep &= ~par_overlap
        if chrom_flank_filter_params is not None:
            quality_keep = np.ones(len(chrom_df), dtype=bool)
            quality_positions = np.flatnonzero(~par_overlap)
            if len(quality_positions) > 0:
                subset_keep, quality_stats = compute_bin_quality_mask(
                    chrom_df.iloc[quality_positions],
                    median_min=chrom_flank_filter_params["median_min"],
                    median_max=chrom_flank_filter_params["median_max"],
                    mad_max=chrom_flank_filter_params["mad_max"],
                    ploidy_map=ploidy_map,
                    par_mask=None,
                )
                quality_keep[quality_positions] = subset_keep
                n_quality = quality_stats["filtered"]
            keep &= quality_keep
        chrom_filtered[chrom] = chrom_df[keep].copy()
        cache_total_bins += len(chrom_filtered[chrom])
        cache_input_bins += len(chrom_df)
        cache_excluded_bins += n_excluded
        cache_par_excluded_bins += n_par_excluded
        cache_quality_filtered_bins += n_quality

    print(
        "  Filtered cache built: "
        f"{cache_total_bins:,}/{cache_input_bins:,} bins retained across "
        f"{len(chrom_filtered)} contigs; exclusion-masked={cache_excluded_bins:,}, "
        f"PAR-excluded={cache_par_excluded_bins:,}, "
        f"quality-filtered={cache_quality_filtered_bins:,}"
    )

    # ------------------------------------------------------------------
    # Helper: create LocusBinMapping entries for a processed locus_df
    # ------------------------------------------------------------------
    def _build_mappings(
        locus_df: pd.DataFrame,
        locus: GDLocus,
        cluster: str,
        flank_regions: List[Tuple[int, int, str]],
        start_idx: int,
    ) -> Tuple[List[LocusBinMapping], int]:
        """Return (mappings_list, next_idx).

        Breakpoint-range bins should already have been removed from
        *locus_df* before this function is called.  If a bin still
        doesn't match any named region it is skipped with a warning.
        """
        all_named_regions = locus.get_intervals() + flank_regions
        mappings: List[LocusBinMapping] = []
        idx_counter = start_idx
        for idx in locus_df.index:
            bin_row = locus_df.loc[idx]
            bin_start, bin_end = bin_row["Start"], bin_row["End"]

            # Assign to the region with the greatest base-pair overlap
            assigned_interval = None
            best_overlap = 0
            for start, end, name in all_named_regions:
                ov = min(bin_end, end) - max(bin_start, start)
                if ov > best_overlap:
                    best_overlap = ov
                    assigned_interval = name

            if assigned_interval is None:
                # Breakpoint-range bins should have been masked upstream;
                # if one slips through, warn and skip rather than silently
                # including unreliable signal.
                print("  WARNING: skipped one unassigned bin outside modeled intervals")
                continue

            mappings.append(LocusBinMapping(
                cluster=cluster,
                locus=locus,
                interval_name=assigned_interval,
                array_idx=idx_counter,
                chrom=bin_row["Chr"],
                start=int(bin_row["Start"]),
                end=int(bin_row["End"]),
            ))
            idx_counter += 1
        return mappings, idx_counter

    # ------------------------------------------------------------------
    # Per-locus processing
    # ------------------------------------------------------------------
    processed_nahr_loci = 0
    skipped_non_nahr_loci = 0
    loci_without_flanks = 0
    exclusion_bypassed_intervals = 0
    highres_rescue_attempts = 0
    for cluster, locus in gd_table.get_all_loci().items():
        # Skip loci outside requested regions
        if regions is not None and not _locus_overlaps_regions(locus, regions):
            continue

        # Only NAHR loci are analyzed
        if not locus.is_nahr:
            skipped_non_nahr_loci += 1
            continue

        processed_nahr_loci += 1
        if _util.VERBOSE:
            print(
                "\nProcessing one NAHR locus "
                f"({len(locus.gd_entries)} GD entries, {locus.n_breakpoints} breakpoints)"
            )

        locus_size = locus.end - locus.start

        # Compute flank coordinates from the full chromosome's filtered bins so
        # that even multi-megabase exclusion deserts don't prevent flank discovery.
        chrom_bins = chrom_filtered.get(locus.chrom, pd.DataFrame())
        chrom_flank_filter_params = get_flank_filter_params(_flank_filter_params, locus.chrom)
        flank_regions = compute_flank_regions_from_bins(
            chrom_bins, locus, locus_size,
            min_flank_bases=min_flank_bases,
            max_flank_bases=max_flank_bases,
            min_flank_bins=min_flank_bins,
            min_flank_coverage=min_flank_coverage,
            filter_params=chrom_flank_filter_params,
            ploidy_map=ploidy_map,
            par_mask=par_mask,
        )
        if flank_regions:
            if _util.VERBOSE:
                print(f"  Flank regions selected: {len(flank_regions)}")
        else:
            loci_without_flanks += 1
            if _util.VERBOSE:
                print("  Warning: no flanking bins found for one locus")

        # Determine extraction bounds: locus body + computed flank extents.
        left_bound = locus.start
        right_bound = locus.end
        for fs, fe, fn in flank_regions:
            if fn == "left_flank":
                left_bound = fs
            elif fn == "right_flank":
                right_bound = fe

        # Compute which body intervals are so heavily overlapped by exclusion
        # regions that masking should be bypassed for bins in those regions.
        exclusion_bypass_regions: List[Tuple[int, int]] = []
        if exclusion_mask is not None:
            for iv_start, iv_end, iv_name in locus.get_intervals():
                iv_overlap = exclusion_mask.get_overlap_fraction(
                    locus.chrom, iv_start, iv_end,
                )
                if iv_overlap >= exclusion_bypass_threshold:
                    exclusion_bypassed_intervals += 1
                    exclusion_bypass_regions.append((iv_start, iv_end))

        # ------------------------------------------------------------------
        # Step 2 (body): extract body bins with exclusion threshold + bypass
        # ------------------------------------------------------------------
        body_df = extract_locus_bins(
            df, locus, exclusion_mask,
            exclusion_threshold=exclusion_threshold,
            padding=locus_padding,
            exclusion_bypass_regions=exclusion_bypass_regions,
        )
        if _util.VERBOSE:
            print(f"  Body bins after filtering: {len(body_df)}")

        # ------------------------------------------------------------------
        # Step 2 (flanks): collect pre-filtered flank bins from cache
        # ------------------------------------------------------------------
        flank_dfs = []
        if len(chrom_bins) > 0:
            _cb_mids = (chrom_bins["Start"] + chrom_bins["End"]) / 2
            for fs, fe, fn in flank_regions:
                side_mask = (chrom_bins["End"] > fs) & (chrom_bins["Start"] < fe)
                if fn == "left_flank":
                    side_mask &= _cb_mids < locus.start
                else:
                    side_mask &= _cb_mids >= locus.end
                matched = chrom_bins[side_mask]
                if flank_exclusion_mask is not None and len(matched) > 0:
                    _fe_ov = flank_exclusion_mask.get_overlap_fractions_batch(
                        locus.chrom, matched["Start"].values, matched["End"].values,
                    )
                    matched = matched[_fe_ov == 0.0]
                if len(matched) > 0:
                    flank_dfs.append(matched)
        flank_df = pd.concat(flank_dfs) if flank_dfs else pd.DataFrame(columns=df.columns)

        # Combine body + flank bins
        parts = [p for p in [body_df, flank_df] if len(p) > 0]
        if parts:
            locus_df = pd.concat(parts).drop_duplicates(
                subset=["Chr", "Start", "End"],
            ).sort_values("Start").reset_index(drop=True)
        else:
            locus_df = pd.DataFrame(columns=df.columns)
        if _util.VERBOSE:
            print(f"  Combined bins (body + flanks): {len(locus_df)}")

        # Rebin to reduce number of bins per interval/flank if requested
        if max_bins_per_interval > 0:
            locus_df_orig = locus_df
            if _util.VERBOSE:
                _sc = get_sample_columns(locus_df_orig)
                _meds = np.median(locus_df_orig[_sc].values, axis=1)
                print(
                    "    [verbose] pre-rebin median depth summary: "
                    f"bins={len(locus_df_orig)}, min={np.min(_meds):.3f}, "
                    f"median={np.median(_meds):.3f}, max={np.max(_meds):.3f}"
                )

            locus_df = rebin_locus_intervals(locus_df, locus, max_bins_per_interval, flank_regions,
                                             min_rebin_coverage=min_rebin_coverage)
            if len(locus_df) < len(locus_df_orig) and _util.VERBOSE:
                print(f"  Bins after rebinning: {len(locus_df)} (reduced from {len(locus_df_orig)})")

            if _util.VERBOSE:
                _sc = get_sample_columns(locus_df)
                _meds = np.median(locus_df[_sc].values, axis=1)
                print(
                    "    [verbose] post-rebin median depth summary: "
                    f"bins={len(locus_df)}, min={np.min(_meds):.3f}, "
                    f"median={np.median(_meds):.3f}, max={np.max(_meds):.3f}"
                )

        # Assign bins to intervals and flanking regions
        interval_bins = assign_bins_to_intervals(locus_df, locus, flank_regions)
        total_assigned = sum(len(v) for v in interval_bins.values())
        if _util.VERBOSE:
            print(f"    total: {total_assigned} bins")

        # Mask breakpoint-range bins: bins inside the breakpoint SD-block
        # ranges (not in any body interval or flank) carry unreliable depth
        # signal and must be excluded from model training.
        bp_range_indices = interval_bins.pop("breakpoint_ranges", [])
        if bp_range_indices:
            if _util.VERBOSE:
                print(f"    Masking {len(bp_range_indices)} breakpoint-range bin(s)")
            locus_df = locus_df.drop(index=bp_range_indices)

        # ------------------------------------------------------------------
        # Per-interval bin count check + optional high-res replacement
        # ------------------------------------------------------------------
        body_intervals = locus.get_intervals()  # [(start, end, name), ...]
        body_interval_names = {name for _, _, name in body_intervals}

        effective_min_bins = min_bins_per_interval

        # Identify under-covered body intervals
        def _undercovered_intervals(
            ivbins: Dict[str, List[int]],
            threshold: int,
        ) -> List[Tuple[str, int]]:
            """Return list of (name, count) for body intervals below *threshold*."""
            return [
                (name, len(ivbins.get(name, [])))
                for name in body_interval_names
                if len(ivbins.get(name, [])) < threshold
            ]

        undercovered = _undercovered_intervals(interval_bins, effective_min_bins)

        used_highres = False
        # Body intervals where even the raw high-res data doesn't have
        # enough bins — the resolution is simply too coarse for this locus
        # size, so we exempt them from the hard-failure threshold later.
        hr_physically_limited: set = set()
        hr_fallback_failed: set = set()
        if (
            undercovered and
            highres_counts_path is not None and
            column_medians is not None and
            lowres_median_bin_size is not None
        ):
            highres_rescue_attempts += 1
            print(f"\n  *** {len(undercovered)} body interval(s) have fewer "
                f"than {effective_min_bins} bins; switching to high-res ***")

            # Query the tabix-indexed high-res file for the active region
            hr_raw_df = query_highres_bins(
                highres_counts_path,
                locus.chrom,
                left_bound,
                right_bound,
                sample_cols,
            )
            print(f"    [highres] queried {len(hr_raw_df)} raw bins for active locus region")

            if len(hr_raw_df) > 0:
                # Before any filtering, count how many raw hi-res bins
                # fall in each body interval.  If a body interval cannot
                # reach min_bins_per_interval even with every raw bin, it
                # is physically limited by the bin resolution and should
                # be exempted from the hard-failure threshold (provided it
                # has > 0 bins after processing).
                hr_raw_mids = (hr_raw_df["Start"] + hr_raw_df["End"]) / 2
                for iv_start, iv_end, iv_name in body_intervals:
                    raw_count = int(((hr_raw_mids >= iv_start) & (hr_raw_mids < iv_end)).sum())
                    if raw_count < min_bins_per_interval:
                        hr_physically_limited.add(iv_name)
                        print(
                            "    [highres] one body interval is resolution-limited; "
                            "will relax threshold"
                        )

                # Normalise using the low-res calibration, scaled by bin size ratio
                hr_norm_df = normalize_highres_bins(
                    hr_raw_df, sample_cols, column_medians, lowres_median_bin_size,
                )

                # Apply the full filter-trim-rebin-assign pipeline
                hr_locus_df, hr_interval_bins = _filter_and_prepare_locus_bins(
                    hr_norm_df,
                    locus,
                    flank_regions,
                    left_bound,
                    right_bound,
                    max_bins_per_interval,
                    exclusion_mask=exclusion_mask,
                    flank_exclusion_mask=flank_exclusion_mask,
                    exclusion_threshold=exclusion_threshold,
                    filter_params=filter_params,
                    exclusion_bypass_regions=exclusion_bypass_regions,
                    min_rebin_coverage=min_rebin_coverage,
                    ploidy_map=ploidy_map,
                    par_mask=par_mask,
                )

                interval_replacements = _select_highres_interval_replacements(
                    undercovered,
                    hr_interval_bins,
                    min_bins_per_interval=min_bins_per_interval,
                )
                normal_replacement_names = {name for name, _, _ in interval_replacements}

                hr_still_undercovered = [
                    (name, len(hr_interval_bins.get(name, [])))
                    for name, _ in undercovered
                    if name not in normal_replacement_names and
                    len(hr_interval_bins.get(name, [])) < min_bins_per_interval
                ]
                fallback_replacements: List[Tuple[str, int, int]] = []
                fallback_locus_df: Optional[pd.DataFrame] = None
                fallback_interval_bins: Dict[str, List[int]] = {}

                if hr_still_undercovered:
                    if max_bins_per_interval > 0:
                        fallback_max_bins_per_interval = max_bins_per_interval * 10
                        print(
                            f"    [highres] {len(hr_still_undercovered)} interval(s) "
                            f"remain under-covered after high-res; retrying rebin with "
                            f"--max-bins-per-interval={fallback_max_bins_per_interval}"
                        )
                        fallback_locus_df, fallback_interval_bins = _filter_and_prepare_locus_bins(
                            hr_norm_df,
                            locus,
                            flank_regions,
                            left_bound,
                            right_bound,
                            fallback_max_bins_per_interval,
                            exclusion_mask=exclusion_mask,
                            flank_exclusion_mask=flank_exclusion_mask,
                            exclusion_threshold=exclusion_threshold,
                            filter_params=filter_params,
                            exclusion_bypass_regions=exclusion_bypass_regions,
                            min_rebin_coverage=min_rebin_coverage,
                            ploidy_map=ploidy_map,
                            par_mask=par_mask,
                        )

                        retry_undercovered = [
                            (name, orig_count)
                            for name, orig_count in undercovered
                            if any(name == retry_name for retry_name, _ in hr_still_undercovered)
                        ]
                        fallback_replacements = _select_highres_interval_replacements(
                            retry_undercovered,
                            fallback_interval_bins,
                            min_bins_per_interval=min_bins_per_interval,
                        )
                        fallback_replacement_names = {name for name, _, _ in fallback_replacements}

                        for name, normal_count in hr_still_undercovered:
                            fallback_count = len(fallback_interval_bins.get(name, []))
                            if name in fallback_replacement_names:
                                print(
                                    "    [highres] one fallback interval reached "
                                    f"{fallback_count} bins; replacing low-res interval"
                                )
                            else:
                                hr_fallback_failed.add(name)
                                print(
                                    "    [highres] one fallback interval still has "
                                    f"{fallback_count} bins (normal high-res had {normal_count}; "
                                    f"need >= {min_bins_per_interval})"
                                )
                    else:
                        for name, cnt in hr_still_undercovered:
                            hr_fallback_failed.add(name)
                            print(
                                "    [highres] one interval remains under-covered "
                                f"after high-res ({cnt} bins), and fallback rebin is not "
                                "available because --max-bins-per-interval=0"
                            )

                if interval_replacements:
                    replacement_names = [name for name, _, _ in interval_replacements]
                    print(
                        f"    [highres] {len(replacement_names)} interval(s) improved; "
                        "replacing low-res intervals"
                    )

                    locus_df, interval_bins = _merge_highres_interval_replacements(
                        locus_df,
                        hr_locus_df,
                        interval_bins,
                        hr_interval_bins,
                        replacement_names,
                        locus,
                        flank_regions,
                    )
                    used_highres = True
                    total_assigned = sum(len(v) for v in interval_bins.values())
                    print("    [highres] ACCEPTED — replaced low-res bins")
                    print(f"      total: {total_assigned} bins")
                if fallback_replacements:
                    fallback_names = [name for name, _, _ in fallback_replacements]
                    if fallback_locus_df is None:
                        raise RuntimeError("High-res fallback replacements exist without fallback bins")
                    locus_df, interval_bins = _merge_highres_interval_replacements(
                        locus_df,
                        fallback_locus_df,
                        interval_bins,
                        fallback_interval_bins,
                        fallback_names,
                        locus,
                        flank_regions,
                    )
                    used_highres = True
                    total_assigned = sum(len(v) for v in interval_bins.values())
                    print("    [highres] FALLBACK ACCEPTED — replaced low-res bins")
                    print(f"      total: {total_assigned} bins")
                else:
                    unreplaced = [
                        (name, orig_count)
                        for name, orig_count in undercovered
                        if name not in normal_replacement_names and name not in hr_fallback_failed
                    ]
                    for name, orig_count in unreplaced:
                        hr_count = len(hr_interval_bins.get(name, []))
                        print(
                            "    [highres] one interval did not improve to criteria "
                            f"({orig_count} -> {hr_count} bins; need >= "
                            f"{min_bins_per_interval}); keeping low-res interval"
                        )
                    if not interval_replacements and not fallback_replacements:
                        print(
                            "    [highres] REJECTED — high-res did not improve any "
                            "under-covered intervals to criteria"
                        )
            else:
                print("    [highres] no bins returned from tabix query")

        # Hard-failure check: each body interval must have at least
        # min_bins_per_interval bins after all processing.
        #
        # Intervals that are "physically limited" (the raw high-res data
        # itself didn't have enough bins) are exempted as long as they
        # have > 0 bins after processing.
        hard_check = _undercovered_intervals(interval_bins, min_bins_per_interval)
        if hard_check:
            hard_failures = [
                (name, cnt) for name, cnt in hard_check
                if name in hr_fallback_failed or name not in hr_physically_limited or cnt == 0
            ]
            soft_warnings = [
                (name, cnt) for name, cnt in hard_check
                if name not in hr_fallback_failed and name in hr_physically_limited and cnt > 0
            ]

            if soft_warnings:
                for name, cnt in soft_warnings:
                    print(f"  NOTE: one body interval has {cnt} bins "
                          f"(< {min_bins_per_interval}), but this is the maximum "
                          f"the bin resolution can provide — proceeding anyway.")

            if hard_failures:
                details = "\n".join(
                    f"    interval: {cnt} bins (need >= {min_bins_per_interval})"
                    for name, cnt in hard_failures
                )
                raise ValueError(
                    f"\n  WARNING: One locus has {len(hard_failures)} body interval(s) with "
                    f"fewer than --min-bins-per-interval={min_bins_per_interval} bins "
                    f"after all processing"
                    f"{' (including high-res fallback)' if hr_fallback_failed else (' (including high-res replacement)' if used_highres else '')}:\n"
                    f"{details}\n"
                    f"  CNV calls spanning these intervals will have reduced "
                    f"statistical power."
                )

        # ------------------------------------------------------------------
        # Create mappings
        # ------------------------------------------------------------------
        new_mappings, current_idx = _build_mappings(
            locus_df, locus, cluster, flank_regions, current_idx,
        )
        all_mappings.extend(new_mappings)

        all_locus_dfs.append(locus_df)
        included_loci[cluster] = locus

    if len(all_locus_dfs) == 0:
        print("\nNo loci with sufficient bins found!")
        return pd.DataFrame(), [], {}

    # Combine all DataFrames
    combined_df = pd.concat(all_locus_dfs, axis=0)

    # Reset index to ensure contiguous integer indices
    combined_df = combined_df.reset_index(drop=True)

    print(f"\n{'=' * 80}")
    print(f"TOTAL: {len(combined_df)} bins across {len(included_loci)} loci")
    print(
        "Locus preprocessing summary: "
        f"processed_NAHR={processed_nahr_loci}, skipped_non_NAHR={skipped_non_nahr_loci}, "
        f"without_flanks={loci_without_flanks}, "
        f"exclusion_bypassed_intervals={exclusion_bypassed_intervals}, "
        f"highres_rescue_attempts={highres_rescue_attempts}"
    )
    print(f"{'=' * 80}\n")

    return combined_df, all_mappings, included_loci


def get_locus_interval_bins(
    mappings: List[LocusBinMapping],
    cluster: str,
) -> Dict[str, List[int]]:
    """
    Get array indices for each interval in a specific locus.

    Args:
        mappings: List of LocusBinMapping objects
        cluster: Cluster name to filter by

    Returns:
        Dict mapping interval name to list of array indices
    """
    interval_bins = {}
    for mapping in mappings:
        if mapping.cluster != cluster:
            continue
        if mapping.interval_name not in interval_bins:
            interval_bins[mapping.interval_name] = []
        interval_bins[mapping.interval_name].append(mapping.array_idx)
    return interval_bins


# =============================================================================
# Preprocessed-data I/O
# =============================================================================


def write_preprocessed_bins(
    combined_df: pd.DataFrame,
    output_dir: str,
) -> str:
    """Write the preprocessed combined bin DataFrame to disk.

    Returns:
        Path to the written file.
    """
    output_path = os.path.join(output_dir, "preprocessed_bins.tsv.gz")
    combined_df.to_csv(output_path, sep="\t", index=False, compression="gzip")
    print("  Saved preprocessed bins table")
    print(f"  Rows: {len(combined_df):,}")
    return output_path


def _build_roi_intervals_from_mappings(
    mappings: List[LocusBinMapping],
) -> Dict[str, List[Tuple[int, int]]]:
    """Build merged genomic ROI intervals from retained locus-bin mappings.

    The returned intervals reflect the exact genomic footprint of the bins
    that survived preprocessing, including flanks and any exclusion-driven
    gaps. These are used to crop SNP/BAF rows to the same modeled regions.
    """
    intervals_by_chrom: Dict[str, List[Tuple[int, int]]] = {}
    for mapping in mappings:
        intervals_by_chrom.setdefault(mapping.chrom, []).append(
            (int(mapping.start), int(mapping.end))
        )

    merged_by_chrom: Dict[str, List[Tuple[int, int]]] = {}
    for chrom, intervals in intervals_by_chrom.items():
        if not intervals:
            merged_by_chrom[chrom] = []
            continue
        intervals = sorted(intervals)
        merged: List[Tuple[int, int]] = []
        cur_start, cur_end = intervals[0]
        for start, end in intervals[1:]:
            if start <= cur_end:
                cur_end = max(cur_end, end)
            else:
                merged.append((cur_start, cur_end))
                cur_start, cur_end = start, end
        merged.append((cur_start, cur_end))
        merged_by_chrom[chrom] = merged

    return merged_by_chrom


def _detect_baf_columns(filepath: str) -> Tuple[Optional[int], List[str]]:
    """Detect whether a BAF file has a header and return its column names.

    Expected logical columns are chromosome, position, BAF value, and sample.
    The input produced by SiteDepthtoBAF is typically headerless, so we also
    support that layout directly.
    """
    opener = gzip.open if filepath.endswith(".gz") else open
    with opener(filepath, "rt") as handle:
        first_line = handle.readline().strip()

    if not first_line:
        return None, ["Chr", "Pos", "BAF", "Sample"]

    parts = first_line.split("\t")
    if len(parts) < 4:
        raise ValueError("BAF table must have at least 4 tab-delimited columns")

    try:
        int(parts[1])
        float(parts[2])
    except ValueError:
        return 0, parts[:4]

    return None, ["Chr", "Pos", "BAF", "Sample"]


def _iter_roi_baf_records(
    baf_path: str,
    roi_intervals: Dict[str, List[Tuple[int, int]]],
    header: Optional[int],
    column_names: List[str],
) -> Iterator[Tuple[str, int, str, float, str, str]]:
    """Yield validated ROI BAF records from a tabix-indexed file."""
    tbx = pysam.TabixFile(baf_path)
    try:
        for chrom, intervals in sorted(roi_intervals.items()):
            for start, end in intervals:
                try:
                    records = tbx.fetch(chrom, max(0, start), end)
                except ValueError:
                    continue

                for line in records:
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) < 4:
                        continue

                    if header == 0 and fields[:4] == column_names[:4]:
                        continue

                    record = dict(zip(column_names, fields[:4]))
                    try:
                        pos = int(record["Pos"])
                        baf = float(record["BAF"])
                    except ValueError:
                        continue
                    if not (0.0 <= baf <= 1.0):
                        continue

                    yield (
                        str(record["Chr"]),
                        pos,
                        str(record["Pos"]),
                        baf,
                        str(record["BAF"]),
                        str(record["Sample"]),
                    )
    finally:
        tbx.close()


def write_preprocessed_baf(
    baf_path: str,
    mappings: List[LocusBinMapping],
    output_dir: str,
) -> str:
    """Filter a genome-wide BAF table down to preprocessed GD regions.

    Rows are retained when the SNP position falls within the genomic
    footprint of the retained preprocessed bins.

    The input BAF file must be bgzipped and tabix-indexed. Records are
    fetched via genomic intervals using :class:`pysam.TabixFile`, so the
    whole file is never loaded or scanned sequentially.
    """
    output_path = os.path.join(output_dir, "preprocessed_baf.tsv.gz")
    summary_path = os.path.join(output_dir, "preprocessed_baf_summary.tsv.gz")
    header, column_names = _detect_baf_columns(baf_path)
    roi_intervals = _build_roi_intervals_from_mappings(mappings)
    written_rows = 0
    sample_baf_stats: Dict[str, Tuple[int, float, float]] = {}

    interval_trees: Dict[str, IntervalTree] = {}
    mapping_by_idx = {int(m.array_idx): m for m in mappings}
    for mapping in mappings:
        tree = interval_trees.setdefault(mapping.chrom, IntervalTree())
        tree.addi(int(mapping.start), int(mapping.end), int(mapping.array_idx))

    for _, _, _, baf, _, sample_id in _iter_roi_baf_records(
        baf_path,
        roi_intervals,
        header,
        column_names,
    ):
        count, min_baf, max_baf = sample_baf_stats.get(sample_id, (0, baf, baf))
        if count == 0:
            sample_baf_stats[sample_id] = (1, baf, baf)
        else:
            sample_baf_stats[sample_id] = (
                count + 1,
                min(min_baf, baf),
                max(max_baf, baf),
            )

    excluded_samples = {
        sample_id
        for sample_id, (count, min_baf, max_baf) in sample_baf_stats.items()
        if count > 1 and min_baf == max_baf
    }
    excluded_rows = sum(sample_baf_stats[sample_id][0] for sample_id in excluded_samples)

    if excluded_samples:
        print(
            f"  Excluding BAF for {len(excluded_samples)} sample(s) with constant "
            f"ROI BAF values ({excluded_rows:,} rows skipped)"
        )

    baf_values_by_bin_sample: Dict[Tuple[int, str], List[float]] = defaultdict(list)

    with gzip.open(output_path, "wt") as out_handle:
        out_handle.write("Chr\tPos\tBAF\tSample\n")
        for chrom, pos, pos_text, baf, baf_text, sample_id in _iter_roi_baf_records(
            baf_path,
            roi_intervals,
            header,
            column_names,
        ):
            if sample_id in excluded_samples:
                continue

            out_handle.write(f"{chrom}\t{pos_text}\t{baf_text}\t{sample_id}\n")
            written_rows += 1

            tree = interval_trees.get(chrom)
            if tree is None:
                continue
            for hit in tree.at(pos):
                baf_values_by_bin_sample[(int(hit.data), sample_id)].append(baf)

    summary_rows: List[dict] = []
    summary_columns = [
        "cluster",
        "interval",
        "chr",
        "start",
        "end",
        "array_idx",
        "sample",
        "baf_median",
        "minor_baf_median",
        "baf_variance",
        "baf_empirical_var",
        "baf_n_sites",
    ]
    prior_site_var = 0.01
    for (array_idx, sample_id), baf_values in baf_values_by_bin_sample.items():
        mapping = mapping_by_idx.get(array_idx)
        if mapping is None:
            continue
        values = np.asarray(baf_values, dtype=float)
        if values.size == 0:
            continue
        minor_values = np.minimum(values, 1.0 - values)
        n_sites = int(values.size)
        raw_median = float(np.median(values))
        minor_median = float(np.median(minor_values))

        empirical_var = float(np.var(minor_values, ddof=1)) if n_sites > 1 else np.nan
        if np.isnan(empirical_var):
            shrunken_site_var = prior_site_var
        else:
            shrunken_site_var = (
                ((n_sites - 1) * empirical_var) + prior_site_var
            ) / n_sites
        median_variance = float((np.pi / (2.0 * n_sites)) * shrunken_site_var)

        summary_rows.append({
            "cluster": mapping.cluster,
            "interval": mapping.interval_name,
            "chr": mapping.chrom,
            "start": mapping.start,
            "end": mapping.end,
            "array_idx": array_idx,
            "sample": sample_id,
            "baf_median": raw_median,
            "minor_baf_median": minor_median,
            "baf_variance": median_variance,
            "baf_empirical_var": empirical_var,
            "baf_n_sites": n_sites,
        })

    summary_df = pd.DataFrame(summary_rows, columns=summary_columns)
    summary_df.to_csv(summary_path, sep="\t", index=False, compression="gzip")

    print("  Saved preprocessed BAF table")
    print(f"  Rows: {written_rows:,} ROI-filtered BAF rows retained")
    print("  Saved preprocessed BAF summary table")
    print(f"  Rows: {len(summary_df):,} bin × sample BAF summaries")
    return output_path


def load_preprocessed_data(
    preprocessed_dir: str,
) -> Tuple[pd.DataFrame, List[LocusBinMapping], Optional[pd.DataFrame]]:
    """Load preprocessed bins and bin-to-interval mappings from disk.

    This is the complement of :func:`write_preprocessed_bins` +
    :func:`~gatk_sv_gd.output.write_locus_metadata`: it reads the files
    produced by the ``preprocess`` subcommand and reconstructs the objects
    needed by :func:`~gatk_sv_gd.infer.run_gd_analysis`.

    Args:
        preprocessed_dir: Directory written by the ``preprocess`` subcommand.

    Returns:
        Tuple of ``(combined_df, mappings, baf_summary_df)`` where
        ``baf_summary_df`` is ``None`` when no BAF summary file exists.
    """
    bins_path = os.path.join(preprocessed_dir, "preprocessed_bins.tsv.gz")
    mappings_path = os.path.join(preprocessed_dir, "bin_mappings.tsv.gz")
    baf_summary_path = os.path.join(preprocessed_dir, "preprocessed_baf_summary.tsv.gz")

    print("  Loading preprocessed bins")
    combined_df = pd.read_csv(bins_path, sep="\t", compression="infer")
    print(f"    {len(combined_df)} bins")

    print("  Loading bin mappings")
    mappings_df = pd.read_csv(mappings_path, sep="\t", compression="infer")
    print(f"    {len(mappings_df)} mapping records")

    mappings: List[LocusBinMapping] = []
    for _, row in mappings_df.iterrows():
        mappings.append(LocusBinMapping(
            cluster=row["cluster"],
            locus=None,
            interval_name=row["interval"],
            array_idx=int(row["array_idx"]),
            chrom=row["chr"],
            start=int(row["start"]),
            end=int(row["end"]),
        ))

    baf_summary_df: Optional[pd.DataFrame] = None
    if os.path.exists(baf_summary_path):
        print("  Loading BAF summaries")
        baf_summary_df = pd.read_csv(baf_summary_path, sep="\t", compression="infer")
        print(f"    {len(baf_summary_df)} BAF summary rows")

    return combined_df, mappings, baf_summary_df


# =============================================================================
# CLI
# =============================================================================


def parse_args():
    """Parse command-line arguments for the preprocess subcommand."""
    parser = argparse.ArgumentParser(
        description="Preprocess read-depth data for GD CNV inference",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Input/Output
    parser.add_argument(
        "-i", "--input", required=True,
        help="Input TSV with normalised read depth (bins × samples)",
    )
    parser.add_argument(
        "-g", "--gd-table", required=True,
        help="GD locus definition table (TSV)",
    )
    parser.add_argument(
        "-e", "--exclusion-intervals", nargs="+", action="append", default=[],
        help="BED file(s) of regions to mask (segdups, centromeres, etc.)",
    )
    parser.add_argument(
        "--flank-exclusion-intervals", nargs="+", action="append", default=[],
        help="BED file(s) of regions to mask only in flanking regions, not GD body intervals",
    )
    parser.add_argument(
        "--par-intervals", nargs="+", action="append", default=[],
        help="BED file(s) of pseudoautosomal intervals to exclude from flanks and ignore during ploidy-aware body filtering",
    )
    parser.add_argument(
        "-o", "--output-dir", required=True,
        help="Output directory for preprocessed files",
    )
    parser.add_argument(
        "--high-res-counts", required=False,
        help="Optional bgzipped tabix-indexed high-res count file "
             "(.tsv.gz + .tbi) for under-covered intervals",
    )
    parser.add_argument(
        "--baf-table", required=False,
        help="Optional BAF table from SiteDepthtoBAF (plain TSV or .gz). "
             "When provided, preprocess writes preprocessed_baf.tsv.gz "
             "filtered to the retained GD regions and analyzed samples.",
    )

    # Locus processing
    parser.add_argument("--locus-padding", type=int, default=10000,
                        help="Padding around locus boundaries (bp)")
    parser.add_argument("--exclusion-threshold", type=float, default=0.1,
                        help="Min overlap fraction with exclusion regions to mask a bin")
    parser.add_argument("--exclusion-bypass-threshold", type=float, default=0.8,
                        help="Body-interval overlap fraction above which masking is skipped")
    parser.add_argument("--min-bins-per-interval", type=int, default=10,
                        help="Hard-failure minimum bins per body interval")
    parser.add_argument("--max-bins-per-interval", type=int, default=20,
                        help="Maximum bins per body interval after rebinning "
                             "(0 = no rebinning)")
    parser.add_argument("--min-rebin-coverage", type=float, default=0.5,
                        help="Min coverage fraction for rebinned bins")
    parser.add_argument("--min-flank-bases", type=int, default=50000,
                        help="Min base pairs each flank must cover")
    parser.add_argument("--max-flank-bases", type=int, default=1000000,
                        help="Max base-pair target per flank (caps locus-size target)")
    parser.add_argument("--min-flank-bins", type=int, default=10,
                        help="Min bins each flank must contain")
    parser.add_argument("--min-flank-coverage", type=float, default=0.5,
                        help="Min fraction of flank bp target that must be covered")

    # Region restriction
    parser.add_argument(
        "--region", dest="regions", action="append", default=None,
        metavar="REGION",
        help=(
            "Restrict processing to GD loci overlapping this region.  "
            "Can be a chromosome (e.g. chr22) or an interval "
            "(e.g. chr1:3000000-4000000).  May be specified multiple times."
        ),
    )

    # Data filtering
    parser.add_argument("--skip-bin-filter", action="store_true", default=False,
                        help="Skip bin quality filtering")
    parser.add_argument("--median-min", type=float, default=1.0,
                        help="Min median depth for bins")
    parser.add_argument("--median-max", type=float, default=3.0,
                        help="Max median depth for bins")
    parser.add_argument("--mad-max", type=float, default=1.0,
                        help="Max MAD for bins")

    # Verbosity
    parser.add_argument(
        "-v", "--verbose", action="store_true", default=False,
        help="Enable detailed per-bin diagnostic logging",
    )

    return parser.parse_args()


def main():
    """Run the preprocessing pipeline (data loading through bin collection)."""
    args = parse_args()
    args.exclusion_intervals = _flatten_multi_args(args.exclusion_intervals)
    args.flank_exclusion_intervals = _flatten_multi_args(args.flank_exclusion_intervals)
    args.par_intervals = _flatten_multi_args(args.par_intervals)
    _util.VERBOSE = args.verbose

    os.makedirs(args.output_dir, exist_ok=True)
    setup_logging(
        args.output_dir,
        filename="preprocess_log.txt",
        verbose=args.verbose,
        command="preprocess",
        args=args,
    )
    print("Output directory configured")

    # Load GD table
    print("\nLoading GD table")
    gd_table = GDTable(args.gd_table)
    validate_gd_table_for_preprocess(gd_table)
    loci_with_breakpoints = sum(1 for locus in gd_table.loci.values() if locus.breakpoints)
    total_breakpoints = sum(locus.n_breakpoints for locus in gd_table.loci.values())
    print(
        f"Loaded {len(gd_table.loci)} loci; "
        f"{loci_with_breakpoints} with breakpoints; "
        f"{total_breakpoints} total breakpoint intervals"
    )

    # Load exclusion mask — all BED files are concatenated first so that
    # cross-file overlapping intervals are merged before the index is built.
    exclusion_mask = None
    if args.exclusion_intervals:
        print(f"\nLoading {len(args.exclusion_intervals)} exclusion interval file(s)")
        exclusion_mask = ExclusionMask(
            args.exclusion_intervals,
            label="exclusion regions",
        )

    flank_exclusion_mask = None
    if args.flank_exclusion_intervals:
        print(f"\nLoading {len(args.flank_exclusion_intervals)} flank exclusion interval file(s)")
        flank_exclusion_mask = ExclusionMask(
            args.flank_exclusion_intervals,
            label="flank exclusion regions",
        )

    par_mask = None
    if args.par_intervals:
        print(f"\nLoading {len(args.par_intervals)} PAR interval file(s)")
        par_mask = ExclusionMask(
            args.par_intervals,
            label="pseudoautosomal intervals",
        )

    # Load and normalise read depth data
    df = read_data(args.input)
    sample_cols = get_sample_columns(df)

    if df["Chr"].map(_is_chr_x).any() and par_mask is None:
        raise ValueError(
            "Ploidy-aware preprocessing requires --par-intervals when chrX bins are present. "
            "Provide PAR BED interval(s) for this reference build so those bins can be "
            "ignored during filtering until we add explicit PAR ploidy support."
        )

    autosome_mask = ~df["Chr"].isin(["chrX", "chrY"])
    if autosome_mask.any():
        column_medians = np.median(df.loc[autosome_mask, sample_cols], axis=0)
    else:
        column_medians = np.median(df[sample_cols], axis=0)
    print(f"Column medians: min={column_medians.min():.3f}, "
          f"max={column_medians.max():.3f}, mean={column_medians.mean():.3f}")

    lowres_bin_sizes = (df["End"] - df["Start"]).values
    lowres_median_bin_size = float(np.median(lowres_bin_sizes))
    print(f"Low-res median bin size: {lowres_median_bin_size:,.0f} bp")

    # Normalise: CN=2 corresponds to depth of 2.0
    df[sample_cols] = 2.0 * df[sample_cols] / column_medians[np.newaxis, :]

    # Estimate ploidy (before quality filtering so the ploidy map is
    # available for ploidy-adjusted median/MAD computation)
    ploidy_df = estimate_ploidy(df, args.output_dir)
    ploidy_map = build_ploidy_map(ploidy_df)

    # Filter low quality bins
    if not args.skip_bin_filter:
        df = filter_low_quality_bins(
            df, median_min=args.median_min,
            median_max=args.median_max, mad_max=args.mad_max,
            ploidy_map=ploidy_map,
            par_mask=par_mask,
        )

    # Build quality-filter params
    filter_params: dict = {
        "median_min": args.median_min,
        "median_max": args.median_max,
        "mad_max": args.mad_max,
    }
    highres_path: Optional[str] = getattr(args, "high_res_counts", None)

    # Parse --region arguments
    parsed_regions = None
    if args.regions:
        parsed_regions = [_parse_region(r) for r in args.regions]
        print(f"\nRestricting to {len(parsed_regions)} region filter(s)")

    # Collect all bins across all loci
    combined_df, mappings, included_loci = collect_all_locus_bins(
        df, gd_table, exclusion_mask,
        flank_exclusion_mask=flank_exclusion_mask,
        par_mask=par_mask,
        exclusion_threshold=args.exclusion_threshold,
        locus_padding=args.locus_padding,
        min_bins_per_interval=args.min_bins_per_interval,
        max_bins_per_interval=args.max_bins_per_interval,
        highres_counts_path=highres_path,
        column_medians=column_medians,
        lowres_median_bin_size=lowres_median_bin_size,
        filter_params=filter_params,
        exclusion_bypass_threshold=args.exclusion_bypass_threshold,
        min_rebin_coverage=args.min_rebin_coverage,
        min_flank_bases=args.min_flank_bases,
        max_flank_bases=args.max_flank_bases,
        min_flank_bins=args.min_flank_bins,
        min_flank_coverage=args.min_flank_coverage,
        regions=parsed_regions,
        ploidy_map=ploidy_map,
    )

    if len(combined_df) == 0:
        raise RuntimeError(
            "No loci survived preprocessing — the output would be empty. "
            "Check that the GD table and --region filters leave at least "
            "one NAHR locus with sufficient bins."
        )

    # Write preprocessed outputs
    print("\n" + "=" * 80)
    print("WRITING PREPROCESSED DATA")
    print("=" * 80)
    write_preprocessed_bins(combined_df, args.output_dir)
    write_locus_metadata(included_loci, mappings, args.output_dir)
    if args.baf_table:
        print("\nFiltering BAF table to retained GD regions...")
        write_preprocessed_baf(
            args.baf_table,
            mappings,
            args.output_dir,
        )

    # Write a filtered GD table containing only the loci that survived
    # region and size filtering, so downstream tools (call, plot, eval)
    # operate on exactly the modeled set.
    #
    # _parse_loci() maps standalone entries (NaN / empty cluster column)
    # to their coordinate-derived standalone locus key, so we must apply
    # the same mapping
    # when filtering the raw DataFrame — otherwise standalone entries
    # would be silently dropped from the output even though their
    # intervals were processed.
    included_clusters = set(included_loci.keys())
    effective_cluster = gd_table.df["cluster"].copy()
    standalone_mask = effective_cluster.isna() | (effective_cluster == "")
    effective_cluster.loc[standalone_mask] = gd_table.df.loc[
        standalone_mask
    ].apply(GDTable._standalone_locus_key, axis=1)
    filtered_gd_df = gd_table.df[effective_cluster.isin(included_clusters)].copy()
    filtered_gd_path = os.path.join(args.output_dir, "gd_table_filtered.tsv")
    filtered_gd_df.to_csv(filtered_gd_path, sep="\t", index=False)
    print("  Saved filtered GD table")
    print(f"  Rows: {len(filtered_gd_df):,} entries "
          f"({len(included_clusters)} loci of {len(gd_table.loci)} total)")

    print("\n" + "=" * 80)
    print("PREPROCESSING COMPLETE")
    print("=" * 80)
    output_count = 6 + (2 if args.baf_table else 0)
    print(f"\nOutput tables written: {output_count}")
    if args.baf_table:
        print("  BAF outputs included")
    print("\nNext: run infer with the preprocessed output directory")
    print("      Use the filtered GD table output for call, plot, and eval")
    print("=" * 80)


if __name__ == "__main__":
    main()


