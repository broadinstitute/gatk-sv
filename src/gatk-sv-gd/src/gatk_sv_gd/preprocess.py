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
import os
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from gatk_sv_gd import _util
from gatk_sv_gd._util import get_sample_columns, setup_logging
from gatk_sv_gd.models import GDLocus, GDTable
from gatk_sv_gd.depth import ExclusionMask
from gatk_sv_gd.bins import (
    LocusBinMapping,
    assign_bins_to_intervals,
    compute_flank_regions_from_bins,
    extract_locus_bins,
    filter_low_quality_bins,
    read_data,
    rebin_locus_intervals,
)
from gatk_sv_gd.highres import normalize_highres_bins, query_highres_bins
from gatk_sv_gd.output import estimate_ploidy, write_locus_metadata


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
                f"Invalid region format '{region_str}': expected chrom:start-end"
            )
        return chrom, int(parts[0]), int(parts[1])
    return region_str, None, None


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
    exclusion_threshold: float = 0.5,
    filter_params: Optional[dict] = None,
    exclusion_bypass_regions: Optional[List[Tuple[int, int]]] = None,
    min_rebin_coverage: float = 0.5,
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
        exclusion_threshold: Overlap fraction threshold for masking.
        filter_params: If provided, a dict with keys ``median_min``,
            ``median_max``, ``mad_max`` to apply
            :func:`filter_low_quality_bins` to *locus_df* before any other
            processing.  Skipped when *None*.
        exclusion_bypass_regions: Genomic ranges where masking is
            skipped (bins whose midpoint falls in a bypass range are kept).

    Returns:
        Tuple of (processed DataFrame, interval_bins dict).
    """
    if len(locus_df) == 0:
        return locus_df, {}

    # Optional exclusion masking (high-res bins need this; low-res already filtered)
    if exclusion_mask is not None:
        bypass = exclusion_bypass_regions or []
        keep_mask = []
        n_bypassed = 0
        for _, row in locus_df.iterrows():
            bin_mid = (row["Start"] + row["End"]) / 2
            in_bypass = any(bs <= bin_mid < be for bs, be in bypass)
            if in_bypass:
                keep_mask.append(True)
                n_bypassed += 1
            else:
                overlap = exclusion_mask.get_overlap_fraction(
                    row["Chr"], row["Start"], row["End"]
                )
                keep_mask.append(overlap < exclusion_threshold)
        n_masked = len(keep_mask) - sum(keep_mask)
        if n_masked > 0:
            print(f"    Masked {n_masked}/{len(locus_df)} bins due to exclusion overlap")
        if n_bypassed > 0:
            print(f"    Bypassed exclusion masking for {n_bypassed} bins in heavily-overlapped intervals")
        locus_df = locus_df[keep_mask].copy()

    # Optional quality filtering (used for high-res bins)
    if filter_params is not None and len(locus_df) > 0:
        sample_cols = get_sample_columns(locus_df)
        depths = locus_df[sample_cols].values
        medians = np.median(depths, axis=1)
        mads = np.median(np.abs(depths - medians[:, np.newaxis]), axis=1)
        keep = (
            (medians >= filter_params["median_min"])
            & (medians <= filter_params["median_max"])
            & (mads <= filter_params["mad_max"])
        )
        n_filt = int((~keep).sum())
        if n_filt > 0:
            print(f"    Quality-filtered {n_filt}/{len(locus_df)} bins")
        if _util.VERBOSE:
            for j in range(len(locus_df)):
                row = locus_df.iloc[j]
                status = "kept" if keep[j] else "FILTERED"
                reason_parts = []
                if medians[j] < filter_params["median_min"]:
                    reason_parts.append(f"median {medians[j]:.3f} < {filter_params['median_min']}")
                if medians[j] > filter_params["median_max"]:
                    reason_parts.append(f"median {medians[j]:.3f} > {filter_params['median_max']}")
                if mads[j] > filter_params["mad_max"]:
                    reason_parts.append(f"MAD {mads[j]:.3f} > {filter_params['mad_max']}")
                reason = "; ".join(reason_parts) if reason_parts else ""
                print(f"      [verbose] bin {row['Chr']}:{row['Start']}-{row['End']}  "
                      f"median={medians[j]:.3f} MAD={mads[j]:.3f}  {status}"
                      f"{' (' + reason + ')' if reason else ''}")
        locus_df = locus_df[keep].copy()

    if len(locus_df) == 0:
        return locus_df, {}

    # Trim to active region
    bin_mids = (locus_df["Start"] + locus_df["End"]) / 2
    locus_df = locus_df[(bin_mids >= left_bound) & (bin_mids < right_bound)].copy()
    print(f"    Bins after trimming to active region [{left_bound:,}, {right_bound:,}): {len(locus_df)}")

    # Rebin
    if max_bins_per_interval > 0 and len(locus_df) > 0:
        n_before = len(locus_df)
        locus_df = rebin_locus_intervals(locus_df, locus, max_bins_per_interval, flank_regions,
                                         min_rebin_coverage=min_rebin_coverage)
        if len(locus_df) < n_before:
            print(f"    Bins after rebinning: {len(locus_df)} (reduced from {n_before})")

    # Assign to intervals
    interval_bins = assign_bins_to_intervals(locus_df, locus, flank_regions) if len(locus_df) > 0 else {}

    # Mask breakpoint-range bins: bins inside the breakpoint SD-block ranges
    # (not in any body interval or flank) carry unreliable depth signal.
    bp_range_indices = interval_bins.pop("breakpoint_ranges", [])
    if bp_range_indices:
        print(f"    Masking {len(bp_range_indices)} breakpoint-range bin(s)")
        locus_df = locus_df.drop(index=bp_range_indices)

    return locus_df, interval_bins


def collect_all_locus_bins(
    df: pd.DataFrame,
    gd_table: GDTable,
    exclusion_mask: Optional[ExclusionMask],
    exclusion_threshold: float = 0.5,
    locus_padding: int = 0,
    min_bins_per_interval: int = 3,
    min_non_nahr_bins_per_interval: int = 100,
    max_bins_per_interval: int = 10,
    non_nahr_max_bins_per_interval: int = 100,
    highres_counts_path: Optional[str] = None,
    column_medians: Optional[np.ndarray] = None,
    lowres_median_bin_size: Optional[float] = None,
    filter_params: Optional[dict] = None,
    exclusion_bypass_threshold: float = 0.9,
    min_rebin_coverage: float = 0.5,
    min_flank_bases: int = 50000,
    min_flank_bins: int = 10,
    min_flank_coverage: float = 0.1,
    regions: Optional[List[Tuple[str, Optional[int], Optional[int]]]] = None,
    min_non_nahr_size: int = 100000,
) -> Tuple[pd.DataFrame, List[LocusBinMapping], Dict[str, GDLocus]]:
    """
    Collect all bins across all GD loci into a single DataFrame.

    Each body interval (region between adjacent breakpoints) must contain
    at least *min_bins_per_interval* bins after all processing (hard
    failure).  Non-NAHR loci use the higher *min_non_nahr_bins_per_interval*
    target to trigger high-res replacement when available.

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
        exclusion_threshold: Minimum overlap fraction with the mask to
            exclude a bin.
        locus_padding: Padding around locus boundaries.
        min_bins_per_interval: Hard-failure minimum: each body interval
            must have at least this many bins after all processing,
            regardless of locus type (NAHR or non-NAHR).
        min_non_nahr_bins_per_interval: Target bin count per body interval
            for non-NAHR loci.  When any body interval of a non-NAHR locus
            has fewer bins than this, high-res replacement is triggered
            (if ``highres_counts_path`` is set).  Default 100.
        max_bins_per_interval: Maximum bins per interval after rebinning
            for NAHR loci (0 = no rebinning).
        non_nahr_max_bins_per_interval: Maximum bins per interval after
            rebinning for non-NAHR loci (default 100).  Higher resolution
            is needed because non-NAHR CNVs can span a small fraction of
            the region.
        highres_counts_path: Optional path to a bgzipped, tabix-indexed
            high-resolution read-count file.  When set, loci with any
            under-covered interval are re-queried at this resolution.
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
        min_flank_bins: Minimum number of bins each flank must contain
            regardless of base-pair thresholds (default 10).
        min_flank_coverage: Minimum fraction of the effective bp target
            that a flank's accumulated bin coverage must reach.  Flanks
            below this threshold are rejected (default 0.1 = 10%).
        regions: Optional list of parsed region tuples ``(chrom, start,
            end)`` to restrict processing to.  Only loci overlapping at
            least one region are included.  ``None`` means process all.
        min_non_nahr_size: Minimum size (bp) for non-NAHR loci.  Non-NAHR
            loci whose span (outermost breakpoints) is smaller than this
            value are skipped.  Set to 0 to disable (default 100 000).

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
    for chrom, chrom_df in df.groupby("Chr"):
        keep = np.ones(len(chrom_df), dtype=bool)
        n_excluded = 0
        n_quality = 0
        if exclusion_mask is not None:
            overlaps = exclusion_mask.get_overlap_fractions_batch(
                chrom, chrom_df["Start"].values, chrom_df["End"].values
            )
            exclusion_keep = overlaps < exclusion_threshold
            n_excluded = int((~exclusion_keep).sum())
            keep &= exclusion_keep
        if _flank_filter_params is not None:
            _sc = get_sample_columns(chrom_df)
            _depths = chrom_df[_sc].values
            _meds = np.median(_depths, axis=1)
            _mads = np.median(np.abs(_depths - _meds[:, np.newaxis]), axis=1)
            quality_keep = (
                (_meds >= _flank_filter_params["median_min"])
                & (_meds <= _flank_filter_params["median_max"])
                & (_mads <= _flank_filter_params["mad_max"])
            )
            n_quality = int((~quality_keep).sum())
            keep &= quality_keep
        chrom_filtered[chrom] = chrom_df[keep].copy()
        parts = [f"{len(chrom_filtered[chrom])} filtered bins (of {len(chrom_df)} total)"]
        if n_excluded > 0:
            parts.append(f"{n_excluded} exclusion-masked")
        if n_quality > 0:
            parts.append(f"{n_quality} quality-filtered")
        print(f"  {chrom}: {', '.join(parts)}")

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
            bin_mid = (bin_row["Start"] + bin_row["End"]) / 2

            assigned_interval = None
            for start, end, name in all_named_regions:
                if start <= bin_mid < end:
                    assigned_interval = name
                    break

            if assigned_interval is None:
                # Breakpoint-range bins should have been masked upstream;
                # if one slips through, warn and skip rather than silently
                # including unreliable signal.
                print(f"  WARNING: bin {bin_row['Chr']}:{int(bin_row['Start'])}-"
                      f"{int(bin_row['End'])} in cluster {cluster} does not match "
                      f"any interval or flank — skipping (likely breakpoint range)")
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
    for cluster, locus in gd_table.get_all_loci().items():
        # Skip loci outside requested regions
        if regions is not None and not _locus_overlaps_regions(locus, regions):
            continue

        # Skip non-NAHR loci below the minimum size threshold
        if not locus.is_nahr and min_non_nahr_size > 0:
            locus_span = locus.end - locus.start
            if locus_span < min_non_nahr_size:
                print(f"\nSkipping non-NAHR locus: {cluster}  "
                      f"({locus.chrom}:{locus.start:,}-{locus.end:,}, "
                      f"{locus_span:,} bp < {min_non_nahr_size:,} bp minimum)")
                continue

        print(f"\nProcessing locus: {cluster}")
        print(f"  Chromosome: {locus.chrom}")
        print(f"  Breakpoints: {locus.breakpoints}")
        print(f"  GD entries: {len(locus.gd_entries)} ({', '.join(locus.svtypes)})")

        locus_size = locus.end - locus.start
        effective_max_bins = max_bins_per_interval if locus.is_nahr else non_nahr_max_bins_per_interval
        if not locus.is_nahr:
            print(f"  Non-NAHR locus: using {effective_max_bins} bins/interval (higher resolution)")

        # Compute flank coordinates from the full chromosome's filtered bins so
        # that even multi-megabase exclusion deserts don't prevent flank discovery.
        chrom_bins = chrom_filtered.get(locus.chrom, pd.DataFrame())
        flank_regions = compute_flank_regions_from_bins(
            chrom_bins, locus, locus_size,
            min_flank_bases=min_flank_bases,
            min_flank_bins=min_flank_bins,
            min_flank_coverage=min_flank_coverage,
            filter_params=_flank_filter_params,
        )
        if flank_regions:
            for fs, fe, fn in flank_regions:
                print(f"  {fn}: {fs:,}-{fe:,} (bin-derived)")
        else:
            print(f"  Warning: no flanking bins found for locus {cluster}")

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
                    print(f"  Interval '{iv_name}' is {iv_overlap:.0%} exclusion-overlapped "
                          f"(>= {exclusion_bypass_threshold:.0%}); bypassing exclusion masking")
                    exclusion_bypass_regions.append((iv_start, iv_end))

        # Extract only the bins within the active region (locus + flanks),
        # applying exclusion masking. locus_padding still applies as a minimum.
        active_padding = max(locus_padding, locus.start - left_bound, right_bound - locus.end)
        locus_df = extract_locus_bins(
            df, locus, exclusion_mask,
            exclusion_threshold=exclusion_threshold,
            padding=active_padding,
            exclusion_bypass_regions=exclusion_bypass_regions,
        )

        print(f"  Bins after filtering: {len(locus_df)}")

        # Trim to active region: drop any bins outside [left_bound, right_bound)
        bin_mids = (locus_df["Start"] + locus_df["End"]) / 2
        locus_df = locus_df[(bin_mids >= left_bound) & (bin_mids < right_bound)].copy()
        print(f"  Bins after trimming to active region [{left_bound:,}, {right_bound:,}): {len(locus_df)}")

        # Rebin to reduce number of bins per interval/flank if requested
        if effective_max_bins > 0:
            locus_df_orig = locus_df
            if _util.VERBOSE:
                _sc = get_sample_columns(locus_df_orig)
                _meds = np.median(locus_df_orig[_sc].values, axis=1)
                print(f"    [verbose] pre-rebin per-bin median depths ({len(locus_df_orig)} bins):")
                for _j in range(len(locus_df_orig)):
                    _r = locus_df_orig.iloc[_j]
                    print(f"      {_r['Chr']}:{_r['Start']}-{_r['End']}  "
                          f"median={_meds[_j]:.3f}")

            locus_df = rebin_locus_intervals(locus_df, locus, effective_max_bins, flank_regions,
                                             min_rebin_coverage=min_rebin_coverage)
            if len(locus_df) < len(locus_df_orig):
                print(f"  Bins after rebinning: {len(locus_df)} (reduced from {len(locus_df_orig)})")

            if _util.VERBOSE:
                _sc = get_sample_columns(locus_df)
                _meds = np.median(locus_df[_sc].values, axis=1)
                print(f"    [verbose] post-rebin per-bin median depths ({len(locus_df)} bins):")
                for _j in range(len(locus_df)):
                    _r = locus_df.iloc[_j]
                    print(f"      {_r['Chr']}:{_r['Start']}-{_r['End']}  "
                          f"median={_meds[_j]:.3f}")

        # Assign bins to intervals and flanking regions
        interval_bins = assign_bins_to_intervals(locus_df, locus, flank_regions)
        total_assigned = sum(len(v) for v in interval_bins.values())
        for region_name, bins in interval_bins.items():
            print(f"    {region_name}: {len(bins)} bins")
        print(f"    total: {total_assigned} bins")

        # Mask breakpoint-range bins: bins inside the breakpoint SD-block
        # ranges (not in any body interval or flank) carry unreliable depth
        # signal and must be excluded from model training.
        bp_range_indices = interval_bins.pop("breakpoint_ranges", [])
        if bp_range_indices:
            print(f"    Masking {len(bp_range_indices)} breakpoint-range bin(s)")
            locus_df = locus_df.drop(index=bp_range_indices)

        # ------------------------------------------------------------------
        # Per-interval bin count check + optional high-res replacement
        # ------------------------------------------------------------------
        body_intervals = locus.get_intervals()  # [(start, end, name), ...]
        body_interval_names = {name for _, _, name in body_intervals}

        # Determine per-interval bin target for this locus type.
        # Non-NAHR loci need finer resolution because the variant can
        # span only a small fraction of the region.
        effective_min_bins = (
            min_non_nahr_bins_per_interval if not locus.is_nahr
            else min_bins_per_interval
        )

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
        if (
            undercovered
            and highres_counts_path is not None
            and column_medians is not None
            and lowres_median_bin_size is not None
        ):
            print(f"\n  *** {len(undercovered)} body interval(s) have fewer "
                  f"than {effective_min_bins} bins; switching to high-res ***")
            for name, cnt in undercovered:
                print(f"      {name}: {cnt} bins")

            # Query the tabix-indexed high-res file for the active region
            hr_raw_df = query_highres_bins(
                highres_counts_path,
                locus.chrom,
                left_bound,
                right_bound,
                sample_cols,
            )
            print(f"    [highres] queried {len(hr_raw_df)} raw bins "
                  f"in {locus.chrom}:{left_bound:,}-{right_bound:,}")

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
                        print(f"    [highres] body interval '{iv_name}' has only "
                              f"{raw_count} raw bins (< {min_bins_per_interval}); "
                              f"resolution-limited, will relax threshold")

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
                    effective_max_bins,
                    exclusion_mask=exclusion_mask,
                    exclusion_threshold=exclusion_threshold,
                    filter_params=filter_params,
                    exclusion_bypass_regions=exclusion_bypass_regions,
                    min_rebin_coverage=min_rebin_coverage,
                )

                hr_undercovered = _undercovered_intervals(hr_interval_bins, effective_min_bins)

                # Accept the high-res data if it improved coverage:
                # (a) fewer under-covered intervals, OR
                # (b) same number of under-covered intervals but more bins
                #     in those intervals (e.g. 5 bins vs 0 is a real
                #     improvement even if both are still below threshold).
                hr_uc_bins = sum(cnt for _, cnt in hr_undercovered)
                orig_uc_bins = sum(cnt for _, cnt in undercovered)
                accept_hr = (
                    len(hr_undercovered) < len(undercovered)
                    or (
                        len(hr_undercovered) <= len(undercovered)
                        and hr_uc_bins > orig_uc_bins
                    )
                )

                if accept_hr:
                    locus_df = hr_locus_df
                    interval_bins = hr_interval_bins
                    used_highres = True
                    total_assigned = sum(len(v) for v in interval_bins.values())
                    print(f"    [highres] ACCEPTED — replaced low-res bins")
                    for region_name, bins in interval_bins.items():
                        print(f"      {region_name}: {len(bins)} bins")
                    print(f"      total: {total_assigned} bins")
                else:
                    print(f"    [highres] REJECTED — high-res did not "
                          f"reduce under-covered intervals")
            else:
                print(f"    [highres] no bins returned from tabix query")

        # Hard-failure check: use the base threshold (min_bins_per_interval)
        # which applies uniformly to all locus types.  Non-NAHR loci may
        # have been checked against a higher target above to trigger
        # high-res replacement, but as long as they meet the base minimum
        # they are kept.
        #
        # Intervals that are "physically limited" (the raw high-res data
        # itself didn't have enough bins) are exempted as long as they
        # have > 0 bins after processing.
        hard_check = _undercovered_intervals(interval_bins, min_bins_per_interval)
        if hard_check:
            hard_failures = [
                (name, cnt) for name, cnt in hard_check
                if name not in hr_physically_limited or cnt == 0
            ]
            soft_warnings = [
                (name, cnt) for name, cnt in hard_check
                if name in hr_physically_limited and cnt > 0
            ]

            if soft_warnings:
                for name, cnt in soft_warnings:
                    print(f"  NOTE: body interval '{name}' has {cnt} bins "
                          f"(< {min_bins_per_interval}), but this is the maximum "
                          f"the bin resolution can provide — proceeding anyway.")

            if hard_failures:
                details = "\n".join(
                    f"    interval '{name}': {cnt} bins (need >= {min_bins_per_interval})"
                    for name, cnt in hard_failures
                )
                raise ValueError(
                    f"\n  WARNING: Locus {cluster} ({locus.chrom}:{locus.start:,}-"
                    f"{locus.end:,}) has {len(hard_failures)} body interval(s) with "
                    f"fewer than --min-bins-per-interval={min_bins_per_interval} bins "
                    f"after all processing"
                    f"{' (including high-res replacement)' if used_highres else ''}:\n"
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
    print(f"  Saved: {output_path}")
    print(f"  Rows: {len(combined_df):,}")
    return output_path


def load_preprocessed_data(
    preprocessed_dir: str,
) -> Tuple[pd.DataFrame, List[LocusBinMapping]]:
    """Load preprocessed bins and bin-to-interval mappings from disk.

    This is the complement of :func:`write_preprocessed_bins` +
    :func:`~gatk_sv_gd.output.write_locus_metadata`: it reads the files
    produced by the ``preprocess`` subcommand and reconstructs the objects
    needed by :func:`~gatk_sv_gd.infer.run_gd_analysis`.

    Args:
        preprocessed_dir: Directory written by the ``preprocess`` subcommand.

    Returns:
        Tuple of (combined_df, mappings).
    """
    bins_path = os.path.join(preprocessed_dir, "preprocessed_bins.tsv.gz")
    mappings_path = os.path.join(preprocessed_dir, "bin_mappings.tsv.gz")

    print(f"  Loading preprocessed bins: {bins_path}")
    combined_df = pd.read_csv(bins_path, sep="\t", compression="infer")
    print(f"    {len(combined_df)} bins")

    print(f"  Loading bin mappings: {mappings_path}")
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

    return combined_df, mappings


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
        "-e", "--exclusion-intervals", nargs="*", default=[],
        help="BED file(s) of regions to mask (segdups, centromeres, etc.)",
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

    # Locus processing
    parser.add_argument("--locus-padding", type=int, default=10000,
                        help="Padding around locus boundaries (bp)")
    parser.add_argument("--exclusion-threshold", type=float, default=1.0,
                        help="Min overlap fraction with exclusion regions to mask a bin")
    parser.add_argument("--exclusion-bypass-threshold", type=float, default=0.6,
                        help="Body-interval overlap fraction above which masking is skipped")
    parser.add_argument("--min-bins-per-interval", type=int, default=10,
                        help="Hard-failure minimum bins per body interval, applied "
                             "uniformly to NAHR and non-NAHR loci")
    parser.add_argument("--min-non-nahr-bins-per-interval", type=int, default=100,
                        help="Target bins per body interval for non-NAHR loci; "
                             "triggers high-res replacement when not met")
    parser.add_argument("--max-bins-per-interval", type=int, default=20,
                        help="Maximum bins per body interval after rebinning for "
                             "NAHR loci (0 = no rebinning)")
    parser.add_argument("--non-nahr-max-bins-per-interval", type=int, default=100,
                        help="Maximum bins per body interval after rebinning for "
                             "non-NAHR loci (higher resolution needed because "
                             "non-NAHR CNVs can span a small fraction of the region)")
    parser.add_argument("--min-rebin-coverage", type=float, default=0.5,
                        help="Min coverage fraction for rebinned bins")
    parser.add_argument("--min-flank-bases", type=int, default=50000,
                        help="Min base pairs each flank must cover")
    parser.add_argument("--min-flank-bins", type=int, default=10,
                        help="Min bins each flank must contain")
    parser.add_argument("--min-flank-coverage", type=float, default=0.5,
                        help="Min fraction of flank bp target that must be covered")
    parser.add_argument("--min-non-nahr-size", type=int, default=1000000,
                        help="Minimum size (bp) for non-NAHR loci; smaller loci are skipped "
                             "(set to 0 to disable)")

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
    parser.add_argument("--mad-max", type=float, default=2.0,
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
    _util.VERBOSE = args.verbose

    os.makedirs(args.output_dir, exist_ok=True)
    setup_logging(args.output_dir, filename="preprocess_log.txt")
    print(f"Output directory: {args.output_dir}")

    # Load GD table
    print(f"\nLoading GD table: {args.gd_table}")
    gd_table = GDTable(args.gd_table)
    print(f"Loaded {len(gd_table.loci)} loci")
    for cluster, locus in gd_table.loci.items():
        if locus.breakpoints:
            bp_start = min(bp[0] for bp in locus.breakpoints)
            bp_end = max(bp[1] for bp in locus.breakpoints)
            print(f"  {cluster}: {locus.chrom}:{bp_start}-{bp_end} "
                  f"({len(locus.gd_entries)} entries, {locus.n_breakpoints} breakpoints)")
        else:
            print(f"  {cluster}: {locus.chrom} - NO BREAKPOINTS DEFINED")

    # Load exclusion mask
    exclusion_mask = None
    for bed_path in args.exclusion_intervals:
        bed_label = os.path.basename(bed_path)
        if exclusion_mask is None:
            print(f"\nLoading exclusion mask: {bed_path}")
            exclusion_mask = ExclusionMask(bed_path, label=bed_label)
        else:
            print(f"\nMerging exclusion intervals: {bed_path}")
            exclusion_mask.merge(bed_path, label=bed_label)

    # Load and normalise read depth data
    df = read_data(args.input)
    sample_cols = get_sample_columns(df)

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

    # Filter low quality bins
    if not args.skip_bin_filter:
        df = filter_low_quality_bins(
            df, median_min=args.median_min,
            median_max=args.median_max, mad_max=args.mad_max,
        )

    # Estimate ploidy
    estimate_ploidy(df, args.output_dir)

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
        print(f"\nRestricting to {len(parsed_regions)} region(s):")
        for chrom, start, end in parsed_regions:
            if start is None:
                print(f"  {chrom}")
            else:
                print(f"  {chrom}:{start:,}-{end:,}")

    # Collect all bins across all loci
    combined_df, mappings, included_loci = collect_all_locus_bins(
        df, gd_table, exclusion_mask,
        exclusion_threshold=args.exclusion_threshold,
        locus_padding=args.locus_padding,
        min_bins_per_interval=args.min_bins_per_interval,
        min_non_nahr_bins_per_interval=args.min_non_nahr_bins_per_interval,
        max_bins_per_interval=args.max_bins_per_interval,
        non_nahr_max_bins_per_interval=args.non_nahr_max_bins_per_interval,
        highres_counts_path=highres_path,
        column_medians=column_medians,
        lowres_median_bin_size=lowres_median_bin_size,
        filter_params=filter_params,
        exclusion_bypass_threshold=args.exclusion_bypass_threshold,
        min_rebin_coverage=args.min_rebin_coverage,
        min_flank_bases=args.min_flank_bases,
        min_flank_bins=args.min_flank_bins,
        min_flank_coverage=args.min_flank_coverage,
        regions=parsed_regions,
        min_non_nahr_size=args.min_non_nahr_size,
    )

    if len(combined_df) == 0:
        raise RuntimeError(
            "No loci survived preprocessing — the output would be empty. "
            "Check that the GD table, --region filters, and "
            "--min-non-nahr-size threshold leave at least one locus with "
            "sufficient bins."
        )

    # Write preprocessed outputs
    print("\n" + "=" * 80)
    print("WRITING PREPROCESSED DATA")
    print("=" * 80)
    write_preprocessed_bins(combined_df, args.output_dir)
    write_locus_metadata(included_loci, mappings, args.output_dir)

    # Write a filtered GD table containing only the loci that survived
    # region and size filtering, so downstream tools (call, plot, eval)
    # operate on exactly the modeled set.
    #
    # _parse_loci() maps standalone entries (NaN / empty cluster column)
    # to their GD_ID as the dict key, so we must apply the same mapping
    # when filtering the raw DataFrame — otherwise standalone entries
    # would be silently dropped from the output even though their
    # intervals were processed.
    included_clusters = set(included_loci.keys())
    effective_cluster = gd_table.df["cluster"].copy()
    standalone_mask = effective_cluster.isna() | (effective_cluster == "")
    effective_cluster[standalone_mask] = gd_table.df.loc[standalone_mask, "GD_ID"]
    filtered_gd_df = gd_table.df[effective_cluster.isin(included_clusters)].copy()
    filtered_gd_path = os.path.join(args.output_dir, "gd_table_filtered.tsv")
    filtered_gd_df.to_csv(filtered_gd_path, sep="\t", index=False)
    print(f"  Saved: {filtered_gd_path}")
    print(f"  Rows: {len(filtered_gd_df):,} entries "
          f"({len(included_clusters)} loci of {len(gd_table.loci)} total)")

    print("\n" + "=" * 80)
    print("PREPROCESSING COMPLETE")
    print("=" * 80)
    print("\nOutput files:")
    print("  - preprocessed_bins.tsv.gz")
    print("  - bin_mappings.tsv.gz")
    print("  - locus_intervals.tsv.gz")
    print("  - gd_entry_intervals.tsv.gz")
    print("  - ploidy_estimates.tsv")
    print("  - gd_table_filtered.tsv")
    print(f"\nNext: run 'gatk-sv-gd infer --preprocessed-dir {args.output_dir}'")
    print(f"       Use --gd-table {args.output_dir}/gd_table_filtered.tsv for call/plot/eval")
    print("=" * 80)


if __name__ == "__main__":
    main()


