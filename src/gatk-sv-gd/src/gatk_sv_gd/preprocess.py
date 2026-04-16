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
from typing import Dict, List, Optional, Tuple

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
                f"Invalid region format '{region_str}': expected chrom:start-end"
            )
        return chrom, int(parts[0]), int(parts[1])
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
        par_mask: Optional pseudoautosomal interval mask. Bins that
            overlap PAR are temporarily exempted from quality filtering.

    Returns:
        Tuple of (processed DataFrame, interval_bins dict).
    """
    if len(locus_df) == 0:
        return locus_df, {}

    # Trim to active region
    locus_df = locus_df[
        (locus_df["End"] > left_bound) & (locus_df["Start"] < right_bound)
    ].copy()
    print(f"    Bins in active region [{left_bound:,}, {right_bound:,}): {len(locus_df)}")

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
            currently exempted from quality filtering while ploidy-aware
            handling remains contig-wide.
        exclusion_threshold: Minimum overlap fraction with the mask to
            exclude a bin.
        locus_padding: Padding around locus boundaries.
        min_bins_per_interval: Hard-failure minimum: each body interval
            must have at least this many bins after all processing.
        max_bins_per_interval: Maximum bins per interval after rebinning
            (0 = no rebinning).
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
    for chrom, chrom_df in df.groupby("Chr"):
        keep = np.ones(len(chrom_df), dtype=bool)
        n_excluded = 0
        n_quality = 0
        n_par_ignored = 0
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
        if chrom_flank_filter_params is not None:
            quality_keep, quality_stats = compute_bin_quality_mask(
                chrom_df,
                median_min=chrom_flank_filter_params["median_min"],
                median_max=chrom_flank_filter_params["median_max"],
                mad_max=chrom_flank_filter_params["mad_max"],
                ploidy_map=ploidy_map,
                par_mask=par_mask,
            )
            n_quality = quality_stats["filtered"]
            n_par_ignored = quality_stats["par_ignored"]
            keep &= quality_keep
        chrom_filtered[chrom] = chrom_df[keep].copy()
        parts = [f"{len(chrom_filtered[chrom])} filtered bins (of {len(chrom_df)} total)"]
        if n_excluded > 0:
            parts.append(f"{n_excluded} exclusion-masked")
        if n_quality > 0:
            parts.append(f"{n_quality} quality-filtered")
        if n_par_ignored > 0:
            parts.append(f"{n_par_ignored} PAR-ignored")
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

        # Only NAHR loci are analyzed
        if not locus.is_nahr:
            print(f"\nSkipping non-NAHR locus: {cluster}  "
                  f"({locus.chrom}:{locus.start:,}-{locus.end:,})")
            continue

        print(f"\nProcessing locus: {cluster}")
        print(f"  Chromosome: {locus.chrom}")
        print(f"  Breakpoints: {locus.breakpoints}")
        print(f"  GD entries: {len(locus.gd_entries)} ({', '.join(locus.svtypes)})")

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

        # ------------------------------------------------------------------
        # Step 2 (body): extract body bins with exclusion threshold + bypass
        # ------------------------------------------------------------------
        body_df = extract_locus_bins(
            df, locus, exclusion_mask,
            exclusion_threshold=exclusion_threshold,
            padding=locus_padding,
            exclusion_bypass_regions=exclusion_bypass_regions,
        )
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
                    print(f"  {fn} bins from filtered cache: {len(matched)}")
        flank_df = pd.concat(flank_dfs) if flank_dfs else pd.DataFrame(columns=df.columns)

        # Combine body + flank bins
        parts = [p for p in [body_df, flank_df] if len(p) > 0]
        if parts:
            locus_df = pd.concat(parts).drop_duplicates(
                subset=["Chr", "Start", "End"],
            ).sort_values("Start").reset_index(drop=True)
        else:
            locus_df = pd.DataFrame(columns=df.columns)
        print(f"  Combined bins (body + flanks): {len(locus_df)}")

        # Rebin to reduce number of bins per interval/flank if requested
        if max_bins_per_interval > 0:
            locus_df_orig = locus_df
            if _util.VERBOSE:
                _sc = get_sample_columns(locus_df_orig)
                _meds = np.median(locus_df_orig[_sc].values, axis=1)
                print(f"    [verbose] pre-rebin per-bin median depths ({len(locus_df_orig)} bins):")
                for _j in range(len(locus_df_orig)):
                    _r = locus_df_orig.iloc[_j]
                    print(f"      {_r['Chr']}:{_r['Start']}-{_r['End']}  "
                          f"median={_meds[_j]:.3f}")

            locus_df = rebin_locus_intervals(locus_df, locus, max_bins_per_interval, flank_regions,
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
        if (
            undercovered and
            highres_counts_path is not None and
            column_medians is not None and
            lowres_median_bin_size is not None
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

                hr_undercovered = _undercovered_intervals(hr_interval_bins, effective_min_bins)

                # Accept the high-res data if it improved coverage:
                # (a) fewer under-covered intervals, OR
                # (b) same number of under-covered intervals but more bins
                #     in those intervals (e.g. 5 bins vs 0 is a real
                #     improvement even if both are still below threshold).
                hr_uc_bins = sum(cnt for _, cnt in hr_undercovered)
                orig_uc_bins = sum(cnt for _, cnt in undercovered)
                accept_hr = (
                    len(hr_undercovered) < len(undercovered) or
                    (
                        len(hr_undercovered) <= len(undercovered) and
                        hr_uc_bins > orig_uc_bins
                    )
                )

                if accept_hr:
                    locus_df = hr_locus_df
                    interval_bins = hr_interval_bins
                    used_highres = True
                    total_assigned = sum(len(v) for v in interval_bins.values())
                    print("    [highres] ACCEPTED — replaced low-res bins")
                    for region_name, bins in interval_bins.items():
                        print(f"      {region_name}: {len(bins)} bins")
                    print(f"      total: {total_assigned} bins")
                else:
                    print(
                        "    [highres] REJECTED — high-res did not "
                        "reduce under-covered intervals"
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
        raise ValueError(
            f"BAF table {filepath} must have at least 4 tab-delimited columns"
        )

    try:
        int(parts[1])
        float(parts[2])
    except ValueError:
        return 0, parts[:4]

    return None, ["Chr", "Pos", "BAF", "Sample"]


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

    interval_trees: Dict[str, IntervalTree] = {}
    mapping_by_idx = {int(m.array_idx): m for m in mappings}
    for mapping in mappings:
        tree = interval_trees.setdefault(mapping.chrom, IntervalTree())
        tree.addi(int(mapping.start), int(mapping.end), int(mapping.array_idx))

    baf_values_by_bin_sample: Dict[Tuple[int, str], List[float]] = defaultdict(list)

    with gzip.open(output_path, "wt") as out_handle:
        out_handle.write("Chr\tPos\tBAF\tSample\n")

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

                        out_handle.write(
                            f"{record['Chr']}\t{record['Pos']}\t"
                            f"{record['BAF']}\t{record['Sample']}\n"
                        )
                        written_rows += 1

                        tree = interval_trees.get(record["Chr"])
                        if tree is None:
                            continue
                        for hit in tree.at(pos):
                            baf_values_by_bin_sample[(int(hit.data), str(record["Sample"]))].append(baf)
        finally:
            tbx.close()

    summary_rows: List[dict] = []
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

    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(summary_path, sep="\t", index=False, compression="gzip")

    print(f"  Saved: {output_path}")
    print(f"  Rows: {written_rows:,} ROI-filtered BAF rows retained")
    print(f"  Saved: {summary_path}")
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

    baf_summary_df: Optional[pd.DataFrame] = None
    if os.path.exists(baf_summary_path):
        print(f"  Loading BAF summaries: {baf_summary_path}")
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
        help="BED file(s) of pseudoautosomal intervals to ignore during ploidy-aware filtering",
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
    parser.add_argument("--exclusion-threshold", type=float, default=0.5,
                        help="Min overlap fraction with exclusion regions to mask a bin")
    parser.add_argument("--exclusion-bypass-threshold", type=float, default=0.6,
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
    setup_logging(args.output_dir, filename="preprocess_log.txt")
    print(f"Output directory: {args.output_dir}")

    # Load GD table
    print(f"\nLoading GD table: {args.gd_table}")
    gd_table = GDTable(args.gd_table)
    validate_gd_table_for_preprocess(gd_table)
    print(f"Loaded {len(gd_table.loci)} loci")
    for cluster, locus in gd_table.loci.items():
        if locus.breakpoints:
            bp_start = min(bp[0] for bp in locus.breakpoints)
            bp_end = max(bp[1] for bp in locus.breakpoints)
            print(f"  {cluster}: {locus.chrom}:{bp_start}-{bp_end} "
                  f"({len(locus.gd_entries)} entries, {locus.n_breakpoints} breakpoints)")
        else:
            print(f"  {cluster}: {locus.chrom} - NO BREAKPOINTS DEFINED")

    # Load exclusion mask — all BED files are concatenated first so that
    # cross-file overlapping intervals are merged before the index is built.
    exclusion_mask = None
    if args.exclusion_intervals:
        print(f"\nLoading {len(args.exclusion_intervals)} exclusion interval file(s):")
        for p in args.exclusion_intervals:
            print(f"  {p}")
        exclusion_mask = ExclusionMask(
            args.exclusion_intervals,
            label="exclusion regions",
        )

    flank_exclusion_mask = None
    if args.flank_exclusion_intervals:
        print(f"\nLoading {len(args.flank_exclusion_intervals)} flank exclusion interval file(s):")
        for p in args.flank_exclusion_intervals:
            print(f"  {p}")
        flank_exclusion_mask = ExclusionMask(
            args.flank_exclusion_intervals,
            label="flank exclusion regions",
        )

    par_mask = None
    if args.par_intervals:
        print(f"\nLoading {len(args.par_intervals)} PAR interval file(s):")
        for p in args.par_intervals:
            print(f"  {p}")
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
        print(f"\nRestricting to {len(parsed_regions)} region(s):")
        for chrom, start, end in parsed_regions:
            if start is None:
                print(f"  {chrom}")
            else:
                print(f"  {chrom}:{start:,}-{end:,}")

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
    if args.baf_table:
        print("  - preprocessed_baf.tsv.gz")
        print("  - preprocessed_baf_summary.tsv.gz")
    print(f"\nNext: run 'gatk-sv-gd infer --preprocessed-dir {args.output_dir}'")
    print(f"       Use --gd-table {args.output_dir}/gd_table_filtered.tsv for call/plot/eval")
    print("=" * 80)


if __name__ == "__main__":
    main()


