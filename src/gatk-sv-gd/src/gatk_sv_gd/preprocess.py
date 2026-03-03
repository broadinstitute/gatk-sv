"""
High-resolution bin processing and locus bin collection.

Handles querying tabix-indexed high-resolution count files, normalising
them to match the low-resolution scale, and assembling all bins across
GD loci into a single combined DataFrame.
"""

from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import pysam

from gatk_sv_gd import _util
from gatk_sv_gd._util import get_sample_columns
from gatk_sv_gd.models import GDLocus, GDTable
from gatk_sv_gd.depth import ExclusionMask
from gatk_sv_gd.bins import (
    LocusBinMapping,
    assign_bins_to_intervals,
    compute_flank_regions_from_bins,
    extract_locus_bins,
    filter_low_quality_bins,
    rebin_locus_intervals,
)

def query_highres_bins(
    highres_path: str,
    chrom: str,
    start: int,
    end: int,
    sample_cols: List[str],
) -> pd.DataFrame:
    """
    Query a tabix-indexed read-count file for bins overlapping a region.

    The file is expected to have a header line starting with ``#Chr`` (or
    ``Chr``) followed by ``Start``, ``End``, and one column per sample with
    raw (un-normalised) read counts.  Only the intersection of columns
    present in the file *and* in *sample_cols* is returned so that the
    result is compatible with the low-resolution DataFrame.

    Args:
        highres_path: Path to a bgzipped, tabix-indexed TSV (.tsv.gz + .tbi).
        chrom: Chromosome name (e.g. ``"chr15"``).
        start: Query start position (0-based, inclusive).
        end: Query end position (0-based, exclusive).
        sample_cols: Sample column names expected in the output.

    Returns:
        DataFrame with columns ``Chr``, ``Start``, ``End``, ``source_file``
        and one column per sample, indexed by a ``Bin`` string.
    """
    tbx = pysam.TabixFile(highres_path)
    try:
        header_line = tbx.header[-1] if tbx.header else ""
        header_cols = header_line.lstrip("#").split("\t")

        rows: List[dict] = []
        for line in tbx.fetch(chrom, max(0, start), end):
            fields = line.split("\t")
            row = dict(zip(header_cols, fields))
            rows.append(row)
    finally:
        tbx.close()

    if len(rows) == 0:
        cols = ["Chr", "Start", "End", "source_file"] + list(sample_cols)
        empty = pd.DataFrame(columns=cols)
        empty["Bin"] = pd.Series(dtype=str)
        return empty.set_index("Bin")

    raw_df = pd.DataFrame(rows)

    # Normalise column naming
    if "#Chr" in raw_df.columns:
        raw_df.rename(columns={"#Chr": "Chr"}, inplace=True)

    raw_df["Start"] = raw_df["Start"].astype(int)
    raw_df["End"] = raw_df["End"].astype(int)
    raw_df["source_file"] = "highres"

    # Keep only sample columns that exist in both files
    common_samples = [s for s in sample_cols if s in raw_df.columns]
    if len(common_samples) == 0:
        raise ValueError(
            "High-resolution counts file shares no sample columns with the "
            "low-resolution file. Check that both files were generated from "
            "the same sample set."
        )
    for s in common_samples:
        raw_df[s] = pd.to_numeric(raw_df[s], errors="coerce")

    # Pad with NaN for any samples absent from the high-res file —
    # add all missing columns at once to avoid DataFrame fragmentation.
    missing = [s for s in sample_cols if s not in raw_df.columns]
    if missing:
        raw_df = pd.concat(
            [raw_df, pd.DataFrame(np.nan, index=raw_df.index, columns=missing)],
            axis=1,
        )

    raw_df["Bin"] = (
        raw_df["Chr"].astype(str) + ":"
        + raw_df["Start"].astype(str) + "-"
        + raw_df["End"].astype(str)
    )
    raw_df = raw_df.set_index("Bin")

    keep_cols = ["Chr", "Start", "End", "source_file"] + list(sample_cols)
    return raw_df[keep_cols]


def normalize_highres_bins(
    highres_df: pd.DataFrame,
    sample_cols: List[str],
    column_medians: np.ndarray,
    lowres_median_bin_size: float,
) -> pd.DataFrame:
    """
    Normalise high-resolution raw counts to the same CN-2 ≈ 2.0 scale used
    for low-resolution bins.

    The low-resolution pipeline normalises each sample by dividing by its
    genome-wide autosomal median count (``column_medians``) and multiplying
    by 2.  Those medians were estimated from bins of a specific size.  For
    high-res bins that are smaller, we must account for the expected
    proportional decrease in raw counts::

        norm_depth = 2.0 * raw_count / (column_median * highres_bin_size / lowres_bin_size)

    This is equivalent to first scaling ``column_medians`` by the bin-size
    ratio, then applying the standard normalisation.

    Args:
        highres_df: DataFrame returned by :func:`query_highres_bins` with
            **raw** (un-normalised) counts.
        sample_cols: Sample column names.
        column_medians: Per-sample autosomal median raw counts estimated
            from the low-resolution file (1-D array, same order as
            *sample_cols*).
        lowres_median_bin_size: Median bin size (bp) in the low-resolution
            file, used to scale the expected counts.

    Returns:
        A copy of *highres_df* with sample columns normalised in-place.
    """
    df = highres_df.copy()
    highres_bin_sizes = (df["End"] - df["Start"]).values
    highres_median_bin_size = float(np.median(highres_bin_sizes))

    bin_size_ratio = highres_median_bin_size / lowres_median_bin_size
    print(f"    [highres] low-res median bin size: {lowres_median_bin_size:,.0f} bp")
    print(f"    [highres] high-res median bin size: {highres_median_bin_size:,.0f} bp")
    print(f"    [highres] bin size ratio: {bin_size_ratio:.4f}")

    # Scale column_medians by the bin-size ratio so the normalisation
    # accounts for the smaller expected raw counts in high-res bins.
    adjusted_medians = column_medians * bin_size_ratio  # shape (n_samples,)

    if _util.VERBOSE:
        raw_vals = df[sample_cols].values
        print(f"    [verbose] high-res pre-normalisation: "
              f"per-bin-median mean={np.nanmean(np.nanmedian(raw_vals, axis=1)):.3f}, "
              f"adjusted_medians range=[{adjusted_medians.min():.3f}, {adjusted_medians.max():.3f}]")

    df[sample_cols] = 2.0 * df[sample_cols].values / adjusted_medians[np.newaxis, :]

    if _util.VERBOSE:
        norm_vals = df[sample_cols].values
        per_bin_medians = np.nanmedian(norm_vals, axis=1)
        print(f"    [verbose] high-res post-normalisation: "
              f"per-bin-median mean={per_bin_medians.mean():.3f}, "
              f"min={per_bin_medians.min():.3f}, max={per_bin_medians.max():.3f}")
        # Log first few bins for spot-checking
        for i in range(min(5, len(df))):
            row = df.iloc[i]
            med = np.nanmedian(row[sample_cols].values.astype(float))
            print(f"      bin {row['Chr']}:{row['Start']}-{row['End']}  "
                  f"norm_median={med:.3f}")
        if len(df) > 5:
            print(f"      ... ({len(df) - 5} more bins)")

    return df


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
    min_bins_per_region: int = 3,
    max_bins_per_interval: int = 10,
    highres_counts_path: Optional[str] = None,
    column_medians: Optional[np.ndarray] = None,
    lowres_median_bin_size: Optional[float] = None,
    filter_params: Optional[dict] = None,
    exclusion_bypass_threshold: float = 0.9,
    min_rebin_coverage: float = 0.5,
    min_flank_bases: int = 50000,
    min_flank_bins: int = 10,
    min_flank_coverage: float = 0.1,
) -> Tuple[pd.DataFrame, List[LocusBinMapping], Dict[str, GDLocus]]:
    """
    Collect all bins across all GD loci into a single DataFrame.

    Each body interval (region between adjacent breakpoints) should contain
    at least *min_bins_per_region* bins after processing.  Intervals that
    fall below this threshold trigger a warning.

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
        min_bins_per_region: Minimum number of bins expected in each
            body interval (region between adjacent breakpoints).  A warning
            is printed if any interval has fewer bins after all processing.
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
        min_flank_bins: Minimum number of bins each flank must contain
            regardless of base-pair thresholds (default 10).
        min_flank_coverage: Minimum fraction of the effective bp target
            that a flank's accumulated bin coverage must reach.  Flanks
            below this threshold are rejected (default 0.1 = 10%).

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
        print(f"\nProcessing locus: {cluster}")
        print(f"  Chromosome: {locus.chrom}")
        print(f"  Breakpoints: {locus.breakpoints}")
        print(f"  GD entries: {len(locus.gd_entries)} ({', '.join(locus.svtypes)})")

        locus_size = locus.end - locus.start

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

        # Identify under-covered body intervals
        def _undercovered_intervals(
            ivbins: Dict[str, List[int]],
        ) -> List[Tuple[str, int]]:
            """Return list of (name, count) for body intervals below threshold."""
            return [
                (name, len(ivbins.get(name, [])))
                for name in body_interval_names
                if len(ivbins.get(name, [])) < min_bins_per_region
            ]

        undercovered = _undercovered_intervals(interval_bins)

        used_highres = False
        # Body intervals where even the raw high-res data doesn't have
        # enough bins — the resolution is simply too coarse for this locus
        # size, so we exempt them from the min-bins threshold later.
        hr_physically_limited: set = set()
        if (
            undercovered
            and highres_counts_path is not None
            and column_medians is not None
            and lowres_median_bin_size is not None
        ):
            print(f"\n  *** {len(undercovered)} body interval(s) have fewer "
                  f"than {min_bins_per_region} bins; switching to high-res ***")
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
                # reach min_bins_per_region even with every raw bin, it is
                # physically limited by the bin resolution and should be
                # exempted from the threshold (provided it has > 0 bins
                # after processing).
                hr_raw_mids = (hr_raw_df["Start"] + hr_raw_df["End"]) / 2
                for iv_start, iv_end, iv_name in body_intervals:
                    raw_count = int(((hr_raw_mids >= iv_start) & (hr_raw_mids < iv_end)).sum())
                    if raw_count < min_bins_per_region:
                        hr_physically_limited.add(iv_name)
                        print(f"    [highres] body interval '{iv_name}' has only "
                              f"{raw_count} raw bins (< {min_bins_per_region}); "
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
                    exclusion_threshold=exclusion_threshold,
                    filter_params=filter_params,
                    exclusion_bypass_regions=exclusion_bypass_regions,
                    min_rebin_coverage=min_rebin_coverage,
                )

                hr_undercovered = _undercovered_intervals(hr_interval_bins)

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
                    undercovered = hr_undercovered
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

        # Warn about under-covered body intervals.  Some GD loci have
        # intervals dominated by exclusion regions so 0 surviving bins
        # is possible even after bypass — downgrade to a warning.
        #
        # Intervals that are "physically limited" (the raw high-res data
        # itself didn't have enough bins) are exempted from the threshold
        # as long as they still have > 0 bins after processing — the bin
        # resolution is simply too coarse for the locus size.
        if undercovered:
            # Partition into truly problematic vs resolution-limited
            hard_failures = [
                (name, cnt) for name, cnt in undercovered
                if name not in hr_physically_limited or cnt == 0
            ]
            soft_warnings = [
                (name, cnt) for name, cnt in undercovered
                if name in hr_physically_limited and cnt > 0
            ]

            if soft_warnings:
                for name, cnt in soft_warnings:
                    print(f"  NOTE: body interval '{name}' has {cnt} bins "
                          f"(< {min_bins_per_region}), but this is the maximum "
                          f"the bin resolution can provide — proceeding anyway.")

            if hard_failures:
                details = "\n".join(
                    f"    interval '{name}': {cnt} bins (need >= {min_bins_per_region})"
                    for name, cnt in hard_failures
                )
                raise ValueError(
                    f"\n  WARNING: Locus {cluster} ({locus.chrom}:{locus.start:,}-"
                    f"{locus.end:,}) has {len(hard_failures)} body interval(s) with "
                    f"fewer than --min-bins-per-region={min_bins_per_region} bins "
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


