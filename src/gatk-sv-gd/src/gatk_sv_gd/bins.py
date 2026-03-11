"""
Bin processing for GD CNV analysis.

Functions for extracting, filtering, rebinning, and assigning genomic bins
to locus intervals. Also includes I/O helpers for reading depth data.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from gatk_sv_gd import _util
from gatk_sv_gd._util import get_sample_columns
from gatk_sv_gd.models import GDLocus
from gatk_sv_gd.depth import ExclusionMask

def extract_locus_bins(
    df: pd.DataFrame,
    locus: GDLocus,
    exclusion_mask: Optional[ExclusionMask] = None,
    exclusion_threshold: float = 0.5,
    padding: int = 0,
    exclusion_bypass_regions: Optional[List[Tuple[int, int]]] = None,
) -> pd.DataFrame:
    """
    Extract bins overlapping a GD locus, optionally masking regions
    that overlap the supplied exclusion mask.

    Args:
        df: DataFrame with all bins
        locus: GDLocus object defining the region
        exclusion_mask: Optional ExclusionMask for filtering.
        exclusion_threshold: Minimum overlap fraction with the mask to
            exclude a bin.
        padding: Extend the region by this many base pairs on each side.
        exclusion_bypass_regions: List of (start, end) genomic ranges where
            masking is skipped.  Bins that have any true interval overlap
            with a bypass range are kept regardless of mask overlap.

    Returns:
        DataFrame with bins for this locus
    """
    chrom = locus.chrom
    # Get overall start and end from breakpoint ranges
    start = min(bp[0] for bp in locus.breakpoints) - padding
    end = max(bp[1] for bp in locus.breakpoints) + padding

    # Filter to bins in this region
    mask = (
        (df["Chr"] == chrom) &
        (df["End"] > start) &
        (df["Start"] < end)
    )
    locus_df = df[mask].copy()

    if len(locus_df) == 0:
        return locus_df

    # Apply exclusion masking if provided
    if exclusion_mask is not None:
        bypass = exclusion_bypass_regions or []
        bin_starts = locus_df["Start"].values
        bin_ends = locus_df["End"].values

        # True interval-overlap bypass check (no midpoint heuristic).
        # A bin is "in bypass" if it has any real overlap with a bypass range:
        #   bin_start < bypass_end  AND  bin_end > bypass_start
        in_bypass = np.zeros(len(locus_df), dtype=bool)
        for bs, be in bypass:
            in_bypass |= (bin_starts < be) & (bin_ends > bs)

        # Batch exclusion overlap only for non-bypass bins
        overlaps = np.zeros(len(locus_df))
        non_bypass = ~in_bypass
        if non_bypass.any():
            overlaps[non_bypass] = exclusion_mask.get_overlap_fractions_batch(
                chrom,
                bin_starts[non_bypass],
                bin_ends[non_bypass],
            )

        keep_mask = in_bypass | (overlaps < exclusion_threshold)
        n_bypassed = int(in_bypass.sum())
        n_masked = int((~keep_mask).sum())

        if n_masked > 0:
            print(f"  Masked {n_masked}/{len(locus_df)} bins due to exclusion overlap")
        if n_bypassed > 0:
            print(f"  Bypassed exclusion masking for {n_bypassed} bins in heavily-overlapped intervals")

        if _util.VERBOSE:
            for i, (_, row) in enumerate(locus_df.iterrows()):
                if in_bypass[i]:
                    print(f"    [verbose] bin {row['Chr']}:{row['Start']}-{row['End']} "
                          f"in exclusion bypass region — kept")
                elif not keep_mask[i]:
                    sample_cols = get_sample_columns(locus_df)
                    med = np.median(row[sample_cols].values.astype(float))
                    print(f"    [verbose] bin {row['Chr']}:{row['Start']}-{row['End']} "
                          f"exclusion overlap={overlaps[i]:.2f} >= {exclusion_threshold} — MASKED "
                          f"(median depth={med:.3f})")

        locus_df = locus_df[keep_mask]

    return locus_df


def compute_flank_regions_from_bins(
    locus_df: pd.DataFrame,
    locus: GDLocus,
    target_size: int,
    min_flank_bases: int = 50000,
    min_flank_bins: int = 10,
    min_flank_coverage: float = 0.1,
    filter_params: Optional[dict] = None,
) -> List[Tuple[int, int, str]]:
    """
    Compute flanking regions based on actual available (filtered) bins.

    Accumulates bins outward from each side of the locus until **all three**
    stopping criteria are satisfied simultaneously:

    1. Accumulated base pairs ≥ *target_size* (typically locus_size).
    2. Accumulated base pairs ≥ *min_flank_bases* (absolute floor, default 50 kb).
    3. Number of accumulated bins ≥ *min_flank_bins* (default 10).

    When *filter_params* is provided, quality filtering (per-bin median
    depth and MAD thresholds) is applied to candidate bins **during**
    accumulation.  Bins that fail the quality check are excluded so that
    only usable bins count toward the coverage targets.  This causes the
    flank to extend further to compensate for quality-failed bins, ensuring
    the resulting flank region contains enough high-quality bins.

    If the accumulated base pairs never reach the target (e.g. because a
    exclusion desert consumed most of the flanking region), the
    flank keeps extending until all available bins on that side of the
    chromosome are consumed.  The flank is always emitted — even if
    accumulated coverage is below ``min_flank_coverage`` × *effective_bp_target*
    — because hitting the chromosome boundary is the only valid reason to
    stop.  A warning is logged when the coverage fraction is not met so
    that downstream consumers can flag these loci.

    The input ``locus_df`` should already be filtered (exclusion-masked
    bins removed) so that only usable bins participate in flank discovery.
    Because masked bins are absent, flanks naturally stop at the edge of
    any large masked region rather than spanning across it.

    This ensures small loci still receive adequately sized flanks for
    establishing a reliable baseline, while large loci keep the existing
    behaviour of matching the body size.

    Args:
        locus_df: DataFrame of filtered bins available for this locus.
            Should have exclusion-masked bins already removed.
        locus: GDLocus object
        target_size: Target cumulative bin coverage in bp for each flank
        min_flank_bases: Minimum cumulative bp each flank must cover
            regardless of *target_size* (default 50 000).
        min_flank_bins: Minimum number of bins each flank must contain
            regardless of the base-pair thresholds (default 10).
        min_flank_coverage: Minimum fraction of *effective_bp_target*
            that a flank's accumulated bin coverage should reach.  When
            coverage falls below this fraction a warning is logged, but
            the flank is still kept (default 0.1 = 10%).
        filter_params: Optional dict with keys ``median_min``,
            ``median_max``, ``mad_max``.  When provided, bins failing
            these quality thresholds are excluded before accumulation
            so that the flank extends until enough *quality* bins are
            found.

    Returns:
        List of (start, end, name) tuples for "left_flank" and/or "right_flank".
        A flank is omitted only if no bins exist on that side of the locus.
    """
    # Apply quality filtering to candidate bins before accumulation
    if filter_params is not None and len(locus_df) > 0:
        sample_cols = get_sample_columns(locus_df)
        depths = locus_df[sample_cols].values
        medians = np.median(depths, axis=1)
        mads = np.median(np.abs(depths - medians[:, np.newaxis]), axis=1)
        quality_keep = (
            (medians >= filter_params["median_min"])
            & (medians <= filter_params["median_max"])
            & (mads <= filter_params["mad_max"])
        )
        n_quality_filtered = int((~quality_keep).sum())
        if n_quality_filtered > 0:
            print(f"    [flanks] quality-filtered {n_quality_filtered}/{len(locus_df)} "
                  f"candidate bins (median [{filter_params['median_min']}, "
                  f"{filter_params['median_max']}], MAD <= {filter_params['mad_max']})")
        locus_df = locus_df[quality_keep].copy()

    flanks = []
    locus_start = locus.start
    locus_end = locus.end

    # Effective base-pair threshold is the larger of target_size and min_flank_bases
    effective_bp_target = max(target_size, min_flank_bases)

    print(f"    [flanks] locus body: {locus_start:,}-{locus_end:,}  "
          f"target_size: {target_size:,} bp  min_flank_bases: {min_flank_bases:,} bp  "
          f"min_flank_bins: {min_flank_bins}  effective_bp_target: {effective_bp_target:,} bp")
    print(f"    [flanks] input bins: {len(locus_df)}  "
          f"range: {int(locus_df['Start'].min()):,}-{int(locus_df['End'].max()):,}" if len(locus_df) > 0
          else f"    [flanks] input bins: 0")

    # Use bin midpoints for candidacy — consistent with how bins are assigned
    # throughout the rest of the codebase. Using strict End/Start comparisons
    # misses bins that straddle the locus boundary (common with zero-width breakpoints).
    bin_mids = (locus_df["Start"] + locus_df["End"]) / 2

    # Left flank: bins whose midpoint is left of locus start, accumulated right-to-left
    left_bins = locus_df[bin_mids < locus_start].sort_values("Start", ascending=False)

    print(f"    [flanks] left candidate bins (mid < {locus_start:,}): {len(left_bins)}", end="")
    if len(left_bins) > 0:
        print(f"  span {int(left_bins['Start'].min()):,}-{int(left_bins['End'].max()):,}  "
              f"total bp: {int((left_bins['End'] - left_bins['Start']).sum()):,}")
        sizes = (left_bins["End"].values - left_bins["Start"].values).astype(int)
        cumsum_bp = np.cumsum(sizes)
        n_bins_arr = np.arange(1, len(sizes) + 1)
        done = (cumsum_bp >= effective_bp_target) & (n_bins_arr >= min_flank_bins)
        idx = int(np.argmax(done)) if done.any() else len(sizes) - 1
        cumulative = int(cumsum_bp[idx])
        n_bins_accumulated = int(n_bins_arr[idx])
        leftmost_start = int(left_bins.iloc[idx]["Start"])
        coverage_frac = cumulative / effective_bp_target if effective_bp_target > 0 else 1.0
        print(f"    [flanks] left_flank: {leftmost_start:,}-{locus_start:,}  "
              f"accumulated {cumulative:,} bp / {n_bins_accumulated} bins "
              f"(targets: {effective_bp_target:,} bp, {min_flank_bins} bins)  "
              f"coverage: {coverage_frac:.1%}")
        if coverage_frac < min_flank_coverage:
            print(f"    [flanks] WARNING: left flank coverage {coverage_frac:.1%} "
                  f"< {min_flank_coverage:.0%} of target — exhausted "
                  f"available bins, clearing {n_bins_accumulated} usable bins")
        else:
            flanks.append((leftmost_start, locus_start, "left_flank"))
    else:
        print()  # newline after the count

    # Right flank: bins whose midpoint is right of locus end, accumulated left-to-right
    right_bins = locus_df[bin_mids >= locus_end].sort_values("Start")

    print(f"    [flanks] right candidate bins (mid >= {locus_end:,}): {len(right_bins)}", end="")
    if len(right_bins) > 0:
        print(f"  span {int(right_bins['Start'].min()):,}-{int(right_bins['End'].max()):,}  "
              f"total bp: {int((right_bins['End'] - right_bins['Start']).sum()):,}")
        sizes = (right_bins["End"].values - right_bins["Start"].values).astype(int)
        cumsum_bp = np.cumsum(sizes)
        n_bins_arr = np.arange(1, len(sizes) + 1)
        done = (cumsum_bp >= effective_bp_target) & (n_bins_arr >= min_flank_bins)
        idx = int(np.argmax(done)) if done.any() else len(sizes) - 1
        cumulative = int(cumsum_bp[idx])
        n_bins_accumulated = int(n_bins_arr[idx])
        rightmost_end = int(right_bins.iloc[idx]["End"])
        coverage_frac = cumulative / effective_bp_target if effective_bp_target > 0 else 1.0
        print(f"    [flanks] right_flank: {locus_end:,}-{rightmost_end:,}  "
              f"accumulated {cumulative:,} bp / {n_bins_accumulated} bins "
              f"(targets: {effective_bp_target:,} bp, {min_flank_bins} bins)  "
              f"coverage: {coverage_frac:.1%}")
        if coverage_frac < min_flank_coverage:
            print(f"    [flanks] WARNING: right flank coverage {coverage_frac:.1%} "
                  f"< {min_flank_coverage:.0%} of target — exhausted "
                  f"available bins, clearing {n_bins_accumulated} usable bins")
        else:
            flanks.append((locus_end, rightmost_end, "right_flank"))
    else:
        print()  # newline after the count

    return flanks


def rebin_locus_intervals(
    df: pd.DataFrame,
    locus: GDLocus,
    max_bins_per_interval: int = 10,
    flank_regions: Optional[List[Tuple[int, int, str]]] = None,
    min_rebin_coverage: float = 0.5,
) -> pd.DataFrame:
    """
    Reduce the number of bins by rebinning each interval and flanking region to have
    at most max_bins_per_interval bins.

    For each region, the genomic range (first bin start → last bin end) is divided
    into ``max_bins_per_interval`` evenly-spaced putative intervals.  Each putative
    interval's depth is the overlap-weighted average of all original bins that
    intersect it, where the weight is the number of base pairs of overlap.  This
    guarantees uniform-width output bins regardless of the input bin size
    distribution, avoiding the problem of a few tiny bins alongside one or two
    very large ones.

    Putative intervals where the fraction of the bin width covered by original
    bins is less than *min_rebin_coverage* are discarded to avoid biased
    estimates from tiny slivers of data.

    Applies to:
      - Intervals between adjacent breakpoints
      - Left and right flanking regions (bin-derived if flank_regions is provided,
        otherwise falls back to geometric 100%-of-locus-size windows)

    Bins inside breakpoint SD-block ranges are NOT rebinned; they are passed
    through as-is and will be masked (dropped) by the caller after
    assign_bins_to_intervals() identifies them as "breakpoint_ranges".

    Args:
        df: DataFrame with bins for this locus
        locus: GDLocus object
        max_bins_per_interval: Maximum number of bins per interval (default: 10)
        min_rebin_coverage: Minimum fraction of each new bin's width that must
            be overlapped by original bins (default: 0.5).  Bins below this
            threshold are dropped.

    Returns:
        DataFrame with rebinned data
    """
    if len(df) == 0:
        return df

    # Get sample columns
    metadata_cols = ["Chr", "Start", "End", "source_file", "Bin"]
    sample_cols = [col for col in df.columns if col not in metadata_cols]

    def rebin_region(region_df: pd.DataFrame) -> List[dict]:
        """Rebin bins into evenly-spaced putative intervals using overlap-weighted averaging.

        The full genomic range of *region_df* is divided into
        *max_bins_per_interval* equal-width putative intervals.  For each
        putative interval the depth values are the overlap-weighted average
        of every original bin that intersects it (weight = bp of overlap).
        Putative intervals with zero overlap are omitted.
        """
        if len(region_df) <= max_bins_per_interval:
            return region_df.to_dict("records")

        n_new_bins = max_bins_per_interval
        region_df = region_df.sort_values("Start")

        # Evenly divide the full range into putative intervals
        range_start = int(region_df["Start"].min())
        range_end = int(region_df["End"].max())
        putative_edges = np.linspace(range_start, range_end, n_new_bins + 1)

        # Original bin coordinates as arrays for vectorised overlap computation
        bin_starts = region_df["Start"].values.astype(float)
        bin_ends = region_df["End"].values.astype(float)

        rows = []
        for i in range(n_new_bins):
            p_start = putative_edges[i]
            p_end = putative_edges[i + 1]

            # Overlap of each original bin with this putative interval (bp)
            overlap_starts = np.maximum(bin_starts, p_start)
            overlap_ends = np.minimum(bin_ends, p_end)
            overlaps = np.maximum(overlap_ends - overlap_starts, 0.0)

            total_overlap = overlaps.sum()
            if total_overlap == 0:
                continue

            # Require that a minimum fraction of the new bin is covered by
            # original bins so that tiny slivers don't produce biased estimates.
            bin_width = p_end - p_start
            if bin_width > 0 and (total_overlap / bin_width) < min_rebin_coverage:
                if _util.VERBOSE:
                    print(f"      [verbose] rebin: dropping putative bin "
                          f"{int(round(p_start))}-{int(round(p_end))} "
                          f"(coverage {total_overlap / bin_width:.1%} "
                          f"< {min_rebin_coverage:.0%})")
                continue

            weights = overlaps / total_overlap

            new_bin = {
                "Chr": region_df.iloc[0]["Chr"],
                "Start": int(round(p_start)),
                "End": int(round(p_end)),
            }
            for col in sample_cols:
                new_bin[col] = (region_df[col].values * weights).sum()
            if "source_file" in region_df.columns:
                new_bin["source_file"] = region_df.iloc[0]["source_file"]
            rows.append(new_bin)

        return rows

    # All rebinnable regions: between-breakpoint intervals + flanking regions
    all_rebin_regions = locus.get_intervals() + (
        flank_regions if flank_regions is not None else locus.get_flanking_regions()
    )

    rebinned_rows = []
    bin_starts = df["Start"].values
    bin_ends = df["End"].values
    processed_mask = np.zeros(len(df), dtype=bool)

    for region_start, region_end, _ in all_rebin_regions:
        # True interval overlap: bin overlaps region if bin_start < region_end
        # and bin_end > region_start
        mask = (bin_starts < region_end) & (bin_ends > region_start)
        if not mask.any():
            continue
        processed_mask |= mask
        rebinned_rows.extend(rebin_region(df[mask].copy()))

    # Bins inside breakpoint SD-block ranges (not covered by any body
    # interval or flank) are dropped — they carry unreliable depth signal
    # and will be masked after assign_bins_to_intervals().  Any remaining
    # unmatched bins (edge cases) are passed through as-is; they will also
    # be identified as "breakpoint_ranges" and masked downstream.
    remaining = df[~processed_mask]
    if len(remaining) > 0:
        rebinned_rows.extend(remaining.to_dict("records"))

    if len(rebinned_rows) == 0:
        return pd.DataFrame(columns=df.columns)

    result_df = pd.DataFrame(rebinned_rows)

    # A bin that straddles two adjacent regions (e.g. a body interval and a
    # flank) is included in both calls to rebin_region().  When a region is
    # small enough to skip rebinning, the original bin rows are returned
    # as-is, producing exact duplicate (Chr, Start, End) rows.  Deduplicate
    # here so that downstream code never sees the same genomic bin twice.
    if len(result_df) > 0:
        result_df = result_df.drop_duplicates(
            subset=["Chr", "Start", "End"],
        ).reset_index(drop=True)

    return result_df


def assign_bins_to_intervals(
    df: pd.DataFrame,
    locus: GDLocus,
    flank_regions: Optional[List[Tuple[int, int, str]]] = None,
) -> Dict[str, List[int]]:
    """
    Assign bins to breakpoint intervals and flanking regions within a locus.

    Each bin is assigned to the named region with which it has the greatest
    base-pair overlap (ties go to the first region in definition order).
    Bins that do not overlap any named region are placed into
    ``"breakpoint_ranges"``.

    Args:
        df: DataFrame with bins for this locus
        locus: GDLocus object
        flank_regions: Pre-computed bin-derived flank regions. If None, falls back
            to geometric windows from locus.get_flanking_regions().

    Returns:
        Dict mapping region name to list of bin indices. Includes between-breakpoint
        interval names (e.g. "1-2"), flanking names ("left_flank", "right_flank"),
        and "breakpoint_ranges" for bins that fall inside the breakpoint SD blocks.
        Callers should pop "breakpoint_ranges" and drop those bins from the
        DataFrame — they carry unreliable depth signal and must not participate
        in model training or CNV calling.
    """
    all_named_regions = locus.get_intervals() + (
        flank_regions if flank_regions is not None else locus.get_flanking_regions()
    )
    interval_bins: Dict[str, List[int]] = {name: [] for _, _, name in all_named_regions}
    interval_bins["breakpoint_ranges"] = []

    for idx, row in df.iterrows():
        bin_start, bin_end = row["Start"], row["End"]

        # Assign to the region with the greatest base-pair overlap
        best_name = None
        best_overlap = 0
        for start, end, name in all_named_regions:
            ov = min(bin_end, end) - max(bin_start, start)
            if ov > best_overlap:
                best_overlap = ov
                best_name = name

        if best_name is not None:
            interval_bins[best_name].append(idx)
        else:
            interval_bins["breakpoint_ranges"].append(idx)

    return interval_bins


def compute_interval_cn_stats(
    cn_posterior: np.ndarray,
    interval_bins: Dict[str, List[int]],
    sample_idx: int,
) -> Dict[str, dict]:
    """
    Compute copy number statistics for each interval.

    Args:
        cn_posterior: Array of shape (n_bins, n_samples, n_states)
        cn_map: Array of shape (n_bins, n_samples) with MAP CN estimates
        interval_bins: Dict mapping interval name to bin indices
        sample_idx: Sample index

    Returns:
        Dict mapping interval name to CN statistics
    """
    result = {}

    for interval_name, bin_indices in interval_bins.items():
        if len(bin_indices) == 0:
            result[interval_name] = {
                "n_bins": 0,
                "cn_probs": np.zeros(6),
            }
            continue

        # Get CN values for this interval
        interval_cn_probs = cn_posterior[bin_indices, sample_idx, :]
        
        # Average posterior probabilities across bins
        mean_probs = interval_cn_probs.mean(axis=0)

        result[interval_name] = {
            "n_bins": len(bin_indices),
            "cn_probs": mean_probs,
        }

    return result


def call_gd_cnv(
    locus: GDLocus,
    interval_stats: Dict[str, dict],
    log_prob_threshold: float = -0.5) -> List[dict]:
    """
    Call GD CNVs based on interval copy number statistics.

    For each GD entry in the locus (each DEL/DUP definition), check if the
    copy number in the corresponding region supports a CNV call.

    Args:
        locus: GDLocus object
        interval_stats: Dict mapping interval name to CN statistics
        log_prob_threshold: Minimum log probability score to call a CNV

    Returns:
        List of CNV call dictionaries
    """
    calls = []

    for entry in locus.gd_entries:
        gd_id = entry["GD_ID"]
        svtype = entry["svtype"]
        gd_start = entry["start_GRCh38"]
        gd_end = entry["end_GRCh38"]
        bp1 = entry["BP1"]
        bp2 = entry["BP2"]

        # Determine which intervals this entry spans using BP1 and BP2.
        # Use the locus's authoritative breakpoint ordering directly
        # instead of reverse-engineering from interval name strings
        # (which breaks when BP names contain hyphens).
        covered_tuples = locus.get_intervals_between(bp1, bp2)
        covered_intervals = [name for _, _, name in covered_tuples]

        if len(covered_intervals) == 0:
            continue

        # Determine if this is a CNV using weighted average of log probabilities
        # Weight by interval size (number of bins)
        log_prob_score = 0.0
        total_weight = 0.0
        
        if svtype == "DEL":
            # For affected intervals: weighted average of log P(CN < 2)
            for interval in covered_intervals:
                if interval in interval_stats and interval_stats[interval]["n_bins"] > 0:
                    probs = interval_stats[interval]["cn_probs"]
                    weight = interval_stats[interval]["n_bins"]
                    p_del = max(probs[0] + probs[1], 1e-3)  # P(CN=0) + P(CN=1)
                    if p_del > 0:
                        log_prob_score += weight * np.log(p_del)
                        total_weight += weight
        
        elif svtype == "DUP":
            # For affected intervals: weighted average of log P(CN > 2)
            for interval in covered_intervals:
                if interval in interval_stats and interval_stats[interval]["n_bins"] > 0:
                    probs = interval_stats[interval]["cn_probs"]
                    weight = interval_stats[interval]["n_bins"]
                    p_dup = max(probs[3:].sum(), 1e-3)  # P(CN >= 3)
                    if p_dup > 0:
                        log_prob_score += weight * np.log(p_dup)
                        total_weight += weight
        
        # Normalize by total weight to get weighted average
        if total_weight > 0:
            log_prob_score = log_prob_score / total_weight
        else:
            log_prob_score = np.nan  # No data to support call
        
        # Determine if this is a CNV based on log probability score
        is_cnv = log_prob_score > log_prob_threshold

        # Check flanking regions for evidence of a large spanning CNV.
        # A true GD only affects the region between its breakpoints; a spanning
        # variant would also show copy number change in the flanking regions.
        flanking_log_prob_score = np.nan
        is_spanning = False
        if is_cnv:
            flank_score_sum = 0.0
            flank_weight_total = 0.0
            for _, _, flank_name in locus.get_flanking_regions():
                if flank_name not in interval_stats:
                    continue
                flank_stat = interval_stats[flank_name]
                if flank_stat["n_bins"] == 0:
                    continue
                flank_probs = flank_stat["cn_probs"]
                flank_weight = flank_stat["n_bins"]
                if svtype == "DEL":
                    p_flank = max(flank_probs[0] + flank_probs[1], 1e-3)
                elif svtype == "DUP":
                    p_flank = max(flank_probs[3:].sum(), 1e-3)
                else:
                    continue
                flank_score_sum += flank_weight * np.log(p_flank)
                flank_weight_total += flank_weight

            if flank_weight_total > 0:
                flanking_log_prob_score = flank_score_sum / flank_weight_total
                # If flanking regions also show a CN change, this is a spanning variant
                if flanking_log_prob_score > log_prob_threshold:
                    is_spanning = True
                    is_cnv = False

        # Count total bins across covered intervals
        n_bins = sum(interval_stats[iv]["n_bins"] for iv in covered_intervals if iv in interval_stats)

        calls.append({
            "GD_ID": gd_id,
            "cluster": locus.cluster,
            "chrom": locus.chrom,
            "start": gd_start,
            "end": gd_end,
            "svtype": svtype,
            "is_terminal": locus.is_terminal,
            "log_prob_score": log_prob_score,
            "flanking_log_prob_score": flanking_log_prob_score,
            "is_carrier": is_cnv,
            "is_spanning": is_spanning,
            "intervals": covered_intervals,
            "n_bins": n_bins
        })

    return calls


def determine_best_breakpoints(
    locus: GDLocus,
    interval_stats: Dict[str, dict],
    calls: List[dict],
) -> Dict[str, Optional[str]]:
    """
    Determine the most likely breakpoint pair for a GD CNV, separately for DEL and DUP.

    For loci with multiple possible breakpoint configurations (e.g., BP1-2, BP1-3),
    determine which best fits the observed data for each svtype.

    Scoring approach:
    - For DEL: log-sum P(CN < 2) in affected intervals + log-sum P(CN >= 2) in unaffected intervals
    - For DUP: log-sum P(CN > 2) in affected intervals + log-sum P(CN <= 2) in unaffected intervals

    Args:
        locus: GDLocus object
        interval_stats: Dict mapping interval name to CN statistics
        calls: List of CNV call dictionaries

    Returns:
        Dict mapping svtype ("DEL", "DUP") to best matching GD_ID or None
    """
    best_by_svtype = {}
    
    # Get all interval names for this locus
    all_intervals = set(name for _, _, name in locus.get_intervals())

    for svtype in ["DEL", "DUP"]:
        # Filter to carrier calls of this svtype
        carrier_calls = [c for c in calls if c["is_carrier"] and c["svtype"] == svtype]

        if len(carrier_calls) == 0:
            best_by_svtype[svtype] = None
            continue

        if len(carrier_calls) == 1:
            best_by_svtype[svtype] = carrier_calls[0]["GD_ID"]
            continue

        # Multiple carrier calls of same svtype - score each based on CN probabilities
        best_score = -np.inf
        best_size = -1
        best_gd_id = None

        for call in carrier_calls:
            covered_intervals = set(call["intervals"])
            uncovered_intervals = all_intervals - covered_intervals

            # Compute weighted average score based on svtype
            score = 0.0
            total_weight = 0.0

            if svtype == "DEL":
                # For affected intervals: weighted average of log P(CN < 2)
                for interval in covered_intervals:
                    if interval in interval_stats and interval_stats[interval]["n_bins"] > 0:
                        probs = interval_stats[interval]["cn_probs"]
                        weight = interval_stats[interval]["n_bins"]
                        p_del = probs[0] + probs[1]  # P(CN=0) + P(CN=1)
                        if p_del > 0:
                            score += weight * np.log(p_del)
                            total_weight += weight

                # For unaffected intervals: weighted average of log P(CN >= 2)
                for interval in uncovered_intervals:
                    if interval in interval_stats and interval_stats[interval]["n_bins"] > 0:
                        probs = interval_stats[interval]["cn_probs"]
                        weight = interval_stats[interval]["n_bins"]
                        p_normal = probs[2:].sum()  # P(CN >= 2)
                        if p_normal > 0:
                            score += weight * np.log(p_normal)
                            total_weight += weight

            else:  # DUP
                # For affected intervals: weighted average of log P(CN > 2)
                for interval in covered_intervals:
                    if interval in interval_stats and interval_stats[interval]["n_bins"] > 0:
                        probs = interval_stats[interval]["cn_probs"]
                        weight = interval_stats[interval]["n_bins"]
                        p_dup = probs[3:].sum()  # P(CN >= 3)
                        if p_dup > 0:
                            score += weight * np.log(p_dup)
                            total_weight += weight

                # For unaffected intervals: weighted average of log P(CN <= 2)
                for interval in uncovered_intervals:
                    if interval in interval_stats and interval_stats[interval]["n_bins"] > 0:
                        probs = interval_stats[interval]["cn_probs"]
                        weight = interval_stats[interval]["n_bins"]
                        p_normal = probs[:3].sum()  # P(CN <= 2)
                        if p_normal > 0:
                            score += weight * np.log(p_normal)
                            total_weight += weight
            
            # Normalize by total weight
            if total_weight > 0:
                score = score / total_weight

            # Update best if this score is higher; tiebreaker: prefer larger variant
            call_size = call["end"] - call["start"]
            if score > best_score or (score == best_score and call_size > best_size):
                best_score = score
                best_size = call_size
                best_gd_id = call["GD_ID"]

        best_by_svtype[svtype] = best_gd_id

    return best_by_svtype


def read_data(file_path: str) -> pd.DataFrame:
    """Load binned read count data from TSV file."""
    print(f"Loading: {file_path}")
    try:
        # Read gzipped TSV file
        df_file = pd.read_csv(file_path, sep="\t", compression="infer")

        # Get source file identifier
        source_file = str(Path(file_path).parent.parent.name)

        # Rename chromosome column if it exists
        if "#Chr" in df_file.columns:
            df_file["Chr"] = df_file["#Chr"]
            df_file = df_file.drop("#Chr", axis=1)

        # Create bin identifier
        df_file["Bin"] = (
            df_file["Chr"].astype(str) +
            ":" +
            df_file["Start"].astype(str) +
            "-" +
            df_file["End"].astype(str)
        )
        df_file = df_file.set_index("Bin")
        df_file["source_file"] = source_file
        return df_file
    except Exception as e:
        print(f"Error loading {file_path}")
        raise e


def filter_low_quality_bins(
    df: pd.DataFrame,
    median_min: float = 1.5,
    median_max: float = 2.5,
    mad_max: float = 0.5,
):
    """
    Filter out low quality bins based on median and MAD thresholds.

    Args:
        df: DataFrame with bins as rows and samples as columns
        median_min: Minimum median depth for bins
        median_max: Maximum median depth for bins
        mad_max: Maximum MAD for bins

    Returns:
        Filtered DataFrame
    """
    # Get sample columns (exclude metadata)
    sample_cols = get_sample_columns(df)

    # Compute median and MAD for each bin across samples
    depths = df[sample_cols].values
    medians = np.median(depths, axis=1)
    mads = np.median(np.abs(depths - medians[:, np.newaxis]), axis=1)

    print(f"\n{'=' * 80}")
    print("FILTERING LOW QUALITY BINS")
    print(f"{'=' * 80}")
    print(f"Starting bins: {len(df)}")

    # Filter based on thresholds
    keep_mask = (
        (medians >= median_min) &
        (medians <= median_max) &
        (mads <= mad_max)
    )

    n_filtered = (~keep_mask).sum()
    print(f"Thresholds: median [{median_min}, {median_max}], MAD <= {mad_max}")
    print(f"Bins filtered: {n_filtered}")
    print(f"Bins remaining: {keep_mask.sum()}")

    if _util.VERBOSE:
        # Break down why bins were filtered
        fail_med_low = (medians < median_min).sum()
        fail_med_high = (medians > median_max).sum()
        fail_mad = (mads > mad_max).sum()
        print(f"  [verbose] Filter breakdown:")
        print(f"    median < {median_min}: {fail_med_low}")
        print(f"    median > {median_max}: {fail_med_high}")
        print(f"    MAD > {mad_max}: {fail_mad}")
        # Distribution of filtered bin statistics
        if n_filtered > 0:
            filt_medians = medians[~keep_mask]
            filt_mads = mads[~keep_mask]
            print(f"    filtered bins median depth: "
                  f"min={filt_medians.min():.3f}, max={filt_medians.max():.3f}, "
                  f"mean={filt_medians.mean():.3f}")
            print(f"    filtered bins MAD: "
                  f"min={filt_mads.min():.3f}, max={filt_mads.max():.3f}, "
                  f"mean={filt_mads.mean():.3f}")
        # Surviving bin distribution
        surv_medians = medians[keep_mask]
        surv_mads = mads[keep_mask]
        if len(surv_medians) > 0:
            print(f"    surviving bins median depth: "
                  f"min={surv_medians.min():.3f}, max={surv_medians.max():.3f}, "
                  f"mean={surv_medians.mean():.3f}")
            print(f"    surviving bins MAD: "
                  f"min={surv_mads.min():.3f}, max={surv_mads.max():.3f}, "
                  f"mean={surv_mads.mean():.3f}")

    print(f"{'=' * 80}\n")

    return df[keep_mask].copy()


@dataclass
class LocusBinMapping:
    """Tracks mapping between array indices and locus/interval assignments."""
    cluster: str
    locus: Optional[GDLocus]
    interval_name: str
    array_idx: int  # Index in the combined depth tensor
    chrom: str  # Chromosome
    start: int  # Bin start position
    end: int  # Bin end position

