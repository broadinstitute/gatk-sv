"""
Preprocess subcommand — read, normalise, and filter depth data.

Reads a raw depth TSV, normalises each sample by its autosomal median so that
diploid depth ≈ 2.0, filters low-quality bins by median/MAD thresholds, and
writes the preprocessed matrix to disk.
"""

from __future__ import annotations

import argparse
import logging
import os
from typing import List, Tuple

import numpy as np
import pandas as pd
from scipy import stats

from gatk_sv_ploidy._util import (
    DEPTH_SPACES,
    get_sample_columns,
    write_observation_type,
)
from gatk_sv_ploidy.data import read_depth_tsv

logger = logging.getLogger(__name__)


def _compute_autosome_medians(df: pd.DataFrame) -> tuple[list[str], np.ndarray]:
    """Compute per-sample autosomal median depth."""
    sample_cols = get_sample_columns(df)
    autosome_mask = ~df["Chr"].isin(["chrX", "chrY"])
    medians = np.median(df.loc[autosome_mask, sample_cols].values, axis=0)
    return sample_cols, medians


def read_sample_list(path: str) -> list[str]:
    """Read a text file of sample IDs to retain.

    Blank lines and comment lines beginning with ``#`` are ignored. Duplicate
    sample IDs are collapsed while preserving first-seen order.

    Args:
        path: Path to a text file with one sample ID per line.

    Returns:
        Ordered list of sample IDs.

    Raises:
        ValueError: If the file does not contain any sample IDs.
    """
    sample_ids: list[str] = []
    seen: set[str] = set()

    with open(path) as handle:
        for line in handle:
            sample_id = line.strip()
            if not sample_id or sample_id.startswith("#"):
                continue
            if sample_id in seen:
                continue
            seen.add(sample_id)
            sample_ids.append(sample_id)

    if not sample_ids:
        raise ValueError(f"No sample IDs found in samples list: {path}")

    logger.info("Loaded sample subset: %d sample(s)", len(sample_ids))
    return sample_ids


def subset_depth_samples(df: pd.DataFrame, sample_ids: list[str]) -> pd.DataFrame:
    """Subset a depth matrix to a requested set of samples.

    Args:
        df: Depth DataFrame with metadata columns plus sample columns.
        sample_ids: Ordered list of sample IDs to retain.

    Returns:
        Copy of *df* containing metadata columns and only the requested sample
        columns, preserving the requested sample order.

    Raises:
        ValueError: If any requested sample is absent from *df*.
    """
    sample_cols = get_sample_columns(df)
    sample_col_set = set(sample_cols)
    missing = [sample_id for sample_id in sample_ids if sample_id not in sample_col_set]
    if missing:
        preview = ", ".join(missing[:5])
        more = "" if len(missing) <= 5 else f" (+{len(missing) - 5} more)"
        raise ValueError(
            "Requested samples not found in depth matrix: "
            f"{preview}{more}"
        )

    metadata_cols = [col for col in df.columns if col not in sample_col_set]
    df_out = df.loc[:, metadata_cols + sample_ids].copy()
    logger.info(
        "Sample subset: retained %d / %d samples",
        len(sample_ids),
        len(sample_cols),
    )
    return df_out


def read_bed_intervals(path: str) -> pd.DataFrame:
    """Read a BED-like interval file with at least 3 columns.

    Args:
        path: Path to BED file containing ``chrom``, ``start``, ``end``.

    Returns:
        DataFrame with columns ``Chr``, ``Start``, ``End`` sorted by
        chromosome and start coordinate.

    Raises:
        ValueError: If any interval has invalid coordinates.
    """
    df = pd.read_csv(
        path,
        sep="\t",
        comment="#",
        header=None,
        usecols=[0, 1, 2],
        names=["Chr", "Start", "End"],
        dtype={"Chr": str, "Start": np.int64, "End": np.int64},
    )
    df = df.dropna(subset=["Chr", "Start", "End"]).copy()

    bad = df["End"] <= df["Start"]
    if bad.any():
        first_bad = df.loc[bad].iloc[0]
        raise ValueError(
            "Invalid poor-region interval with end <= start: "
            f"{first_bad['Chr']}:{int(first_bad['Start'])}-{int(first_bad['End'])}"
        )

    df = df.sort_values(["Chr", "Start", "End"]).reset_index(drop=True)
    logger.info(
        "Loaded poor regions: %d intervals across %d chromosome(s)",
        len(df),
        df["Chr"].nunique(),
    )
    return df


def _merge_intervals_by_chr(
    intervals: pd.DataFrame,
) -> dict[str, list[tuple[int, int]]]:
    """Merge overlapping/adjacent intervals by chromosome."""
    merged: dict[str, list[tuple[int, int]]] = {}

    for chrom, chrom_df in intervals.groupby("Chr", sort=False):
        merged_intervals: list[tuple[int, int]] = []
        for start, end in chrom_df[["Start", "End"]].itertuples(index=False):
            start_i = int(start)
            end_i = int(end)
            if not merged_intervals or start_i > merged_intervals[-1][1]:
                merged_intervals.append((start_i, end_i))
            else:
                prev_start, prev_end = merged_intervals[-1]
                merged_intervals[-1] = (prev_start, max(prev_end, end_i))
        merged[chrom] = merged_intervals

    return merged


def filter_poor_region_bins(
    df: pd.DataFrame,
    poor_regions: pd.DataFrame,
    *,
    min_poor_region_coverage: float = 0.5,
) -> pd.DataFrame:
    """Filter bins substantially covered by poor genomic regions.

    A bin is removed when the fraction of its length overlapped by the union
    of poor-region intervals on the same chromosome is at least
    *min_poor_region_coverage*.

    Args:
        df: Depth DataFrame with ``Chr``, ``Start``, and ``End`` columns.
        poor_regions: BED-like DataFrame of poor intervals.
        min_poor_region_coverage: Minimum overlap fraction required to
            remove a bin.

    Returns:
        Filtered DataFrame.
    """
    if not 0.0 <= min_poor_region_coverage <= 1.0:
        raise ValueError("min_poor_region_coverage must be between 0 and 1")

    if len(df) == 0 or len(poor_regions) == 0:
        logger.info(
            "Poor-region filter: no intervals or no bins to filter; kept %d bins",
            len(df),
        )
        return df.copy()

    merged_by_chr = _merge_intervals_by_chr(poor_regions)
    chr_values = df["Chr"].to_numpy()
    starts_all = df["Start"].to_numpy(dtype=np.int64)
    ends_all = df["End"].to_numpy(dtype=np.int64)
    overlap_fraction = np.zeros(len(df), dtype=np.float64)

    for chrom, merged_intervals in merged_by_chr.items():
        chrom_pos = np.flatnonzero(chr_values == chrom)
        if len(chrom_pos) == 0:
            continue

        chrom_starts = starts_all[chrom_pos]
        order = np.argsort(chrom_starts, kind="mergesort")
        sorted_pos = chrom_pos[order]

        interval_idx = 0
        n_intervals = len(merged_intervals)

        for row_pos in sorted_pos:
            bin_start = int(starts_all[row_pos])
            bin_end = int(ends_all[row_pos])
            bin_len = max(bin_end - bin_start, 1)

            while interval_idx < n_intervals and merged_intervals[interval_idx][1] <= bin_start:
                interval_idx += 1

            overlap_bp = 0
            scan_idx = interval_idx
            while scan_idx < n_intervals and merged_intervals[scan_idx][0] < bin_end:
                overlap_start = max(bin_start, merged_intervals[scan_idx][0])
                overlap_end = min(bin_end, merged_intervals[scan_idx][1])
                if overlap_end > overlap_start:
                    overlap_bp += overlap_end - overlap_start
                if merged_intervals[scan_idx][1] >= bin_end:
                    break
                scan_idx += 1

            overlap_fraction[row_pos] = overlap_bp / bin_len

    poor_bad = overlap_fraction >= min_poor_region_coverage
    n_flagged = int(poor_bad.sum())
    df_out = df.loc[~poor_bad].copy()
    logger.info(
        "Poor-region filter (coverage ≥ %.0f%%): removed %d / %d bins; %d remain",
        min_poor_region_coverage * 100,
        n_flagged,
        len(df),
        len(df_out),
    )
    return df_out


# ── normalisation ───────────────────────────────────────────────────────────


def normalise_depth(df: pd.DataFrame) -> pd.DataFrame:
    """Normalise per-sample depth so that autosomal CN = 2 corresponds to 2.0.

    For each sample the median depth across autosomal bins is computed and all
    bins (including sex chromosomes) are scaled by ``2 / median``.

    Args:
        df: DataFrame with ``Chr`` metadata column and per-sample depth columns.

    Returns:
        Copy of *df* with normalised depth values.
    """
    sample_cols, medians = _compute_autosome_medians(df)

    logger.info(
        "Autosomal medians: min=%.3f, max=%.3f, mean=%.3f",
        medians.min(),
        medians.max(),
        medians.mean(),
    )

    df = df.copy()
    df[sample_cols] = 2.0 * df[sample_cols].values / medians[np.newaxis, :]
    return df


def clamp_depth_ratio(df: pd.DataFrame, *, depth_ratio_clamp: float = 6.0) -> pd.DataFrame:
    """Clamp per-sample depth to a maximum approximate copy number.

    Each sample's depth is converted to an approximate copy number via
    ``2 * depth / autosomal_median``. Values above *depth_ratio_clamp* are
    capped at that threshold, equivalent to capping the raw depth at
    ``autosomal_median * depth_ratio_clamp / 2``.

    Args:
        df: Depth DataFrame with autosomal bins available for median scaling.
        depth_ratio_clamp: Maximum approximate copy number to retain.

    Returns:
        Copy of *df* with per-sample depth values clamped.

    Raises:
        ValueError: If *depth_ratio_clamp* is not positive.
    """
    if depth_ratio_clamp <= 0:
        raise ValueError("depth_ratio_clamp must be > 0")

    sample_cols, medians = _compute_autosome_medians(df)
    clamp_values = medians * (depth_ratio_clamp / 2.0)

    df_out = df.copy()
    original = df_out[sample_cols].to_numpy(dtype=np.float64)
    clamped = np.minimum(original, clamp_values[np.newaxis, :])
    df_out[sample_cols] = clamped

    n_clamped = int(np.count_nonzero(original > clamped))
    logger.info(
        "Depth-ratio clamp (CN ≤ %.1f): clamped %d sample-bin value(s)",
        depth_ratio_clamp,
        n_clamped,
    )
    return df_out


def _infer_sex_depth_groups(df: pd.DataFrame) -> pd.Series:
    """Infer rough XX/XY groups from normalized chrX/chrY depth.

    The grouping is used only for preprocess-time filtering of sex
    chromosomes. Each sample is assigned to the closer of the normalized
    depth prototypes ``XX=(chrX=2, chrY=0)`` and ``XY=(chrX=1, chrY=1)``
    using whichever sex chromosomes are available in *df*.

    Args:
        df: Normalized depth DataFrame.

    Returns:
        Series indexed by sample ID with values ``"XX"`` or ``"XY"``.
        Returns an empty Series when no sex-chromosome bins are present.
    """
    sample_cols = get_sample_columns(df)
    chr_x_mask = df["Chr"] == "chrX"
    chr_y_mask = df["Chr"] == "chrY"

    if len(sample_cols) == 0 or (not chr_x_mask.any() and not chr_y_mask.any()):
        return pd.Series(dtype=object)

    chr_x_depth = np.full(len(sample_cols), np.nan, dtype=np.float64)
    chr_y_depth = np.full(len(sample_cols), np.nan, dtype=np.float64)

    if chr_x_mask.any():
        chr_x_depth = np.median(
            df.loc[chr_x_mask, sample_cols].to_numpy(dtype=np.float64),
            axis=0,
        )
    if chr_y_mask.any():
        chr_y_depth = np.median(
            df.loc[chr_y_mask, sample_cols].to_numpy(dtype=np.float64),
            axis=0,
        )

    group_by_sample: dict[str, str] = {}
    for sample_idx, sample_id in enumerate(sample_cols):
        xx_distance = 0.0
        xy_distance = 0.0
        used_dimension = False

        if np.isfinite(chr_x_depth[sample_idx]):
            xx_distance += (chr_x_depth[sample_idx] - 2.0) ** 2
            xy_distance += (chr_x_depth[sample_idx] - 1.0) ** 2
            used_dimension = True
        if np.isfinite(chr_y_depth[sample_idx]):
            xx_distance += chr_y_depth[sample_idx] ** 2
            xy_distance += (chr_y_depth[sample_idx] - 1.0) ** 2
            used_dimension = True

        if not used_dimension:
            continue
        group_by_sample[sample_id] = "XX" if xx_distance <= xy_distance else "XY"

    sex_groups = pd.Series(group_by_sample, dtype=object)
    if not sex_groups.empty:
        counts = sex_groups.value_counts().to_dict()
        logger.info(
            "Approximate sex groups for preprocess filtering: %s",
            ", ".join(f"{label}={count}" for label, count in sorted(counts.items())),
        )
    return sex_groups


# ── bin quality filtering ───────────────────────────────────────────────────


def filter_low_quality_bins(
    df: pd.DataFrame,
    *,
    autosome_median_min: float = 1.5,
    autosome_median_max: float = 2.5,
    autosome_mad_max: float = 2.0,
    chrX_median_min: float = 0.0,
    chrX_median_max: float = 3.0,
    chrX_mad_max: float = 2.0,
    chrY_median_min: float = 0.0,
    chrY_median_max: float = 0.85,
    chrY_mad_max: float = 2.0,
    cohort_deviation_threshold: float = 0.3,
    cohort_deviation_fraction_max: float = 0.75,
    min_bins_per_chr: int = 10,
) -> pd.DataFrame:
    """Filter bins whose cross-sample statistics fall outside thresholds.

    Autosomes are filtered using pooled cross-sample median/MAD thresholds.
    Sex chromosomes are filtered separately within rough XX and XY depth
    groups inferred from normalized chrX/chrY depth, because pooled chrX/chrY
    depth is intrinsically bimodal in mixed-sex cohorts.

    In addition to per-bin median/MAD thresholds, a **cohort deviation
    fraction** filter removes bins where more than
    *cohort_deviation_fraction_max* of samples have depth deviating from
    the expected ploidy by more than *cohort_deviation_threshold*.  This
    catches systematic mapping biases that shift all samples uniformly,
    which can be misinterpreted as real CNVs.

    Args:
        df: Normalised depth DataFrame.
        autosome_median_min/max: Median depth range for autosomal bins.
        autosome_mad_max: Maximum MAD for autosomal bins.
        chrX_median_min/max: Legacy pooled chrX median range. Used only when
            sex-stratified filtering cannot be applied.
        chrX_mad_max: Maximum within-group MAD for chrX bins.
        chrY_median_min/max: Legacy pooled chrY median range. Used only when
            sex-stratified filtering cannot be applied.
        chrY_mad_max: Maximum within-group MAD for chrY bins.
        cohort_deviation_threshold: Per-sample depth deviation from the
            expected ploidy that counts as "deviant" (default 0.3). The same
            threshold is used for autosome diploid filtering and for
            sex-stratified chrX/chrY filtering.
        cohort_deviation_fraction_max: Maximum fraction of samples that
            may deviate before the bin is removed (default 0.75).
        min_bins_per_chr: Raise if any chromosome has fewer bins after
            filtering.

    Returns:
        Filtered DataFrame.

    Raises:
        ValueError: If a chromosome has fewer than *min_bins_per_chr* bins.
    """
    sample_cols = get_sample_columns(df)
    depths = df[sample_cols].values
    medians = np.median(depths, axis=1)
    mads = stats.median_abs_deviation(depths, axis=1)

    logger.info("Starting bins: %d", len(df))

    keep = np.ones(len(df), dtype=bool)
    n_threshold_removed = 0

    # -- pooled autosome thresholds --
    _thresholds = {
        "autosome": (
            ~df["Chr"].isin(["chrX", "chrY"]),
            autosome_median_min,
            autosome_median_max,
            autosome_mad_max,
        ),
    }

    for label, (mask, med_min, med_max, mad_max) in _thresholds.items():
        if not mask.any():
            continue
        ok = (medians >= med_min) & (medians <= med_max) & (mads <= mad_max)
        n_before = int(mask.sum())
        keep[mask.values] &= ok[mask.values]
        n_after = int((mask & keep).sum())
        logger.info(
            "%s: median [%.1f, %.1f], MAD ≤ %.1f → %d / %d bins kept",
            label,
            med_min,
            med_max,
            mad_max,
            n_after,
            n_before,
        )

    n_threshold_removed = int((~keep).sum())

    # -- sex-stratified chrX / chrY thresholds --
    sex_groups = _infer_sex_depth_groups(df)
    if not sex_groups.empty:
        group_indices = {
            label: np.array(
                [
                    idx
                    for idx, sample_id in enumerate(sample_cols)
                    if sex_groups.get(sample_id) == label
                ],
                dtype=np.intp,
            )
            for label in ("XX", "XY")
        }
        expected_depths = {
            "chrX": {"XX": 2.0, "XY": 1.0},
            "chrY": {"XX": 0.0, "XY": 1.0},
        }
        mad_limits = {"chrX": chrX_mad_max, "chrY": chrY_mad_max}

        for chrom in ("chrX", "chrY"):
            chrom_mask = (df["Chr"] == chrom).values
            if not chrom_mask.any():
                continue

            for label in ("XX", "XY"):
                sample_idx = group_indices[label]
                if len(sample_idx) == 0:
                    continue

                group_depths = depths[:, sample_idx]
                group_mads = stats.median_abs_deviation(group_depths, axis=1)
                group_bad = group_mads > mad_limits[chrom]

                if cohort_deviation_threshold > 0 and cohort_deviation_fraction_max < 1.0:
                    expected = expected_depths[chrom][label]
                    deviant_fraction = (
                        np.abs(group_depths - expected) > cohort_deviation_threshold
                    ).mean(axis=1)
                    group_bad |= deviant_fraction > cohort_deviation_fraction_max

                group_bad &= chrom_mask
                newly_removed = int((group_bad & keep).sum())
                keep &= ~group_bad
                logger.info(
                    "%s %s: MAD ≤ %.1f and |depth−expected|>%.2f in >%.0f%% of samples → %d / %d bins kept (%d new removed)",
                    chrom,
                    label,
                    mad_limits[chrom],
                    cohort_deviation_threshold,
                    cohort_deviation_fraction_max * 100,
                    int(((df["Chr"] == chrom).values & keep).sum()),
                    int(chrom_mask.sum()),
                    newly_removed,
                )
        n_threshold_removed = int((~keep).sum())
    else:
        # Fallback for datasets lacking chrX / chrY context: retain the
        # legacy pooled sex-chromosome thresholds.
        legacy_thresholds = {
            "chrX": (
                df["Chr"] == "chrX",
                chrX_median_min,
                chrX_median_max,
                chrX_mad_max,
            ),
            "chrY": (
                df["Chr"] == "chrY",
                chrY_median_min,
                chrY_median_max,
                chrY_mad_max,
            ),
        }
        for label, (mask, med_min, med_max, mad_max) in legacy_thresholds.items():
            if not mask.any():
                continue
            ok = (medians >= med_min) & (medians <= med_max) & (mads <= mad_max)
            n_before = int(mask.sum())
            keep[mask.values] &= ok[mask.values]
            n_after = int((mask & keep).sum())
            logger.info(
                "%s legacy pooled filter: median [%.1f, %.1f], MAD ≤ %.1f → %d / %d bins kept",
                label,
                med_min,
                med_max,
                mad_max,
                n_after,
                n_before,
            )
        n_threshold_removed = int((~keep).sum())

    # -- cohort deviation fraction filter --
    # For each bin, compute the expected ploidy (2.0 for autosomes) and
    # reject bins where most samples deviate from it uniformly.
    if cohort_deviation_threshold > 0 and cohort_deviation_fraction_max < 1.0:
        is_auto = ~df["Chr"].isin(["chrX", "chrY"]).values
        # Expected diploid depth for autosomes; sex chromosomes are
        # not tested because the expected depth depends on the unknown
        # sex composition of the cohort.
        auto_deviations = np.abs(depths - 2.0)  # (n_bins, n_samples)
        deviant_frac = (
            auto_deviations > cohort_deviation_threshold
        ).mean(axis=1)
        cohort_bad = is_auto & (deviant_frac > cohort_deviation_fraction_max)
        n_cohort = int(cohort_bad.sum())
        # Only remove bins that weren't already removed
        newly_removed = int((cohort_bad & keep).sum())
        keep &= ~cohort_bad
        logger.info(
            "Cohort deviation filter (|depth−expected|>%.2f in >%.0f%% "
            "of samples): flagged %d bins (%d new)",
            cohort_deviation_threshold,
            cohort_deviation_fraction_max * 100,
            n_cohort,
            newly_removed,
        )

    df_out = df[keep].copy()
    logger.info(
        "Bins after filtering: %d (removed %d: %d threshold, %d cohort)",
        len(df_out),
        len(df) - len(df_out),
        n_threshold_removed,
        len(df) - len(df_out) - n_threshold_removed,
    )

    # Verify sufficient bins per chromosome
    counts = df_out.groupby("Chr").size()
    bad = counts[counts < min_bins_per_chr]
    if len(bad) > 0:
        msg = "; ".join(f"{c}={n}" for c, n in bad.items())
        raise ValueError(
            f"Insufficient bins after filtering ({msg}). "
            "Consider --skip-bin-filter or relaxing thresholds."
        )

    return df_out


def collapse_bins_per_contig(
    df: pd.DataFrame,
    *,
    bins_per_contig: int = 30,
    aggregation: str = "weighted_mean",
) -> pd.DataFrame:
    """Collapse contiguous bins within each contig when a contig is too dense.

    For each contig with more than *bins_per_contig* bins, choose the largest
    chunk size *k* such that collapsing contiguous groups of *k* bins still
    yields at least *bins_per_contig* output bins. Depth values are averaged
    using bin-length weights so larger intervals contribute proportionally
    more to the collapsed depth. When ``aggregation='sum'``, per-sample depth
    values are summed instead; this is used for raw-count output so the
    collapsed bins remain integer count totals.

    Args:
        df: Preprocessed depth DataFrame.
        bins_per_contig: Target minimum number of bins to retain per contig.
            Values <= 0 disable collapsing.
        aggregation: ``'weighted_mean'`` for normalized depth or ``'sum'``
            for raw counts.

    Returns:
        DataFrame with collapsed bins and refreshed bin identifiers.
    """
    if aggregation not in {"weighted_mean", "sum"}:
        raise ValueError("aggregation must be 'weighted_mean' or 'sum'")

    if bins_per_contig <= 0 or len(df) == 0:
        return df.copy()

    sample_cols = get_sample_columns(df)
    collapsed_parts: list[pd.DataFrame] = []
    n_rebinned_contigs = 0
    n_bins_before = len(df)

    for chrom, chrom_df in df.groupby("Chr", sort=False):
        chrom_df = chrom_df.sort_values("Start", kind="mergesort").copy()
        if "BinLengthBp" in chrom_df.columns:
            bin_lengths = chrom_df["BinLengthBp"].to_numpy(dtype=np.int64)
        else:
            bin_lengths = np.maximum(
                chrom_df["End"].to_numpy(dtype=np.int64) -
                chrom_df["Start"].to_numpy(dtype=np.int64),
                1,
            )
            chrom_df["BinLengthBp"] = bin_lengths
        n_bins = len(chrom_df)
        if n_bins <= bins_per_contig:
            collapsed_parts.append(chrom_df)
            continue

        chunk_size = 1
        while ((n_bins + chunk_size) // (chunk_size + 1)) >= bins_per_contig:
            chunk_size += 1

        rows: list[dict] = []
        source_values = (
            chrom_df["source_file"].tolist()
            if "source_file" in chrom_df.columns
            else None
        )
        starts = chrom_df["Start"].to_numpy(dtype=np.int64)
        ends = chrom_df["End"].to_numpy(dtype=np.int64)
        depths = chrom_df[sample_cols].to_numpy(dtype=np.float64)

        for start_idx in range(0, n_bins, chunk_size):
            end_idx = min(start_idx + chunk_size, n_bins)
            group_lengths = bin_lengths[start_idx:end_idx]
            group_depths = depths[start_idx:end_idx]
            if aggregation == "sum":
                aggregated_depth = group_depths.sum(axis=0)
            else:
                aggregated_depth = np.average(
                    group_depths,
                    axis=0,
                    weights=group_lengths,
                )

            row = {
                "Chr": chrom,
                "Start": int(starts[start_idx]),
                "End": int(ends[end_idx - 1]),
                "BinLengthBp": int(group_lengths.sum()),
            }
            if source_values is not None:
                row["source_file"] = source_values[start_idx]
            for sample_col, value in zip(sample_cols, aggregated_depth):
                row[sample_col] = float(value)
            rows.append(row)

        collapsed_df = pd.DataFrame(rows)
        collapsed_df.index = [
            f"{chrom}:{int(row.Start)}-{int(row.End)}"
            for row in collapsed_df.itertuples(index=False)
        ]
        collapsed_parts.append(collapsed_df)
        n_rebinned_contigs += 1
        logger.info(
            "Collapsed %s from %d bins to %d using contiguous groups of %d",
            chrom,
            n_bins,
            len(collapsed_df),
            chunk_size,
        )

    collapsed = pd.concat(collapsed_parts, axis=0)
    logger.info(
        "Contig bin collapsing: %d / %d contigs rebinned; %d → %d total bins",
        n_rebinned_contigs,
        df["Chr"].nunique(),
        n_bins_before,
        len(collapsed),
    )
    return collapsed

# ── allele fraction preprocessing (per-site joint model) ─────────────────────


def read_site_depth_tsv(
    path: str, stride: int = 1, position_stride: int = 0,
) -> pd.DataFrame:
    """Read a GATK CollectSVEvidence site-depth file.

    SD files are headerless, tab-delimited, optionally bgzipped, with columns:

    1. contig (str)
    2. position (int, 0-based)
    3. sample (str)
    4. A depth (int)
    5. C depth (int)
    6. G depth (int)
    7. T depth (int)

    Args:
        path: Path to ``.sd.txt.gz`` or plain SD file.
        stride: Keep every *stride*-th row by file-row index (1 = keep all).
            Ignored when *position_stride* > 0.
        position_stride: Keep positions where ``position % position_stride == 0``.
            This is slower than *stride* (the full file must be decompressed)
            but guarantees that every sample keeps the **same** set of genomic
            positions, which is critical for multi-sample allele-fraction
            analysis.  When > 0, *stride* is ignored.

    Returns:
        DataFrame with columns ``contig``, ``position``, ``sample``,
        ``A``, ``C``, ``G``, ``T``.
    """
    _dtypes = {
        "contig": str,
        "position": np.int64,
        "sample": str,
        "A": np.int32,
        "C": np.int32,
        "G": np.int32,
        "T": np.int32,
    }
    _names = ["contig", "position", "sample", "A", "C", "G", "T"]

    position_stride = max(0, int(position_stride))
    stride = max(1, int(stride))

    if position_stride > 1:
        # Position-based striding: read in chunks and keep only rows whose
        # genomic position is divisible by position_stride.  This ensures
        # cross-sample consistency at the cost of a full file scan.
        logger.info("Reading site depth file (position_stride=%d)", position_stride)
        chunks: list[pd.DataFrame] = []
        reader = pd.read_csv(
            path, sep="\t", compression="infer", header=None,
            names=_names, dtype=_dtypes, chunksize=2_000_000,
        )
        for chunk in reader:
            filtered = chunk[chunk["position"] % position_stride == 0]
            if len(filtered) > 0:
                chunks.append(filtered)
        df = (
            pd.concat(chunks, ignore_index=True)
            if chunks
            else pd.DataFrame(columns=_names)
        )
    else:
        logger.info("Reading site depth file (stride=%d)", stride)
        skiprows = (lambda i: i % stride != 0) if stride > 1 else None
        df = pd.read_csv(
            path, sep="\t", compression="infer", header=None,
            names=_names, dtype=_dtypes, skiprows=skiprows,
        )

    logger.info(
        "  %d sites for %d sample(s)",
        len(df),
        df["sample"].nunique(),
    )
    return df


def _infer_cohort_major_base(
    a: np.ndarray,
    c: np.ndarray,
    g: np.ndarray,
    t: np.ndarray,
) -> int:
    """Return the cohort-major base index for a site.

    The returned index follows ``A,C,G,T`` ordering.
    """
    pooled = np.array(
        [
            np.sum(a, dtype=np.int64),
            np.sum(c, dtype=np.int64),
            np.sum(g, dtype=np.int64),
            np.sum(t, dtype=np.int64),
        ],
        dtype=np.int64,
    )
    return int(np.argmax(pooled))


def _site_nonmajor_total(
    a: np.ndarray,
    c: np.ndarray,
    g: np.ndarray,
    t: np.ndarray,
    major_base_idx: int,
) -> Tuple[np.ndarray, np.ndarray]:
    """Compute non-major allele count and total depth from per-base counts.

    The major allele is defined once per site from the pooled cohort counts.
    Each sample's alt count is then all reads not matching that shared allele.
    """
    counts = np.stack([a, c, g, t], axis=-1)
    total = counts.sum(axis=-1)
    major = counts[:, major_base_idx]
    alt = total - major
    return alt, total


def _estimate_site_pop_af(alt: np.ndarray, total: np.ndarray) -> float:
    """Estimate cohort AF from pooled alt and total counts with smoothing."""
    pooled_alt = float(np.sum(alt, dtype=np.float64))
    pooled_total = float(np.sum(total, dtype=np.float64))
    if pooled_total <= 0.0:
        return 0.5
    # Jeffreys smoothing avoids brittle 0/1 AFs at low coverage.
    return float((pooled_alt + 0.5) / (pooled_total + 1.0))


def _process_sd_file(
    sd_path: str,
    sample_to_idx: dict[str, int],
    bin_lookup: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]],
    *,
    stride: int,
    min_site_depth: int,
) -> tuple[str, int, dict[int, dict[int, dict[str, object]]]]:
    """Process one SD file into sparse per-bin site records."""
    sd_df = read_site_depth_tsv(sd_path, stride=stride)

    file_site_data: dict[int, dict[int, dict[str, object]]] = {}
    total_sites_used = 0

    for sample_id in sd_df["sample"].unique():
        if sample_id not in sample_to_idx:
            logger.warning(
                "Encountered site-depth rows for a sample absent from the depth matrix; skipping those rows.",
            )
            continue

        si = sample_to_idx[sample_id]
        smask = sd_df["sample"].values == sample_id
        s_contigs = sd_df["contig"].values[smask]
        s_positions = sd_df["position"].values[smask].astype(np.int64)
        s_a = sd_df["A"].values[smask]
        s_c = sd_df["C"].values[smask]
        s_g = sd_df["G"].values[smask]
        s_t = sd_df["T"].values[smask]

        for chrom, (chrom_bin_idx, starts, ends) in bin_lookup.items():
            cmask = s_contigs == chrom
            if not cmask.any():
                continue

            c_pos = s_positions[cmask]
            c_a = s_a[cmask]
            c_c = s_c[cmask]
            c_g = s_g[cmask]
            c_t = s_t[cmask]

            bi = np.searchsorted(starts, c_pos, side="right") - 1
            bi_safe = np.clip(bi, 0, len(starts) - 1)
            in_bin = ((bi >= 0) & (bi < len(starts)) & (c_pos < ends[bi_safe]))

            # Drop sites that don't fall inside any surviving bin.
            if not in_bin.all():
                n_dropped = int((~in_bin).sum())
                logger.debug(
                    "  site-depth filtering: dropping %d sites outside bin boundaries on %s",
                    n_dropped,
                    chrom,
                )
                c_pos = c_pos[in_bin]
                c_a = c_a[in_bin]
                c_c = c_c[in_bin]
                c_g = c_g[in_bin]
                c_t = c_t[in_bin]
                bi_safe = bi_safe[in_bin]

            g_bi_arr = chrom_bin_idx[bi_safe]
            if len(c_pos) == 0:
                continue

            totals = (c_a + c_c + c_g + c_t).astype(np.int64)

            depth_ok = totals >= min_site_depth
            c_pos = c_pos[depth_ok]
            c_a, c_c, c_g, c_t = (
                c_a[depth_ok], c_c[depth_ok], c_g[depth_ok], c_t[depth_ok],
            )
            g_bi_arr = g_bi_arr[depth_ok]
            totals = totals[depth_ok]
            if len(c_pos) == 0:
                continue

            if len(c_pos) == 0:
                continue

            n_sites_chrom = len(c_pos)
            total_sites_used += n_sites_chrom
            for i in range(n_sites_chrom):
                pos = int(c_pos[i])
                g_bi = int(g_bi_arr[i])
                bin_entry = file_site_data.setdefault(g_bi, {})
                site_entry = bin_entry.setdefault(pos, {"values": []})
                site_entry["values"].append(
                    (si, int(c_a[i]), int(c_c[i]), int(c_g[i]), int(c_t[i]))
                )

            if n_sites_chrom > 0:
                logger.debug(
                    "  site-depth aggregation: %d sites accumulated on %s",
                    n_sites_chrom,
                    chrom,
                )

    return sd_path, total_sites_used, file_site_data


def build_per_site_data(
    sd_paths: List[str],
    bins_df: pd.DataFrame,
    *,
    min_site_depth: int = 10,
    max_sites_per_bin: int = 50,
    sd_stride: int = 1,
    seed: int = 42,
) -> dict:
    """Build per-site allele-count arrays for the joint genotype/CN model.

    Every retained SD position inside a genomic bin is assigned a stable,
    cohort-derived allele identity. For each site, the pooled cohort counts
    define a single major allele, and each sample's alt count is all reads not
    matching that allele. The site population AF is estimated directly from the
    pooled alt and total counts.

    Args:
        sd_paths: Paths to per-sample ``.sd.txt.gz`` files.
        bins_df: Preprocessed depth DataFrame with ``Chr``, ``Start``, ``End``
            metadata columns.
        min_site_depth: Minimum total depth at a site to retain a sample-site
            observation during cohort AF inference.
        max_sites_per_bin: Pad/subsample to this many sites per bin.
        sd_stride: Keep every *sd_stride*-th row when reading SD files
            (1 = keep all, 100 = 100× downsample).
        seed: Random seed for reproducible subsampling.

    Returns:
        Dictionary with the following NumPy arrays (ready for
        :func:`numpy.savez`):

        - ``site_alt``: ``int32 (n_bins, max_sites, n_samples)`` alt counts
        - ``site_total``: ``int32 (n_bins, max_sites, n_samples)`` total counts
        - ``site_pop_af``: ``float32 (n_bins, max_sites)`` population AFs
        - ``site_mask``: ``bool (n_bins, max_sites, n_samples)`` observed entries
        - ``sample_ids``: ``str (n_samples,)``
        - ``bin_chr``: ``str (n_bins,)``
        - ``bin_start``: ``int64 (n_bins,)``
        - ``bin_end``: ``int64 (n_bins,)``
    """
    sample_cols = get_sample_columns(bins_df)
    n_bins = len(bins_df)
    n_samples = len(sample_cols)
    sample_to_idx = {s: i for i, s in enumerate(sample_cols)}

    # Build bin-lookup per chromosome
    bin_lookup: dict = {}
    for chrom in bins_df["Chr"].unique():
        chrom_mask = (bins_df["Chr"] == chrom).values
        chrom_indices = np.where(chrom_mask)[0]
        starts = bins_df["Start"].values[chrom_mask]
        ends = bins_df["End"].values[chrom_mask]
        bin_lookup[chrom] = (chrom_indices, starts, ends)

    # Intermediate: collect per-bin lists of (site_pos, per-sample A/C/G/T counts)
    # For each bin we accumulate site-level data across all SD files
    # site_data[global_bin_idx] = {pos: {"a": [n_samples], ..., "observed": [n_samples]}}
    site_data: List[dict] = [{} for _ in range(n_bins)]

    total_sites_used = 0
    logger.info("Loading %d SD file(s)", len(sd_paths))
    worker_results = [
        _process_sd_file(
            sd_path,
            sample_to_idx,
            bin_lookup,
            stride=sd_stride,
            min_site_depth=min_site_depth,
        )
        for sd_path in sd_paths
    ]

    for sd_path, file_total_sites, file_site_data in worker_results:
        total_sites_used += file_total_sites
        logger.info(
            "Completed site depth ingestion: %d retained site-sample entries",
            file_total_sites,
        )
        for g_bi, bin_sites in file_site_data.items():
            entry = site_data[g_bi]
            for pos, site_record in bin_sites.items():
                if pos not in entry:
                    entry[pos] = {
                        "a": np.zeros(n_samples, dtype=np.int32),
                        "c": np.zeros(n_samples, dtype=np.int32),
                        "g": np.zeros(n_samples, dtype=np.int32),
                        "t": np.zeros(n_samples, dtype=np.int32),
                        "observed": np.zeros(n_samples, dtype=bool),
                    }
                for si, a, c, g, t in site_record["values"]:
                    entry[pos]["a"][si] = a
                    entry[pos]["c"][si] = c
                    entry[pos]["g"][si] = g
                    entry[pos]["t"][si] = t
                    entry[pos]["observed"][si] = True

    # ── pack into padded arrays ──────────────────────────────────────────
    site_alt = np.zeros((n_bins, max_sites_per_bin, n_samples), dtype=np.int32)
    site_total = np.zeros((n_bins, max_sites_per_bin, n_samples), dtype=np.int32)
    site_pop_af = np.zeros((n_bins, max_sites_per_bin), dtype=np.float32)
    site_mask = np.zeros((n_bins, max_sites_per_bin, n_samples), dtype=bool)

    bins_with_data = 0
    for g_bi in range(n_bins):
        entries = site_data[g_bi]
        if not entries:
            continue
        bins_with_data += 1

        positions = sorted(entries.keys())
        if len(positions) > max_sites_per_bin:
            # Prioritise positions with data for the most samples.
            coverage = np.array(
                [int(np.sum(entries[p]["observed"])) for p in positions],
            )
            # Stable argsort descending; among ties, keep original order
            order = np.argsort(-coverage, kind="mergesort")
            positions = sorted([positions[i] for i in order[:max_sites_per_bin]])

        for slot, pos in enumerate(positions):
            e = entries[pos]
            major_base_idx = _infer_cohort_major_base(
                e["a"], e["c"], e["g"], e["t"],
            )
            alt, total = _site_nonmajor_total(
                e["a"], e["c"], e["g"], e["t"], major_base_idx,
            )
            site_alt[g_bi, slot, :] = alt.astype(np.int32)
            site_total[g_bi, slot, :] = total.astype(np.int32)
            site_pop_af[g_bi, slot] = _estimate_site_pop_af(alt, total)
            site_mask[g_bi, slot, :] = e["observed"]

    logger.info(
        "Per-site data: %d total site-sample entries, %d / %d bins with data, "
        "padded to %d sites/bin",
        total_sites_used, bins_with_data, n_bins, max_sites_per_bin,
    )

    # Per-chromosome summary for diagnostics
    chr_labels = bins_df["Chr"].values
    for chrom in sorted(set(chr_labels)):
        chrom_mask_arr = chr_labels == chrom
        n_chrom_bins = int(chrom_mask_arr.sum())
        chrom_sites = site_mask[chrom_mask_arr].any(axis=2).sum()
        n_unique = sum(len(site_data[i]) for i in np.where(chrom_mask_arr)[0])
        logger.info(
            "  %s: %d bins, %d unique positions, %d retained (max %d/bin)",
            chrom, n_chrom_bins, n_unique, int(chrom_sites),
            max_sites_per_bin,
        )

    return {
        "site_alt": site_alt,
        "site_total": site_total,
        "site_pop_af": site_pop_af,
        "site_mask": site_mask,
        "sample_ids": np.array(sample_cols, dtype=object),
        "bin_chr": bins_df["Chr"].values.astype(object),
        "bin_start": bins_df["Start"].values.astype(np.int64),
        "bin_end": bins_df["End"].values.astype(np.int64),
    }

# ── CLI ─────────────────────────────────────────────────────────────────────


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the preprocess subcommand."""
    p = argparse.ArgumentParser(
        description="Preprocess depth data for aneuploidy inference",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "-i", "--input", required=True,
        help="Input TSV file with raw read-depth (bins × samples)",
    )
    p.add_argument(
        "-o", "--output-dir", required=True,
        help="Output directory",
    )
    p.add_argument(
        "--samples-list", default=None,
        help="Text file listing sample IDs to retain (one per line)",
    )
    p.add_argument(
        "--viable-only", action="store_true", default=False,
        help="Subset to chr13, chr18, chr21, chrX, chrY only",
    )
    p.add_argument(
        "--skip-bin-filter", action="store_true", default=False,
        help="Skip bin quality filtering",
    )
    p.add_argument(
        "--poor-regions", default=None,
        help="BED file of poor regions (for example segmental duplications); "
             "bins with sufficient overlap are removed",
    )
    p.add_argument(
        "--min-poor-region-coverage", type=float, default=0.5,
        help="Minimum fraction of bin length overlapped by poor regions to "
             "remove the bin",
    )
    p.add_argument(
        "--bins-per-contig", type=int, default=30,
        help="When a contig has more than this many bins, collapse contiguous bins into larger bins while retaining at least this many bins per contig",
    )
    p.add_argument(
        "--output-space", choices=list(DEPTH_SPACES), default="raw",
        help="Write filtered normalized depth (historical behavior) or filtered raw counts. Bin-quality filters always run on normalized depth.",
    )
    p.add_argument(
        "--depth-ratio-clamp", type=float, default=6.0,
        help="Clamp each sample-bin depth to at most this approximate copy number, where approximate CN is computed as 2 * depth / autosomal_median",
    )

    # Filter thresholds
    g = p.add_argument_group("bin-filter thresholds")
    g.add_argument("--autosome-median-min", type=float, default=1.5)
    g.add_argument("--autosome-median-max", type=float, default=2.5)
    g.add_argument("--autosome-mad-max", type=float, default=2.0)
    g.add_argument("--chrX-median-min", type=float, default=0.0)
    g.add_argument("--chrX-median-max", type=float, default=3.0)
    g.add_argument("--chrX-mad-max", type=float, default=2.0)
    g.add_argument("--chrY-median-min", type=float, default=0.0)
    g.add_argument("--chrY-median-max", type=float, default=0.85)
    g.add_argument("--chrY-mad-max", type=float, default=2.0)
    g.add_argument(
        "--cohort-deviation-threshold", type=float, default=0.3,
        help="Per-sample depth deviation from expected ploidy that "
             "counts as deviant (default 0.3)",
    )
    g.add_argument(
        "--cohort-deviation-fraction-max", type=float, default=0.75,
        help="Max fraction of samples that may deviate before the bin "
             "is removed (default 0.75)",
    )
    g.add_argument(
        "--min-bins-per-chr", type=int, default=10,
        help="Minimum number of bins required per chromosome after filtering",
    )

    # Allele fraction arguments
    a = p.add_argument_group("allele fraction (site depth)")
    a.add_argument(
        "--site-depth-list", default=None,
        help="Text file listing paths to per-sample SD files (one per line)",
    )
    a.add_argument(
        "--site-depth", nargs="*", default=None,
        help="One or more per-sample SD file paths directly on the command line",
    )
    a.add_argument("--min-site-depth", type=int, default=3,
                   help="Minimum total depth at a site to retain a sample-site observation during cohort AF inference")
    a.add_argument("--max-sites-per-bin", type=int, default=50,
                   help="Maximum sites per bin (pad/subsample to this)")
    a.add_argument("--sd-stride", type=int, default=10,
                   help="Keep every Nth row from SD files (100 = 100x downsample)")
    return p.parse_args()


def main() -> None:
    """Entry point for ``gatk-sv-ploidy preprocess``."""
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    # 1. Read
    raw_df = read_depth_tsv(args.input)
    n_input_bins = len(raw_df)

    if args.samples_list:
        sample_ids = read_sample_list(args.samples_list)
        raw_df = subset_depth_samples(raw_df, sample_ids)

    poor_regions_df = None
    if args.poor_regions:
        poor_regions_df = read_bed_intervals(args.poor_regions)

    # 2. Optionally subset to viable trisomy chromosomes
    if args.viable_only:
        viable = {"chr13", "chr18", "chr21", "chrX", "chrY"}
        n_before = len(raw_df)
        raw_df = raw_df[raw_df["Chr"].isin(viable)].copy()
        logger.info(
            "Viable-only filter: %d → %d bins (%s)",
            n_before,
            len(raw_df),
            sorted(raw_df["Chr"].unique()),
        )

    raw_df = clamp_depth_ratio(raw_df, depth_ratio_clamp=args.depth_ratio_clamp)

    # 3. Normalise for filter evaluation.
    normalized_df = normalise_depth(raw_df)

    # 4. Filter
    if args.skip_bin_filter:
        logger.info("Skipping bin quality filtering (--skip-bin-filter)")
    else:
        normalized_df = filter_low_quality_bins(
            normalized_df,
            autosome_median_min=args.autosome_median_min,
            autosome_median_max=args.autosome_median_max,
            autosome_mad_max=args.autosome_mad_max,
            chrX_median_min=args.chrX_median_min,
            chrX_median_max=args.chrX_median_max,
            chrX_mad_max=args.chrX_mad_max,
            chrY_median_min=args.chrY_median_min,
            chrY_median_max=args.chrY_median_max,
            chrY_mad_max=args.chrY_mad_max,
            cohort_deviation_threshold=args.cohort_deviation_threshold,
            cohort_deviation_fraction_max=args.cohort_deviation_fraction_max,
            min_bins_per_chr=args.min_bins_per_chr,
        )

    raw_df = raw_df.loc[normalized_df.index].copy()

    if poor_regions_df is not None:
        normalized_df = filter_poor_region_bins(
            normalized_df,
            poor_regions_df,
            min_poor_region_coverage=args.min_poor_region_coverage,
        )
        raw_df = raw_df.loc[normalized_df.index].copy()

    normalized_df = collapse_bins_per_contig(
        normalized_df,
        bins_per_contig=args.bins_per_contig,
        aggregation="weighted_mean",
    )
    raw_df = collapse_bins_per_contig(
        raw_df,
        bins_per_contig=args.bins_per_contig,
        aggregation="sum",
    )

    df = normalized_df if args.output_space == "normalized" else raw_df

    logger.info(
        "Total bins retained after preprocess filters: %d / %d",
        len(df),
        n_input_bins,
    )
    if args.output_space == "raw":
        logger.info(
            "Writing filtered raw counts (--output-space raw); all quality filters were evaluated on normalized depth."
        )

    # 5. Write preprocessed depth
    out_path = os.path.join(args.output_dir, "preprocessed_depth.tsv")
    df.to_csv(out_path, sep="\t")
    logger.info("Preprocessed depth written.")
    write_observation_type(args.output_dir, args.output_space)
    logger.info("Observation type written.")

    # 6. Build per-site allele data from site-depth files (optional)
    sd_paths: List[str] = []
    if args.site_depth_list:
        with open(args.site_depth_list) as fh:
            sd_paths.extend(line.strip() for line in fh if line.strip())
    if args.site_depth:
        sd_paths.extend(args.site_depth)

    if sd_paths:
        logger.info("Building per-site allele data from %d SD file(s) …", len(sd_paths))
        site_arrays = build_per_site_data(
            sd_paths,
            df,
            min_site_depth=args.min_site_depth,
            max_sites_per_bin=args.max_sites_per_bin,
            sd_stride=args.sd_stride,
        )
        site_data_path = os.path.join(args.output_dir, "site_data.npz")
        np.savez_compressed(site_data_path, **site_arrays)
        logger.info("Per-site allele data written.")
    else:
        logger.info("No site-depth files provided; skipping allele data")


if __name__ == "__main__":
    main()
