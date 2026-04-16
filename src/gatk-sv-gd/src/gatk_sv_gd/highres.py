"""
High-resolution read-count helpers.

Querying tabix-indexed high-resolution count files and normalising them
to the same CN-2 ≈ 2.0 scale used for the low-resolution bins.
"""

from typing import List, Optional

import numpy as np
import pandas as pd
import pysam

from gatk_sv_gd import _util


def query_highres_bins(
    highres_path: str,
    chrom: str,
    start: int,
    end: int,
    sample_cols: List[str],
    max_bins: Optional[int] = None,
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
        max_bins: Optional maximum number of bins to return. When set,
            consecutive raw bins are aggregated into at most this many
            wider bins by summing counts prior to normalisation.

    Returns:
        DataFrame with columns ``Chr``, ``Start``, ``End``, ``source_file``
        and one column per sample, indexed by a ``Bin`` string.
    """
    tbx = pysam.TabixFile(highres_path)
    try:
        header_line = tbx.header[-1] if tbx.header else ""
        header_cols = header_line.lstrip("#").split("\t")
        if "#Chr" in header_cols:
            header_cols = ["Chr" if col == "#Chr" else col for col in header_cols]

        # Keep only sample columns that exist in both files
        common_samples = [s for s in sample_cols if s in header_cols]
        if len(common_samples) == 0:
            raise ValueError(
                "High-resolution counts file shares no sample columns with the "
                "low-resolution file. Check that both files were generated from "
                "the same sample set."
            )

        if max_bins is not None:
            start_idx = header_cols.index("Start")
            end_idx = header_cols.index("End")
            sample_indices = [header_cols.index(sample) for sample in common_samples]
            target_width = max(int(np.ceil((end - start) / max_bins)), 1)

            group_starts: List[int] = []
            group_ends: List[int] = []
            group_sums: List[np.ndarray] = []

            current_group = None
            current_start = None
            current_end = None
            current_sum = None

            for line in tbx.fetch(chrom, max(0, start), end):
                fields = line.split("\t")
                row_start = int(fields[start_idx])
                row_end = int(fields[end_idx])
                group_id = int((row_start - start) // target_width)
                row_values = np.fromiter(
                    (float(fields[idx]) for idx in sample_indices),
                    dtype=float,
                    count=len(sample_indices),
                )

                if current_group != group_id:
                    if current_group is not None:
                        group_starts.append(current_start)
                        group_ends.append(current_end)
                        group_sums.append(current_sum)
                    current_group = group_id
                    current_start = row_start
                    current_end = row_end
                    current_sum = row_values
                else:
                    current_end = row_end
                    current_sum += row_values

            if current_group is not None:
                group_starts.append(current_start)
                group_ends.append(current_end)
                group_sums.append(current_sum)

            if not group_starts:
                cols = ["Chr", "Start", "End", "source_file"] + list(sample_cols)
                empty = pd.DataFrame(columns=cols)
                empty["Bin"] = pd.Series(dtype=str)
                return empty.set_index("Bin")

            coord_df = pd.DataFrame({
                "Chr": [chrom] * len(group_starts),
                "Start": group_starts,
                "End": group_ends,
                "source_file": "highres",
            })
            sample_df = pd.DataFrame(group_sums, columns=common_samples).reindex(
                columns=list(sample_cols)
            )
            result_df = pd.concat([coord_df, sample_df], axis=1)
            result_df.index = (
                result_df["Chr"].astype(str) + ":" +
                result_df["Start"].astype(str) + "-" +
                result_df["End"].astype(str)
            )
            result_df.index.name = "Bin"
            return result_df

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

    # Keep only sample columns that exist in both files
    common_samples = [s for s in sample_cols if s in raw_df.columns]

    if max_bins is not None and len(raw_df) > max_bins:
        target_width = max(int(np.ceil((end - start) / max_bins)), 1)
        raw_df = raw_df.copy()
        raw_df["Start"] = raw_df["Start"].astype(int)
        raw_df["End"] = raw_df["End"].astype(int)
        raw_df["_group"] = ((raw_df["Start"] - start) // target_width).astype(int)
        raw_df[common_samples] = raw_df[common_samples].apply(pd.to_numeric, errors="coerce")
        raw_df = (
            raw_df
            .groupby("_group", sort=True)
            .agg(
                Chr=("Chr", "first"),
                Start=("Start", "min"),
                End=("End", "max"),
                **{sample: (sample, "sum") for sample in common_samples},
            )
            .reset_index(drop=True)
        )

    # Build the output DataFrame in one concat to avoid repeated column
    # insertions that fragment the internal block structure.
    coord_df = pd.DataFrame({
        "Chr": raw_df["Chr"],
        "Start": raw_df["Start"].astype(int),
        "End": raw_df["End"].astype(int),
        "source_file": "highres",
    })

    # Convert all sample columns to numeric at once, then reindex to the full
    # sample list (fills any absent samples with NaN in a single allocation).
    sample_df = (
        raw_df[common_samples]
        .apply(pd.to_numeric, errors="coerce")
        .reindex(columns=list(sample_cols))
    )

    result_df = pd.concat([coord_df, sample_df], axis=1)

    result_df.index = (
        result_df["Chr"].astype(str) + ":" +
        result_df["Start"].astype(str) + "-" +
        result_df["End"].astype(str)
    )
    result_df.index.name = "Bin"

    return result_df


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

    return df.copy()
