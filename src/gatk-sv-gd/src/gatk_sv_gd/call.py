"""
GD CNV Calling from Model Posteriors

Call genomic disorder (GD) copy-number variants from posterior probabilities
produced by gd_cnv_pyro.py using the Viterbi segmentation strategy.

The Viterbi algorithm is run on per-bin CN posteriors using a user-supplied
state transition matrix to produce a smooth copy-number segmentation.  The
resulting path is partitioned into contiguous regions of altered copy number
(DEL: CN < ploidy, DUP: CN > ploidy), and each such segment is tested for
reciprocal overlap with every GD entry of matching svtype.  A carrier call
is made when the reciprocal overlap meets or exceeds
``--reciprocal-overlap-threshold`` (default 90 %).

Usage:
    python gd_cnv_call.py \\
        --cn-posteriors cn_posteriors.tsv.gz \\
        --bin-mappings bin_mappings.tsv.gz \\
        --gd-table gd_table.tsv \\
        --transition-matrix transition_probs.tsv \\
        --output-dir results/
"""

import argparse
import os
import re
import sys
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from gatk_sv_gd import _util
from gatk_sv_gd.models import GDTable
from gatk_sv_gd.viterbi import (
    load_transition_matrix,
    viterbi_call_gd_cnv,
)


# =============================================================================
# Logging
# =============================================================================


def setup_logging(output_dir: str, filename: str = "call_log.txt"):
    """Redirect stdout and stderr to both the console and a log file.

    Also sets :data:`_util._log_fh` so that :func:`_util.vlog` can write
    verbose diagnostics to the log file without flooding stdout.
    """
    log_path = os.path.join(output_dir, filename)
    log_fh = open(log_path, "w")
    _util._log_fh = log_fh
    sys.stdout = _util.TeeStream(sys.__stdout__, log_fh)
    sys.stderr = _util.TeeStream(sys.__stderr__, log_fh)
    return log_fh


# =============================================================================
# CNV Calling Functions
# =============================================================================


def get_locus_interval_bins(
    bin_mappings_df: pd.DataFrame,
    cluster: str,
) -> Dict[str, List[int]]:
    """
    Get bin indices for each interval in a specific locus.

    Args:
        bin_mappings_df: DataFrame with bin-to-interval mappings
        cluster: Cluster name to filter by

    Returns:
        Dict mapping interval name to list of array_idx values
            (positional indices into the per-sample posterior arrays)
    """
    locus_bins = bin_mappings_df[bin_mappings_df["cluster"] == cluster]
    interval_bins = {}
    for interval_name, group in locus_bins.groupby("interval"):
        interval_bins[interval_name] = group["array_idx"].tolist()
    return interval_bins


def determine_best_breakpoints(
    calls: List[dict],
) -> Dict[str, Optional[str]]:
    """Pick the best GD_ID per svtype among carrier calls.

    When multiple GD entries of the same svtype are called as carriers,
    the one with the highest reciprocal overlap is selected.  Ties are
    broken by genomic span (larger wins).

    Args:
        calls: List of call dicts from ``viterbi_call_gd_cnv()``.

    Returns:
        Dict mapping svtype ("DEL", "DUP") to best matching GD_ID or None.
    """
    best_by_svtype: Dict[str, Optional[str]] = {}
    for svtype in ["DEL", "DUP"]:
        carrier_calls = [
            c for c in calls if c["is_carrier"] and c["svtype"] == svtype
        ]
        if not carrier_calls:
            best_by_svtype[svtype] = None
        elif len(carrier_calls) == 1:
            best_by_svtype[svtype] = carrier_calls[0]["GD_ID"]
        else:
            best = max(
                carrier_calls,
                key=lambda c: (
                    c.get("reciprocal_overlap", 0),
                    c["end"] - c["start"],
                ),
            )
            best_by_svtype[svtype] = best["GD_ID"]
    return best_by_svtype


def get_pair_state_columns(
    cn_posteriors_df: pd.DataFrame,
) -> Tuple[List[str], List[Tuple[int, int]]]:
    """Return pair-state posterior columns and their parsed labels."""
    pair_cols = [
        column for column in cn_posteriors_df.columns
        if column.startswith("prob_pair_")
    ]
    if not pair_cols:
        raise ValueError(
            "cn_posteriors.tsv.gz is missing pair-state posterior columns "
            "(expected columns like prob_pair_0_1). Re-run infer first."
        )

    pair_labels: List[Tuple[int, int]] = []
    for column in pair_cols:
        match = re.fullmatch(r"prob_pair_(\d+)_(\d+)", column)
        if match is None:
            raise ValueError(f"Unrecognized pair-state column name: {column}")
        pair_labels.append((int(match.group(1)), int(match.group(2))))

    return pair_cols, pair_labels


# =============================================================================
# Orchestrator
# =============================================================================


def call_cnvs_from_posteriors(
    cn_posteriors_df: pd.DataFrame,
    bin_mappings_df: pd.DataFrame,
    gd_table: GDTable,
    transition_matrix: np.ndarray,
    ploidy_df: Optional[pd.DataFrame] = None,
    verbose: bool = False,
    reciprocal_overlap_threshold: float = 0.90,
    breakpoint_transition_matrix: Optional[np.ndarray] = None,
) -> pd.DataFrame:
    """
    Call CNVs from posterior probabilities using the Viterbi strategy.

    The Viterbi algorithm is run once per sample per locus to produce a
    smooth copy-number segmentation.  The resulting path is partitioned
    into contiguous DEL / DUP segments and matched against GD entries by
    reciprocal overlap.

    Args:
        cn_posteriors_df: DataFrame with CN posterior probabilities per bin/sample
        bin_mappings_df: DataFrame with bin-to-interval mappings
        gd_table: GDTable with locus definitions
        transition_matrix: (n_states, n_states) CN-state transition probability
            matrix.
        ploidy_df: Optional DataFrame with columns (sample, contig, ploidy).
            If None, ploidy=2 is assumed for all sample/contig pairs.
        verbose: If True, print per-sample log probability scores for every
            GD entry at every locus.
        reciprocal_overlap_threshold: Minimum reciprocal overlap between a
            Viterbi segment and a GD entry to call a carrier.  Default 0.90.
        breakpoint_transition_matrix: Optional (n_states, n_states) transition
            matrix applied at known recurrent breakpoint boundaries during
            Viterbi.  Should have lower diagonal values than *transition_matrix*
            to encourage CN state changes at those positions.

    Returns:
        Tuple of:
          - DataFrame with CNV calls for all samples and loci
          - DataFrame with Viterbi category segments (sample, cluster,
            start, end, mean_cn, category)
    """
    print("\n" + "=" * 80)
    print("CALLING CNVs FROM POSTERIORS")
    bp_str = " + breakpoint matrix" if breakpoint_transition_matrix is not None else ""
    print(f"  NAHR strategy: Viterbi segmentation{bp_str}  "
          f"(reciprocal overlap threshold={reciprocal_overlap_threshold:.0%})")
    print("=" * 80)

    all_results = []
    all_path_records = []
    sample_ids = cn_posteriors_df["sample"].unique()

    # Build ploidy lookup: (sample, contig) -> ploidy, default 2
    ploidy_lookup: Dict[Tuple[str, str], int] = {}
    if ploidy_df is not None:
        for _, row in ploidy_df.iterrows():
            ploidy_lookup[(str(row["sample"]), str(row["contig"]))] = int(row["ploidy"])
        print(f"  Loaded ploidy for {len(ploidy_lookup)} sample/contig pairs")
    else:
        print("  No ploidy table provided; assuming diploid (ploidy=2) everywhere")

    # ---- Validate alignment between cn_posteriors and bin_mappings ----
    n_bins = len(bin_mappings_df)
    n_samples = len(sample_ids)
    expected_rows = n_bins * n_samples
    if len(cn_posteriors_df) != expected_rows:
        print(f"  WARNING: cn_posteriors has {len(cn_posteriors_df)} rows, "
              f"expected {expected_rows} ({n_bins} bins × {n_samples} samples)")

    first_sample = sample_ids[0]
    first_sample_rows = cn_posteriors_df[cn_posteriors_df["sample"] == first_sample]
    if len(first_sample_rows) == n_bins:
        post_coords = list(zip(
            first_sample_rows["chr"].values,
            first_sample_rows["start"].values,
            first_sample_rows["end"].values,
        ))
        map_coords = list(zip(
            bin_mappings_df["chr"].values,
            bin_mappings_df["start"].values,
            bin_mappings_df["end"].values,
        ))
        if post_coords != map_coords:
            raise ValueError(
                "Bin coordinates in cn_posteriors do not match bin_mappings. "
                "The files may have been generated by different runs or with "
                "a version that had a bin-ordering bug. Please re-run "
                "gd_cnv_pyro.py to regenerate both files."
            )
        print(f"  Validated: bin coordinates match between posteriors and mappings ({n_bins} bins)")
    else:
        print(f"  WARNING: sample {first_sample} has {len(first_sample_rows)} bins, "
              f"expected {n_bins}; skipping alignment check")

    # Pre-extract numpy arrays for fast access.
    print("  Organizing data for fast access...")
    pair_prob_cols, pair_state_labels = get_pair_state_columns(cn_posteriors_df)
    n_pair_states = len(pair_prob_cols)
    pair_prob_3d = np.empty((n_samples, n_bins, n_pair_states))
    depth_2d = np.empty((n_samples, n_bins))
    for s_idx, sample_id in enumerate(sample_ids):
        mask = cn_posteriors_df["sample"] == sample_id
        sample_rows = cn_posteriors_df.loc[mask]
        pair_prob_3d[s_idx] = sample_rows[pair_prob_cols].values
        depth_2d[s_idx] = sample_rows["depth"].values
    print(f"    Extracted {n_samples} x {n_bins} x {n_pair_states} pair-state probability array")

    # Build a per-bin coordinate lookup used by viterbi_call_gd_cnv to identify
    # which transitions cross known breakpoint boundaries.
    bin_coords_by_idx: Dict[int, Tuple[int, int]] = dict(zip(
        bin_mappings_df["array_idx"].astype(int),
        zip(
            bin_mappings_df["start"].astype(int),
            bin_mappings_df["end"].astype(int),
        ),
    ))

    for cluster, locus in gd_table.loci.items():
        print(f"\nCalling CNVs for locus: {cluster}")

        interval_bins = get_locus_interval_bins(bin_mappings_df, cluster)

        # Breakpoint-range bins (inside SD-block ranges, not in any body
        # interval or flank) carry unreliable depth signal and must be
        # excluded from CNV calling.  They should already be absent from
        # the mappings file but this is a safety guard.
        bp_masked = interval_bins.pop("breakpoint_ranges", [])
        if bp_masked:
            print(f"  Masking {len(bp_masked)} breakpoint-range bin(s)")

        # Skip loci that have no bins (e.g. excluded by --region during
        # preprocessing).
        if not interval_bins:
            print("  Skipping — no bins in mappings for this locus")
            continue

        for interval_name, bins in interval_bins.items():
            print(f"  {interval_name}: {len(bins)} bins")

        # Warn when flanks are absent — the Viterbi trace will not cover
        # the flanking region, which is usually a preprocessing issue.
        for flank_name in ("left_flank", "right_flank"):
            if flank_name not in interval_bins or len(interval_bins[flank_name]) == 0:
                print(f"  WARNING: no {flank_name} bins for {cluster} — "
                      f"Viterbi trace will not cover that flank")

        interval_bin_arrays = {
            name: np.array(bins, dtype=int)
            for name, bins in interval_bins.items()
        }

        # Build the complete set of non-breakpoint bins for the cluster.
        # This ensures the Viterbi runs over ALL bins (including flanks)
        # even if some are not assigned to a named interval in bin_mappings.
        all_cluster_idxs = set(
            bin_mappings_df[
                bin_mappings_df["cluster"] == cluster
            ]["array_idx"].astype(int)
        )
        bp_set = set(int(b) for b in bp_masked)
        all_cluster_bins = sorted(all_cluster_idxs - bp_set)

        for s_idx, sample_id in enumerate(sample_ids):
            sample_ploidy = ploidy_lookup.get((str(sample_id), locus.chrom), 2)

            if verbose:
                strategy = "PAIR-VIT"
                _util.vlog(f"\n  [{strategy}] Sample: {sample_id}  "
                           f"({locus.chrom}, ploidy={sample_ploidy})")

            calls, path_records = viterbi_call_gd_cnv(
                locus,
                pair_prob_3d[s_idx],
                pair_state_labels,
                transition_matrix,
                interval_bin_arrays,
                ploidy=sample_ploidy,
                reciprocal_overlap_threshold=reciprocal_overlap_threshold,
                verbose=verbose,
                sample_id=str(sample_id),
                breakpoint_transition_matrix=breakpoint_transition_matrix,
                bin_coords=bin_coords_by_idx,
                all_cluster_bins=all_cluster_bins,
            )
            for start, end, cn_state, category, haplotype in path_records:
                all_path_records.append({
                    "sample": sample_id,
                    "cluster": cluster,
                    "start": start,
                    "end": end,
                    "cn_state": cn_state,
                    "category": category,
                    "haplotype": haplotype,
                })
            best_by_svtype = determine_best_breakpoints(calls)

            if verbose:
                for call in calls:
                    tag = " [CARRIER]" if call["is_carrier"] else ""
                    ro = call.get("reciprocal_overlap", 0.0)
                    _util.vlog(
                        f"    [VIT] {sample_id:30s}  "
                        f"{call['GD_ID']:25s}  "
                        f"{call['svtype']:4s}  ploidy={sample_ploidy}  "
                        f"score={call['log_prob_score']:+.4f}  "
                        f"RO={ro:.2%}  "
                        f"intervals={','.join(call['intervals'])}"
                        f"{tag}"
                    )

            for call in calls:
                svtype = call["svtype"]
                best_gd_for_svtype = best_by_svtype.get(svtype)

                covered_bin_indices = []
                for interval_name in call["intervals"]:
                    if interval_name in interval_bin_arrays:
                        covered_bin_indices.extend(interval_bin_arrays[interval_name])

                if len(covered_bin_indices) > 0:
                    mean_depth = float(depth_2d[s_idx, covered_bin_indices].mean())
                else:
                    mean_depth = np.nan

                result = {
                    "sample": sample_id,
                    "cluster": cluster,
                    "GD_ID": call["GD_ID"],
                    "chrom": call["chrom"],
                    "start": call["start"],
                    "end": call["end"],
                    "svtype": svtype,
                    "BP1": call["BP1"],
                    "BP2": call["BP2"],
                    "is_terminal": call["is_terminal"],
                    "n_bins": call["n_bins"],
                    "mean_depth": mean_depth,
                    "sample_ploidy": call.get("sample_ploidy", sample_ploidy),
                    "matched_haplotype": call.get("haplotype", np.nan),
                    "hap_cn_state": call.get("hap_cn_state", np.nan),
                    "matched_seg_start": call.get("matched_seg_start", np.nan),
                    "matched_seg_end": call.get("matched_seg_end", np.nan),
                    "matched_seg_n_bins": call.get("matched_seg_n_bins", 0),
                    "reciprocal_overlap": call.get("reciprocal_overlap", np.nan),
                    "is_carrier": call["is_carrier"],
                    "is_best_match": (
                        call["GD_ID"] == best_gd_for_svtype
                        if best_gd_for_svtype else False
                    ),
                    "log_prob_score": call["log_prob_score"],
                }
                all_results.append(result)

    paths_df = pd.DataFrame(all_path_records)
    return pd.DataFrame(all_results), paths_df


# =============================================================================
# Main
# =============================================================================


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Call GD CNVs from model posterior probabilities",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--cn-posteriors",
        required=True,
        help="CN posteriors file (cn_posteriors.tsv.gz) with depth values",
    )
    parser.add_argument(
        "--bin-mappings",
        required=True,
        help="Bin mappings file (bin_mappings.tsv.gz) from gd_cnv_pyro.py",
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
        "--output-dir", "-o",
        required=True,
        help="Output directory for calls",
    )
    parser.add_argument(
        "--transition-matrix",
        required=True,
        help="CN-state transition probability matrix (TSV) for NAHR loci.  "
             "Square K×K matrix where entry (i,j) is the probability of "
             "transitioning from CN state i to CN state j between adjacent "
             "bins.  Rows will be automatically normalized to sum to 1.",
    )
    parser.add_argument(
        "--breakpoint-transition-matrix",
        required=False,
        help="CN-state transition probability matrix (TSV) applied "
             "specifically at known recurrent breakpoint boundaries — both "
             "the start and end coordinate of each breakpoint's SD-block "
             "range.  Should have lower diagonal values than "
             "--transition-matrix to encourage CN state changes at these "
             "positions.  Same TSV format as --transition-matrix.  "
             "Only used when --transition-matrix is also set (Viterbi only).",
    )
    parser.add_argument(
        "--reciprocal-overlap-threshold",
        type=float,
        default=0.90,
        help="Minimum reciprocal overlap between a Viterbi CN segment and a "
             "GD entry to call a carrier.  Reciprocal overlap is "
             "min(overlap/segment_length, overlap/entry_length).  "
             "Default: 0.90 (90%%).",
    )
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Print detailed per-sample log probability scores for all GD "
             "entries at every locus, including flanking region scores.",
    )
    return parser.parse_args()


def main():
    """Main function."""
    args = parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    setup_logging(args.output_dir)

    print(f"Output directory: {args.output_dir}")

    # Load data
    print("\nLoading data...")

    print(f"  Loading CN posteriors: {args.cn_posteriors}")
    cn_posteriors_df = pd.read_csv(args.cn_posteriors, sep="\t", compression="infer")
    print(f"    {len(cn_posteriors_df)} bin-sample records")

    print(f"  Loading bin mappings: {args.bin_mappings}")
    bin_mappings_df = pd.read_csv(args.bin_mappings, sep="\t", compression="infer")
    print(f"    {len(bin_mappings_df)} bin mappings")

    print(f"  Loading GD table: {args.gd_table}")
    gd_table = GDTable(args.gd_table)
    print(f"    {len(gd_table.loci)} loci")

    # Load ploidy table if provided
    ploidy_df = None
    if args.ploidy_table:
        print(f"  Loading ploidy table: {args.ploidy_table}")
        ploidy_df = pd.read_csv(args.ploidy_table, sep="\t")
        print(f"    {len(ploidy_df)} sample/contig ploidy records")

    # Load transition matrix
    print(f"\n  Loading transition matrix: {args.transition_matrix}")
    transition_matrix = load_transition_matrix(args.transition_matrix)

    # Load breakpoint-specific transition matrix if provided
    breakpoint_transition_matrix = None
    if args.breakpoint_transition_matrix:
        print(f"  Loading breakpoint transition matrix: "
              f"{args.breakpoint_transition_matrix}")
        breakpoint_transition_matrix = load_transition_matrix(
            args.breakpoint_transition_matrix)

    # Call CNVs from posteriors
    calls_df, paths_df = call_cnvs_from_posteriors(
        cn_posteriors_df,
        bin_mappings_df,
        gd_table,
        transition_matrix=transition_matrix,
        ploidy_df=ploidy_df,
        verbose=args.verbose,
        reciprocal_overlap_threshold=args.reciprocal_overlap_threshold,
        breakpoint_transition_matrix=breakpoint_transition_matrix,
    )

    # Save calls
    output_file = os.path.join(args.output_dir, "gd_cnv_calls.tsv.gz")
    calls_df.to_csv(output_file, sep="\t", index=False, compression="gzip")
    print(f"\n  Saved calls to: {output_file}")
    print(f"    {len(calls_df)} call records")

    # Save Viterbi paths
    paths_file = os.path.join(args.output_dir, "viterbi_paths.tsv.gz")
    paths_df.to_csv(paths_file, sep="\t", index=False, compression="gzip")
    print(f"  Saved Viterbi paths to: {paths_file}")
    print(f"    {len(paths_df)} path records")

    # Print carrier summary
    if len(calls_df) > 0:
        carriers = calls_df[calls_df["is_carrier"]]
        n_carriers = carriers["sample"].nunique()
        n_sites = carriers["GD_ID"].nunique()
        print(f"\n  {n_carriers} carrier samples across {n_sites} GD sites")
    else:
        print("\n  No calls produced — check that the GD table, bin mappings, "
              "and CN posteriors refer to the same loci.")

    print("\n" + "=" * 80)
    print("Calling complete!")
    print("=" * 80)


if __name__ == "__main__":
    main()
