#!/bin/python

"""
GD CNV Calling from Model Posteriors

Call genomic disorder (GD) copy-number variants from posterior probabilities
produced by gd_cnv_pyro.py.  For each GD locus in the table the script
evaluates every possible breakpoint configuration per sample and emits a
carrier / non-carrier call with a log-probability score.

Usage:
    python gd_cnv_call.py \
        --cn-posteriors cn_posteriors.tsv.gz \
        --bin-mappings bin_mappings.tsv.gz \
        --gd-table gd_table.tsv \
        --output-dir results/
"""

import argparse
import os
import sys
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

# Import GDLocus and GDTable from gd_cnv_pyro.py to avoid duplication
from gd_cnv_pyro import GDLocus, GDTable


# =============================================================================
# Logging
# =============================================================================


class TeeStream:
    """Write to both the original stream and a log file."""

    def __init__(self, original_stream, log_file):
        self.original_stream = original_stream
        self.log_file = log_file

    def write(self, message: str):
        self.original_stream.write(message)
        self.log_file.write(message)

    def flush(self):
        self.original_stream.flush()
        self.log_file.flush()

    def __getattr__(self, name):
        return getattr(self.original_stream, name)


def setup_logging(output_dir: str, filename: str = "call_log.txt"):
    """Redirect stdout and stderr to both the console and a log file."""
    log_path = os.path.join(output_dir, filename)
    log_fh = open(log_path, "w")
    sys.stdout = TeeStream(sys.__stdout__, log_fh)
    sys.stderr = TeeStream(sys.__stderr__, log_fh)
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


def compute_interval_cn_stats(
    cn_posteriors_df: pd.DataFrame,
    interval_bins: Dict[str, List[int]],
    sample_id: str,
) -> Dict[str, dict]:
    """
    Compute copy number statistics for each interval for a specific sample.

    Args:
        cn_posteriors_df: DataFrame with CN posterior probabilities
        interval_bins: Dict mapping interval name to list of array_idx values
            (positional indices into the per-sample posterior rows)
        sample_id: Sample identifier

    Returns:
        Dict mapping interval name to CN statistics
    """
    result = {}

    sample_posteriors = (
        cn_posteriors_df[cn_posteriors_df["sample"] == sample_id]
        .reset_index(drop=True)
    )

    prob_cols = [c for c in sample_posteriors.columns if c.startswith("prob_cn_")]
    n_states = len(prob_cols)

    for interval_name, bin_indices in interval_bins.items():
        interval_data = sample_posteriors.iloc[bin_indices]

        if len(interval_data) == 0:
            result[interval_name] = {
                "n_bins": 0,
                "cn_probs": np.zeros(n_states),
            }
            continue

        mean_probs = interval_data[prob_cols].mean(axis=0).values

        result[interval_name] = {
            "n_bins": len(interval_data),
            "cn_probs": mean_probs,
        }

    return result


def call_gd_cnv(
    locus: GDLocus,
    interval_stats: Dict[str, dict],
    log_prob_threshold: float = -0.5,
    flanking_log_prob_threshold: float = -1.0,
    ploidy: int = 2,
) -> List[dict]:
    """
    Call GD CNVs based on interval copy number statistics.

    For each GD entry in the locus (each DEL/DUP definition), check if the
    copy number in the corresponding region supports a CNV call.

    Args:
        locus: GDLocus object
        interval_stats: Dict mapping interval name to CN statistics
        log_prob_threshold: Minimum log probability score to call a CNV
        flanking_log_prob_threshold: Minimum log probability score in flanking
            regions to classify as spanning
        ploidy: Expected copy number for this sample on this chromosome

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

        covered_tuples = locus.get_intervals_between(bp1, bp2)
        covered_intervals = [name for _, _, name in covered_tuples]

        if len(covered_intervals) == 0:
            continue

        # Compute weighted average log probability score
        log_prob_score = 0.0
        total_weight = 0.0

        if svtype == "DEL":
            for interval in covered_intervals:
                if interval in interval_stats and interval_stats[interval]["n_bins"] > 0:
                    probs = interval_stats[interval]["cn_probs"]
                    weight = interval_stats[interval]["n_bins"]
                    p_del = max(probs[:ploidy].sum(), 1e-5)
                    if p_del > 0:
                        log_prob_score += weight * np.log(p_del)
                        total_weight += weight

        elif svtype == "DUP":
            for interval in covered_intervals:
                if interval in interval_stats and interval_stats[interval]["n_bins"] > 0:
                    probs = interval_stats[interval]["cn_probs"]
                    weight = interval_stats[interval]["n_bins"]
                    p_dup = max(probs[ploidy + 1:].sum(), 1e-5)
                    if p_dup > 0:
                        log_prob_score += weight * np.log(p_dup)
                        total_weight += weight

        if total_weight > 0:
            log_prob_score = log_prob_score / total_weight
        else:
            log_prob_score = np.nan

        is_cnv = log_prob_score > log_prob_threshold

        # Check flanking regions for evidence of a large spanning CNV
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
                    p_flank = max(flank_probs[:ploidy].sum(), 1e-5)
                elif svtype == "DUP":
                    p_flank = max(flank_probs[ploidy + 1:].sum(), 1e-5)
                else:
                    continue
                flank_score_sum += flank_weight * np.log(p_flank)
                flank_weight_total += flank_weight

            if flank_weight_total > 0:
                flanking_log_prob_score = flank_score_sum / flank_weight_total
                if flanking_log_prob_score > flanking_log_prob_threshold:
                    is_spanning = True
                    is_cnv = False

        n_bins = sum(interval_stats[iv]["n_bins"] for iv in covered_intervals if iv in interval_stats)

        calls.append({
            "GD_ID": gd_id,
            "cluster": locus.cluster,
            "chrom": locus.chrom,
            "start": gd_start,
            "end": gd_end,
            "svtype": svtype,
            "BP1": bp1,
            "BP2": bp2,
            "is_nahr": locus.is_nahr,
            "is_terminal": locus.is_terminal,
            "log_prob_score": log_prob_score,
            "flanking_log_prob_score": flanking_log_prob_score,
            "is_carrier": is_cnv,
            "is_spanning": is_spanning,
            "intervals": covered_intervals,
            "n_bins": n_bins,
        })

    return calls


def determine_best_breakpoints(
    locus: GDLocus,
    interval_stats: Dict[str, dict],
    calls: List[dict],
    ploidy: int = 2,
) -> Dict[str, Optional[str]]:
    """
    Determine the most likely breakpoint pair for a GD CNV, separately for
    DEL and DUP.

    Args:
        locus: GDLocus object
        interval_stats: Dict mapping interval name to CN statistics
        calls: List of CNV call dictionaries
        ploidy: Expected copy number for this sample on this chromosome

    Returns:
        Dict mapping svtype ("DEL", "DUP") to best matching GD_ID or None
    """
    best_by_svtype = {}

    all_intervals = set(name for _, _, name in locus.get_intervals())

    for svtype in ["DEL", "DUP"]:
        carrier_calls = [c for c in calls if c["is_carrier"] and c["svtype"] == svtype]

        if len(carrier_calls) == 0:
            best_by_svtype[svtype] = None
            continue

        if len(carrier_calls) == 1:
            best_by_svtype[svtype] = carrier_calls[0]["GD_ID"]
            continue

        best_score = -np.inf
        best_size = -1
        best_gd_id = None

        for call in carrier_calls:
            covered_intervals = set(call["intervals"])
            uncovered_intervals = all_intervals - covered_intervals

            score = 0.0
            total_weight = 0.0

            if svtype == "DEL":
                for interval in covered_intervals:
                    if interval in interval_stats and interval_stats[interval]["n_bins"] > 0:
                        probs = interval_stats[interval]["cn_probs"]
                        weight = 1.
                        p_del = probs[:ploidy].sum()
                        if p_del > 0:
                            score += weight * np.log(p_del)
                            total_weight += weight

                for interval in uncovered_intervals:
                    if interval in interval_stats and interval_stats[interval]["n_bins"] > 0:
                        probs = interval_stats[interval]["cn_probs"]
                        weight = 1.
                        p_normal = probs[ploidy:].sum()
                        if p_normal > 0:
                            score += weight * np.log(p_normal)
                            total_weight += weight

            else:  # DUP
                for interval in covered_intervals:
                    if interval in interval_stats and interval_stats[interval]["n_bins"] > 0:
                        probs = interval_stats[interval]["cn_probs"]
                        weight = 1.
                        p_dup = probs[ploidy + 1:].sum()
                        if p_dup > 0:
                            score += weight * np.log(p_dup)
                            total_weight += weight

                for interval in uncovered_intervals:
                    if interval in interval_stats and interval_stats[interval]["n_bins"] > 0:
                        probs = interval_stats[interval]["cn_probs"]
                        weight = 1.
                        p_normal = probs[:ploidy + 1].sum()
                        if p_normal > 0:
                            score += weight * np.log(p_normal)
                            total_weight += weight

            if total_weight > 0:
                score = score / total_weight

            call_size = call["end"] - call["start"]
            if score > best_score or (score == best_score and call_size > best_size):
                best_score = score
                best_size = call_size
                best_gd_id = call["GD_ID"]

        best_by_svtype[svtype] = best_gd_id

    return best_by_svtype


def call_cnvs_from_posteriors(
    cn_posteriors_df: pd.DataFrame,
    bin_mappings_df: pd.DataFrame,
    gd_table: GDTable,
    log_prob_threshold: float = -0.5,
    flanking_log_prob_threshold: float = -1.0,
    ploidy_df: Optional[pd.DataFrame] = None,
    verbose: bool = False,
) -> pd.DataFrame:
    """
    Call CNVs from posterior probabilities.

    Args:
        cn_posteriors_df: DataFrame with CN posterior probabilities per bin/sample
        bin_mappings_df: DataFrame with bin-to-interval mappings
        gd_table: GDTable with locus definitions
        log_prob_threshold: Minimum log probability score to call a CNV
        flanking_log_prob_threshold: Minimum log probability score in flanking
            regions to classify as spanning
        ploidy_df: Optional DataFrame with columns (sample, contig, ploidy).
            If None, ploidy=2 is assumed for all sample/contig pairs.
        verbose: If True, print per-sample log probability scores for every
            GD entry at every locus, including flanking region scores.

    Returns:
        DataFrame with CNV calls for all samples and loci
    """
    print("\n" + "=" * 80)
    print("CALLING CNVs FROM POSTERIORS")
    print("=" * 80)

    all_results = []
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
    prob_cols = sorted(c for c in cn_posteriors_df.columns if c.startswith("prob_cn_"))
    n_states = len(prob_cols)
    prob_3d = np.empty((n_samples, n_bins, n_states))
    depth_2d = np.empty((n_samples, n_bins))
    for s_idx, sample_id in enumerate(sample_ids):
        mask = cn_posteriors_df["sample"] == sample_id
        sample_rows = cn_posteriors_df.loc[mask]
        prob_3d[s_idx] = sample_rows[prob_cols].values
        depth_2d[s_idx] = sample_rows["depth"].values
    print(f"    Extracted {n_samples} x {n_bins} x {n_states} probability array")

    for cluster, locus in gd_table.loci.items():
        print(f"\nCalling CNVs for locus: {cluster}")

        interval_bins = get_locus_interval_bins(bin_mappings_df, cluster)

        for interval_name, bins in interval_bins.items():
            print(f"  {interval_name}: {len(bins)} bins")

        interval_bin_arrays = {
            name: np.array(bins, dtype=int)
            for name, bins in interval_bins.items()
        }

        for s_idx, sample_id in enumerate(sample_ids):
            interval_stats = {}
            for interval_name, bin_idx in interval_bin_arrays.items():
                if len(bin_idx) == 0:
                    interval_stats[interval_name] = {
                        "n_bins": 0,
                        "cn_probs": np.zeros(n_states),
                    }
                else:
                    mean_probs = prob_3d[s_idx, bin_idx].mean(axis=0)
                    interval_stats[interval_name] = {
                        "n_bins": len(bin_idx),
                        "cn_probs": mean_probs,
                    }

            sample_ploidy = ploidy_lookup.get((str(sample_id), locus.chrom), 2)

            calls = call_gd_cnv(locus, interval_stats, log_prob_threshold,
                                flanking_log_prob_threshold, ploidy=sample_ploidy)

            best_by_svtype = determine_best_breakpoints(locus, interval_stats,
                                                        calls, ploidy=sample_ploidy)

            if verbose:
                for call in calls:
                    tag_parts = []
                    if call["is_carrier"]:
                        tag_parts.append("CARRIER")
                    if call.get("is_spanning"):
                        tag_parts.append("SPANNING")
                    tag = " [" + ",".join(tag_parts) + "]" if tag_parts else ""
                    flank_str = (
                        f"{call['flanking_log_prob_score']:.4f}"
                        if not np.isnan(call.get("flanking_log_prob_score", np.nan))
                        else "NA"
                    )
                    print(
                        f"    {sample_id:30s}  {call['GD_ID']:25s}  "
                        f"{call['svtype']:4s}  ploidy={sample_ploidy}  "
                        f"logP={call['log_prob_score']:+.4f}  "
                        f"flank={flank_str}  "
                        f"intervals={','.join(call['intervals'])}{tag}"
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
                    "is_nahr": call["is_nahr"],
                    "is_terminal": call["is_terminal"],
                    "n_bins": call["n_bins"],
                    "mean_depth": mean_depth,
                    "is_carrier": call["is_carrier"],
                    "is_best_match": call["GD_ID"] == best_gd_for_svtype if best_gd_for_svtype else False,
                    "log_prob_score": call["log_prob_score"],
                }
                all_results.append(result)

    return pd.DataFrame(all_results)


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
        "--log-prob-threshold",
        type=float,
        default=-0.3,
        help="Minimum log probability score to call a CNV",
    )
    parser.add_argument(
        "--flanking-log-prob-threshold",
        type=float,
        default=-1.0,
        help="Minimum log probability score in flanking regions to classify "
             "a call as spanning",
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
    log_fh = setup_logging(args.output_dir)

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

    # Call CNVs from posteriors
    calls_df = call_cnvs_from_posteriors(
        cn_posteriors_df,
        bin_mappings_df,
        gd_table,
        log_prob_threshold=args.log_prob_threshold,
        flanking_log_prob_threshold=args.flanking_log_prob_threshold,
        ploidy_df=ploidy_df,
        verbose=args.verbose,
    )

    # Save calls
    output_file = os.path.join(args.output_dir, "gd_cnv_calls.tsv.gz")
    calls_df.to_csv(output_file, sep="\t", index=False, compression="gzip")
    print(f"\n  Saved calls to: {output_file}")
    print(f"    {len(calls_df)} call records")

    # Print carrier summary
    carriers = calls_df[calls_df["is_carrier"]]
    n_carriers = carriers["sample"].nunique()
    n_sites = carriers["GD_ID"].nunique()
    print(f"\n  {n_carriers} carrier samples across {n_sites} GD sites")

    print("\n" + "=" * 80)
    print("Calling complete!")
    print("=" * 80)


if __name__ == "__main__":
    main()
