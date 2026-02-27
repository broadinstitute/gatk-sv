#!/bin/python

"""
GD CNV Evaluation

Evaluate GD copy-number variant calls produced by gd_cnv_call.py against a
manually curated truth table.  Reports per-site and overall sensitivity and
precision.

Usage:
    python gd_cnv_eval.py \
        --calls gd_cnv_calls.tsv.gz \
        --truth-table truth.tsv \
        --output-dir results/
"""

import argparse
import os
import sys
from typing import Dict, Optional

import numpy as np
import pandas as pd


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


def setup_logging(output_dir: str, filename: str = "eval_log.txt"):
    """Redirect stdout and stderr to both the console and a log file."""
    log_path = os.path.join(output_dir, filename)
    log_fh = open(log_path, "w")
    sys.stdout = TeeStream(sys.__stdout__, log_fh)
    sys.stderr = TeeStream(sys.__stderr__, log_fh)
    return log_fh


# =============================================================================
# Truth Table Loading
# =============================================================================


def load_truth_table(filepath: str) -> pd.DataFrame:
    """
    Load a truth table of manually labelled GD carrier calls.

    Expected columns (tab-separated, header may start with ``#``)::

        #chrom  start  end  name  svtype  samples  NAHR_GD  NAHR_GD_atypical

    For rows where ``NAHR_GD`` is True the *name* column is treated as the
    GD region ID.  Rows where ``NAHR_GD_atypical`` is True (atypical
    breakpoints) are dropped.

    The *samples* column is a comma-separated list of carrier sample IDs.

    Args:
        filepath: Path to TSV truth table.

    Returns:
        DataFrame with one row per GD site, carriers parsed into a set.
        Columns are normalised to the internal schema used by
        :func:`evaluate_against_truth` (``GD_ID``, ``chr``, ``start``,
        ``end``, ``SVTYPE``, ``carrier_set``).
    """
    df = pd.read_csv(filepath, sep="\t", dtype=str, comment=None)
    # Strip whitespace and normalise the leading '#' that BED-style
    # headers sometimes carry (e.g. "#chrom" → "chrom").
    df.columns = [c.strip().lstrip("#") for c in df.columns]

    required = {"chrom", "start", "end", "name", "svtype", "samples",
                "NAHR_GD", "NAHR_GD_atypical"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Truth table missing required columns: {missing}")

    # Parse boolean columns (accept True/False, true/false, TRUE/FALSE)
    def _to_bool(val):
        return str(val).strip().lower() == "true"

    df["NAHR_GD"] = df["NAHR_GD"].apply(_to_bool)
    df["NAHR_GD_atypical"] = df["NAHR_GD_atypical"].apply(_to_bool)

    # Keep only canonical NAHR-mediated GD rows; drop atypical records
    n_before = len(df)
    df = df[df["NAHR_GD"] & ~df["NAHR_GD_atypical"]].copy()
    n_dropped = n_before - len(df)
    if n_dropped > 0:
        print(f"  Dropped {n_dropped} non-NAHR / atypical rows from truth table")

    # Map columns to internal schema
    df.rename(columns={
        "chrom": "chr",
        "name": "GD_ID",
        "svtype": "SVTYPE",
        "samples": "carriers",
    }, inplace=True)
    # No cluster_ID in this format
    df["cluster_ID"] = ""

    # Parse carriers into sets; treat empty / NaN as no carriers
    def _parse_carriers(val):
        if pd.isna(val) or str(val).strip() == "":
            return set()
        return set(s.strip() for s in str(val).split(",") if s.strip())

    df["carrier_set"] = df["carriers"].apply(_parse_carriers)
    return df


# =============================================================================
# Evaluation
# =============================================================================


def evaluate_against_truth(
    calls_df: pd.DataFrame,
    truth_df: pd.DataFrame,
    output_dir: str,
    batch_samples: Optional[set] = None,
) -> pd.DataFrame:
    """
    Cross-reference predicted GD calls against a truth table and report
    sensitivity and precision for every site with at least one truth or
    predicted carrier.

    Matching is performed by *GD_ID*.  For each GD_ID the set of predicted
    carriers (samples where ``is_carrier == True``) is compared to the truth
    carrier set.

    Truth carrier sets are intersected with *batch_samples* (when provided)
    so that sensitivity is not penalised for labelled carriers absent from
    the current batch.

    Args:
        calls_df: Predicted calls DataFrame (from gd_cnv_call.py or loaded
            from ``--calls``).
        truth_df: Truth table DataFrame produced by :func:`load_truth_table`.
        output_dir: Directory to write the report TSV.
        batch_samples: Optional set of sample IDs present in the current
            batch.  If provided, truth carriers not in this set are removed
            before scoring.

    Returns:
        Per-site report DataFrame.
    """
    print("\n" + "=" * 80)
    print("EVALUATING PREDICTIONS AGAINST TRUTH TABLE")
    print("=" * 80)

    if batch_samples is not None:
        print(f"  Batch contains {len(batch_samples)} samples; "
              "truth carriers will be restricted to this set.")

    # Build predicted carrier sets keyed by GD_ID.
    # Only count samples whose call is both a carrier AND the best-match
    # breakpoint configuration for its svtype.
    pred_by_gd: Dict[str, set] = {}
    pred_meta_by_gd: Dict[str, dict] = {}
    for gd_id, grp in calls_df.groupby("GD_ID"):
        gd_id_str = str(gd_id)
        carrier_mask = grp["is_carrier"] == True  # noqa: E712
        if "is_best_match" in grp.columns:
            carrier_mask = carrier_mask & (grp["is_best_match"] == True)  # noqa: E712
        pred_by_gd[gd_id_str] = set(
            grp.loc[carrier_mask, "sample"].unique()
        )
        first = grp.iloc[0]
        pred_meta_by_gd[gd_id_str] = {
            "chr": first.get("chrom", ""),
            "start": first.get("start", ""),
            "end": first.get("end", ""),
            "cluster_ID": first.get("cluster", ""),
            "SVTYPE": first.get("svtype", ""),
        }

    # Build truth carrier sets keyed by GD_ID
    truth_by_gd: Dict[str, dict] = {}
    for _, row in truth_df.iterrows():
        gd_id = str(row["GD_ID"])
        carrier_set = row["carrier_set"]
        if batch_samples is not None:
            carrier_set = carrier_set & batch_samples
        truth_by_gd[gd_id] = {
            "carrier_set": carrier_set,
            "chr": row.get("chr", ""),
            "start": row.get("start", ""),
            "end": row.get("end", ""),
            "cluster_ID": row.get("cluster_ID", ""),
            "SVTYPE": row.get("SVTYPE", ""),
        }

    # Union of all GD_IDs with at least one truth or predicted carrier
    all_gd_ids = sorted(
        {gd for gd, s in pred_by_gd.items() if len(s) > 0}
        | {gd for gd, d in truth_by_gd.items() if len(d["carrier_set"]) > 0}
    )

    rows = []
    total_tp = total_fp = total_fn = 0
    for gd_id in all_gd_ids:
        truth_set = truth_by_gd.get(gd_id, {}).get("carrier_set", set())
        pred_set = pred_by_gd.get(gd_id, set())

        tp = len(truth_set & pred_set)
        fp = len(pred_set - truth_set)
        fn = len(truth_set - pred_set)
        sensitivity = tp / (tp + fn) if (tp + fn) > 0 else float("nan")
        precision = tp / (tp + fp) if (tp + fp) > 0 else float("nan")

        total_tp += tp
        total_fp += fp
        total_fn += fn

        meta = truth_by_gd.get(gd_id, {})
        fallback = pred_meta_by_gd.get(gd_id, {})
        rows.append({
            "GD_ID": gd_id,
            "chr": meta.get("chr", "") or fallback.get("chr", ""),
            "start": meta.get("start", "") or fallback.get("start", ""),
            "end": meta.get("end", "") or fallback.get("end", ""),
            "cluster_ID": meta.get("cluster_ID", "") or fallback.get("cluster_ID", ""),
            "SVTYPE": meta.get("SVTYPE", "") or fallback.get("SVTYPE", ""),
            "n_truth_carriers": len(truth_set),
            "n_pred_carriers": len(pred_set),
            "TP": tp,
            "FP": fp,
            "FN": fn,
            "sensitivity": round(sensitivity, 4) if not np.isnan(sensitivity) else "NA",
            "precision": round(precision, 4) if not np.isnan(precision) else "NA",
            "FP_samples": ",".join(sorted(pred_set - truth_set)) if fp > 0 else "",
            "FN_samples": ",".join(sorted(truth_set - pred_set)) if fn > 0 else "",
        })

    report_df = pd.DataFrame(rows)

    # Print summary
    overall_sens = total_tp / (total_tp + total_fn) if (total_tp + total_fn) > 0 else float("nan")
    overall_prec = total_tp / (total_tp + total_fp) if (total_tp + total_fp) > 0 else float("nan")
    print(f"\n  Sites evaluated: {len(all_gd_ids)}")
    print(f"  Overall TP={total_tp}  FP={total_fp}  FN={total_fn}")
    print(f"  Overall sensitivity: {overall_sens:.4f}")
    print(f"  Overall precision:   {overall_prec:.4f}")

    # Per-site summary to stdout
    print(f"\n  {'GD_ID':30s}  {'SVTYPE':6s}  truth  pred    TP    FP    FN   sens   prec")
    print("  " + "-" * 100)
    for r in rows:
        sens_str = f"{r['sensitivity']:.2f}" if r["sensitivity"] != "NA" else "  NA"
        prec_str = f"{r['precision']:.2f}" if r["precision"] != "NA" else "  NA"
        print(f"  {r['GD_ID']:30s}  {str(r['SVTYPE']):6s}  "
              f"{r['n_truth_carriers']:5d}  {r['n_pred_carriers']:4d}  "
              f"{r['TP']:4d}  {r['FP']:4d}  {r['FN']:4d}  "
              f"{sens_str:>5s}  {prec_str:>5s}")

    # Write report
    output_path = os.path.join(output_dir, "truth_evaluation_report.tsv")
    report_df.to_csv(output_path, sep="\t", index=False)
    print(f"\n  Saved report: {output_path}")
    print("=" * 80 + "\n")

    return report_df


# =============================================================================
# Main
# =============================================================================


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Evaluate GD CNV calls against a truth table",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--calls", "-c",
        required=True,
        help="GD CNV calls file (gd_cnv_calls.tsv.gz) produced by gd_cnv_call.py",
    )
    parser.add_argument(
        "--truth-table", "-t",
        required=True,
        help="Truth table TSV with manually labelled GD carriers. "
             "Columns: #chrom, start, end, name, svtype, samples, NAHR_GD, "
             "NAHR_GD_atypical.  Only rows with NAHR_GD=True (and "
             "NAHR_GD_atypical=False) are used.",
    )
    parser.add_argument(
        "--sample-posteriors",
        required=False,
        help="Sample posteriors table (TSV) with columns: sample, sample_var_map. "
             "Used to identify the set of samples in the current batch. "
             "Truth carriers absent from this sample set are excluded from "
             "sensitivity scoring.",
    )
    parser.add_argument(
        "--output-dir", "-o",
        required=True,
        help="Output directory for evaluation report",
    )
    return parser.parse_args()


def main():
    """Main function."""
    args = parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    log_fh = setup_logging(args.output_dir)

    print(f"Output directory: {args.output_dir}")

    # Load calls
    print(f"\n  Loading calls: {args.calls}")
    calls_df = pd.read_csv(args.calls, sep="\t", compression="infer")
    print(f"    {len(calls_df)} call records")

    # Load sample posteriors to identify batch samples (if provided)
    batch_samples = None
    if args.sample_posteriors:
        print(f"\n  Loading sample posteriors: {args.sample_posteriors}")
        sample_post_df = pd.read_csv(args.sample_posteriors, sep="\t")
        batch_samples = set(sample_post_df["sample"].astype(str).unique())
        print(f"    {len(batch_samples)} samples in batch")

    # Load truth table
    print(f"\n  Loading truth table: {args.truth_table}")
    truth_df = load_truth_table(args.truth_table)
    print(f"    {len(truth_df)} truth entries")

    # Evaluate
    evaluate_against_truth(
        calls_df, truth_df, args.output_dir,
        batch_samples=batch_samples,
    )

    print("\n" + "=" * 80)
    print("Evaluation complete!")
    print("=" * 80)


if __name__ == "__main__":
    main()
