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

from gatk_sv_gd._util import TeeStream, setup_logging


# =============================================================================
# Truth Table Loading
# =============================================================================


def _load_truth_table_bed_format(df: pd.DataFrame) -> pd.DataFrame:
    """Load the classic BED-style truth table (8-column format).

    Expected columns (after ``#`` stripping)::

        chrom  start  end  name  svtype  samples  NAHR_GD  NAHR_GD_atypical

    Returns the internal schema: ``GD_ID``, ``chr``, ``start``, ``end``,
    ``SVTYPE``, ``cluster_ID``, ``carrier_set``.
    """
    required = {"chrom", "start", "end", "name", "svtype", "samples",
                "NAHR_GD", "NAHR_GD_atypical"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Truth table missing required columns: {missing}")

    def _to_bool(val):
        return str(val).strip().lower() == "true"

    df["NAHR_GD"] = df["NAHR_GD"].apply(_to_bool)
    df["NAHR_GD_atypical"] = df["NAHR_GD_atypical"].apply(_to_bool)

    # Drop atypical records
    n_before = len(df)
    df = df[~df["NAHR_GD_atypical"]].copy()
    n_atypical = n_before - len(df)
    if n_atypical > 0:
        print(f"  Dropped {n_atypical} atypical rows from truth table")

    # Drop non-NAHR entries
    n_before_nahr = len(df)
    df = df[df["NAHR_GD"]].copy()
    n_non_nahr = n_before_nahr - len(df)
    if n_non_nahr > 0:
        print(f"  Dropped {n_non_nahr} non-NAHR rows from truth table")
    print(f"  Truth table: {len(df)} NAHR entries (BED format)")

    df.rename(columns={
        "chrom": "chr",
        "name": "GD_ID",
        "svtype": "SVTYPE",
        "samples": "carriers",
    }, inplace=True)
    df["cluster_ID"] = ""

    def _parse_carriers(val):
        if pd.isna(val) or str(val).strip() == "":
            return set()
        return set(s.strip() for s in str(val).split(",") if s.strip())

    df["carrier_set"] = df["carriers"].apply(_parse_carriers)
    return df


def _load_truth_table_synthesize_format(df: pd.DataFrame) -> pd.DataFrame:
    """Load the two-column truth table produced by ``gatk-sv-gd synthesize``.

    Expected columns::

        sample_id  GD_ID

    Each row represents a single sample carrying the given GD.  Rows are
    grouped by ``GD_ID`` to produce one output row per GD with a
    ``carrier_set``.

    Returns the internal schema: ``GD_ID``, ``chr``, ``start``, ``end``,
    ``SVTYPE``, ``cluster_ID``, ``carrier_set`` (genomic coordinate columns
    are empty because the synthesize table does not carry them).
    """
    required = {"sample_id", "GD_ID"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"Synthesize truth table missing required columns: {missing}"
        )

    rows = []
    for gd_id, grp in df.groupby("GD_ID"):
        carrier_set = set(grp["sample_id"].astype(str).str.strip())
        rows.append({
            "GD_ID": str(gd_id),
            "chr": "",
            "start": "",
            "end": "",
            "SVTYPE": "",
            "cluster_ID": "",
            "carrier_set": carrier_set,
        })

    result = pd.DataFrame(rows)
    n_carriers = df["sample_id"].nunique()
    print(f"  Truth table: {len(result)} GD entries, "
          f"{n_carriers} total carriers (synthesize format)")
    return result


def load_truth_table(filepath: str) -> pd.DataFrame:
    """
    Load a truth table of GD carrier calls.

    Two formats are accepted and auto-detected from the column names:

    **BED-style format** (manually curated)::

        #chrom  start  end  name  svtype  samples  NAHR_GD  NAHR_GD_atypical

    - Only NAHR rows (``NAHR_GD=True``, ``NAHR_GD_atypical=False``) are kept.
    - The *name* column is used as the GD region ID.
    - The *samples* column is a comma-separated list of carrier sample IDs.

    **Synthesize format** (produced by ``gatk-sv-gd synthesize``)::

        sample_id  GD_ID

    - Each row is one sample carrying the given GD.
    - All rows are used as-is (no NAHR/atypical filtering).

    Args:
        filepath: Path to TSV truth table.

    Returns:
        DataFrame with one row per truth site.  Columns are normalised to
        ``GD_ID``, ``chr``, ``start``, ``end``, ``SVTYPE``, ``cluster_ID``,
        ``carrier_set``.
    """
    df = pd.read_csv(filepath, sep="\t", dtype=str, comment=None)
    # Strip whitespace and normalise the leading '#' that BED-style
    # headers sometimes carry (e.g. "#chrom" → "chrom").
    df.columns = [c.strip().lstrip("#") for c in df.columns]

    # Auto-detect format based on column names
    if ("sample_id" in df.columns and "GD_ID" in df.columns
            and "chrom" not in df.columns and "NAHR_GD" not in df.columns):
        print("  Detected synthesize truth table format")
        return _load_truth_table_synthesize_format(df)

    return _load_truth_table_bed_format(df)


def _get_confidence_column(calls_df: pd.DataFrame) -> Optional[str]:
    """Return the preferred confidence column present in the calls file."""
    if "confidence_score" in calls_df.columns:
        return "confidence_score"
    if "log_prob_score" in calls_df.columns:
        return "log_prob_score"
    return None


# =============================================================================
# Evaluation
# =============================================================================


def evaluate_against_truth(
    calls_df: pd.DataFrame,
    truth_df: pd.DataFrame,
    output_dir: str,
    batch_samples: Optional[set] = None,
    min_confidence: Optional[float] = None,
) -> pd.DataFrame:
    """
    Cross-reference predicted GD calls against a truth table and report
    sensitivity and precision for every site with at least one truth or
    predicted carrier.

    Truth entries are matched to calls by ``GD_ID``.

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
        min_confidence: Optional minimum confidence required for a predicted
            carrier call to count during evaluation. Uses
            ``confidence_score`` when present, otherwise falls back to
            ``log_prob_score``. If ``None``, no confidence filter is applied.

    Returns:
        Per-site report DataFrame.
    """
    print("\n" + "=" * 80)
    print("EVALUATING PREDICTIONS AGAINST TRUTH TABLE")
    print("=" * 80)

    if batch_samples is not None:
        print(f"  Batch contains {len(batch_samples)} samples; "
              "truth carriers will be restricted to this set.")
    score_column = None
    if min_confidence is not None:
        score_column = _get_confidence_column(calls_df)
        if score_column is None:
            raise ValueError(
                "--min-confidence requires a calls file with either "
                "confidence_score or log_prob_score column"
            )
        print(
            f"  Enforcing call confidence threshold: "
            f"{score_column} >= {min_confidence:.3f}"
        )

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
        if min_confidence is not None:
            carrier_mask = carrier_mask & (grp[score_column] >= min_confidence)
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
    pred_gd_ids = {gd for gd, s in pred_by_gd.items() if len(s) > 0}
    truth_gd_ids = {gd for gd, d in truth_by_gd.items() if len(d["carrier_set"]) > 0}
    all_gd_ids = sorted(
        pred_gd_ids | truth_gd_ids
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
        help="Truth table TSV.  Two formats are accepted: "
             "(1) BED-style manually curated table with columns: "
             "#chrom, start, end, name, svtype, samples, NAHR_GD, "
             "NAHR_GD_atypical (only NAHR_GD=True rows are used); "
             "(2) two-column synthesize table produced by "
             "'gatk-sv-gd synthesize' with columns: sample_id, GD_ID. "
             "The format is auto-detected from the header.",
    )
    parser.add_argument(
        "--gd-table", "-g",
        required=False, default=None,
        help="Filtered GD table (gd_table_filtered.tsv) produced by the "
             "preprocess step.  When provided, truth entries whose GD_ID "
             "does not appear in this table are excluded from evaluation "
             "so that loci removed by region or size filtering do not "
             "produce spurious false negatives.  Recommended: use the "
             "gd_table_filtered.tsv written to the preprocess output "
             "directory.",
    )
    parser.add_argument(
        "--ploidy-table",
        required=True,
        help="Ploidy estimates table (TSV) produced by the preprocess step "
             "(ploidy_estimates.tsv).  Used to identify the set of samples "
             "in the current batch.  Truth carriers absent from this sample "
             "set are excluded from sensitivity scoring.",
    )
    parser.add_argument(
        "--output-dir", "-o",
        required=True,
        help="Output directory for evaluation report",
    )
    parser.add_argument(
        "--min-confidence",
        nargs="?",
        const=-0.3,
        default=None,
        type=float,
        help="Optional minimum confidence required for a predicted carrier to "
             "count in evaluation. Uses confidence_score when present, "
             "otherwise falls back to log_prob_score. If the flag is provided "
             "without a value, uses -0.3. If omitted entirely, no confidence "
             "threshold is enforced.",
    )
    return parser.parse_args()


def main():
    """Main function."""
    args = parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    setup_logging(args.output_dir, filename="eval_log.txt")

    print(f"Output directory: {args.output_dir}")

    # Load calls
    print(f"\n  Loading calls: {args.calls}")
    calls_df = pd.read_csv(args.calls, sep="\t", compression="infer")
    print(f"    {len(calls_df)} call records")

    # Load ploidy table to identify batch samples
    print(f"\n  Loading ploidy table: {args.ploidy_table}")
    ploidy_df = pd.read_csv(args.ploidy_table, sep="\t")
    batch_samples = set(ploidy_df["sample"].astype(str).unique())
    print(f"    {len(batch_samples)} samples in batch")

    # Load truth table
    print(f"\n  Loading truth table: {args.truth_table}")
    truth_df = load_truth_table(args.truth_table)
    print(f"    {len(truth_df)} truth entries")

    # Restrict truth entries to loci that were actually modeled in this run.
    # The calls file contains one row per modeled GD_ID (not just carriers),
    # so it is the safest default scope even when --gd-table is omitted.
    modeled_gd_ids = set(calls_df["GD_ID"].astype(str).unique())
    print(f"\n  Derived {len(modeled_gd_ids)} modeled GD_IDs from calls file")

    if args.gd_table is not None:
        print(f"  Loading filtered GD table: {args.gd_table}")
        gd_df = pd.read_csv(args.gd_table, sep="\t")
        gd_table_ids = set(gd_df["GD_ID"].astype(str).unique())
        print(f"    {len(gd_table_ids)} modeled GD_IDs in filtered table")
        modeled_gd_ids &= gd_table_ids

    n_before = len(truth_df)
    truth_df = truth_df[
        truth_df["GD_ID"].astype(str).isin(modeled_gd_ids)
    ].copy()
    n_dropped = n_before - len(truth_df)
    if n_dropped > 0:
        print(f"    Dropped {n_dropped} truth entries at unmodeled loci")
    print(f"    {len(truth_df)} truth entries after filtering")

    # Evaluate
    evaluate_against_truth(
        calls_df, truth_df, args.output_dir,
        batch_samples=batch_samples,
        min_confidence=args.min_confidence,
    )

    print("\n" + "=" * 80)
    print("Evaluation complete!")
    print("=" * 80)


if __name__ == "__main__":
    main()
