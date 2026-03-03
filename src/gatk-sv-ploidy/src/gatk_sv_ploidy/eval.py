"""
Eval subcommand — evaluate predictions against a truth set.

Loads aneuploidy-type predictions and a truth JSON, then computes a confusion
matrix, classification report, and per-type metrics.
"""

from __future__ import annotations

import argparse
import json
import logging
import os

import pandas as pd
from sklearn.metrics import classification_report, confusion_matrix

logger = logging.getLogger(__name__)


# ── metrics ─────────────────────────────────────────────────────────────────


def calculate_metrics(
    pred_df: pd.DataFrame,
    output_dir: str,
) -> None:
    """Compute and save accuracy, classification report, and confusion matrix.

    Args:
        pred_df: DataFrame with ``true_aneuploidy_type`` and
            ``predicted_aneuploidy_type`` columns.
        output_dir: Directory for the metrics report.
    """
    y_true = pred_df["true_aneuploidy_type"]
    y_pred = pred_df["predicted_aneuploidy_type"]

    correct = int((y_true == y_pred).sum())
    total = len(pred_df)
    accuracy = correct / total if total else 0.0

    logger.info("Overall accuracy: %d / %d = %.4f (%.2f%%)",
                correct, total, accuracy, 100 * accuracy)

    report_str = classification_report(y_true, y_pred, zero_division=0)
    logger.info("\n%s", report_str)

    labels = sorted(y_true.unique())
    cm = confusion_matrix(y_true, y_pred, labels=labels)
    cm_df = pd.DataFrame(cm, index=labels, columns=labels)
    logger.info("Confusion matrix:\n%s", cm_df)

    # Save text report
    report_path = os.path.join(output_dir, "metrics_report.txt")
    with open(report_path, "w") as fh:
        fh.write("=" * 80 + "\n")
        fh.write("ANEUPLOIDY TYPE PREDICTION METRICS\n")
        fh.write("=" * 80 + "\n\n")
        fh.write(
            f"Overall Accuracy: {correct}/{total} = "
            f"{accuracy:.4f} ({100 * accuracy:.2f}%)\n\n"
        )
        fh.write("-" * 80 + "\n")
        fh.write("Classification Report:\n")
        fh.write("-" * 80 + "\n")
        fh.write(report_str)
        fh.write("\n" + "-" * 80 + "\n")
        fh.write("Confusion Matrix:\n")
        fh.write("-" * 80 + "\n")
        fh.write(cm_df.to_string())
        fh.write("\n")

    logger.info("Metrics report saved to %s", report_path)


# ── CLI ─────────────────────────────────────────────────────────────────────


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the eval subcommand."""
    p = argparse.ArgumentParser(
        description="Evaluate aneuploidy predictions against a truth set",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "-p", "--predictions", required=True,
        help="aneuploidy_type_predictions.tsv (from 'call')",
    )
    p.add_argument(
        "-t", "--truth-json", required=True,
        help="JSON mapping sample ID → true aneuploidy type",
    )
    p.add_argument(
        "-o", "--output-dir", required=True,
        help="Output directory for metrics report",
    )
    p.add_argument(
        "-e", "--exclusion-list", default=None,
        help="Text file with sample IDs to exclude (one per line)",
    )
    return p.parse_args()


def main() -> None:
    """Entry point for ``gatk-sv-ploidy eval``."""
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    # Load predictions
    logger.info("Loading predictions: %s", args.predictions)
    pred_df = pd.read_csv(args.predictions, sep="\t")

    # Load truth and merge
    logger.info("Loading truth: %s", args.truth_json)
    with open(args.truth_json) as fh:
        truth = json.load(fh)

    # Update truth column from JSON (overrides any existing values)
    pred_df["true_aneuploidy_type"] = pred_df["sample"].map(
        lambda s: truth.get(str(s), "NORMAL")
    )

    # Exclude samples if requested
    if args.exclusion_list:
        from gatk_sv_ploidy._util import load_exclusion_ids

        excl = set(load_exclusion_ids(args.exclusion_list))
        n_before = len(pred_df)
        pred_df = pred_df[~pred_df["sample"].isin(excl)]
        logger.info("Excluded %d samples", n_before - len(pred_df))

    # Calculate and save metrics
    calculate_metrics(pred_df, args.output_dir)

    # Also save the merged predictions table
    merged_path = os.path.join(args.output_dir, "predictions_with_truth.tsv")
    pred_df.to_csv(merged_path, sep="\t", index=False)
    logger.info("Merged predictions saved to %s", merged_path)

    logger.info("Done.")


if __name__ == "__main__":
    main()
