"""
Eval subcommand — evaluate predictions against a truth set.

Loads aneuploidy-type predictions and a truth JSON, then computes a confusion
matrix, classification report, and per-type metrics.
"""

from __future__ import annotations

import argparse
import json
import os

import pandas as pd
from sklearn.metrics import classification_report, confusion_matrix

from gatk_sv_ploidy._logging import log_output_artifacts, tool_logging_context


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

    labels = sorted(set(y_true) | set(y_pred))
    if len(labels) == 1:
        only_label = labels[0]
        support = int((y_true == only_label).sum())
        predicted = int((y_pred == only_label).sum())
        tp = int(((y_true == only_label) & (y_pred == only_label)).sum())
        precision = tp / predicted if predicted else 0.0
        recall = tp / support if support else 0.0
        f1_score = (
            2 * precision * recall / (precision + recall)
            if (precision + recall) > 0
            else 0.0
        )
        report_str = (
            "              precision    recall  f1-score   support\n\n"
            f"{only_label:>12}       {precision:0.2f}      {recall:0.2f}      {f1_score:0.2f}         {support}\n\n"
            f"    accuracy                           {accuracy:0.2f}         {total}\n"
        )
    else:
        report_str = classification_report(
            y_true,
            y_pred,
            labels=labels,
            zero_division=0,
        )

    if len(labels) == 1:
        cm_df = pd.DataFrame(
            [[int(((y_true == labels[0]) & (y_pred == labels[0])).sum())]],
            index=labels,
            columns=labels,
        )
    else:
        cm = confusion_matrix(y_true, y_pred, labels=labels)
        cm_df = pd.DataFrame(cm, index=labels, columns=labels)

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


def _run_eval(args: argparse.Namespace, logger) -> None:
    """Run evaluation after CLI logging is configured."""

    # Load predictions
    pred_df = pd.read_csv(args.predictions, sep="\t")
    logger.info(
        "Loaded prediction table: rows=%d samples=%d",
        len(pred_df),
        int(pred_df["sample"].nunique()) if "sample" in pred_df.columns else 0,
    )

    # Load truth and merge
    with open(args.truth_json) as fh:
        truth = json.load(fh)
    logger.info("Loaded truth labels: n_labels=%d", len(truth))

    # Update truth column from JSON (overrides any existing values)
    pred_df["true_aneuploidy_type"] = pred_df["sample"].map(
        lambda s: truth.get(str(s), "NORMAL")
    )

    # Exclude samples if requested
    if args.exclusion_list:
        from gatk_sv_ploidy._util import load_exclusion_ids

        excl = set(load_exclusion_ids(args.exclusion_list))
        pred_df = pred_df[~pred_df["sample"].isin(excl)]
        logger.info("Applied exclusion list: excluded_samples=%d", len(excl))

    # Calculate and save metrics
    calculate_metrics(pred_df, args.output_dir)

    # Also save the merged predictions table
    merged_path = os.path.join(args.output_dir, "predictions_with_truth.tsv")
    pred_df.to_csv(merged_path, sep="\t", index=False)
    report_path = os.path.join(args.output_dir, "metrics_report.txt")
    accuracy = (
        float((pred_df["true_aneuploidy_type"] == pred_df["predicted_aneuploidy_type"]).mean())
        if len(pred_df) else 0.0
    )
    logger.info(
        "Evaluation complete: samples=%d accuracy=%.4f predicted_types=%d",
        len(pred_df),
        accuracy,
        int(pred_df["predicted_aneuploidy_type"].nunique()),
    )
    log_output_artifacts(logger, [report_path, merged_path])


def main() -> None:
    """Entry point for ``gatk-sv-ploidy eval``."""
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    with tool_logging_context(
        tool_name="eval",
        output_dir=args.output_dir,
        args=args,
    ) as logger:
        _run_eval(args, logger)


if __name__ == "__main__":
    main()
