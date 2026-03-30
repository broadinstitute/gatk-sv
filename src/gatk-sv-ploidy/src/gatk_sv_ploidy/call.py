"""
Call subcommand — assign sex karyotype and aneuploidy type.

Reads ``chromosome_stats.tsv`` (output of ``infer``), classifies each sample's
sex from chrX/chrY copy numbers, predicts an aneuploidy type, and writes
sex-assignment and aneuploid-sample tables.
"""

from __future__ import annotations

import argparse
import logging
import os
from typing import Dict

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


# ── classification ──────────────────────────────────────────────────────────


def assign_sex_and_aneuploidy_types(
    df: pd.DataFrame,
    truth_dict: Dict[str, str] | None = None,
) -> pd.DataFrame:
    """Classify each sample's sex and aneuploidy type.

    Sex is determined from chrX/chrY integer copy numbers:

    ========== ======== ========
    Sex         chrX_CN  chrY_CN
    ========== ======== ========
    MALE         1        1
    FEMALE       2        0
    TURNER       1        0
    TRIPLE_X     3        0
    KLINEFELTER  2        1
    JACOBS       1        2
    OTHER       *else*   *else*
    ========== ======== ========

    Aneuploidy type is inferred from the set of aneuploid chromosomes.

    Args:
        df: ``chromosome_stats.tsv`` DataFrame with columns ``sample``,
            ``chromosome``, ``copy_number``, ``is_aneuploid``,
            ``mean_cn_probability``, ``median_depth``.
        truth_dict: Optional mapping of sample ID → true aneuploidy type.

    Returns:
        One-row-per-sample DataFrame with columns ``sample``, ``sex``,
        ``chrX_CN``, ``chrY_CN``, ``chrX_depth``, ``chrY_depth``,
        ``true_aneuploidy_type``, ``predicted_aneuploidy_type``, ``score``.
    """
    if truth_dict is None:
        truth_dict = {}

    rows: list[dict] = []

    for sample_id, sdf in df.groupby("sample"):
        true_type = truth_dict.get(str(sample_id), "NORMAL")

        cn_map = dict(zip(sdf["chromosome"], sdf["copy_number"]))
        depth_map = dict(zip(sdf["chromosome"], sdf["median_depth"]))
        aneu_map = dict(zip(sdf["chromosome"], sdf["is_aneuploid"]))

        x_cn = cn_map.get("chrX", 2)
        y_cn = cn_map.get("chrY", 0)
        x_depth = depth_map.get("chrX", 2.0)
        y_depth = depth_map.get("chrY", 0.0)

        # ── sex assignment ──────────────────────────────────────────
        sex = _classify_sex(x_cn, y_cn)

        # ── aneuploidy type ─────────────────────────────────────────
        pred_type = _classify_aneuploidy(aneu_map, cn_map, x_cn, y_cn)

        rows.append(
            {
                "sample": sample_id,
                "sex": sex,
                "chrX_CN": x_cn,
                "chrY_CN": y_cn,
                "chrX_depth": x_depth,
                "chrY_depth": y_depth,
                "true_aneuploidy_type": true_type,
                "predicted_aneuploidy_type": pred_type,
                "score": float(sdf["mean_cn_probability"].min()),
            }
        )

    return pd.DataFrame(rows)


def _classify_sex(x_cn: int, y_cn: int) -> str:
    """Map chrX/chrY copy numbers to a sex-karyotype label."""
    _SEX_TABLE = {
        (1, 1): "MALE",
        (2, 0): "FEMALE",
        (1, 0): "TURNER",
        (3, 0): "TRIPLE_X",
        (2, 1): "KLINEFELTER",
        (1, 2): "JACOBS",
        # Tetraploid patterns
        (4, 0): "TETRAPLOID_FEMALE",   # XXXX
        (2, 2): "TETRAPLOID_MALE",     # XXYY
        (3, 1): "TRIPLE_X_Y",          # XXXY
    }
    return _SEX_TABLE.get((x_cn, y_cn), "OTHER")


def _classify_aneuploidy(
    aneu_map: dict,
    cn_map: dict,
    x_cn: int,
    y_cn: int,
) -> str:
    """Determine predicted aneuploidy type from per-chromosome calls."""
    aneuploid_chrs = [c for c, v in aneu_map.items() if v]

    if not aneuploid_chrs:
        return "NORMAL"

    sex_aneu = [c for c in aneuploid_chrs if c in ("chrX", "chrY")]
    auto_aneu = [c for c in aneuploid_chrs if c not in ("chrX", "chrY")]

    # Check for genome-wide tetraploidy: ≥80% of autosomes at CN=4
    autosome_cns = {
        c: cn_map.get(c, 2) for c in cn_map if c not in ("chrX", "chrY")
    }
    n_tetra_auto = sum(1 for cn in autosome_cns.values() if cn == 4)
    if autosome_cns and n_tetra_auto >= len(autosome_cns) * 0.8:
        return "TETRAPLOID"

    # Multiple autosomal or mixed → MULTIPLE
    if len(auto_aneu) > 1 or (auto_aneu and sex_aneu):
        return "MULTIPLE"

    # Sex-only aneuploidies
    if sex_aneu and not auto_aneu:
        _sex_type = {
            (2, 1): "KLINEFELTER",
            (3, 0): "TRIPLE_X",
            (1, 0): "TURNER",
            (1, 2): "JACOBS",
            (4, 0): "TETRAPLOID_FEMALE",
            (2, 2): "TETRAPLOID_MALE",
            (3, 1): "TRIPLE_X_Y",
        }
        return _sex_type.get((x_cn, y_cn), "OTHER")

    # Single autosomal aneuploidy
    if len(auto_aneu) == 1:
        chrom = auto_aneu[0]
        cn = cn_map[chrom]
        _auto_type = {
            ("chr13", 3): "TRISOMY_13",
            ("chr18", 3): "TRISOMY_18",
            ("chr21", 3): "TRISOMY_21",
            ("chr13", 4): "TETRASOMY_13",
            ("chr18", 4): "TETRASOMY_18",
            ("chr21", 4): "TETRASOMY_21",
        }
        return _auto_type.get((chrom, cn), "OTHER")

    return "OTHER"


# ── output helpers ──────────────────────────────────────────────────────────


def save_sex_assignments(pred_df: pd.DataFrame, output_dir: str) -> None:
    """Write ``sex_assignments.txt.gz``.

    Args:
        pred_df: DataFrame returned by :func:`assign_sex_and_aneuploidy_types`.
        output_dir: Destination directory.
    """
    out = pred_df[["sample", "chrX_CN", "chrY_CN", "sex"]].copy()
    out.columns = ["sample_id", "chrX.CN", "chrY.CN", "Assignment"]
    path = os.path.join(output_dir, "sex_assignments.txt.gz")
    out.to_csv(path, sep="\t", index=False, compression="gzip")
    logger.info("Sex assignments saved to %s", path)


def export_aneuploid_data(df: pd.DataFrame, output_path: str) -> None:
    """Export rows with ``is_aneuploid=True`` to a TSV.

    Converts ``sample_var_map`` to ``inferred_sample_std`` and drops internal
    columns before writing.

    Args:
        df: ``chromosome_stats.tsv`` DataFrame.
        output_path: Destination file path.
    """
    out = df[df["is_aneuploid"]].copy()
    if "sample_var_map" in out.columns:
        out["inferred_sample_std"] = np.sqrt(out["sample_var_map"])
        out = out.drop(columns=["sample_var_map"], errors="ignore")
    out = out.drop(columns=["chr_type"], errors="ignore")
    out.to_csv(output_path, sep="\t", index=False)
    logger.info("Exported %d aneuploid rows to %s", len(out), output_path)


# ── CLI ─────────────────────────────────────────────────────────────────────


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the call subcommand."""
    p = argparse.ArgumentParser(
        description="Assign sex karyotype and aneuploidy type",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "-c", "--chrom-stats", required=True,
        help="Input chromosome_stats.tsv (from 'infer')",
    )
    p.add_argument(
        "-o", "--output-dir", required=True,
        help="Output directory",
    )
    p.add_argument(
        "-t", "--truth-json", default=None,
        help="Optional JSON mapping sample ID → true aneuploidy type",
    )
    p.add_argument(
        "-e", "--exclusion-list", default=None,
        help="Text file with sample IDs to exclude (one per line)",
    )
    return p.parse_args()


def main() -> None:
    """Entry point for ``gatk-sv-ploidy call``."""
    import json

    logging.basicConfig(level=logging.INFO, format="%(message)s")
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    # ── load chromosome stats ───────────────────────────────────────────
    logger.info("Reading chromosome stats: %s", args.chrom_stats)
    df = pd.read_csv(args.chrom_stats, sep="\t")
    logger.info(
        "Loaded %d rows for %d samples", len(df), df["sample"].nunique()
    )

    # ── optional exclusion ──────────────────────────────────────────────
    if args.exclusion_list:
        from gatk_sv_ploidy._util import load_exclusion_ids

        excl = set(load_exclusion_ids(args.exclusion_list))
        n_before = df["sample"].nunique()
        df = df[~df["sample"].isin(excl)]
        logger.info(
            "Excluded %d samples → %d remaining",
            n_before - df["sample"].nunique(),
            df["sample"].nunique(),
        )

    # ── optional truth ──────────────────────────────────────────────────
    truth: Dict[str, str] = {}
    if args.truth_json:
        with open(args.truth_json) as fh:
            truth = json.load(fh)
        logger.info("Loaded %d truth entries", len(truth))

    # ── classify ────────────────────────────────────────────────────────
    pred_df = assign_sex_and_aneuploidy_types(df, truth)
    save_sex_assignments(pred_df, args.output_dir)

    # ── save predictions ────────────────────────────────────────────────
    pred_path = os.path.join(args.output_dir, "aneuploidy_type_predictions.tsv")
    pred_df.to_csv(pred_path, sep="\t", index=False)
    logger.info("Predictions saved to %s", pred_path)

    # ── export aneuploid data ───────────────────────────────────────────
    aneu_path = os.path.join(args.output_dir, "aneuploid_samples.tsv")
    export_aneuploid_data(df, aneu_path)

    logger.info("Done.")


if __name__ == "__main__":
    main()
