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

from gatk_sv_ploidy._util import (
    compute_cnq_from_probabilities,
    summarize_contig_ploidy_from_bin_calls,
)

logger = logging.getLogger(__name__)

_SEX_CHROMS = {"chrX", "chrY"}
_CN_PROB_COLUMNS = [f"cn_prob_{cn}" for cn in range(6)]
_PLOIDY_PROB_COLUMNS = [f"ploidy_prob_{cn}" for cn in range(6)]
_BINQ_FIELD_OPTIONS = ["auto", "BINQ15", "BINQ20", "CALLQ15", "CALLQ20"]


def _resolve_binq_field(
    bin_quality_df: pd.DataFrame,
    requested_field: str,
) -> str:
    """Resolve a quality field, preferring posterior-predictive BINQ."""
    if requested_field != "auto":
        return requested_field
    if "BINQ20" in bin_quality_df.columns:
        return "BINQ20"
    return "CALLQ20"


def _median_absolute_deviation(values: pd.Series) -> float:
    """Compute the unscaled median absolute deviation."""
    arr = values.to_numpy(dtype=float)
    if arr.size == 0:
        return float("nan")
    median = np.median(arr)
    return float(np.median(np.abs(arr - median)))


def _weighted_mean(
    values: pd.Series,
    weights: pd.Series,
) -> float:
    """Compute a weighted mean while ignoring NaN values and zero weights."""
    val_arr = values.to_numpy(dtype=float)
    wt_arr = weights.to_numpy(dtype=float)
    mask = np.isfinite(val_arr) & np.isfinite(wt_arr) & (wt_arr > 0)
    if not np.any(mask):
        return float("nan")
    return float(np.average(val_arr[mask], weights=wt_arr[mask]))


def _safe_min_probability(values: pd.Series) -> float:
    """Return the minimum finite probability in a series, or NaN if absent."""
    finite = values[np.isfinite(values)]
    if finite.empty:
        return float("nan")
    return float(finite.min())


def _aggregate_bin_stats_to_chromosome_stats(bin_df: pd.DataFrame) -> pd.DataFrame:
    """Recompute chromosome-level summaries from retained bin-level statistics."""
    rows: list[dict] = []

    for (sample_id, chrom), sdf in bin_df.groupby(["sample", "chr"], sort=False):
        if "cn_map" in sdf.columns:
            chr_cn_map = sdf["cn_map"].to_numpy(dtype=np.int64)
        else:
            chr_cn_map = np.argmax(
                sdf[_CN_PROB_COLUMNS].to_numpy(dtype=np.float64),
                axis=1,
            )
        if "cnq" in sdf.columns:
            chr_cnq = sdf["cnq"].to_numpy(dtype=np.float64)
        else:
            chr_cnq = compute_cnq_from_probabilities(
                sdf[_CN_PROB_COLUMNS].to_numpy(dtype=np.float64),
            )
        copy_number, mean_cn_probability, ploidy_fractions, plq = summarize_contig_ploidy_from_bin_calls(
            chr_cn_map,
            chr_cnq,
            n_states=len(_CN_PROB_COLUMNS),
        )

        row = {
            "sample": sample_id,
            "chromosome": chrom,
            "copy_number": copy_number,
            "mean_cn_probability": mean_cn_probability,
            "plq": plq,
            "n_bins": int(len(sdf)),
            "mean_depth": float(sdf["observed_depth"].mean()),
            "std_depth": float(sdf["observed_depth"].std(ddof=0)),
            "median_depth": float(sdf["observed_depth"].median()),
            "mad_depth": _median_absolute_deviation(sdf["observed_depth"]),
        }
        for cn_state, ploidy_prob in enumerate(ploidy_fractions):
            row[f"ploidy_prob_{cn_state}"] = float(ploidy_prob)

        if "sample_var" in sdf.columns:
            row["sample_var_map"] = float(sdf["sample_var"].iloc[0])
        if "sample_depth" in sdf.columns:
            row["sample_depth_map"] = float(sdf["sample_depth"].iloc[0])
        if "plot_depth" in sdf.columns:
            row["plot_mean_depth"] = float(sdf["plot_depth"].mean())
            row["plot_std_depth"] = float(sdf["plot_depth"].std(ddof=0))
            row["plot_median_depth"] = float(sdf["plot_depth"].median())
            row["plot_mad_depth"] = _median_absolute_deviation(sdf["plot_depth"])
        if "sample_overdispersion" in sdf.columns:
            row["sample_overdispersion_map"] = float(
                sdf["sample_overdispersion"].iloc[0]
            )
        if "n_sites" in sdf.columns:
            row["n_sites"] = int(sdf["n_sites"].sum())
        if "mean_observed_af" in sdf.columns:
            if "n_sites" in sdf.columns:
                row["mean_observed_af"] = _weighted_mean(
                    sdf["mean_observed_af"],
                    sdf["n_sites"],
                )
            else:
                row["mean_observed_af"] = float(sdf["mean_observed_af"].mean())
        if "n_het_sites" in sdf.columns:
            row["n_het_sites"] = int(sdf["n_het_sites"].sum())
        if "mean_het_af" in sdf.columns:
            if "n_het_sites" in sdf.columns:
                row["mean_het_af"] = _weighted_mean(
                    sdf["mean_het_af"],
                    sdf["n_het_sites"],
                )
            else:
                row["mean_het_af"] = float(sdf["mean_het_af"].mean())
        if "af_log_lik" in sdf.columns:
            row["af_log_lik"] = float(sdf["af_log_lik"].sum())

        rows.append(row)

    return pd.DataFrame(rows)


def _annotate_aneuploidy_flags(
    df: pd.DataFrame,
    prob_threshold: float,
) -> pd.DataFrame:
    """Recompute chromosome-level aneuploid flags from chromosome summaries."""
    out = df.copy()
    out["is_aneuploid"] = False

    for sample_id, sdf in out.groupby("sample"):
        cn_map = dict(zip(sdf["chromosome"], sdf["copy_number"]))
        prob_map = dict(
            zip(
                sdf["chromosome"],
                sdf["mean_cn_probability"].fillna(0.0),
            )
        )

        auto_mask = ~sdf["chromosome"].isin(_SEX_CHROMS)
        auto_aneu = auto_mask & (sdf["copy_number"] != 2)
        auto_aneu &= sdf["mean_cn_probability"].fillna(0.0) > prob_threshold
        out.loc[sdf.index[auto_aneu], "is_aneuploid"] = True

        x_cn = int(cn_map.get("chrX", 2))
        y_cn = int(cn_map.get("chrY", 0))
        is_xx = x_cn == 2 and y_cn == 0
        is_xy = x_cn == 1 and y_cn == 1
        x_ok = ("chrX" not in prob_map) or (prob_map["chrX"] > prob_threshold)
        y_ok = ("chrY" not in prob_map) or (prob_map["chrY"] > prob_threshold)
        if not (is_xx or is_xy) and x_ok and y_ok:
            sex_mask = sdf["chromosome"].isin(_SEX_CHROMS)
            out.loc[sdf.index[sex_mask], "is_aneuploid"] = True

    return out


def _annotate_binq_filter_for_bin_stats(
    bin_stats_df: pd.DataFrame,
    bin_quality_df: pd.DataFrame,
    min_binq: float,
    binq_field: str,
) -> pd.DataFrame:
    """Annotate each sample-bin row with BINQ-based call filtering metadata."""
    binq_field = _resolve_binq_field(bin_quality_df, binq_field)
    required_quality_cols = {"chr", "start", "end", binq_field}
    missing_quality_cols = required_quality_cols - set(bin_quality_df.columns)
    if missing_quality_cols:
        missing_cols = ", ".join(sorted(missing_quality_cols))
        raise ValueError(
            f"PPD bin quality file is missing required columns: {missing_cols}"
        )

    quality_cols = ["chr", "start", "end", binq_field]
    out = bin_stats_df.merge(
        bin_quality_df[quality_cols],
        on=["chr", "start", "end"],
        how="left",
    )
    out["binq_field"] = binq_field
    out["binq_value"] = out[binq_field]
    out["min_binq"] = float(min_binq)
    out["ignored_in_call"] = out["binq_value"].fillna(-np.inf) < float(min_binq)
    out["retained_in_call"] = ~out["ignored_in_call"]
    return out


def summarize_binq_filter_stats(annotated_bin_df: pd.DataFrame) -> pd.DataFrame:
    """Summarize BINQ filtering over unique bins for each contig."""
    unique_bins = (
        annotated_bin_df.loc[:, ["chr", "start", "end", "ignored_in_call"]]
        .drop_duplicates()
        .copy()
    )
    if unique_bins.empty:
        return pd.DataFrame(
            columns=[
                "chromosome",
                "n_bins",
                "n_bins_filtered",
                "pct_coverage_remaining",
                "coverage_remaining_mb",
            ]
        )

    unique_bins["bin_width_bp"] = (
        unique_bins["end"].to_numpy(dtype=float) -
        unique_bins["start"].to_numpy(dtype=float)
    )
    unique_bins["retained_bp"] = np.where(
        unique_bins["ignored_in_call"].to_numpy(dtype=bool),
        0.0,
        unique_bins["bin_width_bp"].to_numpy(dtype=float),
    )

    summary = (
        unique_bins.groupby("chr", sort=False, as_index=False)
        .agg(
            n_bins=("ignored_in_call", "size"),
            n_bins_filtered=("ignored_in_call", "sum"),
            total_bp=("bin_width_bp", "sum"),
            retained_bp=("retained_bp", "sum"),
        )
        .rename(columns={"chr": "chromosome"})
    )
    total_bp = summary["total_bp"].to_numpy(dtype=float)
    retained_bp = summary["retained_bp"].to_numpy(dtype=float)
    summary["pct_coverage_remaining"] = np.divide(
        retained_bp,
        total_bp,
        out=np.zeros(len(summary), dtype=float),
        where=total_bp > 0,
    ) * 100.0
    summary["coverage_remaining_mb"] = retained_bp / 1e6
    return summary[
        [
            "chromosome",
            "n_bins",
            "n_bins_filtered",
            "pct_coverage_remaining",
            "coverage_remaining_mb",
        ]
    ]


def log_binq_filter_stats(annotated_bin_df: pd.DataFrame) -> None:
    """Log per-contig BINQ filtering statistics over unique bins."""
    summary = summarize_binq_filter_stats(annotated_bin_df)
    if summary.empty:
        logger.info("No bins were available for BINQ filtering summary.")
        return

    logger.info(
        "BINQ filtering summary by contig:"
    )
    for row in summary.itertuples(index=False):
        logger.info(
            "%s: bins=%d filtered=%d coverage_remaining=%.2f%% remaining_mb=%.3f",
            row.chromosome,
            int(row.n_bins),
            int(row.n_bins_filtered),
            float(row.pct_coverage_remaining),
            float(row.coverage_remaining_mb),
        )


def _apply_binq_filter_to_annotated_bins(
    chrom_df: pd.DataFrame,
    annotated_bin_df: pd.DataFrame,
    prob_threshold: float = 0.5,
) -> pd.DataFrame:
    """Recompute chromosome summaries from bin stats annotated with call filters."""
    retained_bin_df = annotated_bin_df[annotated_bin_df["retained_in_call"]].copy()
    retained_chr_df = _aggregate_bin_stats_to_chromosome_stats(retained_bin_df)
    retained_chr_df = retained_chr_df.rename(
        columns={
            col: f"{col}_filtered"
            for col in retained_chr_df.columns
            if col not in {"sample", "chromosome"}
        }
    )

    out = chrom_df.copy()
    if "n_bins" in out.columns:
        out["n_bins_original"] = out["n_bins"]
    out = out.merge(retained_chr_df, on=["sample", "chromosome"], how="left")

    retained_mask = out.get(
        "n_bins_filtered",
        pd.Series(index=out.index, dtype=float),
    ).notna()
    sex_conflict_mask = (
        retained_mask &
        out["chromosome"].isin(_SEX_CHROMS) &
        out.get(
            "copy_number_filtered",
            pd.Series(index=out.index, dtype=float),
        ).notna() &
        (
            out["copy_number"].to_numpy(dtype=float) !=
            out["copy_number_filtered"].to_numpy(dtype=float)
        )
    )
    replacement_mask = retained_mask & ~sex_conflict_mask
    replacement_pairs = [
        ("copy_number", "copy_number_filtered"),
        ("mean_cn_probability", "mean_cn_probability_filtered"),
        ("plq", "plq_filtered"),
        ("n_bins", "n_bins_filtered"),
        ("mean_depth", "mean_depth_filtered"),
        ("std_depth", "std_depth_filtered"),
        ("median_depth", "median_depth_filtered"),
        ("mad_depth", "mad_depth_filtered"),
        ("sample_var_map", "sample_var_map_filtered"),
        ("sample_depth_map", "sample_depth_map_filtered"),
        ("plot_mean_depth", "plot_mean_depth_filtered"),
        ("plot_std_depth", "plot_std_depth_filtered"),
        ("plot_median_depth", "plot_median_depth_filtered"),
        ("plot_mad_depth", "plot_mad_depth_filtered"),
        ("sample_overdispersion_map", "sample_overdispersion_map_filtered"),
        ("n_sites", "n_sites_filtered"),
        ("mean_observed_af", "mean_observed_af_filtered"),
        ("n_het_sites", "n_het_sites_filtered"),
        ("mean_het_af", "mean_het_af_filtered"),
        ("af_log_lik", "af_log_lik_filtered"),
    ]
    replacement_pairs.extend(
        (col, f"{col}_filtered") for col in _PLOIDY_PROB_COLUMNS
    )
    for base_col, filtered_col in replacement_pairs:
        if filtered_col in out.columns:
            out.loc[replacement_mask, base_col] = out.loc[replacement_mask, filtered_col]

    out["n_bins_retained"] = out.get(
        "n_bins_filtered",
        pd.Series(np.zeros(len(out)), index=out.index),
    ).fillna(0).astype(int)
    if "n_bins_original" in out.columns:
        denom = out["n_bins_original"].fillna(0).to_numpy(dtype=float)
        out["frac_bins_retained"] = np.divide(
            out["n_bins_retained"].to_numpy(dtype=float),
            denom,
            out=np.zeros(len(out), dtype=float),
            where=denom > 0,
        )

    missing_auto_mask = (
        (out["n_bins_retained"] == 0) &
        ~out["chromosome"].isin(_SEX_CHROMS)
    )
    out.loc[missing_auto_mask, "copy_number"] = 2
    out.loc[missing_auto_mask, "mean_cn_probability"] = np.nan
    if "plq" in out.columns:
        out.loc[missing_auto_mask, "plq"] = np.nan
    for col in _PLOIDY_PROB_COLUMNS:
        if col in out.columns:
            out.loc[missing_auto_mask, col] = np.nan
    out.loc[missing_auto_mask, "is_aneuploid"] = False
    for col in [
        "n_bins",
        "mean_depth",
        "std_depth",
        "median_depth",
        "mad_depth",
        "plot_mean_depth",
        "plot_std_depth",
        "plot_median_depth",
        "plot_mad_depth",
        "n_sites",
        "mean_observed_af",
        "n_het_sites",
        "mean_het_af",
        "af_log_lik",
    ]:
        if col in out.columns:
            out.loc[missing_auto_mask, col] = np.nan

    out = _annotate_aneuploidy_flags(out, prob_threshold=prob_threshold)

    filtered_cols = [col for col in out.columns if col.endswith("_filtered")]
    return out.drop(columns=filtered_cols, errors="ignore").sort_values(
        ["sample", "chromosome"]
    )


def apply_binq_filter_to_chromosome_stats(
    chrom_df: pd.DataFrame,
    bin_stats_df: pd.DataFrame,
    bin_quality_df: pd.DataFrame,
    min_binq: float,
    binq_field: str = "BINQ20",
    prob_threshold: float = 0.5,
) -> pd.DataFrame:
    """Recompute chromosome summaries after removing low-BINQ bins.

    Autosomes with no retained bins are neutralized so they no longer influence
    ploidy calling. Sex chromosomes fall back to the original chromosome-level
    summaries when fully removed, and partial filtered subsets are not allowed
    to overturn the original chrX/chrY copy number on their own.
    """
    annotated_bin_df = _annotate_binq_filter_for_bin_stats(
        bin_stats_df,
        bin_quality_df,
        min_binq=min_binq,
        binq_field=binq_field,
    )
    return _apply_binq_filter_to_annotated_bins(
        chrom_df,
        annotated_bin_df,
        prob_threshold=prob_threshold,
    )


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
                "score": _safe_min_probability(sdf["mean_cn_probability"]),
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


def export_ignored_bins(df: pd.DataFrame, output_path: str) -> None:
    """Export sample-bin rows that were ignored during BINQ-filtered calling."""
    cols = [
        "sample",
        "chr",
        "start",
        "end",
        "ignored_in_call",
        "binq_field",
        "binq_value",
        "min_binq",
    ]
    optional_cols = ["observed_depth", "plot_depth", "cn_map", "cnq"]
    use_cols = [col for col in cols + optional_cols if col in df.columns]
    out = df[df["ignored_in_call"]].copy()
    out = out[use_cols]
    out.to_csv(output_path, sep="\t", index=False, compression="gzip")
    n_unique_bins = len(out[["chr", "start", "end"]].drop_duplicates()) if len(out) else 0
    logger.info(
        "Exported %d ignored sample-bin rows across %d bins to %s",
        len(out),
        n_unique_bins,
        output_path,
    )


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
    p.add_argument(
        "--bin-stats", default=None,
        help="Optional infer/bin_stats.tsv.gz used to rebuild chromosome calls after BINQ filtering.",
    )
    p.add_argument(
        "--ppd-bin-quality", default=None,
        help="Optional ppd/ppd_bin_quality.tsv used for BINQ-based bin filtering.",
    )
    p.add_argument(
        "--binq-field", choices=_BINQ_FIELD_OPTIONS, default="BINQ20",
        help="Per-bin quality field to threshold when --min-binq is provided. "
             "'auto' prefers CALLQ20 when available, otherwise BINQ20.",
    )
    p.add_argument(
        "--min-binq", type=float, default=None,
        help="If set, rebuild chromosome calls after excluding bins with BINQ below this threshold.",
    )
    p.add_argument(
        "--prob-threshold", type=float, default=0.5,
        help="Minimum mean chromosome CN probability used when rebuilding aneuploid flags from filtered bins.",
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

    annotated_bin_df: pd.DataFrame | None = None

    if args.min_binq is not None:
        if not args.bin_stats or not args.ppd_bin_quality:
            raise ValueError(
                "--min-binq requires both --bin-stats and --ppd-bin-quality."
            )
        logger.info(
            "Applying %s >= %.2f filter using %s and rebuilding chromosome calls ...",
            args.binq_field,
            args.min_binq,
            args.ppd_bin_quality,
        )
        bin_stats_df = pd.read_csv(args.bin_stats, sep="\t")
        bin_quality_df = pd.read_csv(args.ppd_bin_quality, sep="\t")
        annotated_bin_df = _annotate_binq_filter_for_bin_stats(
            bin_stats_df,
            bin_quality_df,
            min_binq=args.min_binq,
            binq_field=args.binq_field,
        )
        log_binq_filter_stats(annotated_bin_df)
        df = _apply_binq_filter_to_annotated_bins(
            df,
            annotated_bin_df,
            prob_threshold=args.prob_threshold,
        )
        filtered_path = os.path.join(args.output_dir, "chromosome_stats.filtered.tsv")
        df.to_csv(filtered_path, sep="\t", index=False)
        logger.info("Filtered chromosome stats saved to %s", filtered_path)

    # ── optional exclusion ──────────────────────────────────────────────
    if args.exclusion_list:
        from gatk_sv_ploidy._util import load_exclusion_ids

        excl = set(load_exclusion_ids(args.exclusion_list))
        n_before = df["sample"].nunique()
        df = df[~df["sample"].isin(excl)]
        if annotated_bin_df is not None:
            annotated_bin_df = annotated_bin_df[~annotated_bin_df["sample"].isin(excl)].copy()
        logger.info(
            "Excluded %d samples → %d remaining",
            n_before - df["sample"].nunique(),
            df["sample"].nunique(),
        )

    if annotated_bin_df is not None:
        ignored_path = os.path.join(args.output_dir, "ignored_bins.tsv.gz")
        export_ignored_bins(annotated_bin_df, ignored_path)

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
