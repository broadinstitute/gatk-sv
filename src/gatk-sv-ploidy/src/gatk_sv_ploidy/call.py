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

from gatk_sv_ploidy._logging import log_output_artifacts, tool_logging_context
from gatk_sv_ploidy._util import (
    BINQ_FIELD_OPTIONS,
    baseline_ploidy_label,
    compute_cnq_from_probabilities,
    is_expected_allosome_copy_number_pair,
    resolve_binq_field,
    summarize_contig_ploidy_from_bin_calls,
)

_SEX_CHROMS = {"chrX", "chrY"}
_CN_PROB_COLUMNS = [f"cn_prob_{cn}" for cn in range(6)]
_PLOIDY_PROB_COLUMNS = [f"ploidy_prob_{cn}" for cn in range(6)]
_RESERVED_POLYPLOID_ANEUPLOIDY_TYPES = ("HAPLOID", "TRIPLOID", "TETRAPLOID")
_NO_ANEUPLOIDY = "NONE"
_POLYPLOID_SEX_LABELS_BY_BASELINE = {
    1: {
        (1, 0): "HAPLOID_X",
        (0, 1): "HAPLOID_Y",
    },
    3: {
        (3, 0): "TRIPLOID_FEMALE",
        (2, 1): "TRIPLOID_MALE",
        (1, 2): "TRIPLOID_XYY",
    },
    4: {
        (4, 0): "TETRAPLOID_FEMALE",
        (2, 2): "TETRAPLOID_MALE",
        (3, 1): "TETRAPLOID_XXXY",
        (1, 3): "TETRAPLOID_XYYY",
    },
}
LOGGER = logging.getLogger(__name__)


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


def _first_finite(values: pd.Series) -> float:
    """Return the first finite value in a series, or NaN if none exist."""
    arr = values.to_numpy(dtype=float)
    finite = arr[np.isfinite(arr)]
    if finite.size == 0:
        return float("nan")
    return float(finite[0])


def _compute_sample_depth_ratios(df: pd.DataFrame) -> tuple[dict[str, float], float]:
    """Estimate cohort-relative sample-depth ratios from chromosome summaries."""
    if "sample_depth_map" not in df.columns:
        return {}, float("nan")

    sample_depths = df.groupby("sample", sort=False)["sample_depth_map"].apply(_first_finite)
    finite = sample_depths[np.isfinite(sample_depths) & (sample_depths > 0.0)]
    if finite.empty:
        return {}, float("nan")

    reference_depth = float(np.median(finite.to_numpy(dtype=float)))
    if not np.isfinite(reference_depth) or reference_depth <= 0.0:
        return {}, float("nan")

    depth_ratios = {
        str(sample_id): float(depth / reference_depth)
        for sample_id, depth in sample_depths.items()
        if np.isfinite(depth) and depth > 0.0
    }
    return depth_ratios, reference_depth


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

        if "sample_overdispersion" in sdf.columns:
            row["sample_overdispersion_map"] = float(
                sdf["sample_overdispersion"].iloc[0]
            )
        elif "sample_var" in sdf.columns:
            row["sample_overdispersion_map"] = float(sdf["sample_var"].iloc[0])
        if "sample_depth" in sdf.columns:
            row["sample_depth_map"] = float(sdf["sample_depth"].iloc[0])
        if "plot_depth" in sdf.columns:
            row["plot_mean_depth"] = float(sdf["plot_depth"].mean())
            row["plot_std_depth"] = float(sdf["plot_depth"].std(ddof=0))
            row["plot_median_depth"] = float(sdf["plot_depth"].median())
            row["plot_mad_depth"] = _median_absolute_deviation(sdf["plot_depth"])
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

    if not rows:
        return pd.DataFrame(columns=["sample", "chromosome"])

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
        autosomal_baseline_cn = 2
        if "autosomal_baseline_cn" in sdf.columns:
            baseline_values = sdf["autosomal_baseline_cn"].dropna()
            if not baseline_values.empty:
                autosomal_baseline_cn = int(baseline_values.iloc[0])

        auto_mask = ~sdf["chromosome"].isin(_SEX_CHROMS)
        auto_aneu = auto_mask & (sdf["copy_number"] != autosomal_baseline_cn)
        auto_aneu &= sdf["mean_cn_probability"].fillna(0.0) > prob_threshold
        out.loc[sdf.index[auto_aneu], "is_aneuploid"] = True

        x_cn = int(cn_map.get("chrX", 2))
        y_cn = int(cn_map.get("chrY", 0))
        is_expected_sex_pair = is_expected_allosome_copy_number_pair(
            x_cn,
            y_cn,
            autosomal_baseline_cn,
        )
        x_ok = ("chrX" not in prob_map) or (prob_map["chrX"] > prob_threshold)
        y_ok = ("chrY" not in prob_map) or (prob_map["chrY"] > prob_threshold)
        if not is_expected_sex_pair and x_ok and y_ok:
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
    binq_field = resolve_binq_field(bin_quality_df, binq_field)
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
    summary_df = summarize_binq_filter_stats(annotated_bin_df)
    if summary_df.empty:
        LOGGER.info("BINQ filter summary: no bins available")
        return
    total_bins = int(summary_df["n_bins"].sum())
    filtered_bins = int(summary_df["n_bins_filtered"].sum())
    retained_pct = summary_df["pct_coverage_remaining"].to_numpy(dtype=float)
    LOGGER.info(
        "BINQ filter summary: contigs=%d bins=%d filtered_bins=%d "
        "median_coverage_remaining_pct=%.2f min_coverage_remaining_pct=%.2f",
        int(len(summary_df)),
        total_bins,
        filtered_bins,
        float(np.nanmedian(retained_pct)),
        float(np.nanmin(retained_pct)),
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
    filtered_copy_number = out.get(
        "copy_number_filtered",
        pd.Series(index=out.index, dtype=float),
    )
    sex_conflict_mask = (
        retained_mask &
        out["chromosome"].isin(_SEX_CHROMS) &
        filtered_copy_number.notna() &
        (
            out["copy_number"].to_numpy(dtype=float) !=
            filtered_copy_number.to_numpy(dtype=float)
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

    sample_depth_ratios, depth_reference = _compute_sample_depth_ratios(df)
    rows: list[dict] = []

    for sample_id, sdf in df.groupby("sample"):
        true_type = truth_dict.get(str(sample_id), "NORMAL")

        raw_cn_map = dict(zip(sdf["chromosome"], sdf["copy_number"]))
        depth_map = dict(zip(sdf["chromosome"], sdf["median_depth"]))
        autosomal_baseline_cn = 2
        if "autosomal_baseline_cn" in sdf.columns:
            baseline_values = sdf["autosomal_baseline_cn"].dropna()
            if not baseline_values.empty:
                autosomal_baseline_cn = int(baseline_values.iloc[0])
        sample_depth_ratio = sample_depth_ratios.get(str(sample_id), float("nan"))
        cn_map = dict(raw_cn_map)
        cn_scale_factor = 1.0
        aneu_map = {
            chrom: bool(cn != autosomal_baseline_cn)
            for chrom, cn in cn_map.items()
            if chrom not in _SEX_CHROMS
        }

        x_cn = cn_map.get("chrX", 2)
        y_cn = cn_map.get("chrY", 0)
        raw_x_cn = raw_cn_map.get("chrX", 2)
        raw_y_cn = raw_cn_map.get("chrY", 0)
        x_depth = depth_map.get("chrX", 2.0)
        y_depth = depth_map.get("chrY", 0.0)

        # ── sex assignment ──────────────────────────────────────────
        sex = _classify_sex(
            x_cn,
            y_cn,
            autosomal_baseline_cn=autosomal_baseline_cn,
        )

        baseline_ploidy_type = _classify_baseline_ploidy(
            autosomal_baseline_cn,
        )
        autosomal_aneuploidy_type = _classify_autosomal_aneuploidy(
            aneu_map,
            cn_map,
            autosomal_baseline_cn=autosomal_baseline_cn,
        )
        allosomal_aneuploidy_type = _classify_allosomal_aneuploidy(
            x_cn,
            y_cn,
            autosomal_baseline_cn=autosomal_baseline_cn,
        )
        pred_type = _compose_aneuploidy_type(
            baseline_ploidy_type,
            autosomal_aneuploidy_type,
            allosomal_aneuploidy_type,
        )

        rows.append(
            {
                "sample": sample_id,
                "sex": sex,
                "chrX_CN": x_cn,
                "chrY_CN": y_cn,
                "raw_chrX_CN": raw_x_cn,
                "raw_chrY_CN": raw_y_cn,
                "chrX_depth": x_depth,
                "chrY_depth": y_depth,
                "autosomal_baseline_cn": autosomal_baseline_cn,
                "sample_depth_ratio": sample_depth_ratio,
                "sample_depth_reference": depth_reference,
                "global_cn_scale_factor": cn_scale_factor,
                "true_aneuploidy_type": true_type,
                "baseline_ploidy_type": baseline_ploidy_type,
                "autosomal_aneuploidy_type": autosomal_aneuploidy_type,
                "allosomal_aneuploidy_type": allosomal_aneuploidy_type,
                "predicted_aneuploidy_type": pred_type,
                "score": _safe_min_probability(sdf["mean_cn_probability"]),
            }
        )

    return pd.DataFrame(rows)


def _classify_sex(
    x_cn: int,
    y_cn: int,
    autosomal_baseline_cn: int = 2,
) -> str:
    """Map chrX/chrY copy numbers to a sex-karyotype label."""
    if int(autosomal_baseline_cn) != 2:
        labels = _POLYPLOID_SEX_LABELS_BY_BASELINE.get(
            int(autosomal_baseline_cn),
            {},
        )
        return labels.get((x_cn, y_cn), "OTHER")

    _SEX_TABLE = {
        (1, 1): "MALE",
        (2, 0): "FEMALE",
        (1, 0): "TURNER",
        (3, 0): "TRIPLE_X",
        (2, 1): "KLINEFELTER",
        (1, 2): "JACOBS",
    }
    return _SEX_TABLE.get((x_cn, y_cn), "OTHER")


def _classify_baseline_ploidy(autosomal_baseline_cn: int = 2) -> str:
    """Map autosomal baseline copy number to a baseline ploidy label."""
    return baseline_ploidy_label(autosomal_baseline_cn)


def _autosomal_aneuploid_chromosomes(aneu_map: dict) -> list[str]:
    """Return autosomes marked as aneuploid in input order."""
    return [
        chrom for chrom, is_aneuploid in aneu_map.items()
        if is_aneuploid and chrom not in _SEX_CHROMS
    ]


def _classify_autosomal_aneuploidy(
    aneu_map: dict,
    cn_map: dict,
    autosomal_baseline_cn: int = 2,
) -> str:
    """Classify autosomal aneuploidy relative to the baseline CN."""
    auto_aneu = _autosomal_aneuploid_chromosomes(aneu_map)
    if not auto_aneu:
        return _NO_ANEUPLOIDY
    if len(auto_aneu) > 1:
        return "MULTIPLE_AUTOSOMAL"

    chrom = auto_aneu[0]
    cn = int(cn_map[chrom])
    baseline = int(autosomal_baseline_cn)
    if baseline == 2:
        _auto_type = {
            ("chr13", 3): "TRISOMY_13",
            ("chr18", 3): "TRISOMY_18",
            ("chr21", 3): "TRISOMY_21",
            ("chr13", 4): "TETRASOMY_13",
            ("chr18", 4): "TETRASOMY_18",
            ("chr21", 4): "TETRASOMY_21",
        }
        known_type = _auto_type.get((chrom, cn))
        if known_type is not None:
            return known_type

    if cn > baseline:
        return "AUTOSOMAL_GAIN"
    if cn < baseline:
        return "AUTOSOMAL_LOSS"
    return _NO_ANEUPLOIDY


def _classify_allosomal_aneuploidy(
    x_cn: int,
    y_cn: int,
    autosomal_baseline_cn: int = 2,
) -> str:
    """Classify chrX/chrY aneuploidy relative to the baseline CN."""
    if is_expected_allosome_copy_number_pair(
        x_cn,
        y_cn,
        autosomal_baseline_cn,
    ):
        return _NO_ANEUPLOIDY

    if int(autosomal_baseline_cn) == 2:
        _sex_type = {
            (2, 1): "KLINEFELTER",
            (3, 0): "TRIPLE_X",
            (1, 0): "TURNER",
            (1, 2): "JACOBS",
        }
        return _sex_type.get((x_cn, y_cn), "ALLOSOME_ANEUPLOIDY")

    return "DISCORDANT_ALLOSOME_CN"


def _compose_aneuploidy_type(
    baseline_ploidy_type: str,
    autosomal_aneuploidy_type: str,
    allosomal_aneuploidy_type: str,
) -> str:
    """Compose the summary aneuploidy label."""
    has_auto = autosomal_aneuploidy_type != _NO_ANEUPLOIDY
    has_allosome = allosomal_aneuploidy_type != _NO_ANEUPLOIDY

    if baseline_ploidy_type == "DIPLOID":
        if not has_auto and not has_allosome:
            return "NORMAL"
        if has_auto and has_allosome:
            return "MULTIPLE"
        if has_allosome:
            return allosomal_aneuploidy_type
        if autosomal_aneuploidy_type == "MULTIPLE_AUTOSOMAL":
            return "MULTIPLE"
        if autosomal_aneuploidy_type in {"AUTOSOMAL_GAIN", "AUTOSOMAL_LOSS"}:
            return "OTHER"
        return autosomal_aneuploidy_type

    if not has_auto and not has_allosome:
        return baseline_ploidy_type
    if has_auto and has_allosome:
        return f"{baseline_ploidy_type}_WITH_MULTIPLE_ANEUPLOIDY"
    if has_allosome:
        return f"{baseline_ploidy_type}_WITH_ALLOSOME_ANEUPLOIDY"
    return f"{baseline_ploidy_type}_WITH_AUTOSOMAL_ANEUPLOIDY"


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


def export_aneuploid_data(df: pd.DataFrame, output_path: str) -> None:
    """Export rows with ``is_aneuploid=True`` to a TSV.

    Drops internal columns before writing.

    Args:
        df: ``chromosome_stats.tsv`` DataFrame.
        output_path: Destination file path.
    """
    out = df[df["is_aneuploid"]].copy()
    out = out.drop(columns=["chr_type"], errors="ignore")
    out.to_csv(output_path, sep="\t", index=False)


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
        help="Optional ppd/ppd_bin_quality.tsv used for per-bin quality filtering.",
    )
    p.add_argument(
        "--use-callq20",
        action="store_true",
        help="Shortcut for using CALLQ20 instead of the default BINQ20 when --min-binq is provided.",
    )
    p.add_argument(
        "--binq-field", choices=BINQ_FIELD_OPTIONS, default="BINQ20",
        help="Per-bin quality field to threshold when --min-binq is provided. "
             "Defaults to BINQ20. 'auto' prefers BINQ20 and falls back to CALLQ20.",
    )
    p.add_argument(
        "--min-binq", type=float, default=None,
        help="If set, rebuild chromosome calls after excluding bins with the selected quality field below this threshold. By default this uses BINQ20; use --use-callq20 to filter on CALLQ20 instead.",
    )
    p.add_argument(
        "--prob-threshold", type=float, default=0.5,
        help="Minimum mean chromosome CN probability used when rebuilding aneuploid flags from filtered bins.",
    )
    return p.parse_args()


def _run_call(args: argparse.Namespace, logger) -> None:
    """Run call assignment after CLI logging is configured."""
    import json

    # ── load chromosome stats ───────────────────────────────────────────
    logger.info("Loading chromosome statistics")
    df = pd.read_csv(args.chrom_stats, sep="\t")
    logger.info(
        "Chromosome statistics loaded: rows=%d samples=%d",
        len(df),
        int(df["sample"].nunique()) if "sample" in df.columns else 0,
    )

    annotated_bin_df: pd.DataFrame | None = None
    output_artifacts: list[str] = []

    if args.min_binq is not None:
        if not args.bin_stats or not args.ppd_bin_quality:
            raise ValueError(
                "--min-binq requires both --bin-stats and --ppd-bin-quality."
            )
        bin_stats_df = pd.read_csv(args.bin_stats, sep="\t")
        bin_quality_df = pd.read_csv(args.ppd_bin_quality, sep="\t")
        requested_binq_field = args.binq_field
        if args.use_callq20:
            if requested_binq_field in {"BINQ15", "CALLQ15"}:
                raise ValueError(
                    "--use-callq20 cannot be combined with --binq-field set to a 15-threshold metric."
                )
            requested_binq_field = "CALLQ20"
        resolved_binq_field = resolve_binq_field(
            bin_quality_df,
            requested_binq_field,
        )
        annotated_bin_df = _annotate_binq_filter_for_bin_stats(
            bin_stats_df,
            bin_quality_df,
            min_binq=args.min_binq,
            binq_field=resolved_binq_field,
        )
        log_binq_filter_stats(annotated_bin_df)
        df = _apply_binq_filter_to_annotated_bins(
            df,
            annotated_bin_df,
            prob_threshold=args.prob_threshold,
        )
        filtered_path = os.path.join(args.output_dir, "chromosome_stats.filtered.tsv")
        df.to_csv(filtered_path, sep="\t", index=False)
        output_artifacts.append(filtered_path)
        logger.info("Wrote BINQ-filtered chromosome statistics")

    # ── optional exclusion ──────────────────────────────────────────────
    if args.exclusion_list:
        from gatk_sv_ploidy._util import load_exclusion_ids

        excl = set(load_exclusion_ids(args.exclusion_list))
        df = df[~df["sample"].isin(excl)]
        if annotated_bin_df is not None:
            annotated_bin_df = annotated_bin_df[~annotated_bin_df["sample"].isin(excl)].copy()
        logger.info("Applied exclusion list: excluded_samples=%d", len(excl))

    if annotated_bin_df is not None:
        ignored_path = os.path.join(args.output_dir, "ignored_bins.tsv.gz")
        export_ignored_bins(annotated_bin_df, ignored_path)
        output_artifacts.append(ignored_path)
        logger.info("Wrote ignored-bin table")

    # ── optional truth ──────────────────────────────────────────────────
    truth: Dict[str, str] = {}
    if args.truth_json:
        with open(args.truth_json) as fh:
            truth = json.load(fh)
        logger.info("Loaded truth labels: n_labels=%d", len(truth))

    # ── classify ────────────────────────────────────────────────────────
    pred_df = assign_sex_and_aneuploidy_types(
        df,
        truth,
    )
    save_sex_assignments(pred_df, args.output_dir)
    sex_path = os.path.join(args.output_dir, "sex_assignments.txt.gz")
    output_artifacts.append(sex_path)
    logger.info(
        "Assigned baseline-aware labels: samples=%d predicted_types=%d",
        len(pred_df),
        int(pred_df["predicted_aneuploidy_type"].nunique()),
    )

    # ── save predictions ────────────────────────────────────────────────
    pred_path = os.path.join(args.output_dir, "aneuploidy_type_predictions.tsv")
    pred_df.to_csv(pred_path, sep="\t", index=False)
    output_artifacts.append(pred_path)

    # ── export aneuploid data ───────────────────────────────────────────
    aneu_path = os.path.join(args.output_dir, "aneuploid_samples.tsv")
    export_aneuploid_data(df, aneu_path)
    output_artifacts.append(aneu_path)
    logger.info("Wrote call tables")
    log_output_artifacts(logger, output_artifacts)


def main() -> None:
    """Entry point for ``gatk-sv-ploidy call``."""
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    with tool_logging_context(
        tool_name="call",
        output_dir=args.output_dir,
        args=args,
    ) as logger:
        _run_call(args, logger)


if __name__ == "__main__":
    main()
