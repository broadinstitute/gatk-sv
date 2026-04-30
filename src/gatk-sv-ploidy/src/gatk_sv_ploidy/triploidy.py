"""Triploidy subcommand — classify autosomal baseline CN from pooled AF data."""

from __future__ import annotations

import argparse
import logging
import os
from typing import Optional, Sequence

import numpy as np
import pandas as pd
from scipy.special import expit, logsumexp

from gatk_sv_ploidy._util import get_sample_columns, save_and_close_plot
from gatk_sv_ploidy.models import _marginalized_af_log_lik_numpy

logger = logging.getLogger(__name__)

_TRIPLOID_BASELINE_CN = 3
_DIPLOID_BASELINE_CN = 2
_TRIPLOID_MINOR_AF = 1.0 / 3.0
_DIPLOID_HET_AF = 0.5
_DEFAULT_DIAGNOSTIC_AF_WINDOW = 0.08
_DEFAULT_DIAGNOSTIC_SAMPLE_LIMIT = 12
_DEFAULT_DIAGNOSTIC_DIPLOID_SAMPLE_LIMIT = 12
_DEFAULT_DIAGNOSTIC_MAX_SITES_PER_SAMPLE = 5000
_DIAGNOSTIC_RANDOM_SEED = 17
_DEFAULT_TRIPLOIDY_PRIOR = 1e-4
_DEFAULT_AF_CONCENTRATION_PRIOR_LOG_SD = 1.25
_PRIVACY_SAFE_QUANTILES = (0.01, 0.05, 0.25, 0.50, 0.75, 0.95, 0.99)


def _parse_af_concentration_grid(value: Optional[str]) -> Optional[np.ndarray]:
    if value is None or value.strip() == "":
        return None
    try:
        grid = np.array(
            [float(part.strip()) for part in value.split(",") if part.strip()],
            dtype=np.float64,
        )
    except ValueError as exc:
        raise ValueError(
            "--af-concentration-grid must be a comma-separated list of positive numbers."
        ) from exc
    if grid.size == 0 or np.any(~np.isfinite(grid)) or np.any(grid <= 0.0):
        raise ValueError(
            "--af-concentration-grid must contain at least one positive finite value."
        )
    return np.unique(np.sort(grid))


def _build_af_concentration_grid(
    af_concentration: float,
    af_concentration_grid: Optional[Sequence[float]],
) -> np.ndarray:
    if not np.isfinite(af_concentration) or af_concentration <= 0.0:
        raise ValueError("af_concentration must be a positive finite value.")
    if af_concentration_grid is not None:
        grid = np.asarray(af_concentration_grid, dtype=np.float64)
        if grid.size == 0 or np.any(~np.isfinite(grid)) or np.any(grid <= 0.0):
            raise ValueError(
                "af_concentration_grid must contain positive finite values."
            )
        return np.unique(np.sort(grid))

    log_grid = np.linspace(
        np.log(af_concentration) - 3.0,
        np.log(af_concentration) + 2.0,
        11,
    )
    grid = np.exp(log_grid)
    grid = np.concatenate([grid, np.array([af_concentration], dtype=np.float64)])
    return np.unique(np.sort(grid))


def _log_af_concentration_prior(
    concentration_grid: np.ndarray,
    af_concentration: float,
    af_concentration_prior_log_sd: float,
) -> np.ndarray:
    if (not np.isfinite(af_concentration_prior_log_sd) or
            af_concentration_prior_log_sd <= 0.0):
        raise ValueError("af_concentration_prior_log_sd must be positive.")
    log_grid = np.log(concentration_grid)
    log_median = np.log(af_concentration)
    log_weights = -0.5 * (
        (log_grid - log_median) / af_concentration_prior_log_sd
    ) ** 2
    return log_weights - logsumexp(log_weights)


def _summed_af_log_lik_by_concentration(
    *,
    site_alt: np.ndarray,
    site_total: np.ndarray,
    site_pop_af: np.ndarray,
    site_mask: np.ndarray,
    cn_state: int,
    n_states: int,
    concentration_grid: np.ndarray,
) -> np.ndarray:
    per_concentration: list[np.ndarray] = []
    for concentration in concentration_grid:
        per_bin_ll = _marginalized_af_log_lik_numpy(
            site_alt,
            site_total,
            site_pop_af,
            site_mask,
            cn_state=cn_state,
            n_states=n_states,
            concentration=float(concentration),
        )
        per_concentration.append(per_bin_ll.sum(axis=0))
    return np.vstack(per_concentration)


def _concentration_posterior_summaries(
    log_lik_by_concentration: np.ndarray,
    log_prior: np.ndarray,
    log_marginal_lik: np.ndarray,
    concentration_grid: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    log_posterior = (
        log_lik_by_concentration +
        log_prior[:, np.newaxis] -
        log_marginal_lik[np.newaxis, :]
    )
    posterior = np.exp(log_posterior)
    map_concentration = concentration_grid[np.argmax(log_posterior, axis=0)]
    mean_concentration = np.sum(
        posterior * concentration_grid[:, np.newaxis],
        axis=0,
    )
    return map_concentration, mean_concentration


def _finite_values(values: pd.Series | np.ndarray) -> np.ndarray:
    arr = np.asarray(values, dtype=np.float64).reshape(-1)
    return arr[np.isfinite(arr)]


def _format_quantile_summary(
    name: str,
    values: pd.Series | np.ndarray,
) -> str:
    finite_values = _finite_values(values)
    if finite_values.size == 0:
        return f"{name}: n={np.asarray(values).size}, finite=0"
    quantiles = np.quantile(finite_values, _PRIVACY_SAFE_QUANTILES)
    quantile_text = ", ".join(
        f"q{int(q * 100):02d}={value:.6g}"
        for q, value in zip(_PRIVACY_SAFE_QUANTILES, quantiles)
    )
    return (
        f"{name}: n={np.asarray(values).size}, finite={finite_values.size}, "
        f"mean={np.mean(finite_values):.6g}, "
        f"sd={np.std(finite_values, ddof=0):.6g}, {quantile_text}"
    )


def _format_category_counts(name: str, values: pd.Series) -> str:
    total = int(values.shape[0])
    if total == 0:
        return f"{name}: total=0"
    counts = values.value_counts(dropna=False).sort_index()
    parts = [f"total={total}"]
    for value, count in counts.items():
        parts.append(f"{value}={int(count)} ({100.0 * count / total:.1f}%)")
    return f"{name}: " + ", ".join(parts)


def _safe_correlation(
    x_values: pd.Series | np.ndarray,
    y_values: pd.Series | np.ndarray,
) -> tuple[float, int]:
    x = np.asarray(x_values, dtype=np.float64).reshape(-1)
    y = np.asarray(y_values, dtype=np.float64).reshape(-1)
    keep = np.isfinite(x) & np.isfinite(y)
    if int(keep.sum()) < 3:
        return float("nan"), int(keep.sum())
    x = x[keep]
    y = y[keep]
    if np.std(x, ddof=0) <= 0.0 or np.std(y, ddof=0) <= 0.0:
        return float("nan"), int(keep.sum())
    return float(np.corrcoef(x, y)[0, 1]), int(keep.sum())


def _log_privacy_safe_quantiles_by_baseline(
    results_df: pd.DataFrame,
    column: str,
) -> None:
    if column not in results_df.columns:
        return
    for baseline_cn in sorted(results_df["autosomal_baseline_cn"].dropna().unique()):
        subset = results_df[results_df["autosomal_baseline_cn"] == baseline_cn]
        logger.info(
            "PRIVACY_SAFE_TRIPLOIDY result_quantiles_by_baseline_cn "
            "column=%s baseline_cn=%s %s",
            column,
            baseline_cn,
            _format_quantile_summary(column, subset[column]),
        )


def _log_privacy_safe_threshold_counts(
    results_df: pd.DataFrame,
    pvalue_threshold: float,
    effect_size_threshold: float,
) -> None:
    posterior = results_df["triploidy_posterior_probability"].to_numpy(
        dtype=np.float64,
    )
    log_bf = results_df["triploidy_log_bayes_factor"].to_numpy(dtype=np.float64)
    effect_size = results_df["triploidy_effect_size_per_site"].to_numpy(
        dtype=np.float64,
    )
    posterior_error = results_df[
        "triploidy_posterior_error_probability"
    ].to_numpy(dtype=np.float64)
    finite_posterior = np.isfinite(posterior)
    thresholds = (0.5, 0.9, 0.99, 0.999, 0.9999)
    posterior_parts = [
        f"posterior_ge_{threshold:g}={int(np.sum(posterior >= threshold))}"
        for threshold in thresholds
    ]
    logger.info(
        "PRIVACY_SAFE_TRIPLOIDY posterior_threshold_counts total=%d "
        "finite=%d %s posterior_error_le_call_threshold=%d "
        "log_bayes_factor_positive=%d effect_size_ge_threshold=%d",
        int(results_df.shape[0]),
        int(np.sum(finite_posterior)),
        ", ".join(posterior_parts),
        int(np.sum(posterior_error <= pvalue_threshold)),
        int(np.sum(log_bf > 0.0)),
        int(np.sum(effect_size >= effect_size_threshold)),
    )


def _log_privacy_safe_decision_margin(
    results_df: pd.DataFrame,
    pvalue_threshold: float,
) -> None:
    prior_values = results_df["triploidy_prior"].dropna().unique()
    if len(prior_values) != 1:
        logger.info(
            "PRIVACY_SAFE_TRIPLOIDY decision_margin_summary "
            "skipped=non_unique_triploidy_prior unique_prior_count=%d",
            len(prior_values),
        )
        return
    triploidy_prior = float(prior_values[0])
    posterior_threshold = 1.0 - pvalue_threshold
    if not (0.0 < posterior_threshold < 1.0):
        return
    log_prior_odds = np.log(triploidy_prior) - np.log1p(-triploidy_prior)
    log_posterior_threshold_odds = (
        np.log(posterior_threshold) - np.log1p(-posterior_threshold)
    )
    required_log_bf = log_posterior_threshold_odds - log_prior_odds
    log_bf = results_df["triploidy_log_bayes_factor"].to_numpy(dtype=np.float64)
    margin = log_bf - required_log_bf
    called = results_df["autosomal_baseline_cn"] == _TRIPLOID_BASELINE_CN
    logger.info(
        "PRIVACY_SAFE_TRIPLOIDY decision_margin_summary "
        "triploidy_prior=%.6g posterior_call_threshold=%.6g "
        "required_log_bayes_factor=%.6g called_cn3=%d "
        "called_cn3_within_2_log_bf=%d called_cn3_within_5_log_bf=%d "
        "called_cn3_exceeds_by_20_log_bf=%d %s",
        triploidy_prior,
        posterior_threshold,
        required_log_bf,
        int(called.sum()),
        int(np.sum(called & (margin >= 0.0) & (margin <= 2.0))),
        int(np.sum(called & (margin >= 0.0) & (margin <= 5.0))),
        int(np.sum(called & (margin >= 20.0))),
        _format_quantile_summary("log_bayes_factor_margin", margin),
    )


def _log_privacy_safe_concentration_boundaries(
    results_df: pd.DataFrame,
) -> None:
    if "af_concentration_grid" not in results_df.columns or results_df.empty:
        return
    grid_values = results_df["af_concentration_grid"].dropna().unique()
    if len(grid_values) != 1:
        logger.info(
            "PRIVACY_SAFE_TRIPLOIDY concentration_grid_summary "
            "skipped=non_unique_grid unique_grid_count=%d",
            len(grid_values),
        )
        return
    concentration_grid = _parse_af_concentration_grid(str(grid_values[0]))
    if concentration_grid is None or concentration_grid.size == 0:
        return
    grid_min = float(concentration_grid[0])
    grid_max = float(concentration_grid[-1])
    called = results_df["autosomal_baseline_cn"] == _TRIPLOID_BASELINE_CN
    for label, mask in (("all", np.ones(len(results_df), dtype=bool)),
                        ("called_cn3", called.to_numpy(dtype=bool))):
        subset = results_df.loc[mask]
        if subset.empty:
            continue
        diploid_map = subset["diploid_af_concentration_map"].to_numpy(
            dtype=np.float64,
        )
        triploid_map = subset["triploid_af_concentration_map"].to_numpy(
            dtype=np.float64,
        )
        logger.info(
            "PRIVACY_SAFE_TRIPLOIDY concentration_boundary_counts group=%s "
            "n=%d grid_size=%d grid_min=%.6g grid_max=%.6g "
            "diploid_map_at_min=%d diploid_map_at_max=%d "
            "triploid_map_at_min=%d triploid_map_at_max=%d",
            label,
            int(subset.shape[0]),
            int(concentration_grid.size),
            grid_min,
            grid_max,
            int(np.sum(np.isclose(diploid_map, grid_min))),
            int(np.sum(np.isclose(diploid_map, grid_max))),
            int(np.sum(np.isclose(triploid_map, grid_min))),
            int(np.sum(np.isclose(triploid_map, grid_max))),
        )
    ratio = np.divide(
        results_df["triploid_af_concentration_mean"].to_numpy(dtype=np.float64),
        results_df["diploid_af_concentration_mean"].to_numpy(dtype=np.float64),
        out=np.full(results_df.shape[0], np.nan, dtype=np.float64),
        where=results_df["diploid_af_concentration_mean"].to_numpy(
            dtype=np.float64,
        ) > 0.0,
    )
    logger.info(
        "PRIVACY_SAFE_TRIPLOIDY concentration_mean_ratio_summary %s",
        _format_quantile_summary("triploid_over_diploid_concentration_mean", ratio),
    )


def _log_privacy_safe_af_peak_metrics(
    metrics_df: pd.DataFrame,
) -> None:
    metric_columns = [
        "fraction_near_diploid_half",
        "fraction_near_triploid_thirds",
        "fraction_intermediate_third_to_half",
        "contamination_mixture_score",
        "raw_af_triploid_minus_diploid_fraction",
        "folded_af_median",
        "folded_af_p10",
        "folded_af_p90",
        "raw_af_sd",
        "median_site_depth",
        "diagnostic_informative_sites",
    ]
    for baseline_cn in sorted(metrics_df["autosomal_baseline_cn"].dropna().unique()):
        subset = metrics_df[metrics_df["autosomal_baseline_cn"] == baseline_cn]
        logger.info(
            "PRIVACY_SAFE_TRIPLOIDY af_peak_counts baseline_cn=%s n=%d "
            "triploid_peak_fraction_gt_diploid_peak_fraction=%d "
            "triploid_peak_fraction_le_diploid_peak_fraction=%d "
            "mixture_score_ge_0.25=%d mixture_score_ge_0.50=%d",
            baseline_cn,
            int(subset.shape[0]),
            int(np.sum(
                subset["fraction_near_triploid_thirds"] >
                subset["fraction_near_diploid_half"]
            )),
            int(np.sum(
                subset["fraction_near_triploid_thirds"] <=
                subset["fraction_near_diploid_half"]
            )),
            int(np.sum(subset["contamination_mixture_score"] >= 0.25)),
            int(np.sum(subset["contamination_mixture_score"] >= 0.50)),
        )
        for column in metric_columns:
            logger.info(
                "PRIVACY_SAFE_TRIPLOIDY af_peak_quantiles_by_baseline_cn "
                "column=%s baseline_cn=%s %s",
                column,
                baseline_cn,
                _format_quantile_summary(column, subset[column]),
            )

    for column in [
        "fraction_near_diploid_half",
        "fraction_near_triploid_thirds",
        "contamination_mixture_score",
        "raw_af_triploid_minus_diploid_fraction",
        "folded_af_median",
        "raw_af_sd",
        "median_site_depth",
    ]:
        corr, n_used = _safe_correlation(
            metrics_df[column],
            metrics_df["triploidy_log_bayes_factor"],
        )
        logger.info(
            "PRIVACY_SAFE_TRIPLOIDY af_peak_correlation_with_log_bayes_factor "
            "metric=%s n=%d corr=%.6g",
            column,
            n_used,
            corr,
        )


def log_privacy_safe_triploidy_diagnostics(
    results_df: pd.DataFrame,
    *,
    pvalue_threshold: float,
    effect_size_threshold: float,
    metrics_df: Optional[pd.DataFrame] = None,
) -> None:
    """Log aggregate-only triploidy diagnostics safe for protected cohorts."""
    logger.info(
        "PRIVACY_SAFE_TRIPLOIDY privacy_contract aggregate_only=true "
        "forbidden=filenames,paths,sample_ids,raw_rows,sample_level_values"
    )
    logger.info(
        "PRIVACY_SAFE_TRIPLOIDY hypothesis "
        "cn3_overcalls_from_posterior_bayes_factor_or_dispersion_boundary"
    )
    logger.info(
        "PRIVACY_SAFE_TRIPLOIDY %s",
        _format_category_counts(
            "autosomal_baseline_cn_counts",
            results_df["autosomal_baseline_cn"],
        ),
    )
    logger.info(
        "PRIVACY_SAFE_TRIPLOIDY %s",
        _format_category_counts("triploidy_call_counts", results_df["triploidy_call"]),
    )
    logger.info(
        "PRIVACY_SAFE_TRIPLOIDY %s",
        _format_category_counts(
            "triploidy_reason_counts",
            results_df["triploidy_reason"],
        ),
    )
    _log_privacy_safe_threshold_counts(
        results_df,
        pvalue_threshold,
        effect_size_threshold,
    )
    _log_privacy_safe_decision_margin(results_df, pvalue_threshold)
    for column in [
        "triploidy_log_bayes_factor",
        "triploidy_log_posterior_odds",
        "triploidy_posterior_probability",
        "triploidy_posterior_error_probability",
        "triploidy_effect_size_per_site",
        "n_informative_bins",
        "n_informative_sites",
        "diploid_af_concentration_map",
        "triploid_af_concentration_map",
        "diploid_af_concentration_mean",
        "triploid_af_concentration_mean",
    ]:
        logger.info(
            "PRIVACY_SAFE_TRIPLOIDY result_quantiles_all %s",
            _format_quantile_summary(column, results_df[column]),
        )
        _log_privacy_safe_quantiles_by_baseline(results_df, column)
    _log_privacy_safe_concentration_boundaries(results_df)

    if metrics_df is None:
        logger.info(
            "PRIVACY_SAFE_TRIPLOIDY af_peak_metrics not_computed "
            "rerun_with_privacy_safe_diagnostics=true"
        )
    else:
        _log_privacy_safe_af_peak_metrics(metrics_df)


def _prepare_autosomal_informative_sites(
    *,
    n_samples: int,
    bin_chr: np.ndarray,
    site_alt: np.ndarray,
    site_total: np.ndarray,
    site_pop_af: np.ndarray,
    site_mask: np.ndarray,
    min_diploid_het_prior: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    if site_alt.shape != site_total.shape or site_alt.shape != site_mask.shape:
        raise ValueError("site_alt, site_total, and site_mask must share the same shape.")
    if site_alt.shape[0] != len(bin_chr):
        raise ValueError("bin_chr length must match the number of bins in site arrays.")
    if site_alt.shape[2] != n_samples:
        raise ValueError("sample_ids length must match the sample axis of site arrays.")

    bin_chr_arr = np.asarray(bin_chr, dtype=object)
    autosome_mask = ~np.isin(bin_chr_arr, ["chrX", "chrY"])
    if not np.any(autosome_mask):
        raise ValueError("Triploidy classification requires at least one autosomal bin.")

    site_pop_af_arr = np.asarray(site_pop_af, dtype=np.float64)
    if site_pop_af_arr.ndim == 2:
        diploid_het_prior = 2.0 * site_pop_af_arr * (1.0 - site_pop_af_arr)
        informative_site_mask = (
            site_mask &
            (site_total > 0) &
            (diploid_het_prior[:, :, np.newaxis] >= min_diploid_het_prior)
        )
    elif site_pop_af_arr.ndim == 3:
        diploid_het_prior = 2.0 * site_pop_af_arr * (1.0 - site_pop_af_arr)
        informative_site_mask = (
            site_mask &
            (site_total > 0) &
            (diploid_het_prior >= min_diploid_het_prior)
        )
    else:
        raise ValueError(
            "site_pop_af must have shape (n_bins, max_sites) or "
            "(n_bins, max_sites, n_samples)."
        )

    site_pop_af_auto = site_pop_af_arr[autosome_mask]
    return (
        autosome_mask,
        site_alt[autosome_mask],
        site_total[autosome_mask],
        site_pop_af_auto,
        informative_site_mask[autosome_mask],
        bin_chr_arr[autosome_mask],
    )


def classify_triploidy_from_site_data(
    *,
    sample_ids: list[str],
    bin_chr: np.ndarray,
    site_alt: np.ndarray,
    site_total: np.ndarray,
    site_pop_af: np.ndarray,
    site_mask: np.ndarray,
    af_concentration: float = 50.0,
    af_concentration_grid: Optional[Sequence[float]] = None,
    af_concentration_prior_log_sd: float = _DEFAULT_AF_CONCENTRATION_PRIOR_LOG_SD,
    triploidy_prior: float = _DEFAULT_TRIPLOIDY_PRIOR,
    min_diploid_het_prior: float = 0.1,
    min_informative_bins: int = 10,
    min_informative_sites: int = 100,
    pvalue_threshold: float = 1e-3,
    effect_size_threshold: float = 0.01,
    n_states: int = 6,
) -> pd.DataFrame:
    """Classify per-sample autosomal baseline CN as diploid or triploid.

    The classifier compares two generative models for pooled autosomal AF data:
    diploid baseline CN=2 and triploid baseline CN=3.  Each model marginalizes
    over latent genotypes and over an unknown beta-binomial concentration
    parameter so that overdispersed AF data can be explained by measurement
    scatter instead of forcing a triploid call.
    """
    _, site_alt_auto, site_total_auto, site_pop_af_auto, site_mask_auto, _ = (
        _prepare_autosomal_informative_sites(
            n_samples=len(sample_ids),
            bin_chr=bin_chr,
            site_alt=site_alt,
            site_total=site_total,
            site_pop_af=site_pop_af,
            site_mask=site_mask,
            min_diploid_het_prior=min_diploid_het_prior,
        )
    )

    if not np.isfinite(triploidy_prior) or not (0.0 < triploidy_prior < 1.0):
        raise ValueError("triploidy_prior must be between 0 and 1.")
    if not np.isfinite(pvalue_threshold) or not (0.0 < pvalue_threshold < 1.0):
        raise ValueError("pvalue_threshold must be between 0 and 1.")

    concentration_grid = _build_af_concentration_grid(
        af_concentration,
        af_concentration_grid,
    )
    log_concentration_prior = _log_af_concentration_prior(
        concentration_grid,
        af_concentration,
        af_concentration_prior_log_sd,
    )
    diploid_ll_by_concentration = _summed_af_log_lik_by_concentration(
        site_alt=site_alt_auto,
        site_total=site_total_auto,
        site_pop_af=site_pop_af_auto,
        site_mask=site_mask_auto,
        cn_state=_DIPLOID_BASELINE_CN,
        n_states=n_states,
        concentration_grid=concentration_grid,
    )
    triploid_ll_by_concentration = _summed_af_log_lik_by_concentration(
        site_alt=site_alt_auto,
        site_total=site_total_auto,
        site_pop_af=site_pop_af_auto,
        site_mask=site_mask_auto,
        cn_state=_TRIPLOID_BASELINE_CN,
        n_states=n_states,
        concentration_grid=concentration_grid,
    )
    diploid_log_marginal = logsumexp(
        diploid_ll_by_concentration + log_concentration_prior[:, np.newaxis],
        axis=0,
    )
    triploid_log_marginal = logsumexp(
        triploid_ll_by_concentration + log_concentration_prior[:, np.newaxis],
        axis=0,
    )
    diploid_concentration_map, diploid_concentration_mean = (
        _concentration_posterior_summaries(
            diploid_ll_by_concentration,
            log_concentration_prior,
            diploid_log_marginal,
            concentration_grid,
        )
    )
    triploid_concentration_map, triploid_concentration_mean = (
        _concentration_posterior_summaries(
            triploid_ll_by_concentration,
            log_concentration_prior,
            triploid_log_marginal,
            concentration_grid,
        )
    )

    log_bayes_factor = triploid_log_marginal - diploid_log_marginal
    log_prior_odds = np.log(triploidy_prior) - np.log1p(-triploidy_prior)
    log_posterior_odds = log_prior_odds + log_bayes_factor
    posterior_triploidy = expit(log_posterior_odds)
    posterior_error = 1.0 - posterior_triploidy
    informative_site_counts = site_mask_auto.sum(axis=1).astype(np.int64, copy=False)

    rows: list[dict[str, object]] = []
    for sample_idx, sample_id in enumerate(sample_ids):
        block_sites = informative_site_counts[:, sample_idx]
        informative_blocks = block_sites > 0
        n_blocks = int(informative_blocks.sum())
        n_sites = int(block_sites[informative_blocks].sum())

        total_diploid_ll = float(diploid_log_marginal[sample_idx])
        total_triploid_ll = float(triploid_log_marginal[sample_idx])
        total_llr = float(log_bayes_factor[sample_idx])
        effect_size = total_llr / n_sites if n_sites > 0 else float("nan")
        sample_posterior = float(posterior_triploidy[sample_idx])
        sample_posterior_error = float(posterior_error[sample_idx])
        sample_log_posterior_odds = float(log_posterior_odds[sample_idx])

        if n_blocks < min_informative_bins:
            call = "INSUFFICIENT_DATA"
            reason = "insufficient_informative_bins"
            baseline_cn = _DIPLOID_BASELINE_CN
        elif n_sites < min_informative_sites:
            call = "INSUFFICIENT_DATA"
            reason = "insufficient_informative_sites"
            baseline_cn = _DIPLOID_BASELINE_CN
        elif total_llr <= 0.0:
            call = "DIPLOID"
            reason = "triploid_bayes_factor_not_positive"
            baseline_cn = _DIPLOID_BASELINE_CN
        elif not np.isfinite(effect_size) or effect_size < effect_size_threshold:
            call = "DIPLOID"
            reason = "effect_size_below_threshold"
            baseline_cn = _DIPLOID_BASELINE_CN
        elif (
            not np.isfinite(sample_posterior_error) or
            sample_posterior_error > pvalue_threshold
        ):
            call = "DIPLOID"
            reason = "posterior_error_above_threshold"
            baseline_cn = _DIPLOID_BASELINE_CN
        else:
            call = "TRIPLOID"
            reason = "triploid_posterior_supported"
            baseline_cn = _TRIPLOID_BASELINE_CN

        rows.append(
            {
                "sample": str(sample_id),
                "autosomal_baseline_cn": int(baseline_cn),
                "triploidy_call": call,
                "triploidy_reason": reason,
                "triploidy_t_stat": sample_log_posterior_odds,
                "triploidy_pvalue": sample_posterior_error,
                "triploidy_posterior_probability": sample_posterior,
                "triploidy_posterior_error_probability": sample_posterior_error,
                "triploidy_log_bayes_factor": total_llr,
                "triploidy_log_prior_odds": float(log_prior_odds),
                "triploidy_log_posterior_odds": sample_log_posterior_odds,
                "triploidy_log_lik_ratio": total_llr,
                "triploidy_effect_size_per_site": float(effect_size),
                "diploid_log_lik": total_diploid_ll,
                "triploid_log_lik": total_triploid_ll,
                "diploid_af_concentration_map": float(
                    diploid_concentration_map[sample_idx]
                ),
                "diploid_af_concentration_mean": float(
                    diploid_concentration_mean[sample_idx]
                ),
                "triploid_af_concentration_map": float(
                    triploid_concentration_map[sample_idx]
                ),
                "triploid_af_concentration_mean": float(
                    triploid_concentration_mean[sample_idx]
                ),
                "af_concentration_prior_median": float(af_concentration),
                "af_concentration_prior_log_sd": float(
                    af_concentration_prior_log_sd
                ),
                "af_concentration_grid": ",".join(
                    f"{value:.6g}" for value in concentration_grid
                ),
                "triploidy_prior": float(triploidy_prior),
                "n_informative_bins": n_blocks,
                "n_informative_sites": n_sites,
            }
        )

    return pd.DataFrame(rows)


def build_triploidy_diagnostic_metrics(
    *,
    sample_ids: list[str],
    bin_chr: np.ndarray,
    site_alt: np.ndarray,
    site_total: np.ndarray,
    site_pop_af: np.ndarray,
    site_mask: np.ndarray,
    results_df: pd.DataFrame,
    min_diploid_het_prior: float,
    af_window: float = _DEFAULT_DIAGNOSTIC_AF_WINDOW,
) -> pd.DataFrame:
    """Summarize raw autosomal AF patterns used by the triploidy test."""
    if af_window <= 0.0 or af_window >= 0.25:
        raise ValueError("af_window must be greater than 0 and less than 0.25.")

    _, site_alt_auto, site_total_auto, _, site_mask_auto, _ = (
        _prepare_autosomal_informative_sites(
            n_samples=len(sample_ids),
            bin_chr=bin_chr,
            site_alt=site_alt,
            site_total=site_total,
            site_pop_af=site_pop_af,
            site_mask=site_mask,
            min_diploid_het_prior=min_diploid_het_prior,
        )
    )
    observed_af = np.divide(
        site_alt_auto,
        site_total_auto,
        out=np.full(site_alt_auto.shape, np.nan, dtype=np.float64),
        where=site_total_auto > 0,
    )

    rows: list[dict[str, object]] = []
    for sample_idx, sample_id in enumerate(sample_ids):
        sample_mask = site_mask_auto[:, :, sample_idx]
        af_values = observed_af[:, :, sample_idx][sample_mask]
        depth_values = site_total_auto[:, :, sample_idx][sample_mask]
        finite_mask = np.isfinite(af_values)
        af_values = af_values[finite_mask]
        depth_values = depth_values[finite_mask]

        if af_values.size == 0:
            rows.append(
                {
                    "sample": str(sample_id),
                    "diagnostic_informative_sites": 0,
                    "median_site_depth": np.nan,
                    "raw_af_mean": np.nan,
                    "raw_af_sd": np.nan,
                    "folded_af_median": np.nan,
                    "folded_af_p10": np.nan,
                    "folded_af_p90": np.nan,
                    "fraction_near_diploid_half": np.nan,
                    "fraction_near_triploid_thirds": np.nan,
                    "fraction_intermediate_third_to_half": np.nan,
                    "contamination_mixture_score": np.nan,
                    "raw_af_triploid_minus_diploid_fraction": np.nan,
                }
            )
            continue

        folded_af = np.minimum(af_values, 1.0 - af_values)
        diploid_distance = np.abs(af_values - _DIPLOID_HET_AF)
        triploid_distance = np.minimum(
            np.abs(af_values - _TRIPLOID_MINOR_AF),
            np.abs(af_values - (1.0 - _TRIPLOID_MINOR_AF)),
        )
        near_diploid = diploid_distance <= af_window
        near_triploid = triploid_distance <= af_window
        intermediate = (
            (folded_af >= _TRIPLOID_MINOR_AF + (af_window / 2.0)) &
            (folded_af <= _DIPLOID_HET_AF - (af_window / 2.0))
        )
        fraction_near_diploid = float(np.mean(near_diploid))
        fraction_near_triploid = float(np.mean(near_triploid))

        rows.append(
            {
                "sample": str(sample_id),
                "diagnostic_informative_sites": int(af_values.size),
                "median_site_depth": float(np.median(depth_values)),
                "raw_af_mean": float(np.mean(af_values)),
                "raw_af_sd": float(np.std(af_values, ddof=0)),
                "folded_af_median": float(np.median(folded_af)),
                "folded_af_p10": float(np.quantile(folded_af, 0.10)),
                "folded_af_p90": float(np.quantile(folded_af, 0.90)),
                "fraction_near_diploid_half": fraction_near_diploid,
                "fraction_near_triploid_thirds": fraction_near_triploid,
                "fraction_intermediate_third_to_half": float(np.mean(intermediate)),
                "contamination_mixture_score": float(
                    np.sqrt(fraction_near_diploid * fraction_near_triploid)
                ),
                "raw_af_triploid_minus_diploid_fraction": float(
                    fraction_near_triploid - fraction_near_diploid
                ),
            }
        )

    metrics_df = pd.DataFrame(rows)
    result_columns = [
        "sample",
        "autosomal_baseline_cn",
        "triploidy_call",
        "triploidy_reason",
        "triploidy_pvalue",
        "triploidy_posterior_probability",
        "triploidy_log_bayes_factor",
        "triploidy_log_lik_ratio",
        "triploidy_effect_size_per_site",
        "diploid_af_concentration_map",
        "triploid_af_concentration_map",
        "n_informative_bins",
        "n_informative_sites",
    ]
    return metrics_df.merge(
        results_df[result_columns],
        on="sample",
        how="left",
        validate="one_to_one",
    )


def _select_diagnostic_samples(
    metrics_df: pd.DataFrame,
    sample_limit: int,
) -> list[str]:
    if sample_limit <= 0 or metrics_df.empty:
        return []
    ranked_df = metrics_df.copy()
    ranked_df["_triploid_sort"] = (
        ranked_df["autosomal_baseline_cn"] == _TRIPLOID_BASELINE_CN
    ).astype(int)
    ranked_df["_mix_sort"] = ranked_df["contamination_mixture_score"].fillna(-1.0)
    ranked_df["_effect_sort"] = ranked_df[
        "triploidy_effect_size_per_site"
    ].fillna(-np.inf)
    ranked_df["_site_sort"] = ranked_df["diagnostic_informative_sites"].fillna(0)
    ranked_df = ranked_df.sort_values(
        ["_triploid_sort", "_mix_sort", "_effect_sort", "_site_sort"],
        ascending=[False, False, False, False],
    )
    return [str(sample_id) for sample_id in ranked_df["sample"].head(sample_limit)]


def _select_diploid_diagnostic_samples(
    metrics_df: pd.DataFrame,
    sample_limit: int = _DEFAULT_DIAGNOSTIC_DIPLOID_SAMPLE_LIMIT,
) -> list[str]:
    if sample_limit <= 0 or metrics_df.empty:
        return []
    diploid_df = metrics_df[
        (metrics_df["autosomal_baseline_cn"] == _DIPLOID_BASELINE_CN) &
        (metrics_df["triploidy_call"] == "DIPLOID")
    ].copy()
    if diploid_df.empty:
        return []
    diploid_df["_mix_sort"] = diploid_df[
        "contamination_mixture_score"
    ].fillna(-1.0)
    diploid_df["_site_sort"] = diploid_df[
        "diagnostic_informative_sites"
    ].fillna(0)
    diploid_df = diploid_df.sort_values(
        ["_mix_sort", "_site_sort", "sample"],
        ascending=[False, False, True],
    )
    return [str(sample_id) for sample_id in diploid_df["sample"].head(sample_limit)]


def _collect_diagnostic_site_rows(
    *,
    sample_ids: list[str],
    selected_samples: list[str],
    bin_chr: np.ndarray,
    site_alt: np.ndarray,
    site_total: np.ndarray,
    site_pop_af: np.ndarray,
    site_mask: np.ndarray,
    min_diploid_het_prior: float,
    max_sites_per_sample: int,
) -> pd.DataFrame:
    if max_sites_per_sample <= 0 or not selected_samples:
        return pd.DataFrame()
    _, site_alt_auto, site_total_auto, _, site_mask_auto, bin_chr_auto = (
        _prepare_autosomal_informative_sites(
            n_samples=len(sample_ids),
            bin_chr=bin_chr,
            site_alt=site_alt,
            site_total=site_total,
            site_pop_af=site_pop_af,
            site_mask=site_mask,
            min_diploid_het_prior=min_diploid_het_prior,
        )
    )

    sample_index = {sample_id: idx for idx, sample_id in enumerate(sample_ids)}
    rng = np.random.default_rng(_DIAGNOSTIC_RANDOM_SEED)
    rows: list[dict[str, object]] = []
    for sample_id in selected_samples:
        sample_idx = sample_index.get(sample_id)
        if sample_idx is None:
            continue
        sample_mask = site_mask_auto[:, :, sample_idx]
        coordinates = np.argwhere(sample_mask)
        if coordinates.shape[0] > max_sites_per_sample:
            selected_idx = rng.choice(
                coordinates.shape[0],
                size=max_sites_per_sample,
                replace=False,
            )
            coordinates = coordinates[np.sort(selected_idx)]

        for bin_idx, site_idx in coordinates:
            total = int(site_total_auto[bin_idx, site_idx, sample_idx])
            alt = int(site_alt_auto[bin_idx, site_idx, sample_idx])
            if total <= 0:
                continue
            observed_af = float(alt / total)
            rows.append(
                {
                    "sample": str(sample_id),
                    "chr": str(bin_chr_auto[bin_idx]),
                    "autosomal_bin_index": int(bin_idx),
                    "site_index": int(site_idx),
                    "site_alt": alt,
                    "site_total": total,
                    "observed_af": observed_af,
                    "folded_observed_af": float(
                        min(observed_af, 1.0 - observed_af)
                    ),
                }
            )
    return pd.DataFrame(rows)


def _plot_triploidy_diagnostic_metrics(
    metrics_df: pd.DataFrame,
    output_dir: str,
) -> None:
    if metrics_df.empty:
        return
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    color_by_call = {
        "DIPLOID": "#1f77b4",
        "TRIPLOID": "#d62728",
        "INSUFFICIENT_DATA": "#7f7f7f",
    }

    ax = axes[0]
    for call, sub_df in metrics_df.groupby("triploidy_call", dropna=False):
        call_label = str(call)
        ax.scatter(
            sub_df["fraction_near_diploid_half"],
            sub_df["fraction_near_triploid_thirds"],
            s=np.clip(np.log10(sub_df["diagnostic_informative_sites"] + 1.0), 1.0, 5.0) * 20.0,
            alpha=0.75,
            color=color_by_call.get(call_label, "#9467bd"),
            edgecolor="white",
            linewidth=0.4,
            label=call_label,
        )
    ax.set_xlabel("Fraction of informative sites near AF=0.5")
    ax.set_ylabel("Fraction near AF=1/3 or 2/3")
    ax.set_title("Raw AF peak balance")
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best", fontsize=8)

    ax = axes[1]
    for call, sub_df in metrics_df.groupby("triploidy_call", dropna=False):
        call_label = str(call)
        ax.scatter(
            sub_df["triploidy_effect_size_per_site"],
            sub_df["contamination_mixture_score"],
            s=np.clip(np.log10(sub_df["diagnostic_informative_sites"] + 1.0), 1.0, 5.0) * 20.0,
            alpha=0.75,
            color=color_by_call.get(call_label, "#9467bd"),
            edgecolor="white",
            linewidth=0.4,
            label=call_label,
        )
    ax.axvline(0.0, color="black", linestyle="--", linewidth=1.0, alpha=0.45)
    ax.set_xlabel("Triploid - diploid log-likelihood per site")
    ax.set_ylabel("Mixed 0.5 and 1/3/2/3 peak score")
    ax.set_title("Triploidy evidence vs mixture-like AF pattern")
    ax.grid(True, alpha=0.25)

    label_df = metrics_df.sort_values(
        ["contamination_mixture_score", "triploidy_effect_size_per_site"],
        ascending=[False, False],
    ).head(8)
    for _, row in label_df.iterrows():
        if not np.isfinite(row["contamination_mixture_score"]):
            continue
        ax.annotate(
            str(row["sample"]),
            (
                row["triploidy_effect_size_per_site"],
                row["contamination_mixture_score"],
            ),
            fontsize=6,
            xytext=(3, 3),
            textcoords="offset points",
        )

    save_and_close_plot(
        output_dir,
        "triploidy_diagnostic_metrics.png",
        subdir="diagnostics",
    )


def _plot_triploidy_raw_af_profiles(
    raw_site_df: pd.DataFrame,
    metrics_df: pd.DataFrame,
    output_dir: str,
    filename: str = "triploidy_raw_af_profiles.png",
) -> None:
    if raw_site_df.empty:
        return
    import matplotlib.pyplot as plt

    samples = list(raw_site_df["sample"].drop_duplicates())
    n_cols = min(3, len(samples))
    n_rows = int(np.ceil(len(samples) / n_cols))
    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=(4.2 * n_cols, 2.8 * n_rows),
        squeeze=False,
    )
    metrics_by_sample = metrics_df.set_index("sample")
    bins = np.linspace(0.0, 1.0, 51)
    for ax, sample_id in zip(axes.ravel(), samples):
        sub_df = raw_site_df[raw_site_df["sample"] == sample_id]
        ax.hist(
            sub_df["observed_af"],
            bins=bins,
            color="#4c78a8",
            alpha=0.78,
            edgecolor="none",
        )
        ax.axvline(_TRIPLOID_MINOR_AF, color="#d62728", linestyle="--", linewidth=1.0)
        ax.axvline(_DIPLOID_HET_AF, color="black", linestyle=":", linewidth=1.1)
        ax.axvline(1.0 - _TRIPLOID_MINOR_AF, color="#d62728", linestyle="--", linewidth=1.0)
        ax.set_xlim(0.0, 1.0)
        ax.set_xlabel("Raw observed alt fraction")
        ax.set_ylabel("Sites")
        title = str(sample_id)
        if sample_id in metrics_by_sample.index:
            row = metrics_by_sample.loc[sample_id]
            title = (
                f"{sample_id}\n{row['triploidy_call']}, "
                f"mix={row['contamination_mixture_score']:.2f}"
            )
        ax.set_title(title, fontsize=9)
        ax.grid(True, axis="y", alpha=0.2)

    for ax in axes.ravel()[len(samples):]:
        ax.axis("off")
    save_and_close_plot(
        output_dir,
        filename,
        subdir="diagnostics",
    )


def write_triploidy_diagnostics(
    *,
    output_dir: str,
    sample_ids: list[str],
    bin_chr: np.ndarray,
    site_alt: np.ndarray,
    site_total: np.ndarray,
    site_pop_af: np.ndarray,
    site_mask: np.ndarray,
    results_df: pd.DataFrame,
    min_diploid_het_prior: float,
    af_window: float,
    sample_limit: int,
    max_sites_per_sample: int,
) -> None:
    """Write triploidy diagnostic metrics and plots."""
    metrics_df = build_triploidy_diagnostic_metrics(
        sample_ids=sample_ids,
        bin_chr=bin_chr,
        site_alt=site_alt,
        site_total=site_total,
        site_pop_af=site_pop_af,
        site_mask=site_mask,
        results_df=results_df,
        min_diploid_het_prior=min_diploid_het_prior,
        af_window=af_window,
    )
    diagnostics_dir = os.path.join(output_dir, "diagnostics")
    os.makedirs(diagnostics_dir, exist_ok=True)
    metrics_path = os.path.join(diagnostics_dir, "triploidy_diagnostic_metrics.tsv")
    metrics_df.to_csv(metrics_path, sep="\t", index=False)
    logger.info("Triploidy diagnostic metrics saved.")

    selected_samples = _select_diagnostic_samples(metrics_df, sample_limit)
    raw_site_df = _collect_diagnostic_site_rows(
        sample_ids=sample_ids,
        selected_samples=selected_samples,
        bin_chr=bin_chr,
        site_alt=site_alt,
        site_total=site_total,
        site_pop_af=site_pop_af,
        site_mask=site_mask,
        min_diploid_het_prior=min_diploid_het_prior,
        max_sites_per_sample=max_sites_per_sample,
    )
    raw_site_path = os.path.join(diagnostics_dir, "triploidy_raw_af_sites.tsv.gz")
    raw_site_df.to_csv(raw_site_path, sep="\t", index=False)
    logger.info("Triploidy sampled raw AF sites saved.")

    diploid_samples = _select_diploid_diagnostic_samples(metrics_df)
    diploid_site_df = _collect_diagnostic_site_rows(
        sample_ids=sample_ids,
        selected_samples=diploid_samples,
        bin_chr=bin_chr,
        site_alt=site_alt,
        site_total=site_total,
        site_pop_af=site_pop_af,
        site_mask=site_mask,
        min_diploid_het_prior=min_diploid_het_prior,
        max_sites_per_sample=max_sites_per_sample,
    )

    _plot_triploidy_diagnostic_metrics(metrics_df, output_dir)
    _plot_triploidy_raw_af_profiles(raw_site_df, metrics_df, output_dir)
    _plot_triploidy_raw_af_profiles(
        diploid_site_df,
        metrics_df,
        output_dir,
        filename="triploidy_raw_af_diploid_profiles.png",
    )
    logger.info("Triploidy diagnostics saved.")


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the triploidy subcommand."""
    parser = argparse.ArgumentParser(
        description="Classify triploidy from pooled autosomal allele counts",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Preprocessed depth TSV (output of 'preprocess')",
    )
    parser.add_argument(
        "--site-data",
        required=True,
        help="Per-site allele data .npz (output of 'preprocess')",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        required=True,
        help="Output directory",
    )
    parser.add_argument(
        "--af-concentration",
        type=float,
        default=50.0,
        help="Prior median for beta-binomial AF concentration / overdispersion",
    )
    parser.add_argument(
        "--af-concentration-grid",
        required=False,
        help=(
            "Optional comma-separated positive concentration grid used to "
            "marginalize over AF overdispersion. Defaults to a log-spaced "
            "grid around --af-concentration."
        ),
    )
    parser.add_argument(
        "--af-concentration-prior-log-sd",
        type=float,
        default=_DEFAULT_AF_CONCENTRATION_PRIOR_LOG_SD,
        help=(
            "Log-scale standard deviation for the discrete prior over "
            "beta-binomial AF concentration. Larger values allow more "
            "overdispersed AF data."
        ),
    )
    parser.add_argument(
        "--triploidy-prior",
        type=float,
        default=_DEFAULT_TRIPLOIDY_PRIOR,
        help="Prior probability that any one sample is triploid",
    )
    parser.add_argument(
        "--min-diploid-het-prior",
        type=float,
        default=0.1,
        help="Minimum diploid heterozygous prior mass 2p(1-p) to retain a site",
    )
    parser.add_argument(
        "--min-informative-bins",
        type=int,
        default=10,
        help="Minimum autosomal bins with informative AF sites required to attempt a call",
    )
    parser.add_argument(
        "--min-informative-sites",
        type=int,
        default=100,
        help="Minimum informative autosomal site count required to attempt a call",
    )
    parser.add_argument(
        "--pvalue-threshold",
        type=float,
        default=1e-3,
        help=(
            "Maximum posterior error probability for calling triploidy. "
            "The flag name is retained for compatibility with older runs."
        ),
    )
    parser.add_argument(
        "--effect-size-threshold",
        type=float,
        default=0.01,
        help="Minimum mean log-likelihood ratio per informative site for calling triploidy",
    )
    parser.add_argument(
        "--diagnostics",
        action="store_true",
        help="Write triploidy diagnostic metrics, plots, and sampled raw AF data",
    )
    parser.add_argument(
        "--privacy-safe-diagnostics",
        action="store_true",
        help=(
            "Log aggregate-only triploidy diagnostics without filenames, paths, "
            "sample identifiers, raw rows, or sample-level values."
        ),
    )
    parser.add_argument(
        "--diagnostic-sample-limit",
        type=int,
        default=_DEFAULT_DIAGNOSTIC_SAMPLE_LIMIT,
        help="Maximum number of samples to include in raw AF diagnostic plots",
    )
    parser.add_argument(
        "--diagnostic-max-sites-per-sample",
        type=int,
        default=_DEFAULT_DIAGNOSTIC_MAX_SITES_PER_SAMPLE,
        help="Maximum sampled informative sites per sample in raw AF diagnostics",
    )
    parser.add_argument(
        "--diagnostic-af-window",
        type=float,
        default=_DEFAULT_DIAGNOSTIC_AF_WINDOW,
        help="AF window around 0.5 and 1/3 or 2/3 for diagnostic summaries",
    )
    return parser.parse_args()


def main() -> None:
    """Entry point for ``gatk-sv-ploidy triploidy``."""
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    depth_df = pd.read_csv(args.input, sep="\t", index_col=0)
    sample_ids = [str(sample_id) for sample_id in get_sample_columns(depth_df)]

    site_npz = np.load(args.site_data, allow_pickle=True)
    site_alt = np.asarray(site_npz["site_alt"], dtype=np.int32)
    site_total = np.asarray(site_npz["site_total"], dtype=np.int32)
    site_pop_af = np.asarray(site_npz["site_pop_af"], dtype=np.float64)
    site_mask = np.asarray(site_npz["site_mask"], dtype=bool)

    if "sample_ids" in site_npz:
        npz_sample_ids = [str(sample_id) for sample_id in np.asarray(site_npz["sample_ids"]).tolist()]
        if npz_sample_ids != sample_ids:
            raise ValueError(
                "site_data sample_ids do not match the preprocessed depth sample order."
            )
    if site_alt.shape[2] != len(sample_ids):
        raise ValueError(
            "site_data sample axis does not match the preprocessed depth sample count."
        )

    if "bin_chr" in site_npz:
        bin_chr = np.asarray(site_npz["bin_chr"], dtype=object)
    else:
        bin_chr = depth_df["Chr"].astype(object).to_numpy()
    if site_alt.shape[0] != len(bin_chr):
        raise ValueError("site_data bin axis does not match available chromosome metadata.")

    af_concentration_grid = _parse_af_concentration_grid(
        args.af_concentration_grid,
    )
    results_df = classify_triploidy_from_site_data(
        sample_ids=sample_ids,
        bin_chr=bin_chr,
        site_alt=site_alt,
        site_total=site_total,
        site_pop_af=site_pop_af,
        site_mask=site_mask,
        af_concentration=args.af_concentration,
        af_concentration_grid=af_concentration_grid,
        af_concentration_prior_log_sd=args.af_concentration_prior_log_sd,
        triploidy_prior=args.triploidy_prior,
        min_diploid_het_prior=args.min_diploid_het_prior,
        min_informative_bins=args.min_informative_bins,
        min_informative_sites=args.min_informative_sites,
        pvalue_threshold=args.pvalue_threshold,
        effect_size_threshold=args.effect_size_threshold,
    )

    privacy_safe_metrics_df = None
    if args.privacy_safe_diagnostics:
        privacy_safe_metrics_df = build_triploidy_diagnostic_metrics(
            sample_ids=sample_ids,
            bin_chr=bin_chr,
            site_alt=site_alt,
            site_total=site_total,
            site_pop_af=site_pop_af,
            site_mask=site_mask,
            results_df=results_df,
            min_diploid_het_prior=args.min_diploid_het_prior,
            af_window=args.diagnostic_af_window,
        )
    log_privacy_safe_triploidy_diagnostics(
        results_df,
        pvalue_threshold=args.pvalue_threshold,
        effect_size_threshold=args.effect_size_threshold,
        metrics_df=privacy_safe_metrics_df,
    )

    results_path = os.path.join(args.output_dir, "triploidy_test_results.tsv")
    results_df.to_csv(results_path, sep="\t", index=False)
    logger.info("Triploidy test results saved.")

    manifest_df = results_df[["sample", "autosomal_baseline_cn"]].copy()
    manifest_path = os.path.join(args.output_dir, "sample_autosomal_baseline_cn.tsv")
    manifest_df.to_csv(manifest_path, sep="\t", index=False)
    logger.info("Sample autosomal baseline CN manifest saved.")

    if args.diagnostics:
        write_triploidy_diagnostics(
            output_dir=args.output_dir,
            sample_ids=sample_ids,
            bin_chr=bin_chr,
            site_alt=site_alt,
            site_total=site_total,
            site_pop_af=site_pop_af,
            site_mask=site_mask,
            results_df=results_df,
            min_diploid_het_prior=args.min_diploid_het_prior,
            af_window=args.diagnostic_af_window,
            sample_limit=args.diagnostic_sample_limit,
            max_sites_per_sample=args.diagnostic_max_sites_per_sample,
        )


if __name__ == "__main__":
    main()
