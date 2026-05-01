"""Autosomal baseline CN classifier from pooled allele-fraction data."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import logging
import os
import time
from typing import Optional, Sequence

import numpy as np
import pandas as pd
from scipy.special import gammaln, logsumexp
from tqdm import tqdm

from gatk_sv_ploidy._util import (
    baseline_ploidy_label,
    get_sample_columns,
    save_and_close_plot,
)

logger = logging.getLogger(__name__)

_BASELINE_CN_STATES = (1, 2, 3, 4)
_NON_DIPLOID_CN_STATES = (1, 3, 4)
_DIPLOID_BASELINE_CN = 2
_TRIPLOID_BASELINE_CN = 3
_TETRAPLOID_BASELINE_CN = 4
_PEAKS_BY_CN = {
    1: (0.0, 1.0),
    2: (0.0, 0.5, 1.0),
    3: (0.0, 1.0 / 3.0, 2.0 / 3.0, 1.0),
    4: (0.0, 0.25, 0.5, 0.75, 1.0),
}
_HAPLOID_ENDPOINT_PEAKS = (0.0, 1.0)
_DIPLOID_COMPATIBLE_PEAKS = (0.0, 0.5, 1.0)
_TRIPLOID_UNIQUE_PEAKS = (1.0 / 3.0, 2.0 / 3.0)
_TETRAPLOID_UNIQUE_PEAKS = (0.25, 0.75)
_DEFAULT_HAPLOIDY_PRIOR = 1e-6
_DEFAULT_TRIPLOIDY_PRIOR = 1e-4
_DEFAULT_TETRAPLOIDY_PRIOR = 1e-6
_DEFAULT_AF_CONCENTRATION_PRIOR_LOG_SD = 1.25
_DEFAULT_PLOIDY_PEAK_AF_WINDOW = 0.04
_DEFAULT_DIAGNOSTIC_AF_WINDOW = 0.08
_DEFAULT_PEAK_EVIDENCE_WEIGHT = 0.25
_DEFAULT_PEAK_KERNEL_SD = 0.035
_DEFAULT_PEAK_OUTLIER_PROBABILITY = 0.02
_DEFAULT_MIN_HAPLOID_ENDPOINT_FRACTION = 0.85
_DEFAULT_MAX_HAPLOID_OTHER_PEAK_FRACTION = 0.05
_DEFAULT_MIN_TRIPLOID_PEAK_FRACTION_ADVANTAGE = 0.0
_DEFAULT_MIN_TETRAPLOID_QUARTER_PEAK_FRACTION = 0.10
_DEFAULT_MIN_TETRAPLOID_QUARTER_PEAK_FRACTION_ADVANTAGE = 0.05
_DEFAULT_DIAGNOSTIC_SAMPLE_LIMIT = 12
_DEFAULT_DIAGNOSTIC_DIPLOID_SAMPLE_LIMIT = 12
_DEFAULT_DIAGNOSTIC_MAX_SITES_PER_SAMPLE = 5000
_DIAGNOSTIC_RANDOM_SEED = 17
_PRIVACY_SAFE_QUANTILES = (0.01, 0.05, 0.25, 0.50, 0.75, 0.95, 0.99)
_AF_LIKELIHOOD_CONCENTRATION_BATCH_SIZE = 4


@dataclass(frozen=True)
class _FlattenedAfObservations:
    site_alt: np.ndarray
    site_total: np.ndarray
    site_pop_af: np.ndarray
    sample_index: np.ndarray
    log_choose_total_alt: np.ndarray
    n_samples: int

    @property
    def n_observations(self) -> int:
        return int(self.site_alt.size)


def _progress_iter(
    iterable,
    *,
    desc: str,
    unit: str,
    show_progress: bool,
    total: Optional[int] = None,
):
    return tqdm(
        iterable,
        desc=desc,
        unit=unit,
        total=total,
        dynamic_ncols=True,
        disable=not show_progress,
    )


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
            "--af-concentration-grid must be a comma-separated list of "
            "positive numbers."
        ) from exc
    if grid.size == 0 or np.any(~np.isfinite(grid)) or np.any(grid <= 0.0):
        raise ValueError(
            "--af-concentration-grid must contain at least one positive "
            "finite value."
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


def _build_baseline_cn_priors(
    *,
    haploidy_prior: float,
    triploidy_prior: float,
    tetraploidy_prior: float,
) -> np.ndarray:
    priors = {
        1: float(haploidy_prior),
        3: float(triploidy_prior),
        4: float(tetraploidy_prior),
    }
    for cn_state, prior in priors.items():
        if not np.isfinite(prior) or prior <= 0.0:
            raise ValueError(f"cn{cn_state} prior must be positive and finite.")
    non_diploid_total = sum(priors.values())
    if non_diploid_total >= 1.0:
        raise ValueError("Non-diploid baseline CN priors must sum to less than 1.")
    priors[2] = 1.0 - non_diploid_total
    return np.asarray([priors[cn] for cn in _BASELINE_CN_STATES], dtype=np.float64)


def _flatten_valid_af_observations(
    *,
    site_alt: np.ndarray,
    site_total: np.ndarray,
    site_pop_af: np.ndarray,
    site_mask: np.ndarray,
    n_samples: int,
) -> _FlattenedAfObservations:
    site_alt_arr = np.asarray(site_alt, dtype=np.float64)
    site_total_safe = np.maximum(np.asarray(site_total, dtype=np.float64), 1.0)
    site_mask_arr = np.asarray(site_mask, dtype=bool)
    if site_alt_arr.shape != site_total_safe.shape:
        raise ValueError("site_alt and site_total must have matching shapes.")
    if site_alt_arr.shape != site_mask_arr.shape:
        raise ValueError("site_mask must match site_alt shape.")
    if site_alt_arr.ndim != 3:
        raise ValueError("site_alt must have shape (n_bins, max_sites, n_samples).")
    if site_alt_arr.shape[2] != n_samples:
        raise ValueError("site_alt sample dimension must match n_samples.")

    bin_index, site_index, sample_index = np.nonzero(site_mask_arr)
    alt_flat = site_alt_arr[bin_index, site_index, sample_index]
    total_flat = site_total_safe[bin_index, site_index, sample_index]

    site_pop_af_arr = np.asarray(site_pop_af, dtype=np.float64)
    if site_pop_af_arr.ndim == 2:
        pop_af_flat = site_pop_af_arr[bin_index, site_index]
    elif site_pop_af_arr.ndim == 3:
        pop_af_flat = site_pop_af_arr[bin_index, site_index, sample_index]
    else:
        raise ValueError(
            "site_pop_af must have shape (n_bins, max_sites) or "
            "(n_bins, max_sites, n_samples)."
        )

    eps = 1e-6
    pop_af_flat = np.clip(pop_af_flat, eps, 1.0 - eps)
    log_choose = (
        gammaln(total_flat + 1.0) -
        gammaln(alt_flat + 1.0) -
        gammaln(total_flat - alt_flat + 1.0)
    )
    log_choose = np.nan_to_num(log_choose, nan=0.0, posinf=0.0, neginf=-100.0)

    return _FlattenedAfObservations(
        site_alt=alt_flat,
        site_total=total_flat,
        site_pop_af=pop_af_flat,
        sample_index=sample_index.astype(np.int64, copy=False),
        log_choose_total_alt=log_choose,
        n_samples=n_samples,
    )


def _beta_binomial_logpmf_from_flat_observations(
    observations: _FlattenedAfObservations,
    *,
    alpha: np.ndarray,
    beta: np.ndarray,
) -> np.ndarray:
    alt = observations.site_alt[np.newaxis, :]
    total = observations.site_total[np.newaxis, :]
    alpha_arr = np.asarray(alpha, dtype=np.float64)[:, np.newaxis]
    beta_arr = np.asarray(beta, dtype=np.float64)[:, np.newaxis]
    log_prob = observations.log_choose_total_alt[np.newaxis, :]
    log_prob = log_prob + gammaln(alt + alpha_arr)
    log_prob = log_prob + gammaln(total - alt + beta_arr)
    log_prob = log_prob - gammaln(total + alpha_arr + beta_arr)
    log_prob = log_prob - gammaln(alpha_arr)
    log_prob = log_prob - gammaln(beta_arr)
    log_prob = log_prob + gammaln(alpha_arr + beta_arr)

    invalid_count = (
        (observations.site_alt < 0.0) |
        (observations.site_alt > observations.site_total)
    )
    if np.any(invalid_count):
        log_prob[:, invalid_count] = -np.inf
    return np.nan_to_num(log_prob, nan=0.0, posinf=0.0, neginf=-100.0)


def _summed_af_log_lik_by_concentration_fast(
    *,
    observations: _FlattenedAfObservations,
    cn_state: int,
    concentration_grid: np.ndarray,
    show_progress: bool = False,
    progress_desc: str = "AF likelihood grid",
) -> np.ndarray:
    result = np.zeros(
        (int(concentration_grid.size), observations.n_samples),
        dtype=np.float64,
    )
    if observations.n_observations == 0:
        return result

    eps = 1e-6
    pop_af = observations.site_pop_af
    log_pop_af = np.log(pop_af)
    log_ref_af = np.log1p(-pop_af)
    batch_size = _AF_LIKELIHOOD_CONCENTRATION_BATCH_SIZE
    batch_starts = range(0, int(concentration_grid.size), batch_size)
    batch_iter = _progress_iter(
        batch_starts,
        desc=progress_desc,
        unit="batch",
        total=int(np.ceil(concentration_grid.size / batch_size)),
        show_progress=show_progress,
    )

    for batch_start in batch_iter:
        batch_end = min(batch_start + batch_size, int(concentration_grid.size))
        concentration_batch = concentration_grid[batch_start:batch_end]
        log_marginal = None
        for genotype_dosage in range(cn_state + 1):
            log_genotype_prior = (
                gammaln(cn_state + 1.0) -
                gammaln(genotype_dosage + 1.0) -
                gammaln(cn_state - genotype_dosage + 1.0) +
                genotype_dosage * log_pop_af +
                (cn_state - genotype_dosage) * log_ref_af
            )
            expected_af = genotype_dosage / max(cn_state, 1)
            alpha = concentration_batch * expected_af + eps
            beta = concentration_batch * (1.0 - expected_af) + eps
            log_bb = _beta_binomial_logpmf_from_flat_observations(
                observations,
                alpha=alpha,
                beta=beta,
            )
            log_term = log_bb + log_genotype_prior[np.newaxis, :]
            if log_marginal is None:
                log_marginal = log_term
            else:
                log_marginal = np.logaddexp(log_marginal, log_term)

        if log_marginal is None:
            continue
        for batch_offset in range(batch_end - batch_start):
            result[batch_start + batch_offset] = np.bincount(
                observations.sample_index,
                weights=log_marginal[batch_offset],
                minlength=observations.n_samples,
            )

    return result


def _summed_af_log_lik_by_concentration(
    *,
    site_alt: np.ndarray,
    site_total: np.ndarray,
    site_pop_af: np.ndarray,
    site_mask: np.ndarray,
    cn_state: int,
    n_states: int,
    concentration_grid: np.ndarray,
    show_progress: bool = False,
    progress_desc: str = "AF likelihood grid",
    observations: Optional[_FlattenedAfObservations] = None,
) -> np.ndarray:
    if observations is None:
        observations = _flatten_valid_af_observations(
            site_alt=site_alt,
            site_total=site_total,
            site_pop_af=site_pop_af,
            site_mask=site_mask,
            n_samples=int(site_alt.shape[2]),
        )
    return _summed_af_log_lik_by_concentration_fast(
        observations=observations,
        cn_state=cn_state,
        concentration_grid=concentration_grid,
        show_progress=show_progress,
        progress_desc=progress_desc,
    )


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
    arr = np.asarray(values)
    finite_values = _finite_values(values)
    if finite_values.size == 0:
        return f"{name}: n={arr.size}, finite=0"
    quantiles = np.quantile(finite_values, _PRIVACY_SAFE_QUANTILES)
    quantile_text = ", ".join(
        f"q{int(q * 100):02d}={value:.6g}"
        for q, value in zip(_PRIVACY_SAFE_QUANTILES, quantiles)
    )
    return (
        f"{name}: n={arr.size}, finite={finite_values.size}, "
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


def _distance_to_peaks(values: np.ndarray, peaks: Sequence[float]) -> np.ndarray:
    peak_arr = np.asarray(peaks, dtype=np.float64)
    return np.min(np.abs(values[:, np.newaxis] - peak_arr[np.newaxis, :]), axis=1)


def _fraction_near_peaks(
    values: np.ndarray,
    peaks: Sequence[float],
    af_window: float,
) -> float:
    if values.size == 0:
        return float("nan")
    return float(np.mean(_distance_to_peaks(values, peaks) <= af_window))


def _compute_observed_af_peak_metrics(
    observed_af: np.ndarray,
    sample_mask: np.ndarray,
    af_window: float,
) -> dict[str, float]:
    af_values = observed_af[sample_mask]
    af_values = af_values[np.isfinite(af_values)]
    if af_values.size == 0:
        return {
            "fraction_near_haploid_endpoints": float("nan"),
            "fraction_near_diploid_half": float("nan"),
            "fraction_near_triploid_thirds": float("nan"),
            "fraction_near_tetraploid_quarters": float("nan"),
            "fraction_near_diploid_compatible": float("nan"),
            "fraction_near_any_modeled_peak": float("nan"),
            "triploid_peak_fraction_advantage": float("nan"),
            "tetraploid_quarter_peak_fraction_advantage": float("nan"),
            "haploid_endpoint_peak_fraction_advantage": float("nan"),
        }

    fraction_near_haploid = _fraction_near_peaks(
        af_values,
        _HAPLOID_ENDPOINT_PEAKS,
        af_window,
    )
    fraction_near_diploid = _fraction_near_peaks(af_values, (0.5,), af_window)
    fraction_near_triploid = _fraction_near_peaks(
        af_values,
        _TRIPLOID_UNIQUE_PEAKS,
        af_window,
    )
    fraction_near_tetraploid = _fraction_near_peaks(
        af_values,
        _TETRAPLOID_UNIQUE_PEAKS,
        af_window,
    )
    fraction_near_diploid_compatible = _fraction_near_peaks(
        af_values,
        _DIPLOID_COMPATIBLE_PEAKS,
        af_window,
    )
    all_peaks = sorted(
        set(_PEAKS_BY_CN[1] + _PEAKS_BY_CN[2] + _PEAKS_BY_CN[3] + _PEAKS_BY_CN[4])
    )
    fraction_near_any = _fraction_near_peaks(af_values, all_peaks, af_window)
    other_non_haploid = max(
        fraction_near_diploid,
        fraction_near_triploid,
        fraction_near_tetraploid,
    )
    return {
        "fraction_near_haploid_endpoints": fraction_near_haploid,
        "fraction_near_diploid_half": fraction_near_diploid,
        "fraction_near_triploid_thirds": fraction_near_triploid,
        "fraction_near_tetraploid_quarters": fraction_near_tetraploid,
        "fraction_near_diploid_compatible": fraction_near_diploid_compatible,
        "fraction_near_any_modeled_peak": fraction_near_any,
        "triploid_peak_fraction_advantage": (
            fraction_near_triploid - fraction_near_diploid
        ),
        "tetraploid_quarter_peak_fraction_advantage": (
            fraction_near_tetraploid - fraction_near_diploid
        ),
        "haploid_endpoint_peak_fraction_advantage": (
            fraction_near_haploid - other_non_haploid
        ),
    }


def _log_normal_pdf(values: np.ndarray, means: np.ndarray, sd: float) -> np.ndarray:
    return -0.5 * ((values - means) / sd) ** 2 - np.log(sd * np.sqrt(2.0 * np.pi))


def _state_peak_log_lik(
    observed_af: np.ndarray,
    site_pop_af: np.ndarray,
    cn_state: int,
    peak_kernel_sd: float,
    peak_outlier_probability: float,
) -> float:
    if observed_af.size == 0:
        return 0.0
    p = np.clip(site_pop_af, 1e-6, 1.0 - 1e-6)
    k_values = np.arange(cn_state + 1, dtype=np.float64)
    log_comb = (
        gammaln(cn_state + 1.0) -
        gammaln(k_values + 1.0) -
        gammaln(cn_state - k_values + 1.0)
    )
    log_weights = (
        log_comb[np.newaxis, :] +
        k_values[np.newaxis, :] * np.log(p[:, np.newaxis]) +
        (cn_state - k_values)[np.newaxis, :] * np.log1p(-p[:, np.newaxis])
    )
    peak_means = k_values / float(cn_state)
    log_kernel = _log_normal_pdf(
        observed_af[:, np.newaxis],
        peak_means[np.newaxis, :],
        peak_kernel_sd,
    )
    log_signal = np.log1p(-peak_outlier_probability) + logsumexp(
        log_weights + log_kernel,
        axis=1,
    )
    log_outlier = np.log(peak_outlier_probability)
    return float(np.sum(np.logaddexp(log_signal, log_outlier)))


def _peak_log_likelihood_matrix(
    *,
    observed_af_auto: np.ndarray,
    site_pop_af_auto: np.ndarray,
    site_mask_auto: np.ndarray,
    peak_kernel_sd: float,
    peak_outlier_probability: float,
    show_progress: bool = False,
) -> np.ndarray:
    if not np.isfinite(peak_kernel_sd) or peak_kernel_sd <= 0.0:
        raise ValueError("peak_kernel_sd must be positive and finite.")
    if (not np.isfinite(peak_outlier_probability) or
            peak_outlier_probability <= 0.0 or
            peak_outlier_probability >= 1.0):
        raise ValueError("peak_outlier_probability must be between 0 and 1.")

    n_samples = observed_af_auto.shape[2]
    peak_ll = np.zeros((len(_BASELINE_CN_STATES), n_samples), dtype=np.float64)
    sample_iter = _progress_iter(
        range(n_samples),
        desc="Baseline CN peak evidence",
        unit="sample",
        total=n_samples,
        show_progress=show_progress,
    )
    for sample_idx in sample_iter:
        sample_af = observed_af_auto[:, :, sample_idx]
        sample_mask = site_mask_auto[:, :, sample_idx] & np.isfinite(sample_af)
        if not np.any(sample_mask):
            continue
        af_values = sample_af[sample_mask]
        if site_pop_af_auto.ndim == 2:
            pop_af_values = site_pop_af_auto[sample_mask]
        else:
            pop_af_values = site_pop_af_auto[:, :, sample_idx][sample_mask]
        for state_idx, cn_state in enumerate(_BASELINE_CN_STATES):
            peak_ll[state_idx, sample_idx] = _state_peak_log_lik(
                af_values,
                pop_af_values,
                cn_state,
                peak_kernel_sd,
                peak_outlier_probability,
            )
    return peak_ll


def _has_direct_state_peak_support(
    metrics: dict[str, float],
    cn_state: int,
    *,
    min_haploid_endpoint_fraction: float,
    max_haploid_other_peak_fraction: float,
    min_triploid_peak_fraction_advantage: float,
    min_tetraploid_quarter_peak_fraction: float,
    min_tetraploid_quarter_peak_fraction_advantage: float,
) -> bool:
    if cn_state == 1:
        other_peak_fraction = max(
            metrics["fraction_near_diploid_half"],
            metrics["fraction_near_triploid_thirds"],
            metrics["fraction_near_tetraploid_quarters"],
        )
        endpoint_ok = (
            np.isfinite(metrics["fraction_near_haploid_endpoints"]) and
            metrics["fraction_near_haploid_endpoints"] >=
            min_haploid_endpoint_fraction
        )
        other_peak_ok = (
            np.isfinite(other_peak_fraction) and
            other_peak_fraction <= max_haploid_other_peak_fraction
        )
        return bool(endpoint_ok and other_peak_ok)
    if cn_state == 3:
        triploid_peak_ok = (
            np.isfinite(metrics["triploid_peak_fraction_advantage"]) and
            metrics["triploid_peak_fraction_advantage"] >
            min_triploid_peak_fraction_advantage
        )
        return bool(triploid_peak_ok)
    if cn_state == 4:
        quarter_fraction_ok = (
            np.isfinite(metrics["fraction_near_tetraploid_quarters"]) and
            metrics["fraction_near_tetraploid_quarters"] >=
            min_tetraploid_quarter_peak_fraction
        )
        quarter_advantage_ok = (
            np.isfinite(metrics["tetraploid_quarter_peak_fraction_advantage"]) and
            metrics["tetraploid_quarter_peak_fraction_advantage"] >=
            min_tetraploid_quarter_peak_fraction_advantage
        )
        return bool(quarter_fraction_ok and quarter_advantage_ok)
    return cn_state == 2


def _direct_support_failure_reason(cn_state: int) -> str:
    if cn_state == 1:
        return "haploid_endpoint_support_below_threshold"
    if cn_state == 3:
        return "triploid_peak_support_below_threshold"
    if cn_state == 4:
        return "tetraploid_quarter_support_below_threshold"
    return "diploid_default_or_ambiguous"


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
        raise ValueError(
            "Autosomal baseline CN classification requires at least one "
            "autosomal bin."
        )

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


def classify_polyploidy_from_site_data(
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
    haploidy_prior: float = _DEFAULT_HAPLOIDY_PRIOR,
    triploidy_prior: float = _DEFAULT_TRIPLOIDY_PRIOR,
    tetraploidy_prior: float = _DEFAULT_TETRAPLOIDY_PRIOR,
    peak_evidence_weight: float = _DEFAULT_PEAK_EVIDENCE_WEIGHT,
    peak_kernel_sd: float = _DEFAULT_PEAK_KERNEL_SD,
    peak_outlier_probability: float = _DEFAULT_PEAK_OUTLIER_PROBABILITY,
    min_diploid_het_prior: float = 0.1,
    min_informative_bins: int = 10,
    min_informative_sites: int = 100,
    pvalue_threshold: float = 1e-3,
    effect_size_threshold: float = 0.01,
    peak_af_window: float = _DEFAULT_PLOIDY_PEAK_AF_WINDOW,
    min_haploid_endpoint_fraction: float = _DEFAULT_MIN_HAPLOID_ENDPOINT_FRACTION,
    max_haploid_other_peak_fraction: float = _DEFAULT_MAX_HAPLOID_OTHER_PEAK_FRACTION,
    min_triploid_peak_fraction_advantage: float = (
        _DEFAULT_MIN_TRIPLOID_PEAK_FRACTION_ADVANTAGE
    ),
    min_tetraploid_quarter_peak_fraction: float = (
        _DEFAULT_MIN_TETRAPLOID_QUARTER_PEAK_FRACTION
    ),
    min_tetraploid_quarter_peak_fraction_advantage: float = (
        _DEFAULT_MIN_TETRAPLOID_QUARTER_PEAK_FRACTION_ADVANTAGE
    ),
    n_states: int = 6,
    show_progress: bool = False,
    log_progress: bool = False,
) -> pd.DataFrame:
    """Classify per-sample autosomal baseline CN as CN1, CN2, CN3, or CN4."""
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

    if not np.isfinite(pvalue_threshold) or not (0.0 < pvalue_threshold < 1.0):
        raise ValueError("pvalue_threshold must be between 0 and 1.")
    if not np.isfinite(effect_size_threshold):
        raise ValueError("effect_size_threshold must be finite.")
    if not np.isfinite(peak_evidence_weight) or peak_evidence_weight < 0.0:
        raise ValueError("peak_evidence_weight must be non-negative and finite.")
    if not np.isfinite(peak_af_window) or not (0.0 < peak_af_window < 0.25):
        raise ValueError("peak_af_window must be greater than 0 and less than 0.25.")

    prior_probs = _build_baseline_cn_priors(
        haploidy_prior=haploidy_prior,
        triploidy_prior=triploidy_prior,
        tetraploidy_prior=tetraploidy_prior,
    )
    log_cn_prior = np.log(prior_probs)
    concentration_grid = _build_af_concentration_grid(
        af_concentration,
        af_concentration_grid,
    )
    log_concentration_prior = _log_af_concentration_prior(
        concentration_grid,
        af_concentration,
        af_concentration_prior_log_sd,
    )
    af_observations = _flatten_valid_af_observations(
        site_alt=site_alt_auto,
        site_total=site_total_auto,
        site_pop_af=site_pop_af_auto,
        site_mask=site_mask_auto,
        n_samples=len(sample_ids),
    )

    if log_progress:
        logger.info(
            "Polyploidy progress: prepared aggregate AF matrix "
            "samples=%d autosomal_bins=%d max_sites_per_bin=%d "
            "informative_site_observations=%d concentration_grid_size=%d "
            "baseline_cn_states=1,2,3,4",
            int(len(sample_ids)),
            int(site_alt_auto.shape[0]),
            int(site_alt_auto.shape[1]),
            int(af_observations.n_observations),
            int(concentration_grid.size),
        )

    genotype_log_marginal_by_cn: dict[int, np.ndarray] = {}
    concentration_map_by_cn: dict[int, np.ndarray] = {}
    concentration_mean_by_cn: dict[int, np.ndarray] = {}
    for cn_state in _BASELINE_CN_STATES:
        if log_progress:
            logger.info(
                "Polyploidy progress: evaluating CN%d AF likelihood grid.",
                cn_state,
            )
        cn_timer_start = time.perf_counter()
        ll_by_concentration = _summed_af_log_lik_by_concentration(
            site_alt=site_alt_auto,
            site_total=site_total_auto,
            site_pop_af=site_pop_af_auto,
            site_mask=site_mask_auto,
            cn_state=cn_state,
            n_states=n_states,
            concentration_grid=concentration_grid,
            show_progress=show_progress,
            progress_desc=f"CN{cn_state} AF likelihood",
            observations=af_observations,
        )
        if log_progress:
            logger.info(
                "Polyploidy progress: CN%d AF likelihood grid complete "
                "elapsed_sec=%.3f valid_observations=%d "
                "concentration_grid_size=%d",
                cn_state,
                float(time.perf_counter() - cn_timer_start),
                int(af_observations.n_observations),
                int(concentration_grid.size),
            )
        log_marginal = logsumexp(
            ll_by_concentration + log_concentration_prior[:, np.newaxis],
            axis=0,
        )
        genotype_log_marginal_by_cn[cn_state] = log_marginal
        concentration_map, concentration_mean = _concentration_posterior_summaries(
            ll_by_concentration,
            log_concentration_prior,
            log_marginal,
            concentration_grid,
        )
        concentration_map_by_cn[cn_state] = concentration_map
        concentration_mean_by_cn[cn_state] = concentration_mean

    observed_af_auto = np.divide(
        site_alt_auto,
        site_total_auto,
        out=np.full(site_alt_auto.shape, np.nan, dtype=np.float64),
        where=site_total_auto > 0,
    )
    if log_progress:
        logger.info("Polyploidy progress: evaluating direct AF peak mixtures.")
    peak_log_lik = _peak_log_likelihood_matrix(
        observed_af_auto=observed_af_auto,
        site_pop_af_auto=site_pop_af_auto,
        site_mask_auto=site_mask_auto,
        peak_kernel_sd=peak_kernel_sd,
        peak_outlier_probability=peak_outlier_probability,
        show_progress=show_progress,
    )

    genotype_log_marginal = np.vstack(
        [genotype_log_marginal_by_cn[cn] for cn in _BASELINE_CN_STATES]
    )
    log_evidence = (
        log_cn_prior[:, np.newaxis] +
        genotype_log_marginal +
        float(peak_evidence_weight) * peak_log_lik
    )
    log_posterior_norm = logsumexp(log_evidence, axis=0)
    posterior = np.exp(log_evidence - log_posterior_norm[np.newaxis, :])
    state_index = {cn: idx for idx, cn in enumerate(_BASELINE_CN_STATES)}
    diploid_idx = state_index[_DIPLOID_BASELINE_CN]

    informative_site_counts = site_mask_auto.sum(axis=1).astype(np.int64, copy=False)
    peak_metrics = [
        _compute_observed_af_peak_metrics(
            observed_af_auto[:, :, sample_idx],
            site_mask_auto[:, :, sample_idx],
            peak_af_window,
        )
        for sample_idx in range(len(sample_ids))
    ]

    if log_progress:
        logger.info("Polyploidy progress: classifying samples.")
    rows: list[dict[str, object]] = []
    sample_indices = _progress_iter(
        range(len(sample_ids)),
        desc="Classifying samples",
        unit="sample",
        total=len(sample_ids),
        show_progress=show_progress,
    )
    for sample_idx in sample_indices:
        sample_id = sample_ids[sample_idx]
        block_sites = informative_site_counts[:, sample_idx]
        informative_blocks = block_sites > 0
        n_blocks = int(informative_blocks.sum())
        n_sites = int(block_sites[informative_blocks].sum())
        sample_metrics = peak_metrics[sample_idx]

        sample_posterior = {
            cn: float(posterior[state_index[cn], sample_idx])
            for cn in _BASELINE_CN_STATES
        }
        sample_log_evidence = {
            cn: float(log_evidence[state_index[cn], sample_idx])
            for cn in _BASELINE_CN_STATES
        }
        sample_genotype_ll = {
            cn: float(genotype_log_marginal[state_index[cn], sample_idx])
            for cn in _BASELINE_CN_STATES
        }
        sample_peak_ll = {
            cn: float(peak_log_lik[state_index[cn], sample_idx])
            for cn in _BASELINE_CN_STATES
        }
        log_bf_vs_cn2 = {
            cn: float(
                sample_log_evidence[cn] -
                sample_log_evidence[_DIPLOID_BASELINE_CN]
            )
            for cn in _BASELINE_CN_STATES
        }
        effect_size_vs_cn2 = {
            cn: (
                log_bf_vs_cn2[cn] / n_sites
                if n_sites > 0 else float("nan")
            )
            for cn in _BASELINE_CN_STATES
        }
        direct_support = {
            cn: _has_direct_state_peak_support(
                sample_metrics,
                cn,
                min_haploid_endpoint_fraction=min_haploid_endpoint_fraction,
                max_haploid_other_peak_fraction=max_haploid_other_peak_fraction,
                min_triploid_peak_fraction_advantage=(
                    min_triploid_peak_fraction_advantage
                ),
                min_tetraploid_quarter_peak_fraction=(
                    min_tetraploid_quarter_peak_fraction
                ),
                min_tetraploid_quarter_peak_fraction_advantage=(
                    min_tetraploid_quarter_peak_fraction_advantage
                ),
            )
            for cn in _NON_DIPLOID_CN_STATES
        }
        posterior_supported = {
            cn: bool(
                np.isfinite(sample_posterior[cn]) and
                (1.0 - sample_posterior[cn]) <= pvalue_threshold and
                log_bf_vs_cn2[cn] > 0.0 and
                np.isfinite(effect_size_vs_cn2[cn]) and
                effect_size_vs_cn2[cn] >= effect_size_threshold
            )
            for cn in _NON_DIPLOID_CN_STATES
        }
        call_supported = {
            cn: bool(posterior_supported[cn] and direct_support[cn])
            for cn in _NON_DIPLOID_CN_STATES
        }

        map_cn = _BASELINE_CN_STATES[int(np.argmax(posterior[:, sample_idx]))]
        if n_blocks < min_informative_bins:
            baseline_cn = _DIPLOID_BASELINE_CN
            reason = "insufficient_informative_bins"
            call = "INSUFFICIENT_DATA"
        elif n_sites < min_informative_sites:
            baseline_cn = _DIPLOID_BASELINE_CN
            reason = "insufficient_informative_sites"
            call = "INSUFFICIENT_DATA"
        else:
            supported_states = [cn for cn in _NON_DIPLOID_CN_STATES if call_supported[cn]]
            if supported_states:
                baseline_cn = max(supported_states, key=lambda cn: sample_posterior[cn])
                call = baseline_ploidy_label(baseline_cn)
                reason = f"cn{baseline_cn}_posterior_and_peak_supported"
            else:
                baseline_cn = _DIPLOID_BASELINE_CN
                call = "DIPLOID"
                top_non_diploid = max(
                    _NON_DIPLOID_CN_STATES,
                    key=lambda cn: sample_posterior[cn],
                )
                if map_cn == _DIPLOID_BASELINE_CN:
                    reason = "diploid_posterior_supported"
                elif not posterior_supported[top_non_diploid]:
                    if log_bf_vs_cn2[top_non_diploid] <= 0.0:
                        reason = f"cn{top_non_diploid}_bayes_factor_not_positive"
                    elif (
                        not np.isfinite(effect_size_vs_cn2[top_non_diploid]) or
                        effect_size_vs_cn2[top_non_diploid] < effect_size_threshold
                    ):
                        reason = f"cn{top_non_diploid}_effect_size_below_threshold"
                    else:
                        reason = f"cn{top_non_diploid}_posterior_error_above_threshold"
                else:
                    reason = _direct_support_failure_reason(top_non_diploid)
                    if top_non_diploid == _TETRAPLOID_BASELINE_CN:
                        reason = "tetraploid_not_conclusive_diploid_default"

        final_posterior = sample_posterior[baseline_cn]
        final_posterior_error = 1.0 - final_posterior
        triploid_log_prior_odds = (
            np.log(prior_probs[state_index[_TRIPLOID_BASELINE_CN]]) -
            np.log(prior_probs[diploid_idx])
        )
        triploid_log_posterior_odds = (
            np.log(sample_posterior[_TRIPLOID_BASELINE_CN]) -
            np.log(sample_posterior[_DIPLOID_BASELINE_CN])
        )

        row: dict[str, object] = {
            "sample": str(sample_id),
            "autosomal_baseline_cn": int(baseline_cn),
            "baseline_cn_call": call,
            "baseline_cn_reason": reason,
            "baseline_cn_map": int(map_cn),
            "baseline_cn_map_call": baseline_ploidy_label(int(map_cn)),
            "baseline_cn_posterior_probability": float(final_posterior),
            "baseline_cn_posterior_error_probability": float(final_posterior_error),
            "baseline_cn_map_posterior_probability": float(sample_posterior[map_cn]),
            "baseline_cn_effect_size_per_site": float(
                effect_size_vs_cn2.get(baseline_cn, 0.0)
            ),
            "triploidy_call": call,
            "triploidy_reason": reason,
            "triploidy_t_stat": float(triploid_log_posterior_odds),
            "triploidy_pvalue": float(1.0 - sample_posterior[_TRIPLOID_BASELINE_CN]),
            "triploidy_posterior_probability": float(
                sample_posterior[_TRIPLOID_BASELINE_CN]
            ),
            "triploidy_posterior_error_probability": float(
                1.0 - sample_posterior[_TRIPLOID_BASELINE_CN]
            ),
            "triploidy_log_bayes_factor": float(
                log_bf_vs_cn2[_TRIPLOID_BASELINE_CN]
            ),
            "triploidy_log_prior_odds": float(triploid_log_prior_odds),
            "triploidy_log_posterior_odds": float(triploid_log_posterior_odds),
            "triploidy_log_lik_ratio": float(log_bf_vs_cn2[_TRIPLOID_BASELINE_CN]),
            "triploidy_effect_size_per_site": float(
                effect_size_vs_cn2[_TRIPLOID_BASELINE_CN]
            ),
            "ploidy_peak_af_window": float(peak_af_window),
            "triploidy_peak_af_window": float(peak_af_window),
            "min_haploid_endpoint_fraction": float(min_haploid_endpoint_fraction),
            "max_haploid_other_peak_fraction": float(max_haploid_other_peak_fraction),
            "min_triploid_peak_fraction_advantage": float(
                min_triploid_peak_fraction_advantage
            ),
            "min_tetraploid_quarter_peak_fraction": float(
                min_tetraploid_quarter_peak_fraction
            ),
            "min_tetraploid_quarter_peak_fraction_advantage": float(
                min_tetraploid_quarter_peak_fraction_advantage
            ),
            "peak_evidence_weight": float(peak_evidence_weight),
            "peak_kernel_sd": float(peak_kernel_sd),
            "peak_outlier_probability": float(peak_outlier_probability),
            "af_concentration_prior_median": float(af_concentration),
            "af_concentration_prior_log_sd": float(af_concentration_prior_log_sd),
            "af_concentration_grid": ",".join(
                f"{value:.6g}" for value in concentration_grid
            ),
            "haploidy_prior": float(prior_probs[state_index[1]]),
            "diploidy_prior": float(prior_probs[state_index[2]]),
            "triploidy_prior": float(prior_probs[state_index[3]]),
            "tetraploidy_prior": float(prior_probs[state_index[4]]),
            "n_informative_bins": n_blocks,
            "n_informative_sites": n_sites,
        }
        row.update(sample_metrics)
        for cn_state in _BASELINE_CN_STATES:
            row[f"posterior_cn_{cn_state}"] = sample_posterior[cn_state]
            row[f"log_evidence_cn_{cn_state}"] = sample_log_evidence[cn_state]
            row[f"log_bayes_factor_cn{cn_state}_vs_cn2"] = log_bf_vs_cn2[cn_state]
            row[f"effect_size_cn{cn_state}_vs_cn2"] = effect_size_vs_cn2[cn_state]
            row[f"genotype_log_lik_cn_{cn_state}"] = sample_genotype_ll[cn_state]
            row[f"peak_log_lik_cn_{cn_state}"] = sample_peak_ll[cn_state]
            row[f"cn{cn_state}_af_concentration_map"] = float(
                concentration_map_by_cn[cn_state][sample_idx]
            )
            row[f"cn{cn_state}_af_concentration_mean"] = float(
                concentration_mean_by_cn[cn_state][sample_idx]
            )
        for cn_state in _NON_DIPLOID_CN_STATES:
            row[f"cn{cn_state}_direct_peak_supported"] = bool(direct_support[cn_state])
            row[f"cn{cn_state}_posterior_supported"] = bool(
                posterior_supported[cn_state]
            )
            row[f"cn{cn_state}_posterior_and_peak_supported"] = bool(
                call_supported[cn_state]
            )
        row["diploid_log_lik"] = row["log_evidence_cn_2"]
        row["triploid_log_lik"] = row["log_evidence_cn_3"]
        row["diploid_af_concentration_map"] = row["cn2_af_concentration_map"]
        row["triploid_af_concentration_map"] = row["cn3_af_concentration_map"]
        row["diploid_af_concentration_mean"] = row["cn2_af_concentration_mean"]
        row["triploid_af_concentration_mean"] = row["cn3_af_concentration_mean"]
        rows.append(row)

    results_df = pd.DataFrame(rows)
    if log_progress:
        logger.info(
            "Polyploidy progress: classification complete %s",
            _format_category_counts(
                "autosomal_baseline_cn_counts",
                results_df["autosomal_baseline_cn"],
            ),
        )
        logger.info(
            "Polyploidy progress: classification call summary %s",
            _format_category_counts(
                "baseline_cn_call_counts",
                results_df["baseline_cn_call"],
            ),
        )
    return results_df


def _log_privacy_safe_quantiles_by_baseline(
    results_df: pd.DataFrame,
    column: str,
) -> None:
    if column not in results_df.columns:
        return
    for baseline_cn in sorted(results_df["autosomal_baseline_cn"].dropna().unique()):
        subset = results_df[results_df["autosomal_baseline_cn"] == baseline_cn]
        logger.info(
            "PRIVACY_SAFE_POLYPLOIDY result_quantiles_by_baseline_cn "
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
    thresholds = (0.5, 0.9, 0.99, 0.999, 0.9999)
    for cn_state in _NON_DIPLOID_CN_STATES:
        posterior = results_df[f"posterior_cn_{cn_state}"].to_numpy(dtype=np.float64)
        posterior_error = 1.0 - posterior
        log_bf = results_df[f"log_bayes_factor_cn{cn_state}_vs_cn2"].to_numpy(
            dtype=np.float64,
        )
        effect_size = results_df[f"effect_size_cn{cn_state}_vs_cn2"].to_numpy(
            dtype=np.float64,
        )
        direct = results_df[f"cn{cn_state}_direct_peak_supported"].to_numpy(dtype=bool)
        posterior_supported = results_df[f"cn{cn_state}_posterior_supported"].to_numpy(
            dtype=bool,
        )
        called = results_df["autosomal_baseline_cn"].to_numpy(dtype=np.int64) == cn_state
        posterior_parts = [
            f"posterior_ge_{threshold:g}={int(np.sum(posterior >= threshold))}"
            for threshold in thresholds
        ]
        logger.info(
            "PRIVACY_SAFE_POLYPLOIDY posterior_threshold_counts state=cn%d "
            "total=%d finite=%d %s posterior_error_le_call_threshold=%d "
            "log_bayes_factor_positive=%d effect_size_ge_threshold=%d "
            "direct_peak_supported=%d posterior_and_peak_supported=%d "
            "called=%d high_posterior_direct_failed=%d "
            "direct_supported_posterior_failed=%d",
            cn_state,
            int(results_df.shape[0]),
            int(np.isfinite(posterior).sum()),
            ", ".join(posterior_parts),
            int(np.sum(posterior_error <= pvalue_threshold)),
            int(np.sum(log_bf > 0.0)),
            int(np.sum(effect_size >= effect_size_threshold)),
            int(np.sum(direct)),
            int(np.sum(direct & posterior_supported)),
            int(np.sum(called)),
            int(np.sum(posterior_supported & ~direct)),
            int(np.sum(direct & ~posterior_supported)),
        )

    cn4_posterior = results_df["posterior_cn_4"].to_numpy(dtype=np.float64)
    cn4_direct = results_df["cn4_direct_peak_supported"].to_numpy(dtype=bool)
    diploid_compatible = results_df["fraction_near_diploid_compatible"].to_numpy(
        dtype=np.float64,
    )
    quarter_fraction = results_df["fraction_near_tetraploid_quarters"].to_numpy(
        dtype=np.float64,
    )
    logger.info(
        "PRIVACY_SAFE_POLYPLOIDY tetraploid_ambiguity_counts "
        "posterior_ge_0.5_direct_failed=%d posterior_ge_0.9_direct_failed=%d "
        "diploid_compatible_ge_0.9=%d quarter_fraction_zero_or_nan=%d "
        "cn4_called=%d",
        int(np.sum((cn4_posterior >= 0.5) & ~cn4_direct)),
        int(np.sum((cn4_posterior >= 0.9) & ~cn4_direct)),
        int(np.sum(diploid_compatible >= 0.9)),
        int(np.sum((~np.isfinite(quarter_fraction)) | (quarter_fraction <= 0.0))),
        int(np.sum(results_df["autosomal_baseline_cn"] == _TETRAPLOID_BASELINE_CN)),
    )


def _log_privacy_safe_decision_margins(results_df: pd.DataFrame) -> None:
    for cn_state in _NON_DIPLOID_CN_STATES:
        column = f"log_bayes_factor_cn{cn_state}_vs_cn2"
        called = results_df["autosomal_baseline_cn"] == cn_state
        logger.info(
            "PRIVACY_SAFE_POLYPLOIDY decision_margin_summary state=cn%d "
            "called=%d log_bayes_factor_positive=%d %s",
            cn_state,
            int(called.sum()),
            int(np.sum(results_df[column].to_numpy(dtype=np.float64) > 0.0)),
            _format_quantile_summary(column, results_df[column]),
        )


def _log_privacy_safe_concentration_boundaries(results_df: pd.DataFrame) -> None:
    if "af_concentration_grid" not in results_df.columns or results_df.empty:
        return
    grid_values = results_df["af_concentration_grid"].dropna().unique()
    if len(grid_values) != 1:
        logger.info(
            "PRIVACY_SAFE_POLYPLOIDY concentration_grid_summary "
            "skipped=non_unique_grid unique_grid_count=%d",
            len(grid_values),
        )
        return
    concentration_grid = _parse_af_concentration_grid(str(grid_values[0]))
    if concentration_grid is None or concentration_grid.size == 0:
        return
    grid_min = float(concentration_grid[0])
    grid_max = float(concentration_grid[-1])
    for cn_state in _BASELINE_CN_STATES:
        map_column = f"cn{cn_state}_af_concentration_map"
        mean_column = f"cn{cn_state}_af_concentration_mean"
        for label, mask in (
            ("all", np.ones(len(results_df), dtype=bool)),
            (
                f"called_cn{cn_state}",
                (results_df["autosomal_baseline_cn"] == cn_state).to_numpy(dtype=bool),
            ),
        ):
            subset = results_df.loc[mask]
            if subset.empty:
                continue
            map_values = subset[map_column].to_numpy(dtype=np.float64)
            logger.info(
                "PRIVACY_SAFE_POLYPLOIDY concentration_boundary_counts "
                "state=cn%d group=%s n=%d grid_size=%d grid_min=%.6g "
                "grid_max=%.6g map_at_min=%d map_at_max=%d",
                cn_state,
                label,
                int(subset.shape[0]),
                int(concentration_grid.size),
                grid_min,
                grid_max,
                int(np.sum(np.isclose(map_values, grid_min))),
                int(np.sum(np.isclose(map_values, grid_max))),
            )
        logger.info(
            "PRIVACY_SAFE_POLYPLOIDY concentration_mean_quantiles state=cn%d %s",
            cn_state,
            _format_quantile_summary(mean_column, results_df[mean_column]),
        )


def _log_privacy_safe_af_peak_metrics(metrics_df: pd.DataFrame) -> None:
    metric_columns = [
        "fraction_near_haploid_endpoints",
        "fraction_near_diploid_half",
        "fraction_near_triploid_thirds",
        "fraction_near_tetraploid_quarters",
        "fraction_near_diploid_compatible",
        "fraction_near_any_modeled_peak",
        "fraction_intermediate_third_to_half",
        "contamination_mixture_score",
        "raw_af_triploid_minus_diploid_fraction",
        "raw_af_tetraploid_quarter_minus_diploid_fraction",
        "raw_af_haploid_endpoint_minus_other_peak_fraction",
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
            "PRIVACY_SAFE_POLYPLOIDY af_peak_counts baseline_cn=%s n=%d "
            "haploid_endpoint_supported=%d triploid_peak_supported=%d "
            "tetraploid_quarter_supported=%d diploid_compatible_ge_0.9=%d "
            "mixed_half_and_thirds_ge_0.25=%d",
            baseline_cn,
            int(subset.shape[0]),
            int(np.sum(subset.get("cn1_direct_peak_supported", False))),
            int(np.sum(subset.get("cn3_direct_peak_supported", False))),
            int(np.sum(subset.get("cn4_direct_peak_supported", False))),
            int(np.sum(subset["fraction_near_diploid_compatible"] >= 0.9)),
            int(np.sum(subset["contamination_mixture_score"] >= 0.25)),
        )
        for column in metric_columns:
            if column not in subset.columns:
                continue
            logger.info(
                "PRIVACY_SAFE_POLYPLOIDY af_peak_quantiles_by_baseline_cn "
                "column=%s baseline_cn=%s %s",
                column,
                baseline_cn,
                _format_quantile_summary(column, subset[column]),
            )

    for cn_state in _NON_DIPLOID_CN_STATES:
        bf_column = f"log_bayes_factor_cn{cn_state}_vs_cn2"
        if bf_column not in metrics_df.columns:
            continue
        for column in [
            "fraction_near_haploid_endpoints",
            "fraction_near_diploid_half",
            "fraction_near_triploid_thirds",
            "fraction_near_tetraploid_quarters",
            "fraction_near_diploid_compatible",
            "raw_af_sd",
            "median_site_depth",
        ]:
            corr, n_used = _safe_correlation(metrics_df[column], metrics_df[bf_column])
            logger.info(
                "PRIVACY_SAFE_POLYPLOIDY af_peak_correlation_with_log_bayes_factor "
                "state=cn%d metric=%s n=%d corr=%.6g",
                cn_state,
                column,
                n_used,
                corr,
            )


def log_privacy_safe_polyploidy_diagnostics(
    results_df: pd.DataFrame,
    *,
    pvalue_threshold: float,
    effect_size_threshold: float,
    metrics_df: Optional[pd.DataFrame] = None,
) -> None:
    """Log aggregate-only baseline CN diagnostics safe for protected cohorts."""
    logger.info(
        "PRIVACY_SAFE_POLYPLOIDY privacy_contract aggregate_only=true "
        "forbidden=filenames,paths,sample_ids,raw_rows,sample_level_values"
    )
    logger.info(
        "PRIVACY_SAFE_POLYPLOIDY hypothesis "
        "multi_state_cn1_cn2_cn3_cn4_bayesian_peak_mixture"
    )
    logger.info(
        "PRIVACY_SAFE_POLYPLOIDY %s",
        _format_category_counts(
            "autosomal_baseline_cn_counts",
            results_df["autosomal_baseline_cn"],
        ),
    )
    logger.info(
        "PRIVACY_SAFE_POLYPLOIDY %s",
        _format_category_counts(
            "baseline_cn_call_counts",
            results_df["baseline_cn_call"],
        ),
    )
    logger.info(
        "PRIVACY_SAFE_POLYPLOIDY %s",
        _format_category_counts(
            "baseline_cn_reason_counts",
            results_df["baseline_cn_reason"],
        ),
    )
    logger.info(
        "PRIVACY_SAFE_POLYPLOIDY %s",
        _format_category_counts(
            "baseline_cn_map_counts",
            results_df["baseline_cn_map"],
        ),
    )
    _log_privacy_safe_threshold_counts(
        results_df,
        pvalue_threshold,
        effect_size_threshold,
    )
    _log_privacy_safe_decision_margins(results_df)
    for column in [
        "baseline_cn_posterior_probability",
        "baseline_cn_posterior_error_probability",
        "baseline_cn_effect_size_per_site",
        "posterior_cn_1",
        "posterior_cn_2",
        "posterior_cn_3",
        "posterior_cn_4",
        "log_bayes_factor_cn1_vs_cn2",
        "log_bayes_factor_cn3_vs_cn2",
        "log_bayes_factor_cn4_vs_cn2",
        "effect_size_cn1_vs_cn2",
        "effect_size_cn3_vs_cn2",
        "effect_size_cn4_vs_cn2",
        "fraction_near_haploid_endpoints",
        "fraction_near_diploid_half",
        "fraction_near_triploid_thirds",
        "fraction_near_tetraploid_quarters",
        "fraction_near_diploid_compatible",
        "triploid_peak_fraction_advantage",
        "tetraploid_quarter_peak_fraction_advantage",
        "haploid_endpoint_peak_fraction_advantage",
        "n_informative_bins",
        "n_informative_sites",
    ]:
        if column not in results_df.columns:
            continue
        logger.info(
            "PRIVACY_SAFE_POLYPLOIDY result_quantiles_all %s",
            _format_quantile_summary(column, results_df[column]),
        )
        _log_privacy_safe_quantiles_by_baseline(results_df, column)
    _log_privacy_safe_concentration_boundaries(results_df)

    if metrics_df is None:
        logger.info(
            "PRIVACY_SAFE_POLYPLOIDY af_peak_metrics not_computed "
            "rerun_with_privacy_safe_diagnostics=true"
        )
    else:
        _log_privacy_safe_af_peak_metrics(metrics_df)


def build_polyploidy_diagnostic_metrics(
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
    show_progress: bool = False,
) -> pd.DataFrame:
    """Summarize raw autosomal AF patterns used by the baseline CN model."""
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
    sample_indices = _progress_iter(
        range(len(sample_ids)),
        desc="Baseline CN diagnostic metrics",
        unit="sample",
        total=len(sample_ids),
        show_progress=show_progress,
    )
    for sample_idx in sample_indices:
        sample_id = sample_ids[sample_idx]
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
                    "fraction_near_haploid_endpoints": np.nan,
                    "fraction_near_diploid_half": np.nan,
                    "fraction_near_triploid_thirds": np.nan,
                    "fraction_near_tetraploid_quarters": np.nan,
                    "fraction_near_diploid_compatible": np.nan,
                    "fraction_near_any_modeled_peak": np.nan,
                    "fraction_intermediate_third_to_half": np.nan,
                    "contamination_mixture_score": np.nan,
                    "raw_af_triploid_minus_diploid_fraction": np.nan,
                    "raw_af_tetraploid_quarter_minus_diploid_fraction": np.nan,
                    "raw_af_haploid_endpoint_minus_other_peak_fraction": np.nan,
                }
            )
            continue

        folded_af = np.minimum(af_values, 1.0 - af_values)
        peak_metrics = _compute_observed_af_peak_metrics(
            observed_af[:, :, sample_idx],
            sample_mask,
            af_window,
        )
        intermediate = (
            (folded_af >= (1.0 / 3.0) + (af_window / 2.0)) &
            (folded_af <= 0.5 - (af_window / 2.0))
        )
        other_peak_fraction = max(
            peak_metrics["fraction_near_diploid_half"],
            peak_metrics["fraction_near_triploid_thirds"],
            peak_metrics["fraction_near_tetraploid_quarters"],
        )

        row = {
            "sample": str(sample_id),
            "diagnostic_informative_sites": int(af_values.size),
            "median_site_depth": float(np.median(depth_values)),
            "raw_af_mean": float(np.mean(af_values)),
            "raw_af_sd": float(np.std(af_values, ddof=0)),
            "folded_af_median": float(np.median(folded_af)),
            "folded_af_p10": float(np.quantile(folded_af, 0.10)),
            "folded_af_p90": float(np.quantile(folded_af, 0.90)),
            "fraction_intermediate_third_to_half": float(np.mean(intermediate)),
            "contamination_mixture_score": float(
                np.sqrt(
                    peak_metrics["fraction_near_diploid_half"] *
                    peak_metrics["fraction_near_triploid_thirds"]
                )
            ),
            "raw_af_triploid_minus_diploid_fraction": float(
                peak_metrics["fraction_near_triploid_thirds"] -
                peak_metrics["fraction_near_diploid_half"]
            ),
            "raw_af_tetraploid_quarter_minus_diploid_fraction": float(
                peak_metrics["fraction_near_tetraploid_quarters"] -
                peak_metrics["fraction_near_diploid_half"]
            ),
            "raw_af_haploid_endpoint_minus_other_peak_fraction": float(
                peak_metrics["fraction_near_haploid_endpoints"] - other_peak_fraction
            ),
        }
        row.update(peak_metrics)
        rows.append(row)

    metrics_df = pd.DataFrame(rows)
    result_columns = [
        "sample",
        "autosomal_baseline_cn",
        "baseline_cn_call",
        "baseline_cn_reason",
        "baseline_cn_posterior_probability",
        "baseline_cn_posterior_error_probability",
        "baseline_cn_effect_size_per_site",
        "triploidy_call",
        "triploidy_reason",
        "triploidy_pvalue",
        "triploidy_posterior_probability",
        "triploidy_log_bayes_factor",
        "triploidy_log_lik_ratio",
        "triploidy_effect_size_per_site",
        "posterior_cn_1",
        "posterior_cn_2",
        "posterior_cn_3",
        "posterior_cn_4",
        "log_bayes_factor_cn1_vs_cn2",
        "log_bayes_factor_cn3_vs_cn2",
        "log_bayes_factor_cn4_vs_cn2",
        "cn1_direct_peak_supported",
        "cn3_direct_peak_supported",
        "cn4_direct_peak_supported",
        "n_informative_bins",
        "n_informative_sites",
    ]
    existing_result_columns = [column for column in result_columns if column in results_df]
    return metrics_df.merge(
        results_df[existing_result_columns],
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
    posterior_cols = [f"posterior_cn_{cn}" for cn in _NON_DIPLOID_CN_STATES]
    ranked_df["_non_diploid_sort"] = (
        ranked_df["autosomal_baseline_cn"] != _DIPLOID_BASELINE_CN
    ).astype(int)
    ranked_df["_direct_support_sort"] = (
        ranked_df[[
            "cn1_direct_peak_supported",
            "cn3_direct_peak_supported",
            "cn4_direct_peak_supported",
        ]]
        .fillna(False)
        .astype(int)
        .sum(axis=1)
    )
    ranked_df["_max_non_diploid_posterior"] = ranked_df[posterior_cols].max(axis=1)
    ranked_df["_site_sort"] = ranked_df["diagnostic_informative_sites"].fillna(0)
    ranked_df = ranked_df.sort_values(
        [
            "_non_diploid_sort",
            "_direct_support_sort",
            "_max_non_diploid_posterior",
            "_site_sort",
        ],
        ascending=[False, False, False, False],
    )
    return [str(sample_id) for sample_id in ranked_df["sample"].head(sample_limit)]


def _select_ploidy_specific_diagnostic_samples(
    metrics_df: pd.DataFrame,
    baseline_cn: int,
    sample_limit: int,
) -> list[str]:
    if sample_limit <= 0 or metrics_df.empty:
        return []

    subset = metrics_df[
        (metrics_df["autosomal_baseline_cn"] == baseline_cn) &
        (metrics_df["baseline_cn_call"] == baseline_ploidy_label(baseline_cn))
    ].copy()
    if subset.empty:
        return []

    subset["_posterior_sort"] = subset["baseline_cn_posterior_probability"].fillna(0.0)
    subset["_site_sort"] = subset["diagnostic_informative_sites"].fillna(0)
    subset = subset.sort_values(
        ["_posterior_sort", "_site_sort", "sample"],
        ascending=[False, False, True],
    )
    return [str(sample_id) for sample_id in subset["sample"].head(sample_limit)]


def _select_diploid_diagnostic_samples(
    metrics_df: pd.DataFrame,
    sample_limit: int = _DEFAULT_DIAGNOSTIC_DIPLOID_SAMPLE_LIMIT,
) -> list[str]:
    return _select_ploidy_specific_diagnostic_samples(
        metrics_df,
        _DIPLOID_BASELINE_CN,
        sample_limit,
    )


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
    show_progress: bool = False,
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
    selected_iter = _progress_iter(
        selected_samples,
        desc="Sampling diagnostic AF sites",
        unit="sample",
        total=len(selected_samples),
        show_progress=show_progress,
    )
    for sample_id in selected_iter:
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
                    "folded_observed_af": float(min(observed_af, 1.0 - observed_af)),
                }
            )
    return pd.DataFrame(rows)


def _plot_polyploidy_diagnostic_metrics(
    metrics_df: pd.DataFrame,
    output_dir: str,
) -> None:
    if metrics_df.empty:
        return
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(13, 10))
    axes_arr = axes.ravel()
    color_by_call = {
        "HAPLOID": "#9467bd",
        "DIPLOID": "#1f77b4",
        "TRIPLOID": "#d62728",
        "TETRAPLOID": "#2ca02c",
        "INSUFFICIENT_DATA": "#7f7f7f",
    }

    def _scatter_by_call(ax, x_col, y_col, xlabel, ylabel, title):
        for call, sub_df in metrics_df.groupby("baseline_cn_call", dropna=False):
            call_label = str(call)
            ax.scatter(
                sub_df[x_col],
                sub_df[y_col],
                s=np.clip(
                    np.log10(sub_df["diagnostic_informative_sites"] + 1.0),
                    1.0,
                    5.0,
                ) * 20.0,
                alpha=0.75,
                color=color_by_call.get(call_label, "#8c564b"),
                edgecolor="white",
                linewidth=0.4,
                label=call_label,
            )
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.set_xlim(-0.02, 1.02)
        ax.set_ylim(-0.02, 1.02)
        ax.grid(True, alpha=0.25)

    _scatter_by_call(
        axes_arr[0],
        "fraction_near_diploid_half",
        "fraction_near_triploid_thirds",
        "Fraction near AF=0.5",
        "Fraction near AF=1/3 or 2/3",
        "Triploid peak evidence",
    )
    _scatter_by_call(
        axes_arr[1],
        "fraction_near_diploid_half",
        "fraction_near_tetraploid_quarters",
        "Fraction near AF=0.5",
        "Fraction near AF=1/4 or 3/4",
        "Tetraploid quarter-peak evidence",
    )
    _scatter_by_call(
        axes_arr[2],
        "fraction_near_diploid_compatible",
        "fraction_near_tetraploid_quarters",
        "Fraction near AF=0, 0.5, or 1",
        "Fraction near AF=1/4 or 3/4",
        "CN4 ambiguity check",
    )
    _scatter_by_call(
        axes_arr[3],
        "fraction_near_haploid_endpoints",
        "fraction_near_any_modeled_peak",
        "Fraction near AF=0 or 1",
        "Fraction near any modeled peak",
        "Haploid endpoint evidence",
    )
    axes_arr[0].legend(loc="best", fontsize=8)
    save_and_close_plot(
        output_dir,
        "polyploidy_diagnostic_metrics.png",
        subdir="diagnostics",
    )


def _plot_polyploidy_raw_af_profiles(
    raw_site_df: pd.DataFrame,
    metrics_df: pd.DataFrame,
    output_dir: str,
    *,
    baseline_cn: int,
    filename: str,
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
    peak_lines = [
        (0.25, "#2ca02c", "--", 1.0),
        (1.0 / 3.0, "#d62728", "--", 1.0),
        (0.5, "black", ":", 1.1),
        (2.0 / 3.0, "#d62728", "--", 1.0),
        (0.75, "#2ca02c", "--", 1.0),
    ]
    for ax, sample_id in zip(axes.ravel(), samples):
        sub_df = raw_site_df[raw_site_df["sample"] == sample_id]
        ax.hist(
            sub_df["observed_af"],
            bins=bins,
            color="#4c78a8",
            alpha=0.78,
            edgecolor="none",
        )
        for xpos, color, linestyle, linewidth in peak_lines:
            ax.axvline(xpos, color=color, linestyle=linestyle, linewidth=linewidth)
        ax.set_xlim(0.0, 1.0)
        ax.set_xlabel("Raw observed alt fraction")
        ax.set_ylabel("Sites")
        title = str(sample_id)
        if sample_id in metrics_by_sample.index:
            row = metrics_by_sample.loc[sample_id]
            title = (
                f"{sample_id}\n{row['baseline_cn_call']}, "
                f"P={row['baseline_cn_posterior_probability']:.2f}"
            )
        ax.set_title(title, fontsize=9)
        ax.grid(True, axis="y", alpha=0.2)

    label = {
        1: "Monoploid",
        2: "Diploid",
        3: "Triploid",
        4: "Tetraploid",
    }.get(baseline_cn, f"CN{baseline_cn}")
    fig.suptitle(f"{label} Raw AF Profiles", fontsize=12)
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.96))

    for ax in axes.ravel()[len(samples):]:
        ax.axis("off")
    save_and_close_plot(
        output_dir,
        filename,
        subdir="diagnostics",
    )


def write_polyploidy_diagnostics(
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
    show_progress: bool = False,
) -> None:
    """Write baseline CN diagnostic metrics and plots."""
    logger.info("Polyploidy progress: computing diagnostic AF metrics.")
    metrics_df = build_polyploidy_diagnostic_metrics(
        sample_ids=sample_ids,
        bin_chr=bin_chr,
        site_alt=site_alt,
        site_total=site_total,
        site_pop_af=site_pop_af,
        site_mask=site_mask,
        results_df=results_df,
        min_diploid_het_prior=min_diploid_het_prior,
        af_window=af_window,
        show_progress=show_progress,
    )
    diagnostics_dir = os.path.join(output_dir, "diagnostics")
    os.makedirs(diagnostics_dir, exist_ok=True)
    metrics_path = os.path.join(diagnostics_dir, "polyploidy_diagnostic_metrics.tsv")
    metrics_df.to_csv(metrics_path, sep="\t", index=False)
    logger.info("Polyploidy diagnostic metrics saved.")

    selected_samples = _select_diagnostic_samples(metrics_df, sample_limit)
    logger.info(
        "Polyploidy progress: sampling diagnostic AF sites selected_samples=%d.",
        int(len(selected_samples)),
    )
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
        show_progress=show_progress,
    )
    raw_site_path = os.path.join(diagnostics_dir, "polyploidy_raw_af_sites.tsv.gz")
    raw_site_df.to_csv(raw_site_path, sep="\t", index=False)
    logger.info("Polyploidy sampled raw AF sites saved.")

    raw_site_frames_by_cn: dict[int, pd.DataFrame] = {}
    for baseline_cn in _BASELINE_CN_STATES:
        ploidy_samples = _select_ploidy_specific_diagnostic_samples(
            metrics_df,
            baseline_cn,
            sample_limit,
        )
        logger.info(
            "Polyploidy progress: sampling CN%d AF profile sites selected_samples=%d.",
            baseline_cn,
            int(len(ploidy_samples)),
        )
        raw_site_frames_by_cn[baseline_cn] = _collect_diagnostic_site_rows(
            sample_ids=sample_ids,
            selected_samples=ploidy_samples,
            bin_chr=bin_chr,
            site_alt=site_alt,
            site_total=site_total,
            site_pop_af=site_pop_af,
            site_mask=site_mask,
            min_diploid_het_prior=min_diploid_het_prior,
            max_sites_per_sample=max_sites_per_sample,
            show_progress=show_progress,
        )

    logger.info("Polyploidy progress: rendering diagnostic plots.")
    _plot_polyploidy_diagnostic_metrics(metrics_df, output_dir)
    for baseline_cn, filename in (
        (1, "polyploidy_raw_af_monoploid_profiles.png"),
        (2, "polyploidy_raw_af_diploid_profiles.png"),
        (3, "polyploidy_raw_af_triploid_profiles.png"),
        (4, "polyploidy_raw_af_tetraploid_profiles.png"),
    ):
        ploidy_site_df = raw_site_frames_by_cn[baseline_cn]
        if ploidy_site_df.empty:
            continue
        _plot_polyploidy_raw_af_profiles(
            ploidy_site_df,
            metrics_df,
            output_dir,
            baseline_cn=baseline_cn,
            filename=filename,
        )
    logger.info("Polyploidy diagnostics saved.")


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the polyploidy subcommand."""
    parser = argparse.ArgumentParser(
        description="Classify autosomal baseline CN from pooled allele counts",
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
    parser.add_argument("-o", "--output-dir", required=True, help="Output directory")
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
        help="Log-scale standard deviation for the AF concentration prior",
    )
    parser.add_argument(
        "--haploidy-prior",
        type=float,
        default=_DEFAULT_HAPLOIDY_PRIOR,
        help="Prior probability that any one sample has autosomal baseline CN=1",
    )
    parser.add_argument(
        "--triploidy-prior",
        type=float,
        default=_DEFAULT_TRIPLOIDY_PRIOR,
        help="Prior probability that any one sample has autosomal baseline CN=3",
    )
    parser.add_argument(
        "--tetraploidy-prior",
        type=float,
        default=_DEFAULT_TETRAPLOIDY_PRIOR,
        help="Prior probability that any one sample has autosomal baseline CN=4",
    )
    parser.add_argument(
        "--peak-evidence-weight",
        type=float,
        default=_DEFAULT_PEAK_EVIDENCE_WEIGHT,
        help="Composite-likelihood weight for direct AF peak mixture evidence",
    )
    parser.add_argument(
        "--peak-kernel-sd",
        type=float,
        default=_DEFAULT_PEAK_KERNEL_SD,
        help="Standard deviation of direct AF peak mixture kernels",
    )
    parser.add_argument(
        "--peak-outlier-probability",
        type=float,
        default=_DEFAULT_PEAK_OUTLIER_PROBABILITY,
        help="Uniform outlier mixture weight in the direct AF peak model",
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
        help="Maximum posterior error probability for non-diploid calls",
    )
    parser.add_argument(
        "--effect-size-threshold",
        type=float,
        default=0.01,
        help="Minimum mean log Bayes factor per informative site for non-diploid calls",
    )
    parser.add_argument(
        "--ploidy-peak-af-window",
        "--polyploidy-peak-af-window",
        "--triploidy-peak-af-window",
        dest="ploidy_peak_af_window",
        type=float,
        default=_DEFAULT_PLOIDY_PEAK_AF_WINDOW,
        help="AF window around canonical CN1/CN2/CN3/CN4 peaks used by call guards",
    )
    parser.add_argument(
        "--min-haploid-endpoint-fraction",
        type=float,
        default=_DEFAULT_MIN_HAPLOID_ENDPOINT_FRACTION,
        help="Minimum observed AF fraction near 0 or 1 for calling CN=1",
    )
    parser.add_argument(
        "--max-haploid-other-peak-fraction",
        type=float,
        default=_DEFAULT_MAX_HAPLOID_OTHER_PEAK_FRACTION,
        help="Maximum observed AF fraction near non-endpoint peaks for calling CN=1",
    )
    parser.add_argument(
        "--min-triploid-peak-fraction-advantage",
        type=float,
        default=_DEFAULT_MIN_TRIPLOID_PEAK_FRACTION_ADVANTAGE,
        help="Minimum excess of thirds-peak fraction over half-peak fraction for CN=3",
    )
    parser.add_argument(
        "--min-tetraploid-quarter-peak-fraction",
        type=float,
        default=_DEFAULT_MIN_TETRAPLOID_QUARTER_PEAK_FRACTION,
        help="Minimum observed AF fraction near 1/4 or 3/4 for calling CN=4",
    )
    parser.add_argument(
        "--min-tetraploid-quarter-peak-fraction-advantage",
        type=float,
        default=_DEFAULT_MIN_TETRAPLOID_QUARTER_PEAK_FRACTION_ADVANTAGE,
        help="Minimum excess of quarter-peak fraction over half-peak fraction for CN=4",
    )
    parser.add_argument(
        "--diagnostics",
        action="store_true",
        help="Write baseline CN diagnostic metrics, plots, and sampled raw AF data",
    )
    parser.add_argument(
        "--privacy-safe-diagnostics",
        action="store_true",
        help="Log aggregate-only baseline CN diagnostics safe for protected cohorts",
    )
    parser.add_argument(
        "--diagnostic-sample-limit",
        type=int,
        default=_DEFAULT_DIAGNOSTIC_SAMPLE_LIMIT,
        help="Maximum number of samples per ploidy-specific raw AF diagnostic plot",
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
        help="AF window around canonical CN1/CN2/CN3/CN4 peaks for diagnostics",
    )
    parser.add_argument(
        "--no-progress",
        action="store_true",
        help="Disable progress bars while keeping aggregate stage logs",
    )
    return parser.parse_args()


def main() -> None:
    """Entry point for ``gatk-sv-ploidy polyploidy``."""
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    show_progress = not args.no_progress

    logger.info("Polyploidy progress: loading depth table and site data.")
    depth_df = pd.read_csv(args.input, sep="\t", index_col=0)
    sample_ids = [str(sample_id) for sample_id in get_sample_columns(depth_df)]

    site_npz = np.load(args.site_data, allow_pickle=True)
    site_alt = np.asarray(site_npz["site_alt"], dtype=np.int32)
    site_total = np.asarray(site_npz["site_total"], dtype=np.int32)
    site_pop_af = np.asarray(site_npz["site_pop_af"], dtype=np.float64)
    site_mask = np.asarray(site_npz["site_mask"], dtype=bool)
    logger.info(
        "Polyploidy progress: loaded aggregate inputs samples=%d bins=%d "
        "max_sites_per_bin=%d valid_site_observations=%d.",
        int(len(sample_ids)),
        int(site_alt.shape[0]),
        int(site_alt.shape[1]),
        int(site_mask.sum()),
    )

    if "sample_ids" in site_npz:
        npz_sample_ids = [
            str(sample_id) for sample_id in np.asarray(site_npz["sample_ids"]).tolist()
        ]
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

    af_concentration_grid = _parse_af_concentration_grid(args.af_concentration_grid)
    logger.info("Polyploidy progress: starting autosomal baseline CN classification.")
    results_df = classify_polyploidy_from_site_data(
        sample_ids=sample_ids,
        bin_chr=bin_chr,
        site_alt=site_alt,
        site_total=site_total,
        site_pop_af=site_pop_af,
        site_mask=site_mask,
        af_concentration=args.af_concentration,
        af_concentration_grid=af_concentration_grid,
        af_concentration_prior_log_sd=args.af_concentration_prior_log_sd,
        haploidy_prior=args.haploidy_prior,
        triploidy_prior=args.triploidy_prior,
        tetraploidy_prior=args.tetraploidy_prior,
        peak_evidence_weight=args.peak_evidence_weight,
        peak_kernel_sd=args.peak_kernel_sd,
        peak_outlier_probability=args.peak_outlier_probability,
        min_diploid_het_prior=args.min_diploid_het_prior,
        min_informative_bins=args.min_informative_bins,
        min_informative_sites=args.min_informative_sites,
        pvalue_threshold=args.pvalue_threshold,
        effect_size_threshold=args.effect_size_threshold,
        peak_af_window=args.ploidy_peak_af_window,
        min_haploid_endpoint_fraction=args.min_haploid_endpoint_fraction,
        max_haploid_other_peak_fraction=args.max_haploid_other_peak_fraction,
        min_triploid_peak_fraction_advantage=(
            args.min_triploid_peak_fraction_advantage
        ),
        min_tetraploid_quarter_peak_fraction=(
            args.min_tetraploid_quarter_peak_fraction
        ),
        min_tetraploid_quarter_peak_fraction_advantage=(
            args.min_tetraploid_quarter_peak_fraction_advantage
        ),
        show_progress=show_progress,
        log_progress=True,
    )

    privacy_safe_metrics_df = None
    if args.privacy_safe_diagnostics:
        logger.info("Polyploidy progress: computing privacy-safe AF diagnostics.")
        privacy_safe_metrics_df = build_polyploidy_diagnostic_metrics(
            sample_ids=sample_ids,
            bin_chr=bin_chr,
            site_alt=site_alt,
            site_total=site_total,
            site_pop_af=site_pop_af,
            site_mask=site_mask,
            results_df=results_df,
            min_diploid_het_prior=args.min_diploid_het_prior,
            af_window=args.diagnostic_af_window,
            show_progress=show_progress,
        )
    logger.info("Polyploidy progress: logging aggregate privacy-safe summary.")
    log_privacy_safe_polyploidy_diagnostics(
        results_df,
        pvalue_threshold=args.pvalue_threshold,
        effect_size_threshold=args.effect_size_threshold,
        metrics_df=privacy_safe_metrics_df,
    )

    logger.info("Polyploidy progress: writing result tables.")
    results_path = os.path.join(args.output_dir, "polyploidy_test_results.tsv")
    results_df.to_csv(results_path, sep="\t", index=False)
    logger.info("Polyploidy test results saved.")

    manifest_df = results_df[["sample", "autosomal_baseline_cn"]].copy()
    manifest_path = os.path.join(args.output_dir, "sample_autosomal_baseline_cn.tsv")
    manifest_df.to_csv(manifest_path, sep="\t", index=False)
    logger.info("Sample autosomal baseline CN manifest saved.")

    if args.diagnostics:
        write_polyploidy_diagnostics(
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
            show_progress=show_progress,
        )


if __name__ == "__main__":
    main()
