"""
Infer subcommand — train model and run copy-number inference.

Loads preprocessed depth data, trains the Pyro CNV model, obtains MAP
estimates and exact discrete posteriors, then writes per-bin and
per-chromosome summary statistics.
"""

from __future__ import annotations

import argparse
import logging
import os
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import pyro
import torch
from scipy import stats

from gatk_sv_ploidy._util import (
    compose_additive_background_matrix,
    DEFAULT_AF_WEIGHT,
    DEPTH_SPACES,
    compute_cnq_from_probabilities,
    default_obs_likelihood_for_depth_space,
    format_count_summary,
    format_numeric_summary,
    read_observation_type,
    summarize_contig_ploidy_from_bin_calls,
    validate_depth_space,
)
from gatk_sv_ploidy.data import DepthData, load_site_data
from gatk_sv_ploidy.models import (
    AUTOSOME_PRIOR_MODES,
    DEFAULT_BACKGROUND_FACTORS,
    DEFAULT_EPSILON_CONCENTRATION,
    DEFAULT_EPSILON_MEAN,
    OBS_LIKELIHOODS,
    CNVModel,
    _resolve_fixed_site_pop_af_numpy,
    _site_level_marginalized_af_log_lik_numpy,
)

logger = logging.getLogger(__name__)
SITE_AF_ESTIMATORS = ("off", "auto", "naive-bayes")


def _gini_coefficient(values: np.ndarray) -> float | None:
    """Return the Gini coefficient for non-negative finite values."""
    arr = np.asarray(values, dtype=np.float64).reshape(-1)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return None
    arr = np.maximum(arr, 0.0)
    total = float(arr.sum())
    if total <= 0.0:
        return 0.0
    sorted_arr = np.sort(arr)
    n = sorted_arr.size
    index = np.arange(1, n + 1, dtype=np.float64)
    gini = (2.0 * np.sum(index * sorted_arr) / (n * total)) - ((n + 1.0) / n)
    return float(gini)


def _format_concentration_summary(
    label: str,
    values: np.ndarray,
    top_ks: tuple[int, ...] = (1, 3, 5),
) -> str:
    """Summarize whether an aggregate burden is diffuse or sample-concentrated."""
    arr = np.asarray(values, dtype=np.float64).reshape(-1)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return f"{label}: unavailable"

    arr = np.maximum(arr, 0.0)
    total = float(arr.sum())
    if total <= 0.0:
        return f"{label}: total=0"

    sorted_arr = np.sort(arr)[::-1]
    parts = []
    for top_k in top_ks:
        top_share = float(sorted_arr[: min(top_k, sorted_arr.size)].sum() / total)
        parts.append(f"top{top_k}_share={top_share:.1%}")

    gini = _gini_coefficient(arr)
    if gini is not None:
        parts.append(f"gini={gini:.3f}")
    return f"{label}: total={total:.4f}, " + ", ".join(parts)


def _format_transition_count_summary(
    label: str,
    from_states: np.ndarray,
    to_states: np.ndarray,
    state_names: list[str],
    top_n: int = 6,
) -> str:
    """Format aggregate state-transition counts without exposing records."""
    from_arr = np.asarray(from_states, dtype=np.int64).reshape(-1)
    to_arr = np.asarray(to_states, dtype=np.int64).reshape(-1)
    mask = np.isfinite(from_arr) & np.isfinite(to_arr)
    from_arr = from_arr[mask]
    to_arr = to_arr[mask]
    if from_arr.size == 0:
        return f"{label}: unavailable"

    n_states = len(state_names)
    counts = np.zeros((n_states, n_states), dtype=np.int64)
    valid = (
        (from_arr >= 0) & (from_arr < n_states) &
        (to_arr >= 0) & (to_arr < n_states)
    )
    np.add.at(counts, (from_arr[valid], to_arr[valid]), 1)

    total = int(counts.sum())
    shifted = int(total - np.trace(counts))
    if shifted == 0:
        return f"{label}: shifted=0/{total} (0.0%)"

    transitions: list[tuple[int, int, int]] = []
    for from_state in range(n_states):
        for to_state in range(n_states):
            if from_state == to_state or counts[from_state, to_state] == 0:
                continue
            transitions.append((from_state, to_state, int(counts[from_state, to_state])))
    transitions.sort(key=lambda item: item[2], reverse=True)
    top_terms = ", ".join(
        f"{state_names[from_state]}->{state_names[to_state]}={count} ({count / shifted:.1%})"
        for from_state, to_state, count in transitions[: min(top_n, len(transitions))]
    )
    return f"{label}: shifted={shifted}/{total} ({shifted / total:.1%}), top={top_terms}"


def _safe_correlation(x: np.ndarray, y: np.ndarray) -> float | None:
    """Return a stable Pearson correlation when both inputs vary."""
    x_arr = np.asarray(x, dtype=np.float64).reshape(-1)
    y_arr = np.asarray(y, dtype=np.float64).reshape(-1)
    mask = np.isfinite(x_arr) & np.isfinite(y_arr)
    if mask.sum() < 3:
        return None
    x_finite = x_arr[mask]
    y_finite = y_arr[mask]
    if np.std(x_finite) < 1e-12 or np.std(y_finite) < 1e-12:
        return None
    return float(np.corrcoef(x_finite, y_finite)[0, 1])


def _coerce_bin_state_matrix(values: np.ndarray, n_bins: int) -> np.ndarray | None:
    """Best-effort reshape of bin-level CN probabilities to ``(n_bins, n_states)``."""
    matrix = np.asarray(values, dtype=np.float64).squeeze()
    if matrix.ndim == 1:
        matrix = matrix.reshape(1, -1)
    elif matrix.ndim > 2:
        matrix = matrix.reshape(matrix.shape[0], -1, matrix.shape[-1]).mean(axis=1)

    if matrix.ndim != 2 or matrix.shape[0] != n_bins:
        return None
    return matrix


def _coerce_background_factor_matrices(
    map_estimates: Dict[str, np.ndarray],
    n_bins: int,
    n_samples: int,
) -> tuple[np.ndarray | None, np.ndarray | None]:
    """Return normalized bin loadings and sample amplitudes when present."""
    if "background_bin_factors" not in map_estimates or "background_sample_factors" not in map_estimates:
        return None, None

    bin_factors = np.asarray(
        map_estimates["background_bin_factors"],
        dtype=np.float64,
    ).squeeze()
    if bin_factors.ndim == 1:
        bin_factors = bin_factors[:, np.newaxis]

    sample_factors = np.asarray(
        map_estimates["background_sample_factors"],
        dtype=np.float64,
    ).squeeze()
    if sample_factors.ndim == 1:
        sample_factors = sample_factors[np.newaxis, :]

    if bin_factors.ndim != 2 or sample_factors.ndim != 2:
        return None, None
    if bin_factors.shape[0] != n_bins:
        return None, None

    if sample_factors.shape[0] != bin_factors.shape[1]:
        if sample_factors.shape[1] == bin_factors.shape[1] and sample_factors.shape[0] == n_samples:
            sample_factors = sample_factors.T
        else:
            return None, None
    if sample_factors.shape[1] != n_samples:
        return None, None

    normalized_bin_factors = bin_factors / np.maximum(
        bin_factors.mean(axis=0, keepdims=True),
        1e-8,
    )
    return normalized_bin_factors, sample_factors


def _format_threshold_exceedance_summary(
    label: str,
    values: np.ndarray,
    thresholds: tuple[float, ...],
) -> str:
    """Format threshold exceedance counts for cohort-level diagnostics."""
    values_arr = np.asarray(values, dtype=np.float64).reshape(-1)
    finite = values_arr[np.isfinite(values_arr)]
    if finite.size == 0:
        return f"{label}: n=0"

    parts = []
    for threshold in thresholds:
        count = int((finite >= threshold).sum())
        parts.append(f">={threshold:.0%}:{count}/{finite.size} ({count / finite.size:.1%})")
    return f"{label}: " + ", ".join(parts)


def _format_numeric_threshold_summary(
    label: str,
    values: np.ndarray,
    thresholds: tuple[float, ...],
    precision: int = 2,
) -> str:
    """Format numeric threshold exceedances for non-fractional values."""
    values_arr = np.asarray(values, dtype=np.float64).reshape(-1)
    finite = values_arr[np.isfinite(values_arr)]
    if finite.size == 0:
        return f"{label}: n=0"

    parts = []
    for threshold in thresholds:
        count = int((finite >= threshold).sum())
        parts.append(
            f">={threshold:.{precision}f}:{count}/{finite.size} ({count / finite.size:.1%})"
        )
    return f"{label}: " + ", ".join(parts)


def _format_numeric_below_threshold_summary(
    label: str,
    values: np.ndarray,
    thresholds: tuple[float, ...],
    precision: int = 2,
) -> str:
    """Format numeric lower-tail threshold counts for privacy-safe diagnostics."""
    values_arr = np.asarray(values, dtype=np.float64).reshape(-1)
    finite = values_arr[np.isfinite(values_arr)]
    if finite.size == 0:
        return f"{label}: n=0"

    parts = []
    for threshold in thresholds:
        count = int((finite <= threshold).sum())
        parts.append(
            f"<={threshold:.{precision}f}:{count}/{finite.size} ({count / finite.size:.1%})"
        )
    return f"{label}: " + ", ".join(parts)


def _nearest_af_grid_distance(
    observed_af: np.ndarray,
    grid: tuple[float, ...],
) -> np.ndarray:
    """Return the distance from each AF value to the nearest genotype grid point."""
    grid_arr = np.asarray(grid, dtype=np.float64)
    return np.min(np.abs(np.asarray(observed_af, dtype=np.float64)[..., np.newaxis] - grid_arr), axis=-1)


def _format_boolean_fraction_summary(label: str, mask: np.ndarray) -> str:
    """Format a boolean mask as a cohort-level count and fraction."""
    mask_arr = np.asarray(mask, dtype=bool).reshape(-1)
    if mask_arr.size == 0:
        return f"{label}: n=0"
    count = int(mask_arr.sum())
    return f"{label}: {count}/{mask_arr.size} ({count / mask_arr.size:.1%})"


def _format_chromosome_rate_summary(
    label: str,
    chromosomes: np.ndarray,
    rates: np.ndarray,
    top_n: int = 5,
) -> str:
    """Format chromosome-level rate summaries using only cohort aggregates."""
    chrom_arr = np.asarray(chromosomes, dtype=object).reshape(-1)
    rate_arr = np.asarray(rates, dtype=np.float64).reshape(-1)
    mask = np.isfinite(rate_arr)
    if mask.sum() == 0:
        return f"{label}: unavailable"

    chrom_finite = chrom_arr[mask]
    rate_finite = rate_arr[mask]
    order = np.argsort(rate_finite)[::-1]
    top_terms = ", ".join(
        f"{chrom_finite[idx]}={rate_finite[idx]:.2%}"
        for idx in order[:min(top_n, rate_finite.size)]
    )
    return (
        f"{label}: median={np.quantile(rate_finite, 0.50):.2%}, "
        f"p95={np.quantile(rate_finite, 0.95):.2%}, "
        f"top={top_terms}"
    )


def build_safe_inference_diagnostic_messages(
    data: DepthData,
    map_estimates: Dict[str, np.ndarray],
    cn_posterior: Dict[str, np.ndarray],
    depth_only_cn_posterior: Dict[str, np.ndarray] | None = None,
    af_table: np.ndarray | None = None,
    af_concentration: float = 50.0,
) -> list[str]:
    """Build privacy-safe cohort diagnostics for troubleshooting model fits."""
    messages: list[str] = [
        "Safe logging enabled: cohort-level summaries only; sample identifiers and per-sample records are suppressed.",
    ]

    chr_type = data.chr_type.detach().cpu().numpy()
    autosome_mask = chr_type == 0
    chr_x_bins = int((chr_type == 1).sum())
    chr_y_bins = int((chr_type == 2).sum())
    messages.append(
        "Cohort dimensions: "
        f"bins={data.n_bins}, samples={data.n_samples}, "
        f"autosomal_bins={int(autosome_mask.sum())}, chrX_bins={chr_x_bins}, chrY_bins={chr_y_bins}"
    )

    if data.depth_space == "raw":
        messages.append(
            format_numeric_summary(
                "Empirical autosomal raw counts/kb across samples",
                _plot_depth_baseline_per_sample(data),
            )
        )

    for key, label in (
        ("sample_var", "Sample variance latent across samples"),
        ("sample_depth", "Sample depth latent across samples"),
        ("bin_bias", "Bin bias latent across bins"),
        ("bin_var", "Bin variance latent across bins"),
    ):
        if key in map_estimates:
            messages.append(format_numeric_summary(label, map_estimates[key]))

    background_keys = {
        "bin_epsilon",
        "background_bin_factors",
        "background_sample_factors",
    }
    bin_epsilon_matrix = None
    background_bin_loadings = None
    background_sample_factors = None
    reported_background_factors: range = range(0)
    if any(key in map_estimates for key in background_keys):
        bin_epsilon_matrix = compose_additive_background_matrix(
            map_estimates.get("bin_epsilon"),
            data.n_bins,
            data.n_samples,
            contig_index=data.contig_index.detach().cpu().numpy(),
            n_contigs=data.n_contigs,
            background_bin_factors=map_estimates.get("background_bin_factors"),
            background_sample_factors=map_estimates.get("background_sample_factors"),
        )
        messages.append(
            format_numeric_summary(
                "Bin epsilon latent across bin-sample pairs",
                bin_epsilon_matrix,
            )
        )
        background_bin_loadings, background_sample_factors = _coerce_background_factor_matrices(
            map_estimates,
            data.n_bins,
            data.n_samples,
        )
        if background_bin_loadings is not None and background_sample_factors is not None:
            max_report_factors = min(background_bin_loadings.shape[1], 6)
            reported_background_factors = range(max_report_factors)
            if background_bin_loadings.shape[1] > max_report_factors:
                messages.append(
                    "Background factor summaries truncated to first "
                    f"{max_report_factors} of {background_bin_loadings.shape[1]} factors."
                )
            messages.append(
                format_numeric_summary(
                    "Background sample-factor amplitudes across factor-sample pairs",
                    background_sample_factors,
                    precision=4,
                )
            )
            for factor_idx in reported_background_factors:
                messages.append(
                    format_numeric_summary(
                        f"Background factor {factor_idx + 1} normalized bin loadings across bins",
                        background_bin_loadings[:, factor_idx],
                        precision=4,
                    )
                )

    depth_ratio = None
    if all((data.depth_space == "raw", "sample_depth" in map_estimates)):
        empirical_baseline = _plot_depth_baseline_per_sample(data)
        sample_depth = np.asarray(map_estimates["sample_depth"], dtype=np.float64)
        depth_ratio = np.divide(
            sample_depth,
            empirical_baseline,
            out=np.full_like(sample_depth, np.nan, dtype=np.float64),
            where=np.asarray(empirical_baseline) > 0,
        )
        messages.append(
            format_numeric_summary(
                "Sample depth / empirical autosomal baseline across samples",
                depth_ratio,
            )
        )

    cn_probs = np.asarray(cn_posterior["cn_posterior"], dtype=np.float64)
    cn_map = np.argmax(cn_probs, axis=2)
    cn_state_names = [f"CN{state}" for state in range(cn_probs.shape[2])]

    for ct_value, chrom_label in ((1, "chrX"), (2, "chrY")):
        chrom_mask = chr_type == ct_value
        if not np.any(chrom_mask):
            continue

        if "bin_var" in map_estimates:
            messages.append(
                format_numeric_summary(
                    f"Bin variance latent across {chrom_label} bins",
                    np.asarray(map_estimates["bin_var"], dtype=np.float64).reshape(-1)[chrom_mask],
                )
            )
        if "bin_bias" in map_estimates:
            messages.append(
                format_numeric_summary(
                    f"Bin bias latent across {chrom_label} bins",
                    np.asarray(map_estimates["bin_bias"], dtype=np.float64).reshape(-1)[chrom_mask],
                )
            )
        if bin_epsilon_matrix is not None:
            messages.append(
                format_numeric_summary(
                    f"Bin epsilon latent across {chrom_label} bin-sample pairs",
                    bin_epsilon_matrix[chrom_mask, :],
                )
            )
        if background_bin_loadings is not None:
            for factor_idx in reported_background_factors:
                messages.append(
                    format_numeric_summary(
                        f"Background factor {factor_idx + 1} normalized loadings across {chrom_label} bins",
                        background_bin_loadings[chrom_mask, factor_idx],
                        precision=4,
                    )
                )

    cn_prob_matrix = None
    if "cn_probs" in map_estimates:
        cn_prob_matrix = _coerce_bin_state_matrix(
            map_estimates["cn_probs"],
            data.n_bins,
        )
        if cn_prob_matrix is not None:
            for ct_value, chrom_label, states in (
                (1, "chrX", (1, 2)),
                (2, "chrY", (0, 1)),
            ):
                chrom_mask = chr_type == ct_value
                if not np.any(chrom_mask):
                    continue
                for state in states:
                    if cn_prob_matrix.shape[1] <= state:
                        continue
                    messages.append(
                        format_numeric_summary(
                            f"Per-bin CN{state} prior mass across {chrom_label} bins",
                            cn_prob_matrix[chrom_mask, state],
                            precision=4,
                        )
                    )

    if data.depth_space == "raw":
        observed_plot_depth = _plot_normalized_depth(
            data.depth.detach().cpu().numpy(),
            data.bin_length_kb.detach().cpu().numpy()[:, np.newaxis],
            _plot_depth_baseline_per_sample(data)[np.newaxis, :],
        )
    else:
        observed_plot_depth = data.depth.detach().cpu().numpy()

    if "sex_posterior" in cn_posterior:
        sex_post_for_allosomes = np.asarray(
            cn_posterior["sex_posterior"],
            dtype=np.float64,
        )
        if (
            sex_post_for_allosomes.ndim == 2 and
            sex_post_for_allosomes.shape[0] == data.n_samples and
            sex_post_for_allosomes.shape[1] >= 2
        ):
            sex_map_for_allosomes = np.argmax(sex_post_for_allosomes, axis=1)
            expected_cn_by_chrom_and_sex = {
                "chrX": {0: 2, 1: 1},
                "chrY": {0: 0, 1: 1},
            }
            bin_bias_arr = None
            if "bin_bias" in map_estimates:
                bin_bias_arr = np.asarray(
                    map_estimates["bin_bias"],
                    dtype=np.float64,
                ).reshape(-1)

            depth_ratio_for_expected = None
            if all((data.depth_space == "raw", "sample_depth" in map_estimates)):
                empirical_baseline = _plot_depth_baseline_per_sample(data)
                sample_depth = np.asarray(map_estimates["sample_depth"], dtype=np.float64)
                depth_ratio_for_expected = np.divide(
                    sample_depth,
                    empirical_baseline,
                    out=np.full_like(sample_depth, np.nan, dtype=np.float64),
                    where=np.asarray(empirical_baseline) > 0,
                )

            for ct_value, chrom_label in ((1, "chrX"), (2, "chrY")):
                bin_idx = np.flatnonzero(chr_type == ct_value)
                if bin_idx.size == 0:
                    continue
                for sex_state, sex_label in ((0, "XX"), (1, "XY")):
                    sample_idx = np.flatnonzero(sex_map_for_allosomes == sex_state)
                    if sample_idx.size == 0:
                        continue

                    cn_subset = cn_map[np.ix_(bin_idx, sample_idx)]
                    messages.append(
                        format_count_summary(
                            f"{chrom_label} dominant CN across {sex_label}-assigned bin-sample pairs",
                            np.bincount(
                                cn_subset.reshape(-1),
                                minlength=cn_probs.shape[2],
                            ),
                            cn_state_names,
                        )
                    )

                    expected_cn = expected_cn_by_chrom_and_sex[chrom_label][sex_state]
                    unexpected_mask = cn_subset != expected_cn
                    messages.append(
                        _format_boolean_fraction_summary(
                            f"{chrom_label} unexpected dominant CN among {sex_label}-assigned bin-sample pairs",
                            unexpected_mask,
                        )
                    )
                    unexpected_by_bin = unexpected_mask.mean(axis=1)
                    messages.append(
                        _format_threshold_exceedance_summary(
                            f"{chrom_label} bins with unexpected dominant CN sample fraction among {sex_label}-assigned samples thresholds",
                            unexpected_by_bin,
                            (0.10, 0.50, 0.90),
                        )
                    )

                    if "bin_var" in map_estimates:
                        sex_chrom_bin_var = np.asarray(
                            map_estimates["bin_var"],
                            dtype=np.float64,
                        ).reshape(-1)[bin_idx]
                        corr = _safe_correlation(sex_chrom_bin_var, unexpected_by_bin)
                        if corr is not None:
                            messages.append(
                                f"Correlation(bin_var, {chrom_label} unexpected dominant CN sample fraction among {sex_label}-assigned samples by bin)="
                                f"{corr:.3f}"
                            )

                    plot_subset = observed_plot_depth[np.ix_(bin_idx, sample_idx)]
                    mean_plot_depth_by_bin = plot_subset.mean(axis=1)
                    messages.append(
                        format_numeric_summary(
                            f"{chrom_label} plot depth across {sex_label}-assigned bin-sample pairs",
                            plot_subset,
                            precision=4,
                        )
                    )
                    if background_bin_loadings is not None:
                        for factor_idx in reported_background_factors:
                            corr = _safe_correlation(
                                background_bin_loadings[bin_idx, factor_idx],
                                mean_plot_depth_by_bin,
                            )
                            if corr is not None:
                                messages.append(
                                    f"Correlation(background_factor_{factor_idx + 1}_loading, {chrom_label} mean plot depth among {sex_label}-assigned samples by bin)="
                                    f"{corr:.3f}"
                                )

                    if bin_bias_arr is not None:
                        if bin_epsilon_matrix is None:
                            epsilon_subset = np.zeros(
                                (bin_idx.size, sample_idx.size),
                                dtype=np.float64,
                            )
                        else:
                            epsilon_subset = bin_epsilon_matrix[np.ix_(bin_idx, sample_idx)]
                        expected_plot_depth = (
                            expected_cn * bin_bias_arr[bin_idx, np.newaxis] +
                            epsilon_subset
                        )
                        if depth_ratio_for_expected is not None:
                            expected_plot_depth = (
                                expected_plot_depth *
                                depth_ratio_for_expected[np.newaxis, sample_idx]
                            )
                        residual_subset = plot_subset - expected_plot_depth
                        messages.append(
                            format_numeric_summary(
                                f"{chrom_label} observed-minus-expected plot depth for expected {sex_label} copy number",
                                residual_subset,
                                precision=4,
                            )
                        )
                        mean_abs_residual_by_bin = np.abs(residual_subset).mean(axis=1)
                        messages.append(
                            format_numeric_summary(
                                f"{chrom_label} mean absolute observed-minus-expected plot depth among {sex_label}-assigned samples by bin",
                                mean_abs_residual_by_bin,
                                precision=4,
                            )
                        )
                        if "bin_var" in map_estimates:
                            corr = _safe_correlation(
                                sex_chrom_bin_var,
                                mean_abs_residual_by_bin,
                            )
                            if corr is not None:
                                messages.append(
                                    f"Correlation(bin_var, {chrom_label} mean absolute observed-minus-expected plot depth among {sex_label}-assigned samples by bin)="
                                    f"{corr:.3f}"
                                )
                        corr = _safe_correlation(
                            bin_bias_arr[bin_idx],
                            mean_abs_residual_by_bin,
                        )
                        if corr is not None:
                            messages.append(
                                f"Correlation(bin_bias, {chrom_label} mean absolute observed-minus-expected plot depth among {sex_label}-assigned samples by bin)="
                                f"{corr:.3f}"
                            )
                        if background_bin_loadings is not None:
                            for factor_idx in reported_background_factors:
                                corr = _safe_correlation(
                                    background_bin_loadings[bin_idx, factor_idx],
                                    mean_abs_residual_by_bin,
                                )
                                if corr is not None:
                                    messages.append(
                                        f"Correlation(background_factor_{factor_idx + 1}_loading, {chrom_label} mean absolute observed-minus-expected plot depth among {sex_label}-assigned samples by bin)="
                                        f"{corr:.3f}"
                                    )

    if np.any(autosome_mask):
        autosome_cn_probs = cn_probs[autosome_mask, :, :]
        autosome_cn_map = cn_map[autosome_mask, :]
        af_shift_sample_burden: np.ndarray | None = None
        af_shift_bin_fraction: np.ndarray | None = None
        messages.append(
            format_count_summary(
                "Autosomal dominant CN across bin-sample pairs",
                np.bincount(
                    autosome_cn_map.reshape(-1),
                    minlength=cn_probs.shape[2],
                ),
                cn_state_names,
            )
        )

        autosome_cn3_burden = (autosome_cn_map >= 3).mean(axis=0)
        autosome_cn4_burden = (autosome_cn_map >= 4).mean(axis=0)
        autosome_cn4_bin_fraction = (autosome_cn_map >= 4).mean(axis=1)
        messages.append(
            format_numeric_summary(
                "Autosomal dominant CN>=3 burden across samples",
                autosome_cn3_burden,
                precision=4,
            )
        )
        messages.append(
            format_numeric_summary(
                "Autosomal dominant CN>=4 burden across samples",
                autosome_cn4_burden,
                precision=4,
            )
        )
        messages.append(
            format_numeric_summary(
                "Autosomal dominant CN>=4 sample fraction across bins",
                autosome_cn4_bin_fraction,
                precision=4,
            )
        )
        messages.append(
            _format_threshold_exceedance_summary(
                "Autosomal bins with dominant CN>=4 sample fraction thresholds",
                autosome_cn4_bin_fraction,
                (0.10, 0.50, 0.90),
            )
        )

        if data.depth_space == "raw":
            empirical_plot_depth = _plot_normalized_depth(
                data.depth.detach().cpu().numpy(),
                data.bin_length_kb.detach().cpu().numpy()[:, np.newaxis],
                _plot_depth_baseline_per_sample(data)[np.newaxis, :],
            )
        else:
            empirical_plot_depth = data.depth.detach().cpu().numpy()
        autosome_mean_plot_depth = empirical_plot_depth[autosome_mask, :].mean(axis=1)
        messages.append(
            format_numeric_summary(
                "Autosomal mean plot depth across bins",
                autosome_mean_plot_depth,
                precision=4,
            )
        )
        corr = _safe_correlation(
            autosome_mean_plot_depth,
            autosome_cn4_bin_fraction,
        )
        if corr is not None:
            messages.append(
                "Correlation(empirical mean autosomal plot depth by bin, autosomal dominant CN>=4 sample fraction by bin)="
                f"{corr:.3f}"
            )
        recurrent_cn4_bins = autosome_cn4_bin_fraction >= 0.50
        if np.any(recurrent_cn4_bins):
            messages.append(
                format_numeric_summary(
                    "Autosomal mean plot depth for bins with dominant CN>=4 sample fraction >=50%",
                    autosome_mean_plot_depth[recurrent_cn4_bins],
                    precision=4,
                )
            )

        if "bin_var" in map_estimates:
            autosome_bin_var = np.asarray(
                map_estimates["bin_var"],
                dtype=np.float64,
            ).reshape(-1)[autosome_mask]
            corr = _safe_correlation(autosome_bin_var, autosome_cn4_bin_fraction)
            if corr is not None:
                messages.append(
                    "Correlation(bin_var, autosomal dominant CN>=4 sample fraction by bin)="
                    f"{corr:.3f}"
                )

        if bin_epsilon_matrix is not None:
            autosome_bin_epsilon = bin_epsilon_matrix[autosome_mask, :].mean(axis=1)
            corr = _safe_correlation(autosome_bin_epsilon, autosome_cn4_bin_fraction)
            if corr is not None:
                messages.append(
                    "Correlation(bin_epsilon, autosomal dominant CN>=4 sample fraction by bin)="
                    f"{corr:.3f}"
                )

        if "cn_probs" in map_estimates and cn_probs.shape[2] > 4:
            cn_prob_matrix = _coerce_bin_state_matrix(
                map_estimates["cn_probs"],
                data.n_bins,
            )
            if cn_prob_matrix is not None and cn_prob_matrix.shape[1] > 4:
                autosome_high_cn_prior = cn_prob_matrix[autosome_mask, 4:].sum(axis=1)
                corr = _safe_correlation(
                    autosome_high_cn_prior,
                    autosome_cn4_bin_fraction,
                )
                if corr is not None:
                    messages.append(
                        "Correlation(per-bin CN>=4 prior mass, autosomal dominant CN>=4 sample fraction by bin)="
                        f"{corr:.3f}"
                    )

        if cn_probs.shape[2] > 4:
            autosome_cn2_posterior = autosome_cn_probs[:, :, 2]
            autosome_cn4_posterior = autosome_cn_probs[:, :, 4]
            autosome_cn24_posterior = autosome_cn2_posterior + autosome_cn4_posterior
            messages.append(
                format_numeric_summary(
                    "Autosomal CN2+CN4 posterior mass across bin-sample pairs",
                    autosome_cn24_posterior,
                    precision=4,
                )
            )
            messages.append(
                _format_threshold_exceedance_summary(
                    "Autosomal bin-sample pairs with CN2+CN4 posterior mass thresholds",
                    autosome_cn24_posterior,
                    (0.25, 0.50, 0.80),
                )
            )

            autosome_top_two_states = np.argsort(autosome_cn_probs, axis=2)[..., -2:]
            cn24_top_two_mask = (
                ((autosome_top_two_states[..., 0] == 2) & (autosome_top_two_states[..., 1] == 4)) |
                ((autosome_top_two_states[..., 0] == 4) & (autosome_top_two_states[..., 1] == 2))
            )
            messages.append(
                _format_boolean_fraction_summary(
                    "Autosomal top-two posterior pair is CN2/CN4",
                    cn24_top_two_mask,
                )
            )
            if np.any(cn24_top_two_mask):
                cn24_top_two_mass = autosome_cn24_posterior[cn24_top_two_mask]
                cn4_share_within_cn24 = np.divide(
                    autosome_cn4_posterior[cn24_top_two_mask],
                    cn24_top_two_mass,
                    out=np.full(cn24_top_two_mass.shape, np.nan, dtype=np.float64),
                    where=cn24_top_two_mass > 0,
                )
                messages.append(
                    format_numeric_summary(
                        "Autosomal combined CN2+CN4 posterior mass for CN2/CN4 top-two pairs",
                        cn24_top_two_mass,
                        precision=4,
                    )
                )
                messages.append(
                    format_numeric_summary(
                        "Autosomal CN4 share within CN2/CN4 posterior mass for CN2/CN4 top-two pairs",
                        cn4_share_within_cn24,
                        precision=4,
                    )
                )

            cn4_map_mask = autosome_cn_map == 4
            if np.any(cn4_map_mask):
                messages.append(
                    format_numeric_summary(
                        "Autosomal CN2 posterior on CN4-MAP bin-sample pairs",
                        autosome_cn2_posterior[cn4_map_mask],
                        precision=4,
                    )
                )
            cn2_map_mask = autosome_cn_map == 2
            if np.any(cn2_map_mask):
                messages.append(
                    format_numeric_summary(
                        "Autosomal CN4 posterior on CN2-MAP bin-sample pairs",
                        autosome_cn4_posterior[cn2_map_mask],
                        precision=4,
                    )
                )

            if "bin_bias" in map_estimates:
                autosome_bin_bias = np.asarray(
                    map_estimates["bin_bias"],
                    dtype=np.float64,
                ).reshape(-1)[autosome_mask]
                corr = _safe_correlation(
                    autosome_bin_bias,
                    autosome_cn4_posterior.mean(axis=1),
                )
                if corr is not None:
                    messages.append(
                        "Correlation(bin_bias, mean autosomal CN4 posterior mass by bin)="
                        f"{corr:.3f}"
                    )
                corr = _safe_correlation(
                    autosome_bin_bias,
                    autosome_cn24_posterior.mean(axis=1),
                )
                if corr is not None:
                    messages.append(
                        "Correlation(bin_bias, mean autosomal CN2+CN4 posterior mass by bin)="
                        f"{corr:.3f}"
                    )

        messages.append(
            _format_threshold_exceedance_summary(
                "Samples exceeding autosomal dominant CN>=3 burden thresholds",
                autosome_cn3_burden,
                (0.01, 0.05, 0.10),
            )
        )

        autosomal_chromosomes: list[str] = []
        autosomal_cn3_rates: list[float] = []
        chr_values = np.asarray(data.chr)
        for chromosome in np.unique(chr_values[autosome_mask]):
            chrom_mask = chr_values == chromosome
            autosomal_chromosomes.append(str(chromosome))
            autosomal_cn3_rates.append(float((cn_map[chrom_mask, :] >= 3).mean()))
        messages.append(
            _format_chromosome_rate_summary(
                "Autosomal dominant CN>=3 rate by chromosome",
                np.asarray(autosomal_chromosomes, dtype=object),
                np.asarray(autosomal_cn3_rates, dtype=np.float64),
            )
        )

        if "cnq" in cn_posterior:
            messages.append(
                format_numeric_summary(
                    "Autosomal CNQ across bin-sample pairs",
                    np.asarray(cn_posterior["cnq"], dtype=np.float64)[autosome_mask, :],
                )
            )
        if "cn_map_stability" in cn_posterior:
            stability = np.asarray(cn_posterior["cn_map_stability"], dtype=np.float64)
            messages.append(
                format_numeric_summary(
                    "Mean autosomal CN MAP stability across samples",
                    stability[autosome_mask, :].mean(axis=0),
                    precision=4,
                )
            )
        if "sample_var" in map_estimates:
            corr = _safe_correlation(map_estimates["sample_var"], autosome_cn3_burden)
            if corr is not None:
                messages.append(
                    f"Correlation(sample_var, autosomal dominant CN>=3 burden)={corr:.3f}"
                )
            corr = _safe_correlation(map_estimates["sample_var"], autosome_cn4_burden)
            if corr is not None:
                messages.append(
                    f"Correlation(sample_var, autosomal dominant CN>=4 burden)={corr:.3f}"
                )
        if depth_ratio is not None:
            corr = _safe_correlation(depth_ratio, autosome_cn3_burden)
            if corr is not None:
                messages.append(
                    "Correlation(sample_depth / empirical autosomal baseline, "
                    f"autosomal dominant CN>=3 burden)={corr:.3f}"
                )
            corr = _safe_correlation(depth_ratio, autosome_cn4_burden)
            if corr is not None:
                messages.append(
                    "Correlation(sample_depth / empirical autosomal baseline, "
                    f"autosomal dominant CN>=4 burden)={corr:.3f}"
                )

        if "af_temperature" in map_estimates:
            af_temperature = float(np.asarray(map_estimates["af_temperature"]).item())
            if np.isfinite(af_temperature):
                messages.append(
                    f"Learned AF temperature scalar={af_temperature:.6f}"
                )

        if (
            af_table is not None and
            af_table.ndim == 3 and
            af_table.shape[0] >= 6 and
            af_table.shape[1] == data.n_bins and
            af_table.shape[2] == data.n_samples
        ):
            autosome_af_table = np.asarray(af_table[:, autosome_mask, :], dtype=np.float64)
            cn4_vs_cn2_margin = autosome_af_table[4] - autosome_af_table[2]
            cn5_vs_cn2_margin = autosome_af_table[5] - autosome_af_table[2]

            messages.append(
                format_numeric_summary(
                    "Autosomal AF CN4-CN2 log-likelihood margin across bin-sample pairs",
                    cn4_vs_cn2_margin,
                    precision=4,
                )
            )
            messages.append(
                _format_numeric_threshold_summary(
                    "Autosomal bin-sample pairs with AF CN4-CN2 margin thresholds",
                    cn4_vs_cn2_margin,
                    (0.0, 1.0, 2.0),
                )
            )
            messages.append(
                format_numeric_summary(
                    "Autosomal AF CN5-CN2 log-likelihood margin across bin-sample pairs",
                    cn5_vs_cn2_margin,
                    precision=4,
                )
            )
            messages.append(
                _format_numeric_threshold_summary(
                    "Autosomal bin-sample pairs with AF CN5-CN2 margin thresholds",
                    cn5_vs_cn2_margin,
                    (0.0, 1.0, 2.0),
                )
            )

            cn4_prefer_sample_burden = (cn4_vs_cn2_margin > 0).mean(axis=0)
            cn5_prefer_sample_burden = (cn5_vs_cn2_margin > 0).mean(axis=0)
            cn4_prefer_bin_fraction = (cn4_vs_cn2_margin > 0).mean(axis=1)
            cn5_prefer_bin_fraction = (cn5_vs_cn2_margin > 0).mean(axis=1)

            messages.append(
                format_numeric_summary(
                    "Autosomal AF preference for CN4 over CN2 across samples",
                    cn4_prefer_sample_burden,
                    precision=4,
                )
            )
            messages.append(
                format_numeric_summary(
                    "Autosomal AF preference for CN5 over CN2 across samples",
                    cn5_prefer_sample_burden,
                    precision=4,
                )
            )
            messages.append(
                format_numeric_summary(
                    "Autosomal AF preference for CN4 over CN2 across bins",
                    cn4_prefer_bin_fraction,
                    precision=4,
                )
            )
            messages.append(
                format_numeric_summary(
                    "Autosomal AF preference for CN5 over CN2 across bins",
                    cn5_prefer_bin_fraction,
                    precision=4,
                )
            )

            if "sample_var" in map_estimates:
                corr = _safe_correlation(map_estimates["sample_var"], cn4_prefer_sample_burden)
                if corr is not None:
                    messages.append(
                        f"Correlation(sample_var, autosomal AF preference for CN4 over CN2)={corr:.3f}"
                    )
                corr = _safe_correlation(map_estimates["sample_var"], cn5_prefer_sample_burden)
                if corr is not None:
                    messages.append(
                        f"Correlation(sample_var, autosomal AF preference for CN5 over CN2)={corr:.3f}"
                    )
            if depth_ratio is not None:
                corr = _safe_correlation(depth_ratio, cn4_prefer_sample_burden)
                if corr is not None:
                    messages.append(
                        "Correlation(sample_depth / empirical autosomal baseline, "
                        f"autosomal AF preference for CN4 over CN2)={corr:.3f}"
                    )
                corr = _safe_correlation(depth_ratio, cn5_prefer_sample_burden)
                if corr is not None:
                    messages.append(
                        "Correlation(sample_depth / empirical autosomal baseline, "
                        f"autosomal AF preference for CN5 over CN2)={corr:.3f}"
                    )
            if "bin_var" in map_estimates:
                autosome_bin_var = np.asarray(
                    map_estimates["bin_var"],
                    dtype=np.float64,
                ).reshape(-1)[autosome_mask]
                corr = _safe_correlation(autosome_bin_var, cn4_prefer_bin_fraction)
                if corr is not None:
                    messages.append(
                        "Correlation(bin_var, autosomal AF preference for CN4 over CN2 by bin)="
                        f"{corr:.3f}"
                    )
                corr = _safe_correlation(autosome_bin_var, cn5_prefer_bin_fraction)
                if corr is not None:
                    messages.append(
                        "Correlation(bin_var, autosomal AF preference for CN5 over CN2 by bin)="
                        f"{corr:.3f}"
                    )

        if depth_only_cn_posterior is not None:
            depth_only_probs = np.asarray(
                depth_only_cn_posterior["cn_posterior"],
                dtype=np.float64,
            )
            if depth_only_probs.shape == cn_probs.shape:
                depth_only_cn_map = np.argmax(depth_only_probs, axis=2)[autosome_mask, :]
                af_shift_mask = autosome_cn_map != depth_only_cn_map
                af_shift_sample_burden = af_shift_mask.mean(axis=0)
                af_shift_bin_fraction = af_shift_mask.mean(axis=1)
                messages.append(
                    _format_transition_count_summary(
                        "Autosomal depth-only to AF-enabled CN MAP transition counts across bin-sample pairs",
                        depth_only_cn_map,
                        autosome_cn_map,
                        cn_state_names,
                    )
                )
                messages.append(
                    format_numeric_summary(
                        "Autosomal AF-induced CN MAP shift burden across samples",
                        af_shift_sample_burden,
                        precision=4,
                    )
                )
                messages.append(
                    _format_threshold_exceedance_summary(
                        "Samples exceeding autosomal AF-induced CN MAP shift burden thresholds",
                        af_shift_sample_burden,
                        (0.01, 0.05, 0.10, 0.20),
                    )
                )
                messages.append(
                    _format_concentration_summary(
                        "Concentration of autosomal AF-induced CN MAP shift burden across samples",
                        af_shift_sample_burden,
                    )
                )
                messages.append(
                    format_numeric_summary(
                        "Autosomal AF-induced CN MAP shift sample fraction across bins",
                        af_shift_bin_fraction,
                        precision=4,
                    )
                )
                messages.append(
                    _format_threshold_exceedance_summary(
                        "Autosomal bins with AF-induced CN MAP shift sample fraction thresholds",
                        af_shift_bin_fraction,
                        (0.10, 0.50, 0.90),
                    )
                )

                depth_only_conf = depth_only_probs[autosome_mask, :, :].max(axis=2)
                af_enabled_conf = autosome_cn_probs.max(axis=2)
                messages.append(
                    format_numeric_summary(
                        "Autosomal AF-enabled minus depth-only max posterior confidence across bin-sample pairs",
                        af_enabled_conf - depth_only_conf,
                        precision=4,
                    )
                )

                if "sample_var" in map_estimates:
                    corr = _safe_correlation(map_estimates["sample_var"], af_shift_sample_burden)
                    if corr is not None:
                        messages.append(
                            f"Correlation(sample_var, autosomal AF-induced CN MAP shift burden)={corr:.3f}"
                        )
                if depth_ratio is not None:
                    corr = _safe_correlation(depth_ratio, af_shift_sample_burden)
                    if corr is not None:
                        messages.append(
                            "Correlation(sample_depth / empirical autosomal baseline, "
                            f"autosomal AF-induced CN MAP shift burden)={corr:.3f}"
                        )
                if "bin_var" in map_estimates:
                    autosome_bin_var = np.asarray(
                        map_estimates["bin_var"],
                        dtype=np.float64,
                    ).reshape(-1)[autosome_mask]
                    corr = _safe_correlation(autosome_bin_var, af_shift_bin_fraction)
                    if corr is not None:
                        messages.append(
                            "Correlation(bin_var, autosomal AF-induced CN MAP shift sample fraction by bin)="
                            f"{corr:.3f}"
                        )
                if "bin_bias" in map_estimates:
                    autosome_bin_bias = np.asarray(
                        map_estimates["bin_bias"],
                        dtype=np.float64,
                    ).reshape(-1)[autosome_mask]
                    corr = _safe_correlation(autosome_bin_bias, af_shift_bin_fraction)
                    if corr is not None:
                        messages.append(
                            "Correlation(bin_bias, autosomal AF-induced CN MAP shift sample fraction by bin)="
                            f"{corr:.3f}"
                        )

                if data.site_mask is not None:
                    observed_site_count = data.site_mask.detach().cpu().numpy()[autosome_mask, :, :].sum(axis=1)
                    mean_observed_sites_per_sample = observed_site_count.mean(axis=0)
                    mean_observed_sites_per_bin = observed_site_count.mean(axis=1)
                    messages.append(
                        format_numeric_summary(
                            "Mean autosomal observed AF site count per bin across samples",
                            mean_observed_sites_per_sample,
                            precision=2,
                        )
                    )
                    corr = _safe_correlation(mean_observed_sites_per_sample, af_shift_sample_burden)
                    if corr is not None:
                        messages.append(
                            "Correlation(mean autosomal observed AF site count per bin across samples, "
                            f"autosomal AF-induced CN MAP shift burden)={corr:.3f}"
                        )
                    corr = _safe_correlation(mean_observed_sites_per_bin, af_shift_bin_fraction)
                    if corr is not None:
                        messages.append(
                            "Correlation(mean autosomal observed AF site count across samples by bin, "
                            f"autosomal AF-induced CN MAP shift sample fraction by bin)={corr:.3f}"
                        )

                autosomal_chromosomes_shift: list[str] = []
                autosomal_shift_rates: list[float] = []
                chr_values = np.asarray(data.chr)
                autosome_chr_values = chr_values[autosome_mask]
                for chromosome in np.unique(autosome_chr_values):
                    chromosome_mask = autosome_chr_values == chromosome
                    autosomal_chromosomes_shift.append(str(chromosome))
                    autosomal_shift_rates.append(float(af_shift_mask[chromosome_mask, :].mean()))
                messages.append(
                    _format_chromosome_rate_summary(
                        "Autosomal AF-induced CN MAP shift rate by chromosome",
                        np.asarray(autosomal_chromosomes_shift, dtype=object),
                        np.asarray(autosomal_shift_rates, dtype=np.float64),
                    )
                )

        if all((data.site_alt is not None, data.site_total is not None, data.site_mask is not None)):
            site_alt = np.asarray(
                data.site_alt.detach().cpu().numpy()[autosome_mask, :, :],
                dtype=np.float64,
            )
            site_total = np.asarray(
                data.site_total.detach().cpu().numpy()[autosome_mask, :, :],
                dtype=np.float64,
            )
            site_mask = np.asarray(
                data.site_mask.detach().cpu().numpy()[autosome_mask, :, :],
                dtype=bool,
            ) & (site_total > 0)
            if np.any(site_mask):
                if data.site_pop_af is not None:
                    input_site_pop_af = np.asarray(
                        data.site_pop_af.detach().cpu().numpy()[autosome_mask, ...],
                        dtype=np.float64,
                    )
                    effective_site_pop_af, used_leave_one_out = _resolve_fixed_site_pop_af_numpy(
                        site_alt,
                        site_total,
                        input_site_pop_af,
                        site_mask,
                    )
                    if effective_site_pop_af.ndim == 2:
                        effective_site_pop_af = np.broadcast_to(
                            effective_site_pop_af[:, :, np.newaxis],
                            site_total.shape,
                        )
                    effective_site_pop_af = np.clip(
                        np.asarray(effective_site_pop_af, dtype=np.float64),
                        1e-6,
                        1.0 - 1e-6,
                    )
                    diploid_het_prior_mass = 2.0 * effective_site_pop_af * (1.0 - effective_site_pop_af)
                    observed_site_sample_count = int(site_mask.sum())
                    observed_site_count_per_bin = site_mask.sum(axis=(1, 2))
                    mean_effective_site_pop_af_by_bin = np.divide(
                        (effective_site_pop_af * site_mask).sum(axis=(1, 2)),
                        observed_site_count_per_bin,
                        out=np.full(observed_site_count_per_bin.shape, np.nan, dtype=np.float64),
                        where=observed_site_count_per_bin > 0,
                    )
                    mean_diploid_het_prior_mass_by_bin = np.divide(
                        (diploid_het_prior_mass * site_mask).sum(axis=(1, 2)),
                        observed_site_count_per_bin,
                        out=np.full(observed_site_count_per_bin.shape, np.nan, dtype=np.float64),
                        where=observed_site_count_per_bin > 0,
                    )

                    messages.append(
                        f"Autosomal effective site_pop_af uses leave-one-out self-pooled resolution={used_leave_one_out}"
                    )
                    messages.append(
                        format_numeric_summary(
                            "Autosomal effective site_pop_af across observed site-sample values",
                            effective_site_pop_af[site_mask],
                            precision=4,
                        )
                    )
                    messages.append(
                        format_numeric_summary(
                            "Autosomal diploid heterozygous prior mass 2p(1-p) across observed site-sample values",
                            diploid_het_prior_mass[site_mask],
                            precision=4,
                        )
                    )
                    messages.append(
                        _format_numeric_below_threshold_summary(
                            "Autosomal observed site-sample pairs with diploid heterozygous prior mass thresholds",
                            diploid_het_prior_mass[site_mask],
                            (0.05, 0.10, 0.20),
                        )
                    )
                    messages.append(
                        format_numeric_summary(
                            "Mean autosomal diploid heterozygous prior mass across bins",
                            mean_diploid_het_prior_mass_by_bin,
                            precision=4,
                        )
                    )

                    if (
                        af_table is not None and
                        af_table.ndim == 3 and
                        af_table.shape[0] >= 6 and
                        af_table.shape[1] == data.n_bins and
                        af_table.shape[2] == data.n_samples and
                        observed_site_sample_count > 0
                    ):
                        autosome_af_table = np.asarray(af_table[:, autosome_mask, :], dtype=np.float64)
                        mean_cn4_vs_cn2_margin_by_bin = (autosome_af_table[4] - autosome_af_table[2]).mean(axis=1)
                        mean_cn5_vs_cn2_margin_by_bin = (autosome_af_table[5] - autosome_af_table[2]).mean(axis=1)

                        corr = _safe_correlation(
                            mean_effective_site_pop_af_by_bin,
                            mean_cn4_vs_cn2_margin_by_bin,
                        )
                        if corr is not None:
                            messages.append(
                                "Correlation(mean autosomal effective site_pop_af by bin, mean autosomal AF CN4-CN2 margin by bin)="
                                f"{corr:.3f}"
                            )
                        corr = _safe_correlation(
                            mean_effective_site_pop_af_by_bin,
                            mean_cn5_vs_cn2_margin_by_bin,
                        )
                        if corr is not None:
                            messages.append(
                                "Correlation(mean autosomal effective site_pop_af by bin, mean autosomal AF CN5-CN2 margin by bin)="
                                f"{corr:.3f}"
                            )
                        corr = _safe_correlation(
                            mean_diploid_het_prior_mass_by_bin,
                            mean_cn4_vs_cn2_margin_by_bin,
                        )
                        if corr is not None:
                            messages.append(
                                "Correlation(mean autosomal diploid heterozygous prior mass by bin, mean autosomal AF CN4-CN2 margin by bin)="
                                f"{corr:.3f}"
                            )
                        corr = _safe_correlation(
                            mean_diploid_het_prior_mass_by_bin,
                            mean_cn5_vs_cn2_margin_by_bin,
                        )
                        if corr is not None:
                            messages.append(
                                "Correlation(mean autosomal diploid heterozygous prior mass by bin, mean autosomal AF CN5-CN2 margin by bin)="
                                f"{corr:.3f}"
                            )

                    site_level_cn2 = _site_level_marginalized_af_log_lik_numpy(
                        site_alt,
                        site_total,
                        effective_site_pop_af,
                        site_mask,
                        cn_state=2,
                        n_states=cn_probs.shape[2],
                        concentration=af_concentration,
                    )
                    site_level_cn4 = _site_level_marginalized_af_log_lik_numpy(
                        site_alt,
                        site_total,
                        effective_site_pop_af,
                        site_mask,
                        cn_state=4,
                        n_states=cn_probs.shape[2],
                        concentration=af_concentration,
                    )
                    site_level_cn5 = _site_level_marginalized_af_log_lik_numpy(
                        site_alt,
                        site_total,
                        effective_site_pop_af,
                        site_mask,
                        cn_state=5,
                        n_states=cn_probs.shape[2],
                        concentration=af_concentration,
                    )
                    site_level_cn4_vs_cn2_margin = site_level_cn4 - site_level_cn2
                    site_level_cn5_vs_cn2_margin = site_level_cn5 - site_level_cn2
                    observed_site_total = site_total[site_mask]
                    cn4_margin_observed = site_level_cn4_vs_cn2_margin[site_mask]
                    cn5_margin_observed = site_level_cn5_vs_cn2_margin[site_mask]
                    total_site_weight = float(observed_site_total.sum())

                    messages.append(
                        format_numeric_summary(
                            "Autosomal site-level AF CN4-CN2 margin across observed site-sample pairs",
                            cn4_margin_observed,
                            precision=4,
                        )
                    )
                    messages.append(
                        _format_numeric_threshold_summary(
                            "Autosomal observed site-sample pairs with site-level AF CN4-CN2 margin thresholds",
                            cn4_margin_observed,
                            (0.0, 0.5, 1.0),
                        )
                    )
                    messages.append(
                        format_numeric_summary(
                            "Autosomal site-level AF CN5-CN2 margin across observed site-sample pairs",
                            cn5_margin_observed,
                            precision=4,
                        )
                    )
                    messages.append(
                        _format_numeric_threshold_summary(
                            "Autosomal observed site-sample pairs with site-level AF CN5-CN2 margin thresholds",
                            cn5_margin_observed,
                            (0.0, 0.5, 1.0),
                        )
                    )
                    messages.append(
                        format_numeric_summary(
                            "Autosomal observed site_total across observed site-sample pairs",
                            observed_site_total,
                            precision=2,
                        )
                    )
                    if total_site_weight > 0:
                        messages.append(
                            f"Site-total-weighted autosomal site-level AF CN4-CN2 margin={np.average(cn4_margin_observed, weights=observed_site_total):.4f}"
                        )
                        messages.append(
                            f"Site-total-weighted autosomal site-level AF CN5-CN2 margin={np.average(cn5_margin_observed, weights=observed_site_total):.4f}"
                        )
                    corr = _safe_correlation(observed_site_total, cn4_margin_observed)
                    if corr is not None:
                        messages.append(
                            f"Correlation(observed site_total, autosomal site-level AF CN4-CN2 margin)={corr:.3f}"
                        )
                    corr = _safe_correlation(observed_site_total, cn5_margin_observed)
                    if corr is not None:
                        messages.append(
                            f"Correlation(observed site_total, autosomal site-level AF CN5-CN2 margin)={corr:.3f}"
                        )
                    strong_cn4_mask = cn4_margin_observed >= 0.5
                    if np.any(strong_cn4_mask):
                        messages.append(
                            format_numeric_summary(
                                "Autosomal observed site_total on site-sample pairs with site-level AF CN4-CN2 margin >=0.5",
                                observed_site_total[strong_cn4_mask],
                                precision=2,
                            )
                        )
                    strong_cn5_mask = cn5_margin_observed >= 0.5
                    if np.any(strong_cn5_mask):
                        messages.append(
                            format_numeric_summary(
                                "Autosomal observed site_total on site-sample pairs with site-level AF CN5-CN2 margin >=0.5",
                                observed_site_total[strong_cn5_mask],
                                precision=2,
                            )
                        )

                observed_af = np.divide(
                    site_alt,
                    site_total,
                    out=np.full(site_alt.shape, np.nan, dtype=np.float64),
                    where=site_mask,
                )
                diploid_grid_distance = _nearest_af_grid_distance(
                    observed_af,
                    (0.0, 0.5, 1.0),
                )
                cn4_grid_distance = _nearest_af_grid_distance(
                    observed_af,
                    (0.0, 0.25, 0.5, 0.75, 1.0),
                )
                cn5_grid_distance = _nearest_af_grid_distance(
                    observed_af,
                    (0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
                )
                cn4_grid_advantage = diploid_grid_distance - cn4_grid_distance
                cn5_grid_advantage = diploid_grid_distance - cn5_grid_distance

                messages.append(
                    format_numeric_summary(
                        "Autosomal raw AF CN4-vs-diploid grid advantage across observed sites",
                        cn4_grid_advantage[site_mask],
                        precision=4,
                    )
                )
                messages.append(
                    _format_numeric_threshold_summary(
                        "Autosomal observed sites with raw AF CN4-vs-diploid grid advantage thresholds",
                        cn4_grid_advantage[site_mask],
                        (0.0, 0.05, 0.10),
                    )
                )
                messages.append(
                    format_numeric_summary(
                        "Autosomal raw AF CN5-vs-diploid grid advantage across observed sites",
                        cn5_grid_advantage[site_mask],
                        precision=4,
                    )
                )
                messages.append(
                    _format_numeric_threshold_summary(
                        "Autosomal observed sites with raw AF CN5-vs-diploid grid advantage thresholds",
                        cn5_grid_advantage[site_mask],
                        (0.0, 0.05, 0.10),
                    )
                )

                observed_site_count_per_sample = site_mask.sum(axis=(0, 1))
                observed_site_count_per_bin = site_mask.sum(axis=(1, 2))
                cn4_grid_preference_sample = np.divide(
                    (site_mask & (cn4_grid_advantage > 0)).sum(axis=(0, 1)),
                    observed_site_count_per_sample,
                    out=np.full(data.n_samples, np.nan, dtype=np.float64),
                    where=observed_site_count_per_sample > 0,
                )
                cn5_grid_preference_sample = np.divide(
                    (site_mask & (cn5_grid_advantage > 0)).sum(axis=(0, 1)),
                    observed_site_count_per_sample,
                    out=np.full(data.n_samples, np.nan, dtype=np.float64),
                    where=observed_site_count_per_sample > 0,
                )
                cn4_grid_preference_bin = np.divide(
                    (site_mask & (cn4_grid_advantage > 0)).sum(axis=(1, 2)),
                    observed_site_count_per_bin,
                    out=np.full(observed_site_count_per_bin.shape, np.nan, dtype=np.float64),
                    where=observed_site_count_per_bin > 0,
                )
                cn5_grid_preference_bin = np.divide(
                    (site_mask & (cn5_grid_advantage > 0)).sum(axis=(1, 2)),
                    observed_site_count_per_bin,
                    out=np.full(observed_site_count_per_bin.shape, np.nan, dtype=np.float64),
                    where=observed_site_count_per_bin > 0,
                )

                messages.append(
                    format_numeric_summary(
                        "Autosomal raw AF preference for CN4 genotype grid over diploid grid across samples",
                        cn4_grid_preference_sample,
                        precision=4,
                    )
                )
                messages.append(
                    format_numeric_summary(
                        "Autosomal raw AF preference for CN5 genotype grid over diploid grid across samples",
                        cn5_grid_preference_sample,
                        precision=4,
                    )
                )
                messages.append(
                    format_numeric_summary(
                        "Autosomal raw AF preference for CN4 genotype grid over diploid grid across bins",
                        cn4_grid_preference_bin,
                        precision=4,
                    )
                )
                messages.append(
                    format_numeric_summary(
                        "Autosomal raw AF preference for CN5 genotype grid over diploid grid across bins",
                        cn5_grid_preference_bin,
                        precision=4,
                    )
                )

                if af_shift_sample_burden is not None:
                    corr = _safe_correlation(cn4_grid_preference_sample, af_shift_sample_burden)
                    if corr is not None:
                        messages.append(
                            "Correlation(autosomal raw AF preference for CN4 genotype grid over diploid grid, "
                            f"autosomal AF-induced CN MAP shift burden)={corr:.3f}"
                        )
                    corr = _safe_correlation(cn5_grid_preference_sample, af_shift_sample_burden)
                    if corr is not None:
                        messages.append(
                            "Correlation(autosomal raw AF preference for CN5 genotype grid over diploid grid, "
                            f"autosomal AF-induced CN MAP shift burden)={corr:.3f}"
                        )
                if af_shift_bin_fraction is not None:
                    corr = _safe_correlation(cn4_grid_preference_bin, af_shift_bin_fraction)
                    if corr is not None:
                        messages.append(
                            "Correlation(autosomal raw AF preference for CN4 genotype grid over diploid grid by bin, "
                            f"autosomal AF-induced CN MAP shift sample fraction by bin)={corr:.3f}"
                        )
                    corr = _safe_correlation(cn5_grid_preference_bin, af_shift_bin_fraction)
                    if corr is not None:
                        messages.append(
                            "Correlation(autosomal raw AF preference for CN5 genotype grid over diploid grid by bin, "
                            f"autosomal AF-induced CN MAP shift sample fraction by bin)={corr:.3f}"
                        )

    if "sex_posterior" in cn_posterior:
        sex_post = np.asarray(cn_posterior["sex_posterior"], dtype=np.float64)
        sex_map = np.argmax(sex_post, axis=1)
        messages.append(
            format_count_summary(
                "Sex karyotype MAP across samples",
                np.bincount(sex_map, minlength=2),
                ("XX", "XY"),
            )
        )
        messages.append(
            format_numeric_summary(
                "Sex karyotype posterior confidence across samples",
                sex_post[np.arange(sex_post.shape[0]), sex_map],
                precision=4,
            )
        )

    return messages


def _plot_normalized_depth(
    depth: np.ndarray | float,
    bin_length_kb: np.ndarray | float,
    diploid_depth_per_kb: np.ndarray | float,
) -> np.ndarray:
    """Convert raw counts to diploid-normalized depth for plotting.

    Plotting uses the observed per-sample diploid baseline counts per kilobase,
    estimated as the median autosomal depth-per-kb, so a normalized depth of 2
    corresponds to a neutral diploid bin.
    """
    depth_arr = np.asarray(depth, dtype=np.float64)
    denom = np.asarray(bin_length_kb, dtype=np.float64) * np.asarray(
        diploid_depth_per_kb,
        dtype=np.float64,
    )
    return np.divide(
        2.0 * depth_arr,
        denom,
        out=np.zeros_like(depth_arr, dtype=np.float64),
        where=denom > 0,
    )


def _plot_depth_baseline_per_sample(data: DepthData) -> np.ndarray:
    """Estimate per-sample diploid baseline counts per kb for plotting."""
    if data.depth_space != "raw":
        return np.ones(data.n_samples, dtype=np.float64)

    autosome_mask = np.isin(data.chr, [f"chr{i}" for i in range(1, 23)])
    if not np.any(autosome_mask):
        autosome_mask = np.ones(data.n_bins, dtype=bool)

    depth_np = data.depth.detach().cpu().numpy()[autosome_mask, :]
    bin_length_kb = data.bin_length_kb.detach().cpu().numpy()[autosome_mask]
    depth_per_kb = np.divide(
        depth_np,
        bin_length_kb[:, np.newaxis],
        out=np.zeros_like(depth_np, dtype=np.float64),
        where=bin_length_kb[:, np.newaxis] > 0,
    )
    baseline = np.median(depth_per_kb, axis=0)
    return np.maximum(baseline, 1e-8)


def _clamp_threshold_for_depth_space(depth_space: str) -> float | None:
    """Use clamping only for the normalized-depth observation families."""
    return 5.0 if depth_space == "normalized" else None


def _resolve_observation_model(
    requested_depth_space: str,
    requested_obs_likelihood: str,
    input_path: str,
) -> Tuple[str, str]:
    """Resolve depth space and observation family using preprocess metadata."""
    marker_depth_space = read_observation_type(input_path)
    effective_depth_space = (
        marker_depth_space
        if marker_depth_space is not None and requested_depth_space == "auto"
        else requested_depth_space
    )

    if requested_obs_likelihood == "auto":
        if marker_depth_space is not None:
            resolved_obs_likelihood = default_obs_likelihood_for_depth_space(
                marker_depth_space,
            )
        elif effective_depth_space != "auto":
            resolved_obs_likelihood = default_obs_likelihood_for_depth_space(
                effective_depth_space,
            )
        else:
            resolved_obs_likelihood = default_obs_likelihood_for_depth_space(
                "normalized",
            )
    else:
        resolved_obs_likelihood = requested_obs_likelihood

    resolved_depth_space = validate_depth_space(
        effective_depth_space,
        resolved_obs_likelihood,
    )
    return resolved_depth_space, resolved_obs_likelihood


def _inference_tensor_dtype() -> torch.dtype:
    """Use float64 for all inference runs.

    The extra branching around observation family and AF evidence is not worth
    the complexity. Infer always runs in float64 for a single, stable numeric
    policy.
    """
    return torch.float64


def estimate_site_pop_af_naive_bayes(
    site_alt: np.ndarray,
    site_total: np.ndarray,
    site_mask: np.ndarray,
    prior_alpha: float = 1.0,
    prior_beta: float = 1.0,
) -> Dict[str, np.ndarray]:
    """Estimate per-site population AFs from pooled cohort read counts.

    This treats the observed alt reads across samples as conditionally
    independent Bernoulli draws from the site-level alt-allele frequency and
    uses a fixed Beta prior. The posterior mean is then used as the naive-
    Bayes population AF estimate for each site slot.
    """
    if prior_alpha <= 0.0 or prior_beta <= 0.0:
        raise ValueError("site AF prior alpha and beta must be positive.")

    observed_sites = np.any(site_mask, axis=2)
    valid_counts = np.asarray(site_mask, dtype=bool) & (site_total > 0)

    sum_alt = (site_alt * valid_counts).sum(axis=2, dtype=np.float64)
    sum_total = (site_total * valid_counts).sum(axis=2, dtype=np.float64)
    n_observed_samples = np.asarray(site_mask, dtype=np.int32).sum(axis=2)

    posterior_mean = np.zeros(sum_total.shape, dtype=np.float64)
    posterior_sd = np.zeros(sum_total.shape, dtype=np.float64)
    pooled_observed_af = np.full(sum_total.shape, np.nan, dtype=np.float64)

    if np.any(observed_sites):
        posterior_alpha = prior_alpha + sum_alt
        posterior_beta = prior_beta + np.maximum(sum_total - sum_alt, 0.0)
        posterior_denom = posterior_alpha + posterior_beta

        posterior_mean[observed_sites] = (
            posterior_alpha[observed_sites] / posterior_denom[observed_sites]
        )
        posterior_denom_observed = posterior_denom[observed_sites]
        posterior_var_denom = (
            posterior_denom_observed ** 2
        ) * (posterior_denom_observed + 1.0)
        posterior_var = (
            posterior_alpha[observed_sites] * posterior_beta[observed_sites]
        ) / posterior_var_denom
        posterior_sd[observed_sites] = np.sqrt(np.maximum(posterior_var, 0.0))
        np.divide(
            sum_alt,
            sum_total,
            out=pooled_observed_af,
            where=sum_total > 0,
        )

    return {
        "posterior_mean": posterior_mean.astype(np.float32),
        "posterior_sd": posterior_sd.astype(np.float32),
        "pooled_observed_af": pooled_observed_af.astype(np.float32),
        "sum_alt": sum_alt.astype(np.int64),
        "sum_total": sum_total.astype(np.int64),
        "n_observed_samples": n_observed_samples.astype(np.int32),
        "observed_sites": observed_sites,
    }


def should_apply_site_af_estimator(
    estimator: str,
    current_site_pop_af: np.ndarray,
    site_mask: np.ndarray,
) -> bool:
    """Return whether naive-Bayes AF estimates should replace input AFs."""
    if estimator not in SITE_AF_ESTIMATORS:
        raise ValueError(
            f"Unknown site AF estimator: {estimator!r}. "
            f"Choose one of {SITE_AF_ESTIMATORS}."
        )
    if estimator == "off":
        return False
    if estimator == "naive-bayes":
        return True

    # Auto mode is intentionally conservative. We only enable AF estimation
    # when the representation feeding the AF model is known to be coherent,
    # and legacy fallback site_data is not.
    return False


def is_fallback_minor_site_data(
    current_site_pop_af: np.ndarray,
    site_mask: np.ndarray,
) -> bool:
    """Return whether site_data matches preprocess fallback minor-count mode."""

    observed_sites = np.any(site_mask, axis=2)
    if not np.any(observed_sites):
        return False
    observed_pop_af = current_site_pop_af[observed_sites]
    return bool(np.allclose(observed_pop_af, 0.5, atol=1e-6, rtol=0.0))


def resolve_site_af_estimator_application(
    estimator: str,
    current_site_pop_af: np.ndarray,
    site_mask: np.ndarray,
) -> bool:
    """Resolve whether site AF estimation is valid for the current inputs."""
    if is_fallback_minor_site_data(current_site_pop_af, site_mask):
        if estimator == "naive-bayes":
            raise ValueError(
                "--site-af-estimator naive-bayes is incompatible with legacy fallback "
                "site_data that stores per-sample minor counts with site_pop_af=0.5. "
                "Disable the estimator or regenerate site_data with a coherent "
                "site-level alt-allele identity."
            )
        return False

    return should_apply_site_af_estimator(
        estimator,
        current_site_pop_af,
        site_mask,
    )


def build_site_af_estimates(
    data: DepthData,
    input_site_pop_af: np.ndarray,
    naive_bayes_site_pop_af: np.ndarray,
    effective_site_pop_af: np.ndarray,
    posterior_sd: np.ndarray,
    pooled_observed_af: np.ndarray,
    sum_alt: np.ndarray,
    sum_total: np.ndarray,
    n_observed_samples: np.ndarray,
) -> pd.DataFrame:
    """Build a per-site AF estimates table for diagnostics and reuse."""
    if data.site_mask is None:
        return pd.DataFrame()

    observed_sites = data.site_mask.detach().cpu().numpy().any(axis=2)
    rows: list[dict] = []

    for bin_idx in range(data.n_bins):
        for site_idx in np.where(observed_sites[bin_idx])[0]:
            input_af = float(input_site_pop_af[bin_idx, site_idx])
            naive_af = float(naive_bayes_site_pop_af[bin_idx, site_idx])
            effective_af = float(effective_site_pop_af[bin_idx, site_idx])
            rows.append(
                {
                    "chr": data.chr[bin_idx],
                    "bin_start": int(data.start[bin_idx]),
                    "bin_end": int(data.end[bin_idx]),
                    "site_index": int(site_idx),
                    "n_observed_samples": int(n_observed_samples[bin_idx, site_idx]),
                    "sum_alt": int(sum_alt[bin_idx, site_idx]),
                    "sum_total": int(sum_total[bin_idx, site_idx]),
                    "pooled_observed_af": float(pooled_observed_af[bin_idx, site_idx]),
                    "input_site_pop_af": input_af,
                    "naive_bayes_site_pop_af": naive_af,
                    "effective_site_pop_af": effective_af,
                    "posterior_sd": float(posterior_sd[bin_idx, site_idx]),
                    "site_pop_af_changed": not np.isclose(
                        input_af,
                        effective_af,
                        atol=1e-8,
                        rtol=0.0,
                    ),
                }
            )

    return pd.DataFrame(rows)


# ── artifact I/O ────────────────────────────────────────────────────────────


def load_inference_artifacts(
    path: str,
) -> Tuple[Dict[str, np.ndarray], Dict[str, np.ndarray]]:
    """Load MAP estimates and CN posterior from a saved ``.npz`` archive.

    Args:
        path: Path to ``inference_artifacts.npz`` (produced by ``infer``).

    Returns:
        Tuple of ``(map_estimates, cn_posterior)`` dictionaries with
        NumPy arrays.
    """
    npz = np.load(path, allow_pickle=True)
    map_estimates: Dict[str, np.ndarray] = {}
    cn_posterior: Dict[str, np.ndarray] = {}
    for key in npz.files:
        if key.startswith("cn_post_"):
            cn_posterior[key[len("cn_post_"):]] = npz[key]
        else:
            map_estimates[key] = npz[key]
    return map_estimates, cn_posterior


# ── summary statistics ──────────────────────────────────────────────────────


def _get_chr_info(
    data: DepthData,
    cn_probs: np.ndarray,
    chr_name: str,
    sample_idx: int,
    n_states: int = 6,
) -> Tuple[int, float, int, np.ndarray, int]:
    """Return majority-vote contig-ploidy summary for one chromosome/sample."""
    mask = data.chr == chr_name
    n_bins = int(mask.sum())
    if n_bins == 0:
        empty = np.zeros(n_states, dtype=np.float32)
        return 0, 0.0, 0, empty, 0

    chr_cn_map = np.argmax(cn_probs[mask, sample_idx, :], axis=1)
    chr_cnq = compute_cnq_from_probabilities(
        cn_probs[mask, sample_idx, :],
    )
    best_cn, best_prob, ploidy_fractions, plq = summarize_contig_ploidy_from_bin_calls(
        chr_cn_map,
        chr_cnq,
        n_states=n_states,
    )
    return best_cn, best_prob, n_bins, ploidy_fractions, plq


def detect_aneuploidies(
    data: DepthData,
    cn_posterior: Dict[str, np.ndarray],
    prob_threshold: float = 0.5,
) -> Dict[int, List[Tuple[str, int, float]]]:
    """Detect per-chromosome aneuploidies for every sample.

    Args:
        data: :class:`DepthData` used for inference.
        cn_posterior: Posterior from :meth:`CNVModel.run_discrete_inference`.
        prob_threshold: Minimum mean CN probability to call an aneuploidy.

    Returns:
        Dictionary mapping sample index → list of
        ``(chr_name, cn_state, mean_prob)`` tuples for aneuploid chromosomes.
    """
    cn_probs = cn_posterior["cn_posterior"]
    unique_chrs = np.unique(data.chr)

    sex_chrs = {"chrX", "chrY"}
    autosomes = [c for c in unique_chrs if c not in sex_chrs]

    aneuploid: Dict[int, List[Tuple[str, int, float]]] = {
        i: [] for i in range(data.n_samples)
    }

    # Autosomes: aneuploidy when CN ≠ 2 and high confidence
    for chr_name in autosomes:
        for si in range(data.n_samples):
            best_cn, mean_prob, _, _, _ = _get_chr_info(
                data,
                cn_probs,
                chr_name,
                si,
            )
            if best_cn != 2 and mean_prob > prob_threshold:
                aneuploid[si].append((chr_name, best_cn, mean_prob))

    # Sex chromosomes: aneuploidy when karyotype is not XX or XY
    for si in range(data.n_samples):
        x_cn, x_prob, x_bins, _, _ = _get_chr_info(
            data,
            cn_probs,
            "chrX",
            si,
        )
        y_cn, y_prob, y_bins, _, _ = _get_chr_info(
            data,
            cn_probs,
            "chrY",
            si,
        )

        if x_cn is None and y_cn is None:
            continue

        x_ok = x_prob > prob_threshold if x_bins > 0 else True
        y_ok = y_prob > prob_threshold if y_bins > 0 else True
        is_XX = x_cn == 2 and y_cn == 0
        is_XY = x_cn == 1 and y_cn == 1

        if not (is_XX or is_XY) and x_ok and y_ok:
            if x_bins > 0:
                aneuploid[si].append(("chrX", x_cn, x_prob))
            if y_bins > 0:
                aneuploid[si].append(("chrY", y_cn, y_prob))

    return aneuploid


# ── result assembly ─────────────────────────────────────────────────────────


def build_bin_stats(
    data: DepthData,
    map_estimates: Dict[str, np.ndarray],
    cn_posterior: Dict[str, np.ndarray],
    af_table: Optional[np.ndarray] = None,
    min_het_alt: int = 3,
    min_het_af: float = 0.2,
    max_het_af: float = 0.8,
) -> pd.DataFrame:
    """Build a per-bin, per-sample results DataFrame.

    Args:
        data: :class:`DepthData` instance.
        map_estimates: MAP estimates dict.
        cn_posterior: Discrete posterior dict.
        af_table: Optional precomputed AF log-likelihood table of shape
            ``(n_states, n_bins, n_samples)``.
        min_het_alt: Minimum alt-allele read count to call a site heterozygous.
        min_het_af: Minimum allele fraction to classify a site as
            heterozygous (default 0.2).
        max_het_af: Maximum allele fraction to classify a site as
            heterozygous (default 0.8).

    Returns:
        DataFrame with columns for chr, start, end, sample, observed depth,
        MAP CN, per-state probabilities, model parameters, and (when site
        data is available) per-bin allele fraction summaries.
    """
    cn_probs = cn_posterior["cn_posterior"]
    cnq = np.asarray(
        cn_posterior.get("cnq", compute_cnq_from_probabilities(cn_probs))
    )
    cn_map_stability = cn_posterior.get("cn_map_stability")

    # Precompute per-bin AF summaries if site data available
    has_af = data.site_alt is not None
    if has_af:
        sa = data.site_alt.cpu().numpy()
        st = data.site_total.cpu().numpy()
        sp = data.site_pop_af.cpu().numpy()
        sm = data.site_mask.cpu().numpy()
        site_af = sa / np.maximum(st, 1)
        # All-sites stats
        n_sites_arr = sm.sum(axis=1).astype(int)   # (n_bins, n_samples)
        af_sum_arr = (site_af * sm).sum(axis=1)    # (n_bins, n_samples)
        site_pop_af_sum_arr = (sp[:, :, np.newaxis] * sm).sum(axis=1)
        mean_af_arr = np.full(n_sites_arr.shape, np.nan, dtype=np.float64)
        mean_site_pop_af_arr = np.full(n_sites_arr.shape, np.nan, dtype=np.float64)
        np.divide(
            af_sum_arr,
            n_sites_arr,
            out=mean_af_arr,
            where=n_sites_arr > 0,
        )
        np.divide(
            site_pop_af_sum_arr,
            n_sites_arr,
            out=mean_site_pop_af_arr,
            where=n_sites_arr > 0,
        )
        # Het-only stats: alt count above threshold AND AF in het range
        het_mask = sm & (sa >= min_het_alt) & (site_af >= min_het_af) & (site_af <= max_het_af)
        n_het_arr = het_mask.sum(axis=1).astype(int)
        het_af_sum = (site_af * het_mask).sum(axis=1)
        mean_het_af_arr = np.full(n_het_arr.shape, np.nan, dtype=np.float64)
        np.divide(
            het_af_sum,
            n_het_arr,
            out=mean_het_af_arr,
            where=n_het_arr > 0,
        )

    rows: list[dict] = []
    plot_baseline = _plot_depth_baseline_per_sample(data)
    background_keys = {
        "bin_epsilon",
        "background_bin_factors",
        "background_sample_factors",
    }
    bin_epsilon_matrix = None
    if any(key in map_estimates for key in background_keys):
        bin_epsilon_matrix = compose_additive_background_matrix(
            map_estimates.get("bin_epsilon"),
            data.n_bins,
            data.n_samples,
            contig_index=data.contig_index.detach().cpu().numpy(),
            n_contigs=data.n_contigs,
            background_bin_factors=map_estimates.get("background_bin_factors"),
            background_sample_factors=map_estimates.get("background_sample_factors"),
        )

    for i in range(data.n_bins):
        for j in range(data.n_samples):
            prob = cn_probs[i, j, :]
            cn_map_val = int(np.argmax(prob))
            row = {
                "chr": data.chr[i],
                "start": int(data.start[i]),
                "end": int(data.end[i]),
                "sample": data.sample_ids[j],
                "observed_depth": float(data.depth[i, j].cpu().numpy()),
                "cn_map": cn_map_val,
                "cnq": int(cnq[i, j]),
                "cn_prob_0": prob[0],
                "cn_prob_1": prob[1],
                "cn_prob_2": prob[2],
                "cn_prob_3": prob[3],
                "cn_prob_4": prob[4],
                "cn_prob_5": prob[5],
                "max_prob": float(prob.max()),
                "bin_bias": float(map_estimates["bin_bias"].flatten()[i]),
                "bin_var": float(map_estimates["bin_var"].flatten()[i]),
                "sample_var": float(map_estimates["sample_var"].flatten()[j]),
            }
            if cn_map_stability is not None:
                row["cn_map_stability"] = float(cn_map_stability[i, j])
            if "sample_depth" in map_estimates:
                sample_depth = float(map_estimates["sample_depth"].flatten()[j])
                bin_length_kb = float(data.bin_length_kb[i].cpu().numpy())
                row["sample_depth"] = float(
                    sample_depth
                )
                row["bin_length_kb"] = bin_length_kb
                row["plot_depth"] = float(
                    _plot_normalized_depth(
                        row["observed_depth"],
                        bin_length_kb,
                        plot_baseline[j],
                    )
                )
                row["bin_overdispersion"] = row["bin_var"]
                row["sample_overdispersion"] = row["sample_var"]
            if bin_epsilon_matrix is not None:
                row["bin_epsilon"] = float(bin_epsilon_matrix[i, j])
            if has_af:
                row["n_sites"] = int(n_sites_arr[i, j])
                row["mean_observed_af"] = float(mean_af_arr[i, j])
                row["mean_site_pop_af"] = float(mean_site_pop_af_arr[i, j])
                row["n_het_sites"] = int(n_het_arr[i, j])
                row["mean_het_af"] = float(mean_het_af_arr[i, j])
            if af_table is not None:
                row["af_log_lik"] = float(af_table[cn_map_val, i, j])
            rows.append(row)

    return pd.DataFrame(rows)


def build_chromosome_stats(
    data: DepthData,
    map_estimates: Dict[str, np.ndarray],
    cn_posterior: Dict[str, np.ndarray],
    aneuploid_map: Dict[int, List[Tuple[str, int, float]]],
    af_table: Optional[np.ndarray] = None,
    min_het_alt: int = 3,
    min_het_af: float = 0.2,
    max_het_af: float = 0.8,
) -> pd.DataFrame:
    """Build per-chromosome, per-sample summary statistics.

    Args:
        data: :class:`DepthData` instance.
        map_estimates: MAP estimates dict.
        cn_posterior: Discrete posterior dict.
        aneuploid_map: Output of :func:`detect_aneuploidies`.
        af_table: Optional precomputed AF log-likelihood table of shape
            ``(n_states, n_bins, n_samples)``.
        min_het_alt: Minimum alt-allele read count to call a site heterozygous.
        min_het_af: Minimum allele fraction to classify a site as
            heterozygous (default 0.2).
        max_het_af: Maximum allele fraction to classify a site as
            heterozygous (default 0.8).

    Returns:
        DataFrame with one row per (sample, chromosome).
    """
    cn_probs = cn_posterior["cn_posterior"]
    cnq = cn_posterior.get("cnq", compute_cnq_from_probabilities(cn_probs))
    unique_chrs = np.unique(data.chr)

    # Precompute per-bin AF summaries if site data available
    has_af = data.site_alt is not None
    if has_af:
        sa = data.site_alt.cpu().numpy()
        st = data.site_total.cpu().numpy()
        sp = data.site_pop_af.cpu().numpy()
        sm = data.site_mask.cpu().numpy()
        site_af = sa / np.maximum(st, 1)
        n_sites_per_bin = sm.sum(axis=1).astype(int)   # (n_bins, n_samples)
        af_sum_per_bin = (site_af * sm).sum(axis=1)    # (n_bins, n_samples)
        site_pop_af_sum_per_bin = (sp[:, :, np.newaxis] * sm).sum(axis=1)
        het_mask = sm & (sa >= min_het_alt) & (site_af >= min_het_af) & (site_af <= max_het_af)
        n_het_per_bin = het_mask.sum(axis=1).astype(int)
        het_af_sum_per_bin = (site_af * het_mask).sum(axis=1)

    rows: list[dict] = []
    plot_baseline = _plot_depth_baseline_per_sample(data)
    n_states = cn_probs.shape[-1]

    for si in range(data.n_samples):
        aneu_set = {c for c, _, _ in aneuploid_map.get(si, [])}
        for chr_name in unique_chrs:
            mask = data.chr == chr_name
            n_bins = int(mask.sum())
            if n_bins == 0:
                continue

            chr_cn_map = np.argmax(cn_probs[mask, si, :], axis=1)
            best_cn, mean_prob, ploidy_fractions, plq = summarize_contig_ploidy_from_bin_calls(
                chr_cn_map,
                cnq[mask, si],
                n_states=n_states,
            )
            depths = data.depth[mask, si].cpu().numpy()

            row = {
                "sample": data.sample_ids[si],
                "chromosome": chr_name,
                "copy_number": best_cn,
                "mean_cn_probability": mean_prob,
                "plq": plq,
                "is_aneuploid": chr_name in aneu_set,
                "n_bins": n_bins,
                "mean_depth": float(np.mean(depths)),
                "std_depth": float(np.std(depths)),
                "median_depth": float(np.median(depths)),
                "mad_depth": float(stats.median_abs_deviation(depths)),
                "sample_var_map": float(
                    map_estimates["sample_var"].flatten()[si]
                ),
            }
            for cn_state in range(n_states):
                row[f"ploidy_prob_{cn_state}"] = float(ploidy_fractions[cn_state])
            if "sample_depth" in map_estimates:
                sample_depth = float(map_estimates["sample_depth"].flatten()[si])
                plot_depths = _plot_normalized_depth(
                    depths,
                    data.bin_length_kb[mask].cpu().numpy(),
                    plot_baseline[si],
                )
                row["sample_depth_map"] = float(
                    sample_depth
                )
                row["plot_mean_depth"] = float(np.mean(plot_depths))
                row["plot_std_depth"] = float(np.std(plot_depths))
                row["plot_median_depth"] = float(np.median(plot_depths))
                row["plot_mad_depth"] = float(
                    stats.median_abs_deviation(plot_depths)
                )
                row["sample_overdispersion_map"] = row["sample_var_map"]
            if has_af:
                chr_n_sites = int(n_sites_per_bin[mask, si].sum())
                chr_af_sum = float(af_sum_per_bin[mask, si].sum())
                row["n_sites"] = chr_n_sites
                row["mean_observed_af"] = (
                    chr_af_sum / chr_n_sites if chr_n_sites > 0
                    else float('nan')
                )
                chr_site_pop_af_sum = float(site_pop_af_sum_per_bin[mask, si].sum())
                row["mean_site_pop_af"] = (
                    chr_site_pop_af_sum / chr_n_sites if chr_n_sites > 0
                    else float("nan")
                )
                chr_n_het = int(n_het_per_bin[mask, si].sum())
                chr_het_af_sum = float(het_af_sum_per_bin[mask, si].sum())
                row["n_het_sites"] = chr_n_het
                row["mean_het_af"] = (
                    chr_het_af_sum / chr_n_het if chr_n_het > 0
                    else float('nan')
                )
            if af_table is not None:
                bin_indices = np.where(mask)[0]
                row["af_log_lik"] = float(
                    af_table[best_cn, bin_indices, si].sum()
                )
            rows.append(row)

    df = pd.DataFrame(rows).sort_values(["sample", "chromosome"])
    return df


# ── CLI ─────────────────────────────────────────────────────────────────────


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the infer subcommand."""
    p = argparse.ArgumentParser(
        description="Train Bayesian model and run CN inference",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "-i", "--input", required=True,
        help="Preprocessed depth TSV (output of 'preprocess')",
    )
    p.add_argument(
        "-o", "--output-dir", required=True,
        help="Output directory",
    )
    p.add_argument(
        "--depth-space", choices=["auto", *DEPTH_SPACES], default="auto",
        help="Interpret the input matrix as normalized depth or raw counts. 'auto' first consults preprocess observation_type.txt and otherwise falls back to the observation likelihood.",
    )

    # Model priors
    g = p.add_argument_group("model priors")
    g.add_argument(
        "--autosome-prior-mode",
        choices=list(AUTOSOME_PRIOR_MODES),
        default="dirichlet",
        help="Autosomal CN prior family",
    )
    g.add_argument("--alpha-ref", type=float, default=50,
                   help="Dirichlet concentration for CN=2 on autosomes")
    g.add_argument("--alpha-non-ref", type=float, default=1.0,
                   help="Dirichlet concentration for other CN states")
    g.add_argument(
        "--autosome-nonref-mean-alpha",
        type=float,
        default=1.0,
        help="Beta prior alpha for the cohort-wide autosomal non-reference mean in shrinkage mode",
    )
    g.add_argument(
        "--autosome-nonref-mean-beta",
        type=float,
        default=19.0,
        help="Beta prior beta for the cohort-wide autosomal non-reference mean in shrinkage mode",
    )
    g.add_argument(
        "--autosome-nonref-concentration",
        type=float,
        default=20.0,
        help="Per-bin shrinkage concentration toward the learned autosomal non-reference mean in shrinkage mode",
    )
    g.add_argument("--var-bias-bin", type=float, default=0.01,
                   help="LogNormal scale for per-bin bias")
    g.add_argument("--var-sample", type=float, default=0.01,
                   help="Exponential mean for per-sample variance")
    g.add_argument("--var-bin", type=float, default=0.01,
                   help="Exponential mean for per-bin variance")
    g.add_argument("--epsilon-mean", type=float, default=DEFAULT_EPSILON_MEAN,
                   help="Mean of the Gamma prior on per-contig additive background depth; "
                        "increase to absorb low-level zero-copy background depth, set to 0 to disable")
    g.add_argument(
        "--epsilon-concentration",
        type=float,
        default=DEFAULT_EPSILON_CONCENTRATION,
        help=(
            "Gamma concentration for per-contig additive background depth. "
            "Values below 1 keep most entries near zero while retaining a heavier "
            "positive tail; 1.0 matches the old Exponential prior."
        ),
    )
    g.add_argument(
        "--background-factors",
        type=int,
        default=DEFAULT_BACKGROUND_FACTORS,
        help=(
            "Number of low-rank additive background factors. "
            "Set to 0 to disable the structured background term."
        ),
    )
    g.add_argument("--guide-type", choices=["delta", "diagonal", "lowrank"], default="delta",
                   help="Variational guide type")
    g.add_argument("--obs-likelihood", choices=["auto", *list(OBS_LIKELIHOODS)], default="auto",
                   help="Observation likelihood for the depth/count matrix. 'auto' resolves from preprocess observation_type.txt when present.")
    g.add_argument("--obs-df", type=float, default=3.5,
                   help="Student-t degrees of freedom when --obs-likelihood=studentt")
    g.add_argument("--sample-depth-max", type=float, default=10000.0,
                   help="Upper clip used when deriving the empirical-Bayes LogNormal prior center for per-sample depth scale when --obs-likelihood=negative_binomial")
    g.add_argument("--freeze-bin-bias", action="store_true",
                   help="Hold per-bin bin_bias fixed at 1 during training. Useful for diagnosing whether low fitted bin_bias drives non-neutral CN calls.")
    g.add_argument("--freeze-cn-prior", action="store_true",
                   help="Hold the learned per-bin CN prior fixed at the default chromosome-type prior during training. Useful for diagnosing whether recurrent high-copy bins are driven by learned per-bin prior entrenchment.")
    g.add_argument("--freeze-sample-var", action="store_true",
                   help="Hold per-sample sample_var fixed at 0 during training. Useful for diagnosing whether a noisy-sample tail is driven by sample-specific overdispersion.")
    g.add_argument("--freeze-sample-depth", action="store_true",
                   help="Hold per-sample sample_depth fixed at the autosomal-median counts/kb "
                        "anchor during training (raw-count runs only). Useful for diagnosing "
                        "sample_depth collapse.")

    # Sex chromosome priors
    g = p.add_argument_group("sex chromosome priors")
    g.add_argument("--alpha-sex-ref", type=float, default=1.0,
                   help="Dirichlet concentration for CN=2 on sex chromosomes "
                        "(flat by default; sex-CN coupling handles "
                        "sex-dependent CN)")
    g.add_argument("--alpha-sex-non-ref", type=float, default=1.0,
                   help="Dirichlet concentration for other CN states on "
                        "sex chromosomes")
    g.add_argument("--sex-cn-weight", type=float, default=3.0,
                   help="Weight of the sex-CN coupling factor "
                        "(0 to disable)")
    g.add_argument("--sex-prior", type=float, nargs=2,
                   default=[0.5, 0.5], metavar=("P_XX", "P_XY"),
                   help="Prior probabilities for XX and XY karyotypes")

    # Allele fraction (per-site model)
    g = p.add_argument_group("allele fraction")
    g.add_argument("--site-data", default=None,
                   help="Per-site allele data .npz (output of 'preprocess')")
    g.add_argument(
        "--site-af-estimator",
        choices=list(SITE_AF_ESTIMATORS),
        default="auto",
        help=(
            "How to derive the site_pop_af values used during infer. "
            "'auto' is conservative and keeps the input AFs for current "
            "preprocess outputs, 'naive-bayes' explicitly replaces them when "
            "the site encoding is coherent, and 'off' always keeps the input values. "
            "When --learn-site-af is enabled, these values are used as prior centers."
        ),
    )
    g.add_argument(
        "--learn-site-af",
        action="store_true",
        default=False,
        help=(
            "Infer site_pop_af directly inside the model using a Delta guide for site AFs. "
            "This disables AF-table precomputation and uses the current site_pop_af values as Beta prior means."
        ),
    )
    g.add_argument(
        "--site-af-prior-strength",
        type=float,
        default=20.0,
        help="Strength of the Beta prior used when --learn-site-af is enabled",
    )
    g.add_argument(
        "--site-af-prior-alpha",
        type=float,
        default=1.0,
        help="Beta prior alpha for naive-Bayes site AF estimation",
    )
    g.add_argument(
        "--site-af-prior-beta",
        type=float,
        default=1.0,
        help="Beta prior beta for naive-Bayes site AF estimation",
    )
    g.add_argument("--af-concentration", type=float, default=50.0,
                   help="BetaBinomial concentration for allele fraction model")
    g.add_argument("--af-weight", type=float, default=DEFAULT_AF_WEIGHT,
                   help="Fixed AF scale when --fixed-af-temperature is used; otherwise the prior median for the learned global AF temperature (0 to disable)")
    g.add_argument(
        "--af-evidence-mode",
        choices=["absolute", "relative"],
        default="relative",
        help=(
            "How allele-fraction evidence enters the model: 'absolute' uses the "
            "summed AF log-likelihood directly, while 'relative' centers it "
            "against a fixed CN reference mixture before scaling."
        ),
    )
    p.set_defaults(learn_af_temperature=True)
    g.add_argument("--learn-af-temperature", dest="learn_af_temperature", action="store_true",
                   help="Learn a single global AF temperature instead of keeping --af-weight fixed (default)")
    g.add_argument("--fixed-af-temperature", dest="learn_af_temperature", action="store_false",
                   help="Keep --af-weight fixed instead of learning a single global AF temperature")
    g.add_argument("--af-temperature-prior-scale", type=float, default=0.5,
                   help="LogNormal prior scale when --learn-af-temperature is enabled")
    g.add_argument("--min-het-alt", type=int, default=3,
                   help="Minimum alt-allele read count to classify a site as "
                        "heterozygous in summary statistics")
    g.add_argument("--min-het-af", type=float, default=0.2,
                   help="Minimum allele fraction to classify a site as "
                        "heterozygous in summary statistics")
    g.add_argument("--max-het-af", type=float, default=0.8,
                   help="Maximum allele fraction to classify a site as "
                        "heterozygous in summary statistics")

    # Training
    g = p.add_argument_group("training")
    g.add_argument("--max-iter", type=int, default=5000)
    g.add_argument("--lr-init", type=float, default=0.02)
    g.add_argument("--lr-min", type=float, default=0.01)
    g.add_argument("--lr-decay", type=float, default=500)
    g.add_argument(
        "--grad-clip-norm",
        type=float,
        default=10.0,
        help=(
            "Clip the global SVI gradient norm to this value during training; "
            "set to 0 or a negative value to disable clipping"
        ),
    )
    g.add_argument(
        "--svi-init-restarts",
        type=int,
        default=100,
        help=(
            "Number of candidate autoguide initializations to score before "
            "gradient descent. The first candidate uses the anchored default "
            "init and the rest sample random prior starts."
        ),
    )
    g.add_argument("--log-freq", type=int, default=50)
    g.add_argument("--jit", action="store_true", default=False)
    g.add_argument("--early-stopping", action="store_true", default=True)
    g.add_argument("--no-early-stopping", dest="early_stopping",
                   action="store_false")
    g.add_argument("--patience", type=int, default=50)
    g.add_argument(
        "--elbo-window",
        type=int,
        default=50,
        help="Iterations per rolling ELBO window for early stopping",
    )
    g.add_argument(
        "--elbo-rtol",
        type=float,
        default=1e-3,
        help="Relative tolerance between successive rolling ELBO windows",
    )

    # Inference
    g = p.add_argument_group("discrete inference")
    g.add_argument("--prob-threshold", type=float, default=0.5,
                   help="Min mean CN probability for aneuploidy call")
    g.add_argument(
        "--cn-inference-method",
        choices=["current", "median", "multi-draw"],
        default="multi-draw",
        help=(
            "How to handle continuous latent uncertainty before CN inference: "
            "'current' uses a single guide draw (historical behavior), "
            "'median' plugs in guide medians, and 'multi-draw' averages CN "
            "posteriors over repeated guide draws."
        ),
    )
    g.add_argument(
        "--cn-inference-draws",
        type=int,
        default=100,
        help="Number of guide draws to average when --cn-inference-method=multi-draw",
    )

    # Device
    p.add_argument("--device", choices=["cpu", "cuda"], default="cpu")

    return p.parse_args()


def main() -> None:
    """Entry point for ``gatk-sv-ploidy infer``."""
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    # ── reproducibility ─────────────────────────────────────────────────
    pyro.enable_validation(True)
    pyro.distributions.enable_validation(True)
    pyro.set_rng_seed(42)
    torch.manual_seed(42)
    np.random.seed(42)

    # ── load data ───────────────────────────────────────────────────────
    logger.info("Loading preprocessed depth.")
    df = pd.read_csv(args.input, sep="\t", index_col=0)
    depth_space, obs_likelihood = _resolve_observation_model(
        args.depth_space,
        args.obs_likelihood,
        args.input,
    )

    # ── optional per-site allele data ─────────────────────────────────
    sd = None
    if args.site_data:
        sd = load_site_data(args.site_data)

    tensor_dtype = _inference_tensor_dtype()
    logger.info("Using float64 tensors for inference.")

    device = args.device
    data = DepthData(
        df, device=device, dtype=tensor_dtype,
        clamp_threshold=_clamp_threshold_for_depth_space(depth_space),
        depth_space=depth_space,
        site_data=sd,
    )
    logger.info(
        "Safe logging enabled: reporting cohort-level summaries only."
    )
    logger.info("Using %s input with %s observation model.", depth_space, obs_likelihood)
    if args.learn_site_af and data.site_alt is None:
        raise ValueError("--learn-site-af requires --site-data with per-site allele counts.")

    input_site_pop_af_np = None
    naive_bayes_site_pop_af_np = None
    effective_site_pop_af_np = None
    site_af_summary = None
    site_af_estimation_applied = False
    if data.site_alt is not None:
        input_site_pop_af_np = data.site_pop_af.detach().cpu().numpy().copy()
        site_alt_np = data.site_alt.detach().cpu().numpy()
        site_total_np = data.site_total.detach().cpu().numpy()
        site_mask_np = data.site_mask.detach().cpu().numpy()
        site_af_summary = estimate_site_pop_af_naive_bayes(
            site_alt_np,
            site_total_np,
            site_mask_np,
            prior_alpha=args.site_af_prior_alpha,
            prior_beta=args.site_af_prior_beta,
        )
        naive_bayes_site_pop_af_np = site_af_summary["posterior_mean"]
        fallback_minor_site_data = is_fallback_minor_site_data(
            input_site_pop_af_np,
            site_mask_np,
        )
        if fallback_minor_site_data and args.site_af_estimator != "off":
            logger.warning(
                "Detected legacy fallback site_data with uniform site_pop_af=0.5. "
                "These inputs store per-sample minor counts, so infer will not "
                "apply site AF estimation automatically."
            )
        site_af_estimation_applied = resolve_site_af_estimator_application(
            args.site_af_estimator,
            input_site_pop_af_np,
            site_mask_np,
        )
        effective_site_pop_af_np = (
            naive_bayes_site_pop_af_np
            if site_af_estimation_applied
            else input_site_pop_af_np
        )
        if site_af_estimation_applied:
            data.site_pop_af = torch.tensor(
                effective_site_pop_af_np,
                dtype=tensor_dtype,
                device=device,
            )
        observed_sites = site_af_summary["observed_sites"]
        observed_input = input_site_pop_af_np[observed_sites]
        observed_naive = naive_bayes_site_pop_af_np[observed_sites]
        logger.info(
            "Site AF estimation: mode=%s, applied=%s, sites=%d, input_median=%.3f, naive_bayes_median=%.3f",
            args.site_af_estimator,
            site_af_estimation_applied,
            int(observed_sites.sum()),
            float(np.median(observed_input)) if observed_input.size else float("nan"),
            float(np.median(observed_naive)) if observed_naive.size else float("nan"),
        )

    af_enabled = data.site_alt is not None and args.af_weight > 0
    learn_af_temperature = args.learn_af_temperature and af_enabled
    learn_site_af = args.learn_site_af and af_enabled

    if af_enabled:
        logger.info(
            "Allele-fraction evidence enabled (af_weight=%.2f, af_concentration=%.1f, mode=%s, summed over observed sites/bin, learn_af_temperature=%s)",
            args.af_weight,
            args.af_concentration,
            args.af_evidence_mode,
            learn_af_temperature,
        )
    elif data.site_alt is not None:
        logger.info(
            "Allele-fraction evidence available but disabled because af_weight=0."
        )
    if args.cn_inference_draws < 1:
        raise ValueError("--cn-inference-draws must be at least 1.")
    if args.guide_type == "delta" and args.cn_inference_method != "current":
        logger.info(
            "Guide type 'delta' is deterministic; %s will reduce to the current point-estimate behavior.",
            args.cn_inference_method,
        )
    if learn_site_af:
        logger.info(
            "Joint site AF inference enabled: using a Delta guide for site AFs with prior strength %.1f.",
            args.site_af_prior_strength,
        )

    # ── build & train model ─────────────────────────────────────────────
    model = CNVModel(
        n_states=6,
        autosome_prior_mode=args.autosome_prior_mode,
        alpha_ref=args.alpha_ref,
        alpha_non_ref=args.alpha_non_ref,
        autosome_nonref_mean_alpha=args.autosome_nonref_mean_alpha,
        autosome_nonref_mean_beta=args.autosome_nonref_mean_beta,
        autosome_nonref_concentration=args.autosome_nonref_concentration,
        var_bias_bin=args.var_bias_bin,
        var_sample=args.var_sample,
        var_bin=args.var_bin,
        epsilon_mean=args.epsilon_mean,
        epsilon_concentration=args.epsilon_concentration,
        background_factors=args.background_factors,
        device=device,
        dtype=tensor_dtype,
        guide_type=args.guide_type,
        af_concentration=args.af_concentration,
        af_weight=args.af_weight if data.site_alt is not None else 0.0,
        af_evidence_mode=args.af_evidence_mode,
        learn_af_temperature=learn_af_temperature,
        learn_site_pop_af=learn_site_af,
        site_af_prior_strength=args.site_af_prior_strength,
        af_temperature_prior_scale=args.af_temperature_prior_scale,
        alpha_sex_ref=args.alpha_sex_ref,
        alpha_sex_non_ref=args.alpha_sex_non_ref,
        sex_prior=tuple(args.sex_prior),
        sex_cn_weight=args.sex_cn_weight,
        obs_likelihood=obs_likelihood,
        obs_df=args.obs_df,
        sample_depth_max=args.sample_depth_max,
        freeze_bin_bias=args.freeze_bin_bias,
        freeze_cn_prior=args.freeze_cn_prior,
        freeze_sample_var=args.freeze_sample_var,
        freeze_sample_depth=args.freeze_sample_depth,
    )

    model.train(
        data,
        max_iter=args.max_iter,
        lr_init=args.lr_init,
        lr_min=args.lr_min,
        lr_decay=args.lr_decay,
        grad_clip_norm=args.grad_clip_norm,
        init_restarts=args.svi_init_restarts,
        log_freq=args.log_freq,
        jit=args.jit,
        early_stopping=args.early_stopping,
        patience=args.patience,
        convergence_window=args.elbo_window,
        convergence_rtol=args.elbo_rtol,
    )

    # ── save training loss ──────────────────────────────────────────────
    loss_df = pd.DataFrame(model.loss_history)
    loss_path = os.path.join(args.output_dir, "training_loss.tsv")
    loss_df.to_csv(loss_path, sep="\t", index=False)
    logger.info("Training loss saved.")

    # ── point estimates + exact discrete inference ──────────────────────
    point_estimate_method = args.cn_inference_method
    if args.cn_inference_method == "multi-draw":
        point_estimate_method = "median"

    map_est = model.get_map_estimates(
        data,
        estimate_method=point_estimate_method,
    )
    if learn_site_af and "site_pop_af_latent" in map_est:
        effective_site_pop_af_np = np.asarray(
            map_est["site_pop_af_latent"],
            dtype=np.float32,
        )
        data.site_pop_af = torch.tensor(
            effective_site_pop_af_np,
            dtype=tensor_dtype,
            device=device,
        )
        observed_sites = site_af_summary["observed_sites"]
        observed_learned = effective_site_pop_af_np[observed_sites]
        logger.info(
            "Learned site AF MAP: sites=%d, median=%.3f",
            int(observed_sites.sum()),
            float(np.median(observed_learned)) if observed_learned.size else float("nan"),
        )

    # ── precompute AF table once for inference and output stats ─────────
    af_table_np = None
    if data.site_alt is not None and model.af_weight > 0:
        logger.info("Computing AF table for analytical inference and output statistics ...")
        with torch.no_grad():
            af_table_np = model._prepare_af_table_torch(
                data.site_alt, data.site_total,
                data.site_pop_af, data.site_mask,
                data.chr_type if hasattr(data, "chr_type") else None,
            ).cpu().numpy()

    posterior_draws: Dict[str, List[np.ndarray]] = {}
    if args.cn_inference_method == "multi-draw":
        cn_post = model.run_discrete_inference_multi_draw(
            data,
            n_draws=args.cn_inference_draws,
            af_table=af_table_np,
            draw_estimate_collector=posterior_draws,
        )
        if "sex_posterior" in cn_post:
            map_est["sex_karyotype"] = np.argmax(
                cn_post["sex_posterior"],
                axis=1,
            ).astype(np.int64)
    else:
        cn_post = model.run_discrete_inference(
            data,
            map_estimates=map_est,
            af_table=af_table_np,
        )
    cn_post["cnq"] = compute_cnq_from_probabilities(cn_post["cn_posterior"])
    if "cn_map_stability" not in cn_post:
        cn_post["cn_map_stability"] = np.ones(
            cn_post["cn_posterior"].shape[:2],
            dtype=np.float32,
        )

    depth_only_cn_post = None
    if af_table_np is not None and model.af_weight > 0:
        logger.info("Computing depth-only discrete posterior under fixed fitted latents for safe AF diagnostics ...")
        depth_only_cn_post = model.run_discrete_inference(
            data,
            map_estimates=map_est,
            af_table=np.zeros_like(af_table_np),
        )
        depth_only_cn_post["cnq"] = compute_cnq_from_probabilities(
            depth_only_cn_post["cn_posterior"]
        )

    safe_diagnostic_messages = build_safe_inference_diagnostic_messages(
        data,
        map_est,
        cn_post,
        depth_only_cn_posterior=depth_only_cn_post,
        af_table=af_table_np,
        af_concentration=args.af_concentration,
    )
    for message in safe_diagnostic_messages:
        logger.info(message)

    safe_diagnostics_path = os.path.join(args.output_dir, "safe_inference_diagnostics.txt")
    with open(safe_diagnostics_path, "w", encoding="utf-8") as handle:
        for message in safe_diagnostic_messages:
            handle.write(message)
            handle.write("\n")
    logger.info("Safe inference diagnostics saved.")

    if site_af_summary is not None:
        site_af_estimates_df = build_site_af_estimates(
            data,
            input_site_pop_af=input_site_pop_af_np,
            naive_bayes_site_pop_af=naive_bayes_site_pop_af_np,
            effective_site_pop_af=effective_site_pop_af_np,
            posterior_sd=site_af_summary["posterior_sd"],
            pooled_observed_af=site_af_summary["pooled_observed_af"],
            sum_alt=site_af_summary["sum_alt"],
            sum_total=site_af_summary["sum_total"],
            n_observed_samples=site_af_summary["n_observed_samples"],
        )
        site_af_path = os.path.join(args.output_dir, "site_af_estimates.tsv.gz")
        site_af_estimates_df.to_csv(
            site_af_path,
            sep="\t",
            index=False,
            compression="gzip",
        )
        logger.info("Site AF estimates saved.")

    # ── save inference artifacts for downstream tools (ppd, plot) ──────
    artifacts_path = os.path.join(args.output_dir, "inference_artifacts.npz")
    artifact_dict = {k: v for k, v in map_est.items()}
    artifact_dict["obs_likelihood"] = np.asarray(model.obs_likelihood)
    artifact_dict["obs_df"] = np.asarray(model.obs_df)
    artifact_dict["depth_space"] = np.asarray(depth_space)
    artifact_dict["model_n_states"] = np.asarray(model.n_states)
    artifact_dict["model_autosome_prior_mode"] = np.asarray(model.autosome_prior_mode)
    artifact_dict["model_alpha_ref"] = np.asarray(model.alpha_ref)
    artifact_dict["model_alpha_non_ref"] = np.asarray(model.alpha_non_ref)
    artifact_dict["model_autosome_nonref_mean_alpha"] = np.asarray(
        model.autosome_nonref_mean_alpha
    )
    artifact_dict["model_autosome_nonref_mean_beta"] = np.asarray(
        model.autosome_nonref_mean_beta
    )
    artifact_dict["model_autosome_nonref_concentration"] = np.asarray(
        model.autosome_nonref_concentration
    )
    artifact_dict["model_var_bias_bin"] = np.asarray(model.var_bias_bin)
    artifact_dict["model_var_sample"] = np.asarray(model.var_sample)
    artifact_dict["model_var_bin"] = np.asarray(model.var_bin)
    artifact_dict["model_af_concentration"] = np.asarray(model.af_concentration)
    artifact_dict["model_af_weight"] = np.asarray(model.af_weight)
    artifact_dict["model_af_evidence_mode"] = np.asarray(model.af_evidence_mode)
    artifact_dict["model_learn_af_temperature"] = np.asarray(model.learn_af_temperature)
    artifact_dict["model_learn_site_pop_af"] = np.asarray(model.learn_site_pop_af)
    artifact_dict["model_site_af_prior_strength"] = np.asarray(
        model.site_af_prior_strength
    )
    artifact_dict["model_af_temperature_prior_scale"] = np.asarray(
        model.af_temperature_prior_scale
    )
    artifact_dict["model_alpha_sex_ref"] = np.asarray(model.alpha_sex_ref)
    artifact_dict["model_alpha_sex_non_ref"] = np.asarray(model.alpha_sex_non_ref)
    artifact_dict["model_sex_prior"] = np.asarray(model.sex_prior)
    artifact_dict["model_sex_cn_weight"] = np.asarray(model.sex_cn_weight)
    artifact_dict["model_guide_type"] = np.asarray(model.guide_type)
    artifact_dict["epsilon_mean"] = np.asarray(model.epsilon_mean)
    artifact_dict["model_epsilon_concentration"] = np.asarray(
        model.epsilon_concentration
    )
    artifact_dict["model_background_factors"] = np.asarray(model.background_factors)
    artifact_dict["model_background_sample_scale"] = np.asarray(
        model.background_sample_scale
    )
    artifact_dict["model_background_bin_scale"] = np.asarray(
        model.background_bin_scale
    )
    artifact_dict["sample_depth_max"] = np.asarray(model.sample_depth_max)
    artifact_dict["freeze_bin_bias"] = np.asarray(model.freeze_bin_bias)
    artifact_dict["freeze_cn_prior"] = np.asarray(model.freeze_cn_prior)
    artifact_dict["freeze_sample_var"] = np.asarray(model.freeze_sample_var)
    artifact_dict["freeze_sample_depth"] = np.asarray(model.freeze_sample_depth)
    artifact_dict["continuous_estimate_method"] = np.asarray(point_estimate_method)
    artifact_dict["cn_inference_method"] = np.asarray(args.cn_inference_method)
    artifact_dict["cn_inference_draws"] = np.asarray(
        args.cn_inference_draws if args.cn_inference_method == "multi-draw" else 1
    )
    if input_site_pop_af_np is not None:
        artifact_dict["site_pop_af_input"] = input_site_pop_af_np.astype(np.float32)
        artifact_dict["site_pop_af_naive_bayes"] = naive_bayes_site_pop_af_np.astype(
            np.float32
        )
        artifact_dict["site_pop_af_effective"] = effective_site_pop_af_np.astype(
            np.float32
        )
        artifact_dict["site_af_estimator"] = np.asarray(args.site_af_estimator)
        artifact_dict["site_af_estimation_applied"] = np.asarray(
            site_af_estimation_applied
        )
        artifact_dict["site_af_learned_in_model"] = np.asarray(learn_site_af)
        artifact_dict["site_af_prior_alpha"] = np.asarray(args.site_af_prior_alpha)
        artifact_dict["site_af_prior_beta"] = np.asarray(args.site_af_prior_beta)
        artifact_dict["site_af_prior_strength"] = np.asarray(args.site_af_prior_strength)
    for site, draws in posterior_draws.items():
        if draws:
            artifact_dict[f"posterior_draws_{site}"] = np.stack(draws, axis=0)
    for k, v in cn_post.items():
        artifact_dict[f"cn_post_{k}"] = v
    np.savez_compressed(artifacts_path, **artifact_dict)
    logger.info("Inference artifacts saved.")

    # ── detect aneuploidies ─────────────────────────────────────────────
    aneuploid_map = detect_aneuploidies(
        data, cn_post, prob_threshold=args.prob_threshold
    )

    # ── write bin stats ─────────────────────────────────────────────────
    bin_df = build_bin_stats(
        data, map_est, cn_post,
        af_table=af_table_np,
        min_het_alt=args.min_het_alt,
        min_het_af=args.min_het_af,
        max_het_af=args.max_het_af,
    )
    bin_path = os.path.join(args.output_dir, "bin_stats.tsv.gz")
    bin_df.to_csv(bin_path, sep="\t", index=False, compression="gzip")
    logger.info("Bin statistics saved.")

    # ── write chromosome stats ──────────────────────────────────────────
    chr_df = build_chromosome_stats(
        data, map_est, cn_post, aneuploid_map,
        af_table=af_table_np,
        min_het_alt=args.min_het_alt,
        min_het_af=args.min_het_af,
        max_het_af=args.max_het_af,
    )
    chr_path = os.path.join(args.output_dir, "chromosome_stats.tsv")
    chr_df.to_csv(chr_path, sep="\t", index=False)
    logger.info("Chromosome statistics saved.")


if __name__ == "__main__":
    main()
