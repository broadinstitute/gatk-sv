"""
Infer subcommand — train model and run copy-number inference.

Loads preprocessed depth data, trains the Pyro CNV model, obtains MAP
estimates and exact discrete posteriors, then writes per-bin and
per-chromosome summary statistics.
"""

from __future__ import annotations

import argparse
import os
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import pyro
import torch
from scipy import stats

from gatk_sv_ploidy._logging import log_output_artifacts, tool_logging_context
from gatk_sv_ploidy._util import (
    compose_additive_background_matrix,
    DEFAULT_AF_WEIGHT,
    compute_cnq_from_probabilities,
    is_expected_allosome_copy_number_pair,
    format_count_summary,
    format_numeric_summary,
    get_sample_columns,
    summarize_contig_ploidy_from_bin_calls,
)
from gatk_sv_ploidy.data import DepthData, load_site_data
from gatk_sv_ploidy.models import (
    AUTOSOME_PRIOR_MODES,
    DEFAULT_AF_BACKGROUND_CONCENTRATION,
    DEFAULT_AF_OUTLIER_WEIGHT,
    DEFAULT_EPSILON_CONCENTRATION,
    DEFAULT_EPSILON_MEAN,
    DEFAULT_RAW_VARIANCE_POWER,
    CNVModel,
    _compose_effective_bin_bias_numpy,
    _resolve_fixed_site_pop_af_numpy,
    _site_level_marginalized_af_log_lik_numpy,
    estimate_af_background_concentration_numpy,
)
SITE_AF_ESTIMATORS = ("off", "auto", "naive-bayes")
_UNDETERMINED_CALL = "UNDETERMINED"


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


def _coerce_bin_bias_matrix(
    map_estimates: Dict[str, np.ndarray],
    n_bins: int,
    n_samples: int,
) -> np.ndarray | None:
    """Return an effective bin-bias matrix when present or derivable."""
    if "bin_bias_matrix" in map_estimates:
        return _compose_effective_bin_bias_numpy(
            n_bins,
            n_samples,
            fixed_bias=map_estimates["bin_bias_matrix"],
        )
    if "bin_bias" in map_estimates:
        return _compose_effective_bin_bias_numpy(
            n_bins,
            n_samples,
            fixed_bias=map_estimates["bin_bias"],
        )
    return None


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
    af_concentration: float | np.ndarray = 50.0,
    af_outlier_weight: float = 0.0,
    af_background_concentration: float | None = None,
    autosomal_baseline_cn: np.ndarray | None = None,
    af_informative_weight: float = 1.0,
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

    background_keys = {"bin_epsilon"}
    bin_epsilon_matrix = None
    if any(key in map_estimates for key in background_keys):
        bin_epsilon_matrix = compose_additive_background_matrix(
            map_estimates.get("bin_epsilon"),
            data.n_bins,
            data.n_samples,
        )
        messages.append(
            format_numeric_summary(
                "Bin epsilon latent across bin-sample pairs",
                bin_epsilon_matrix,
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
    bin_bias_matrix = _coerce_bin_bias_matrix(
        map_estimates,
        data.n_bins,
        data.n_samples,
    )

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
                    messages.append(
                        format_numeric_summary(
                            f"{chrom_label} plot depth across {sex_label}-assigned bin-sample pairs",
                            plot_subset,
                            precision=4,
                        )
                    )
                    if bin_bias_matrix is not None:
                        if bin_epsilon_matrix is None:
                            epsilon_subset = np.zeros(
                                (bin_idx.size, sample_idx.size),
                                dtype=np.float64,
                            )
                        else:
                            epsilon_subset = bin_epsilon_matrix[np.ix_(bin_idx, sample_idx)]
                        bias_subset = bin_bias_matrix[np.ix_(bin_idx, sample_idx)]
                        expected_plot_depth = (
                            expected_cn * bias_subset +
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
                            bias_subset.mean(axis=1),
                            mean_abs_residual_by_bin,
                        )
                        if corr is not None:
                            messages.append(
                                f"Correlation(bin_bias, {chrom_label} mean absolute observed-minus-expected plot depth among {sex_label}-assigned samples by bin)="
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
                        outlier_weight=af_outlier_weight,
                        background_concentration=af_background_concentration,
                        background_baseline_cn=autosomal_baseline_cn,
                        informative_weight=af_informative_weight,
                    )
                    site_level_cn4 = _site_level_marginalized_af_log_lik_numpy(
                        site_alt,
                        site_total,
                        effective_site_pop_af,
                        site_mask,
                        cn_state=4,
                        n_states=cn_probs.shape[2],
                        concentration=af_concentration,
                        outlier_weight=af_outlier_weight,
                        background_concentration=af_background_concentration,
                        background_baseline_cn=autosomal_baseline_cn,
                        informative_weight=af_informative_weight,
                    )
                    site_level_cn5 = _site_level_marginalized_af_log_lik_numpy(
                        site_alt,
                        site_total,
                        effective_site_pop_af,
                        site_mask,
                        cn_state=5,
                        n_states=cn_probs.shape[2],
                        concentration=af_concentration,
                        outlier_weight=af_outlier_weight,
                        background_concentration=af_background_concentration,
                        background_baseline_cn=autosomal_baseline_cn,
                        informative_weight=af_informative_weight,
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


def _load_autosomal_baseline_cn(
    path: str | None,
    sample_ids: Sequence[str],
) -> np.ndarray:
    """Load per-sample autosomal baseline CN values or default to diploid."""
    if path is None:
        return np.full(len(sample_ids), 2, dtype=np.int64)

    baseline_df = pd.read_csv(path, sep="\t")
    required_columns = {"sample", "autosomal_baseline_cn"}
    missing_columns = required_columns - set(baseline_df.columns)
    if missing_columns:
        raise ValueError(
            "Autosomal baseline CN TSV is missing required columns: " +
            ", ".join(sorted(missing_columns))
        )

    if baseline_df["sample"].duplicated().any():
        duplicates = baseline_df.loc[
            baseline_df["sample"].duplicated(),
            "sample",
        ].astype(str).tolist()
        raise ValueError(
            "Autosomal baseline CN TSV contains duplicate sample rows: " +
            ", ".join(duplicates[:5])
        )

    baseline_series = pd.Series(
        baseline_df["autosomal_baseline_cn"].to_numpy(),
        index=baseline_df["sample"].astype(str),
    )
    missing_samples = [sample_id for sample_id in sample_ids if sample_id not in baseline_series.index]
    if missing_samples:
        preview = ", ".join(missing_samples[:5])
        more = "" if len(missing_samples) <= 5 else f" (+{len(missing_samples) - 5} more)"
        raise ValueError(
            "Autosomal baseline CN TSV is missing samples from the depth matrix: "
            f"{preview}{more}"
        )

    baseline = np.asarray(
        [baseline_series.loc[str(sample_id)] for sample_id in sample_ids],
        dtype=np.int64,
    )
    if np.any((baseline < 1) | (baseline > 5)):
        raise ValueError(
            "autosomal_baseline_cn values must be integers between 1 and 5."
        )
    return baseline


def _coerce_manifest_include_values(values: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(values):
        return values.fillna(True).astype(bool)
    normalized = values.astype(str).str.strip().str.lower()
    false_values = {"0", "false", "f", "no", "n", "exclude", "excluded"}
    true_values = {"1", "true", "t", "yes", "y", "include", "included"}
    include = pd.Series(True, index=values.index)
    include[normalized.isin(false_values)] = False
    include[normalized.isin(true_values)] = True
    include[values.isna()] = True
    return include.astype(bool)


def _subset_site_data_samples(
    site_data: Dict[str, np.ndarray] | None,
    selected_indices: Sequence[int],
    selected_sample_ids: Sequence[str],
) -> Dict[str, np.ndarray] | None:
    if site_data is None:
        return None
    index = np.asarray(selected_indices, dtype=np.int64)
    out = dict(site_data)
    for key in ("site_alt", "site_total", "site_mask"):
        out[key] = out[key][:, :, index]
    if "sample_ids" in out:
        out["sample_ids"] = np.asarray(list(selected_sample_ids), dtype=object)
    return out


def _filter_inputs_by_baseline_manifest(
    df: pd.DataFrame,
    site_data: Dict[str, np.ndarray] | None,
    baseline_path: str | None,
) -> tuple[pd.DataFrame, Dict[str, np.ndarray] | None]:
    if baseline_path is None:
        return df, site_data

    baseline_df = pd.read_csv(baseline_path, sep="\t")
    required_columns = {"sample", "autosomal_baseline_cn"}
    missing_columns = required_columns - set(baseline_df.columns)
    if missing_columns:
        raise ValueError(
            "Autosomal baseline CN TSV is missing required columns: " +
            ", ".join(sorted(missing_columns))
        )

    if baseline_df["sample"].duplicated().any():
        duplicates = baseline_df.loc[
            baseline_df["sample"].duplicated(),
            "sample",
        ].astype(str).tolist()
        raise ValueError(
            "Autosomal baseline CN TSV contains duplicate sample rows: " +
            ", ".join(duplicates[:5])
        )

    sample_cols = get_sample_columns(df)
    manifest = baseline_df.copy()
    manifest["sample"] = manifest["sample"].astype(str)
    manifest = manifest.set_index("sample", drop=False)
    missing_samples = [sample_id for sample_id in sample_cols if sample_id not in manifest.index]
    if missing_samples:
        preview = ", ".join(missing_samples[:5])
        more = "" if len(missing_samples) <= 5 else f" (+{len(missing_samples) - 5} more)"
        raise ValueError(
            "Autosomal baseline CN TSV is missing samples from the depth matrix: "
            f"{preview}{more}"
        )

    include = pd.Series(True, index=manifest.index)
    if "include_in_infer" in manifest.columns:
        include &= _coerce_manifest_include_values(manifest["include_in_infer"])
    if "baseline_cn_call" in manifest.columns:
        include &= (
            manifest["baseline_cn_call"].astype(str).str.upper() !=
            _UNDETERMINED_CALL
        )
    baseline_cn = pd.to_numeric(
        manifest["autosomal_baseline_cn"],
        errors="coerce",
    )
    include &= baseline_cn != 0

    selected_sample_ids = [sample_id for sample_id in sample_cols if bool(include.loc[sample_id])]
    if len(selected_sample_ids) == len(sample_cols):
        return df, site_data
    if not selected_sample_ids:
        raise ValueError(
            "All samples were excluded from infer by the autosomal baseline CN TSV."
        )

    selected_indices = [sample_cols.index(sample_id) for sample_id in selected_sample_ids]
    if site_data is not None and "sample_ids" in site_data:
        site_sample_ids = [
            str(sample_id) for sample_id in np.asarray(site_data["sample_ids"]).tolist()
        ]
        if site_sample_ids != sample_cols:
            raise ValueError(
                "site_data sample_ids do not match the preprocessed depth sample order."
            )

    metadata_cols = [column for column in df.columns if column not in sample_cols]
    filtered_df = df.loc[:, metadata_cols + selected_sample_ids].copy()
    filtered_site_data = _subset_site_data_samples(
        site_data,
        selected_indices,
        selected_sample_ids,
    )
    return filtered_df, filtered_site_data


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
    # when the representation feeding the AF model is known to be coherent.
    return False


def resolve_site_af_estimator_application(
    estimator: str,
    current_site_pop_af: np.ndarray,
    site_mask: np.ndarray,
) -> bool:
    """Resolve whether site AF estimation is valid for the current inputs."""
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
    autosomal_baseline_cn: np.ndarray | None = None,
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
    if autosomal_baseline_cn is None:
        autosomal_baseline_cn = np.full(data.n_samples, 2, dtype=np.int64)
    autosomal_baseline_cn = np.asarray(
        autosomal_baseline_cn,
        dtype=np.int64,
    ).reshape(-1)
    if autosomal_baseline_cn.shape[0] != data.n_samples:
        raise ValueError(
            "autosomal_baseline_cn must provide one baseline copy number per sample."
        )

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
            if best_cn != int(autosomal_baseline_cn[si]) and mean_prob > prob_threshold:
                aneuploid[si].append((chr_name, best_cn, mean_prob))

    # Sex chromosomes: aneuploidy when the X/Y pair is not consistent with
    # the sample's autosomal baseline CN.
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
        baseline_cn = int(autosomal_baseline_cn[si])
        is_expected_sex_pair = is_expected_allosome_copy_number_pair(
            x_cn,
            y_cn,
            baseline_cn,
        )

        if not is_expected_sex_pair and x_ok and y_ok:
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
    background_keys = {"bin_epsilon"}
    bin_epsilon_matrix = None
    if any(key in map_estimates for key in background_keys):
        bin_epsilon_matrix = compose_additive_background_matrix(
            map_estimates.get("bin_epsilon"),
            data.n_bins,
            data.n_samples,
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
    autosomal_baseline_cn = np.asarray(
        map_estimates.get(
            "autosomal_baseline_cn",
            np.full(data.n_samples, 2, dtype=np.int64),
        ),
        dtype=np.int64,
    ).reshape(-1)

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
                "sample_overdispersion_map": float(
                    map_estimates["sample_var"].flatten()[si]
                ),
                "autosomal_baseline_cn": int(autosomal_baseline_cn[si]),
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
    io = p.add_argument_group("inputs and outputs")
    io.add_argument(
        "-i",
        "--input",
        required=True,
        help="Preprocessed depth TSV (output of 'preprocess')",
    )
    io.add_argument(
        "-o",
        "--output-dir",
        required=True,
        help="Output directory",
    )
    io.add_argument(
        "--autosomal-baseline-cn-tsv",
        default=None,
        help=(
            "Optional TSV with columns 'sample' and 'autosomal_baseline_cn' "
            "that fixes each sample's neutral autosomal baseline CN. "
            "Samples default to CN=2 when this is not provided."
        ),
    )

    model_priors = p.add_argument_group("copy-number priors")
    model_priors.add_argument(
        "--autosome-prior-mode",
        choices=list(AUTOSOME_PRIOR_MODES),
        default="dirichlet",
        help="Autosomal CN prior family",
    )
    model_priors.add_argument(
        "--alpha-ref",
        type=float,
        default=50,
        help="Dirichlet concentration for CN=2 on autosomes",
    )
    model_priors.add_argument(
        "--alpha-non-ref",
        type=float,
        default=1.0,
        help="Dirichlet concentration for other CN states",
    )
    model_priors.add_argument(
        "--autosome-nonref-mean-alpha",
        type=float,
        default=1.0,
        help="Beta prior alpha for the cohort-wide autosomal non-reference mean in shrinkage mode",
    )
    model_priors.add_argument(
        "--autosome-nonref-mean-beta",
        type=float,
        default=19.0,
        help="Beta prior beta for the cohort-wide autosomal non-reference mean in shrinkage mode",
    )
    model_priors.add_argument(
        "--autosome-nonref-concentration",
        type=float,
        default=20.0,
        help="Per-bin shrinkage concentration toward the learned autosomal non-reference mean in shrinkage mode",
    )
    model_priors.add_argument(
        "--alpha-sex-ref",
        type=float,
        default=1.0,
        help=(
            "Dirichlet concentration for CN=2 on sex chromosomes "
            "(flat by default; sex-CN coupling handles sex-dependent CN)"
        ),
    )
    model_priors.add_argument(
        "--alpha-sex-non-ref",
        type=float,
        default=1.0,
        help="Dirichlet concentration for other CN states on sex chromosomes",
    )
    model_priors.add_argument(
        "--sex-prior",
        type=float,
        nargs=2,
        default=[0.5, 0.5],
        metavar=("P_XX", "P_XY"),
        help="Prior probabilities for XX and XY karyotypes",
    )
    model_priors.add_argument(
        "--sex-cn-weight",
        type=float,
        default=3.0,
        help="Weight of the sex-CN coupling factor (0 to disable)",
    )
    depth_model = p.add_argument_group("depth model")
    depth_model.add_argument(
        "--epsilon-mean",
        type=float,
        default=DEFAULT_EPSILON_MEAN,
        help=(
            "Mean of the Gamma prior on per-bin additive background depth; "
            "increase to absorb low-level zero-copy background depth, set to 0 to disable"
        ),
    )
    depth_model.add_argument(
        "--epsilon-concentration",
        type=float,
        default=DEFAULT_EPSILON_CONCENTRATION,
        help=(
            "Gamma concentration for per-bin additive background depth. "
            "Values below 1 keep most entries near zero while retaining a heavier "
            "positive tail; 1.0 matches the old Exponential prior."
        ),
    )
    depth_model.add_argument(
        "--var-sample",
        type=float,
        default=0.01,
        help="Exponential mean for per-sample variance",
    )
    depth_model.add_argument(
        "--raw-variance-power",
        type=float,
        default=DEFAULT_RAW_VARIANCE_POWER,
        help=(
            "Power-law exponent for raw-count extra-Poisson variance in the "
            "negative-binomial model. 2.0 is standard NB2 variance; lower "
            "values make residual scale grow sub-linearly with raw depth."
        ),
    )
    depth_model.add_argument(
        "--sample-depth-max",
        type=float,
        default=10000.0,
        help="Upper clip used when deriving the empirical-Bayes LogNormal prior center for per-sample depth scale",
    )
    allele_fraction = p.add_argument_group("allele fraction")
    allele_fraction.add_argument(
        "--site-data",
        default=None,
        help="Per-site allele data .npz (output of 'preprocess')",
    )
    allele_fraction.add_argument(
        "--site-af-estimator",
        choices=list(SITE_AF_ESTIMATORS),
        default="auto",
        help=(
            "How to derive the site_pop_af values used during infer. "
            "'auto' is conservative and keeps the input AFs for current "
            "preprocess outputs, 'naive-bayes' explicitly replaces them when "
            "the site encoding is coherent, and 'off' always keeps the input values."
        ),
    )
    allele_fraction.add_argument(
        "--site-af-prior-alpha",
        type=float,
        default=1.0,
        help="Beta prior alpha for naive-Bayes site AF estimation",
    )
    allele_fraction.add_argument(
        "--site-af-prior-beta",
        type=float,
        default=1.0,
        help="Beta prior beta for naive-Bayes site AF estimation",
    )
    allele_fraction.add_argument(
        "--af-concentration",
        type=float,
        default=50.0,
        help="BetaBinomial concentration for allele fraction model",
    )
    allele_fraction.add_argument(
        "--af-weight",
        type=float,
        default=DEFAULT_AF_WEIGHT,
        help=(
            "Prior median for the learned global genotype-informative allele-fraction "
            "mixture probability (0 to disable)"
        ),
    )
    allele_fraction.add_argument(
        "--af-outlier-weight",
        type=float,
        default=DEFAULT_AF_OUTLIER_WEIGHT,
        help=(
            "Fixed mixture weight for a uniform Beta-Binomial AF outlier component "
            "used to absorb contamination or mismatched sites"
        ),
    )
    allele_fraction.add_argument(
        "--af-background-concentration",
        type=float,
        default=None,
        help=(
            "Beta-Binomial concentration for the baseline-copy-number allele-fraction "
            "background. Defaults to empirical-Bayes estimation from the run."
        ),
    )
    allele_fraction.add_argument(
        "--af-temperature-prior-scale",
        type=float,
        default=0.5,
        help="LogitNormal prior scale when AF temperature learning is enabled (default)",
    )
    allele_fraction.add_argument(
        "--min-het-alt",
        type=int,
        default=3,
        help="Minimum alt-allele read count to classify a site as heterozygous in summary statistics",
    )
    allele_fraction.add_argument(
        "--min-het-af",
        type=float,
        default=0.2,
        help="Minimum allele fraction to classify a site as heterozygous in summary statistics",
    )
    allele_fraction.add_argument(
        "--max-het-af",
        type=float,
        default=0.8,
        help="Maximum allele fraction to classify a site as heterozygous in summary statistics",
    )

    guide = p.add_argument_group("variational guide")
    guide.add_argument(
        "--svi-init-restarts",
        type=int,
        default=10,
        help=(
            "Number of candidate autoguide initializations to score before "
            "gradient descent. The first candidate uses the anchored default "
            "init and the rest sample random prior starts."
        ),
    )

    training = p.add_argument_group("training")
    training.add_argument("--max-iter", type=int, default=5000)
    training.add_argument("--lr-init", type=float, default=0.02)
    training.add_argument("--lr-min", type=float, default=0.01)
    training.add_argument("--lr-decay", type=float, default=500)
    training.add_argument(
        "--grad-clip-norm",
        type=float,
        default=10.0,
        help=(
            "Clip the global SVI gradient norm to this value during training; "
            "set to 0 or a negative value to disable clipping"
        ),
    )
    training.add_argument("--log-freq", type=int, default=50)
    training.add_argument("--jit", action="store_true", default=False)
    training.add_argument("--early-stopping", action="store_true", default=True)
    training.add_argument(
        "--no-early-stopping",
        dest="early_stopping",
        action="store_false",
    )
    training.add_argument("--patience", type=int, default=50)
    training.add_argument(
        "--elbo-window",
        type=int,
        default=50,
        help="Iterations per rolling ELBO window for early stopping",
    )
    training.add_argument(
        "--elbo-rtol",
        type=float,
        default=1e-3,
        help="Relative tolerance between successive rolling ELBO windows",
    )

    discrete = p.add_argument_group("discrete inference")
    discrete.add_argument(
        "--prob-threshold",
        type=float,
        default=0.5,
        help="Min mean CN probability for aneuploidy call",
    )
    discrete.add_argument(
        "--cn-inference-method",
        choices=["single", "median", "multi-draw"],
        default="multi-draw",
        help=(
            "How to handle continuous latent uncertainty before CN inference: "
            "'single' uses a single guide draw, "
            "'median' plugs in guide medians, and 'multi-draw' averages CN "
            "posteriors over repeated guide draws."
        ),
    )
    discrete.add_argument(
        "--cn-inference-draws",
        type=int,
        default=100,
        help="Number of guide draws to average when --cn-inference-method=multi-draw",
    )

    runtime = p.add_argument_group("runtime")
    runtime.add_argument("--device", choices=["cpu", "cuda"], default="cpu")

    return p.parse_args()


def _build_infer_model(
    args: argparse.Namespace,
    *,
    device: str,
    tensor_dtype: torch.dtype,
    autosomal_baseline_cn: np.ndarray,
    af_concentration: float | np.ndarray,
    af_weight: float,
    af_outlier_weight: float,
    af_background_concentration: float,
    learn_af_temperature: bool,
) -> CNVModel:
    """Construct a CNVModel for infer with explicit AF overrides."""
    return CNVModel(
        n_states=6,
        autosome_prior_mode=args.autosome_prior_mode,
        alpha_ref=args.alpha_ref,
        alpha_non_ref=args.alpha_non_ref,
        autosome_nonref_mean_alpha=args.autosome_nonref_mean_alpha,
        autosome_nonref_mean_beta=args.autosome_nonref_mean_beta,
        autosome_nonref_concentration=args.autosome_nonref_concentration,
        var_sample=args.var_sample,
        raw_variance_power=args.raw_variance_power,
        epsilon_mean=args.epsilon_mean,
        epsilon_concentration=args.epsilon_concentration,
        device=device,
        dtype=tensor_dtype,
        af_concentration=af_concentration,
        af_weight=af_weight,
        af_outlier_weight=af_outlier_weight,
        af_background_concentration=af_background_concentration,
        learn_af_temperature=learn_af_temperature,
        af_temperature_prior_scale=args.af_temperature_prior_scale,
        alpha_sex_ref=args.alpha_sex_ref,
        alpha_sex_non_ref=args.alpha_sex_non_ref,
        sex_prior=tuple(args.sex_prior),
        sex_cn_weight=args.sex_cn_weight,
        sample_depth_max=args.sample_depth_max,
        autosomal_baseline_cn=autosomal_baseline_cn,
    )


def _train_infer_model(
    model: CNVModel,
    data: DepthData,
    args: argparse.Namespace,
) -> None:
    """Train a CNVModel using the infer CLI training controls."""
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


def _run_infer(args: argparse.Namespace, logger) -> None:
    """Run model fitting and copy-number inference after logging is configured."""

    # ── reproducibility ─────────────────────────────────────────────────
    pyro.enable_validation(True)
    pyro.distributions.enable_validation(True)
    pyro.set_rng_seed(42)
    torch.manual_seed(42)
    np.random.seed(42)

    # ── load data ───────────────────────────────────────────────────────
    logger.info("Loading preprocessed depth matrix")
    df = pd.read_csv(args.input, sep="\t", index_col=0)
    depth_space = "raw"
    logger.info(
        "Loaded depth matrix for infer: bins=%d samples=%d",
        len(df),
        len(get_sample_columns(df)),
    )

    # ── optional per-site allele data ─────────────────────────────────
    sd = None
    if args.site_data:
        sd = load_site_data(args.site_data)
        logger.info("Loaded per-site allele tensors for infer")

    df, sd = _filter_inputs_by_baseline_manifest(
        df,
        sd,
        args.autosomal_baseline_cn_tsv,
    )

    tensor_dtype = _inference_tensor_dtype()

    device = args.device
    data = DepthData(
        df, device=device, dtype=tensor_dtype,
        clamp_threshold=None,
        depth_space=depth_space,
        site_data=sd,
    )
    logger.info(
        "Prepared inference tensors: bins=%d samples=%d site_data=%s",
        data.n_bins,
        data.n_samples,
        bool(data.site_alt is not None),
    )
    autosomal_baseline_cn = _load_autosomal_baseline_cn(
        args.autosomal_baseline_cn_tsv,
        data.sample_ids,
    )
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

    af_enabled = data.site_alt is not None and args.af_weight > 0
    learn_af_temperature = af_enabled
    logger.info(
        "Allele-fraction evidence: enabled=%s learn_temperature=%s",
        af_enabled,
        learn_af_temperature,
    )
    if args.cn_inference_draws < 1:
        raise ValueError("--cn-inference-draws must be at least 1.")
    model_af_background_concentration = DEFAULT_AF_BACKGROUND_CONCENTRATION
    if af_enabled:
        if args.af_background_concentration is None:
            site_pop_for_background, used_leave_one_out = _resolve_fixed_site_pop_af_numpy(
                data.site_alt.detach().cpu().numpy(),
                data.site_total.detach().cpu().numpy(),
                effective_site_pop_af_np,
                data.site_mask.detach().cpu().numpy(),
            )
            model_af_background_concentration = estimate_af_background_concentration_numpy(
                data.site_alt.detach().cpu().numpy(),
                data.site_total.detach().cpu().numpy(),
                site_pop_for_background,
                data.site_mask.detach().cpu().numpy(),
                background_baseline_cn=autosomal_baseline_cn,
                n_states=6,
            )
        else:
            model_af_background_concentration = float(args.af_background_concentration)
            if not np.isfinite(model_af_background_concentration):
                raise ValueError("--af-background-concentration must be finite.")
            if model_af_background_concentration <= 0.0:
                raise ValueError("--af-background-concentration must be positive.")

    point_estimate_method = (
        "current" if args.cn_inference_method == "single" else args.cn_inference_method
    )
    if args.cn_inference_method == "multi-draw":
        point_estimate_method = "median"

    model_af_concentration: float | np.ndarray = args.af_concentration

    # ── build & train model ─────────────────────────────────────────────
    model = _build_infer_model(
        args,
        device=device,
        tensor_dtype=tensor_dtype,
        autosomal_baseline_cn=autosomal_baseline_cn,
        af_concentration=model_af_concentration,
        af_weight=args.af_weight if data.site_alt is not None else 0.0,
        af_outlier_weight=args.af_outlier_weight,
        af_background_concentration=model_af_background_concentration,
        learn_af_temperature=learn_af_temperature,
    )
    logger.info("Training copy-number model")
    _train_infer_model(model, data, args)

    # ── save training loss ──────────────────────────────────────────────
    loss_df = pd.DataFrame(model.loss_history)
    loss_path = os.path.join(args.output_dir, "training_loss.tsv")
    loss_df.to_csv(loss_path, sep="\t", index=False)
    output_artifacts = [loss_path]
    logger.info("Wrote training loss history: epochs=%d", len(loss_df))

    # ── point estimates + exact discrete inference ──────────────────────
    logger.info("Computing point estimates and discrete copy-number posteriors")
    map_est = model.get_map_estimates(
        data,
        estimate_method=point_estimate_method,
    )
    map_est["autosomal_baseline_cn"] = autosomal_baseline_cn.astype(np.int64, copy=False)
    # ── precompute AF table once for inference and output stats ─────────
    af_table_np = None
    if data.site_alt is not None and model.af_weight > 0:
        af_informative_weight = model._af_scale_numpy(map_est)
        with torch.no_grad():
            af_table_np = model._prepare_af_table_torch(
                data.site_alt, data.site_total,
                data.site_pop_af, data.site_mask,
                data.chr_type if hasattr(data, "chr_type") else None,
                informative_weight=torch.tensor(
                    af_informative_weight,
                    dtype=tensor_dtype,
                    device=device,
                ),
            ).cpu().numpy()

    if args.cn_inference_method == "multi-draw":
        cn_post = model.run_discrete_inference_multi_draw(
            data,
            n_draws=args.cn_inference_draws,
            af_table=af_table_np,
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
        af_concentration=model_af_concentration,
        af_outlier_weight=model.af_outlier_weight,
        af_background_concentration=model_af_background_concentration,
        autosomal_baseline_cn=autosomal_baseline_cn,
        af_informative_weight=float(np.asarray(map_est.get("af_temperature", 1.0))),
    )

    safe_diagnostics_path = os.path.join(args.output_dir, "safe_inference_diagnostics.txt")
    with open(safe_diagnostics_path, "w", encoding="utf-8") as handle:
        for message in safe_diagnostic_messages:
            handle.write(message)
            handle.write("\n")
    output_artifacts.append(safe_diagnostics_path)
    logger.info(
        "Wrote privacy-safe inference diagnostics: n_messages=%d",
        len(safe_diagnostic_messages),
    )

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
        output_artifacts.append(site_af_path)
        logger.info(
            "Wrote site allele-frequency estimates: rows=%d",
            len(site_af_estimates_df),
        )

    # ── save inference artifacts for downstream tools (ppd, plot) ──────
    artifacts_path = os.path.join(args.output_dir, "inference_artifacts.npz")
    artifact_dict = {k: v for k, v in map_est.items()}
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
    artifact_dict["model_var_sample"] = np.asarray(model.var_sample)
    artifact_dict["model_raw_variance_power"] = np.asarray(model.raw_variance_power)
    artifact_dict["model_af_concentration"] = np.asarray(model.af_concentration)
    artifact_dict["model_af_weight"] = np.asarray(model.af_weight)
    artifact_dict["model_af_outlier_weight"] = np.asarray(model.af_outlier_weight)
    artifact_dict["model_af_background_concentration"] = np.asarray(
        model.af_background_concentration
    )
    artifact_dict["model_af_background_mode"] = np.asarray("baseline-copy-number")
    artifact_dict["model_learn_af_temperature"] = np.asarray(model.learn_af_temperature)
    artifact_dict["model_af_temperature_prior_scale"] = np.asarray(
        model.af_temperature_prior_scale
    )
    artifact_dict["model_alpha_sex_ref"] = np.asarray(model.alpha_sex_ref)
    artifact_dict["model_alpha_sex_non_ref"] = np.asarray(model.alpha_sex_non_ref)
    artifact_dict["model_sex_prior"] = np.asarray(model.sex_prior)
    artifact_dict["model_sex_cn_weight"] = np.asarray(model.sex_cn_weight)
    artifact_dict["epsilon_mean"] = np.asarray(model.epsilon_mean)
    artifact_dict["model_epsilon_concentration"] = np.asarray(
        model.epsilon_concentration
    )
    artifact_dict["sample_depth_max"] = np.asarray(model.sample_depth_max)
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
        artifact_dict["site_af_prior_alpha"] = np.asarray(args.site_af_prior_alpha)
        artifact_dict["site_af_prior_beta"] = np.asarray(args.site_af_prior_beta)
    for k, v in cn_post.items():
        artifact_dict[f"cn_post_{k}"] = v
    np.savez_compressed(artifacts_path, **artifact_dict)
    output_artifacts.append(artifacts_path)

    baseline_df = pd.DataFrame(
        {
            "sample": data.sample_ids,
            "autosomal_baseline_cn": autosomal_baseline_cn.astype(np.int64, copy=False),
        }
    )
    baseline_path = os.path.join(args.output_dir, "sample_autosomal_baseline_cn.tsv")
    baseline_df.to_csv(baseline_path, sep="\t", index=False)
    output_artifacts.append(baseline_path)

    # ── detect aneuploidies ─────────────────────────────────────────────
    aneuploid_map = detect_aneuploidies(
        data,
        cn_post,
        prob_threshold=args.prob_threshold,
        autosomal_baseline_cn=autosomal_baseline_cn,
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
    output_artifacts.append(bin_path)

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
    output_artifacts.append(chr_path)
    logger.info(
        "Wrote inference summaries: bin_rows=%d chromosome_rows=%d",
        len(bin_df),
        len(chr_df),
    )
    log_output_artifacts(logger, output_artifacts)


def main() -> None:
    """Entry point for ``gatk-sv-ploidy infer``."""
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    with tool_logging_context(
        tool_name="infer",
        output_dir=args.output_dir,
        args=args,
        random_seeds={"numpy": 42, "pyro": 42, "torch": 42},
    ) as logger:
        _run_infer(args, logger)


if __name__ == "__main__":
    main()
