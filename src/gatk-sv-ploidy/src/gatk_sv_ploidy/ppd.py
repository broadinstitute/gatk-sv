r"""
Posterior predictive check (ppd) subcommand.

Generates posterior predictive draws from the trained model and computes
summary statistics comparing the predictive distribution to observed data.
This is a model-checking diagnostic: if the model fits well, the observed
data should look like a typical sample from the posterior predictive
distribution.

For each bin × sample pair the predictive distribution is a finite mixture
over the copy-number posterior.

For the continuous observation families, the predictive depth is a mixture of
location-scale residual distributions:

.. math::

    p(\tilde{d} | x) = \sum_{c=0}^{5} p(c | x)
        \cdot \mathcal{N}(\tilde{d} \mid c \cdot b_i + \epsilon_{ij},
                           \sigma^2_{ij})

where *b_{ij}* is the MAP effective multiplicative bias surface, *\epsilon_{ij}* is the MAP additive
background depth, and *p(c | x)* is the analytical CN posterior.

For ``negative_binomial`` fits on raw counts, the predictive mean becomes

.. math::

    \mu_{ij}(c) = \ell_i \cdot s_j \cdot (c \cdot b_{ij} + I(c=0)\epsilon_{ij}) / 2

with *\ell_i* the bin length in kilobases, *s_j* the MAP sample-depth
latent. The predictive variance follows

.. math::

    \mu_{ij}(c) + (\alpha_i + \beta_j + \gamma_{t(i)}) \mu_{ij}(c)^p

where \gamma_{t(i)} is the learned chrX/chrY excess-overdispersion term
and is anchored to zero on autosomes, and *p* is the raw-count variance
power saved by ``infer``.
"""

from __future__ import annotations

import argparse
import os
from typing import Dict

import numpy as np
import pandas as pd
import torch
from scipy import stats as sp_stats

from gatk_sv_ploidy._logging import log_output_artifacts, tool_logging_context
from gatk_sv_ploidy._util import (
    compose_additive_background_matrix,
    get_sample_columns,
)
from gatk_sv_ploidy.data import DepthData, load_site_data
from gatk_sv_ploidy.infer import (
    _filter_inputs_by_baseline_manifest,
    load_inference_artifacts,
)
from gatk_sv_ploidy.models import (
    CNVModel,
    _compose_effective_bin_bias_numpy,
    _effective_negative_binomial_overdispersion_numpy,
    _raw_expected_depth_units,
)


def apply_effective_site_pop_af(
    site_data: Dict[str, np.ndarray] | None,
    map_estimates: Dict[str, np.ndarray],
) -> Dict[str, np.ndarray] | None:
    """Replace site-data AFs with the effective AFs persisted by infer."""
    if site_data is None or "site_pop_af_effective" not in map_estimates:
        return site_data

    effective_site_pop_af = np.asarray(map_estimates["site_pop_af_effective"])
    if effective_site_pop_af.shape != site_data["site_pop_af"].shape:
        raise ValueError(
            "Saved effective site_pop_af shape does not match the supplied "
            "site_data.npz. Regenerate both artifacts from the same infer run."
        )

    updated = dict(site_data)
    updated["site_pop_af"] = effective_site_pop_af.astype(
        site_data["site_pop_af"].dtype,
        copy=False,
    )
    return updated


def _align_ppd_inputs_to_sample_ids(
    df: pd.DataFrame,
    site_data: Dict[str, np.ndarray] | None,
    fitted_sample_ids: np.ndarray,
) -> tuple[pd.DataFrame, Dict[str, np.ndarray] | None]:
    """Subset and reorder PPD inputs to the fitted sample order."""
    sample_cols = get_sample_columns(df)
    desired_sample_ids = [
        str(sample_id) for sample_id in np.asarray(fitted_sample_ids).tolist()
    ]
    missing_depth_samples = [
        sample_id for sample_id in desired_sample_ids if sample_id not in sample_cols
    ]
    if missing_depth_samples:
        preview = ", ".join(missing_depth_samples[:5])
        more = "" if len(missing_depth_samples) <= 5 else f" (+{len(missing_depth_samples) - 5} more)"
        raise ValueError(
            "Inference artifacts refer to samples missing from the PPD depth input: "
            f"{preview}{more}"
        )

    metadata_cols = [column for column in df.columns if column not in sample_cols]
    aligned_df = df.loc[:, metadata_cols + desired_sample_ids].copy()

    if site_data is None:
        return aligned_df, None

    if "sample_ids" in site_data:
        site_sample_ids = [
            str(sample_id) for sample_id in np.asarray(site_data["sample_ids"]).tolist()
        ]
    else:
        site_sample_ids = sample_cols
        if int(site_data["site_alt"].shape[2]) != len(sample_cols):
            raise ValueError(
                "site_data is missing sample_ids metadata and its sample axis does not "
                "match the PPD depth input sample count. Regenerate site_data.npz."
            )

    missing_site_samples = [
        sample_id for sample_id in desired_sample_ids if sample_id not in site_sample_ids
    ]
    if missing_site_samples:
        preview = ", ".join(missing_site_samples[:5])
        more = "" if len(missing_site_samples) <= 5 else f" (+{len(missing_site_samples) - 5} more)"
        raise ValueError(
            "Inference artifacts refer to samples missing from the PPD site-data input: "
            f"{preview}{more}"
        )

    selected_indices = [site_sample_ids.index(sample_id) for sample_id in desired_sample_ids]
    aligned_site_data = dict(site_data)
    for key in ("site_alt", "site_total", "site_mask"):
        if key in aligned_site_data:
            aligned_site_data[key] = aligned_site_data[key][:, :, selected_indices]
    aligned_site_data["sample_ids"] = np.asarray(desired_sample_ids, dtype=object)
    return aligned_df, aligned_site_data


def _resolve_ppd_baseline_manifest(
    artifacts_path: str,
    explicit_path: str | None,
) -> str | None:
    """Find the baseline manifest used to exclude samples before infer."""
    if explicit_path is not None:
        return explicit_path

    sibling_path = os.path.join(
        os.path.dirname(os.path.abspath(artifacts_path)),
        "sample_autosomal_baseline_cn.tsv",
    )
    if os.path.exists(sibling_path):
        return sibling_path
    return None


def _align_ppd_inputs_from_baseline_manifest(
    df: pd.DataFrame,
    site_data: Dict[str, np.ndarray] | None,
    baseline_path: str,
    expected_n_samples: int,
) -> tuple[pd.DataFrame, Dict[str, np.ndarray] | None]:
    """Align PPD inputs using either full-manifest or subset-manifest TSVs."""
    baseline_df = pd.read_csv(baseline_path, sep="\t")
    if "sample" not in baseline_df.columns:
        raise ValueError(
            "Autosomal baseline CN TSV is missing required columns: sample"
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
    manifest_sample_ids = [
        str(sample_id) for sample_id in baseline_df["sample"].astype(str).tolist()
    ]
    manifest_sample_set = set(manifest_sample_ids)
    sample_col_set = set(sample_cols)

    if sample_col_set.issubset(manifest_sample_set):
        filtered_df, filtered_site_data = _filter_inputs_by_baseline_manifest(
            df,
            site_data,
            baseline_path,
        )
        return filtered_df, filtered_site_data

    missing_from_depth = [
        sample_id for sample_id in manifest_sample_ids if sample_id not in sample_col_set
    ]
    if missing_from_depth:
        preview = ", ".join(missing_from_depth[:5])
        more = "" if len(missing_from_depth) <= 5 else f" (+{len(missing_from_depth) - 5} more)"
        raise ValueError(
            "Autosomal baseline CN TSV contains samples missing from the PPD depth input: "
            f"{preview}{more}"
        )

    if len(manifest_sample_ids) != expected_n_samples:
        raise ValueError(
            "Auto-discovered autosomal baseline CN TSV appears to contain only the "
            "infer-retained sample subset, but its sample count does not match the "
            f"inference artifacts. Manifest has {len(manifest_sample_ids)} samples, "
            f"while cn_posterior expects {expected_n_samples}."
        )
    return _align_ppd_inputs_to_sample_ids(df, site_data, np.asarray(manifest_sample_ids, dtype=object))


def _prepare_ppd_inputs(
    df: pd.DataFrame,
    site_data: Dict[str, np.ndarray] | None,
    map_estimates: Dict[str, np.ndarray],
    cn_posterior: Dict[str, np.ndarray],
    artifacts_path: str,
    baseline_path: str | None,
) -> tuple[pd.DataFrame, Dict[str, np.ndarray] | None]:
    """Align PPD inputs to the sample subset used by infer."""
    cn_probs = np.asarray(cn_posterior["cn_posterior"])
    if cn_probs.ndim != 3:
        raise ValueError("cn_posterior must have shape (n_bins, n_samples, n_states).")

    fitted_sample_ids = map_estimates.get("sample_ids")
    if fitted_sample_ids is not None:
        fitted_sample_ids = np.asarray(fitted_sample_ids)
        if int(fitted_sample_ids.shape[0]) != int(cn_probs.shape[1]):
            raise ValueError(
                "Inference artifacts contain sample_ids with length "
                f"{int(fitted_sample_ids.shape[0])}, but cn_posterior expects "
                f"{int(cn_probs.shape[1])} samples."
            )
        return _align_ppd_inputs_to_sample_ids(df, site_data, fitted_sample_ids)

    input_n_samples = len(get_sample_columns(df))
    expected_n_samples = int(cn_probs.shape[1])
    if input_n_samples == expected_n_samples:
        return df, site_data

    resolved_baseline_path = _resolve_ppd_baseline_manifest(artifacts_path, baseline_path)
    if resolved_baseline_path is not None:
        filtered_df, filtered_site_data = _align_ppd_inputs_from_baseline_manifest(
            df,
            site_data,
            resolved_baseline_path,
            expected_n_samples,
        )
        filtered_n_samples = len(get_sample_columns(filtered_df))
        if filtered_n_samples == expected_n_samples:
            return filtered_df, filtered_site_data
        raise ValueError(
            "PPD input sample count still does not match inference artifacts after "
            f"filtering with {resolved_baseline_path}. Depth input has {filtered_n_samples} "
            f"retained samples, but cn_posterior expects {expected_n_samples}."
        )

    raise ValueError(
        "PPD input sample count does not match inference artifacts: depth input has "
        f"{input_n_samples} samples, but cn_posterior expects {expected_n_samples}. "
        "This infer run likely excluded samples before fitting. Provide "
        "--autosomal-baseline-cn-tsv or rerun infer with artifact sample_ids support."
    )


def _phred_scale_error_probability(
    error_prob: float,
    max_quality: float = 99.0,
) -> float:
    """Convert an error probability into a Phred-scaled quality score."""
    return float(min(max_quality, -10.0 * np.log10(max(error_prob, 1e-300))))


def _build_model_from_artifacts(
    map_estimates: Dict[str, np.ndarray],
    device: str,
) -> CNVModel:
    """Reconstruct a model instance for analytical inference in PPD."""
    af_concentration = np.asarray(
        map_estimates.get("model_af_concentration", 50.0),
        dtype=np.float64,
    )
    if af_concentration.ndim == 0:
        af_concentration_value: float | np.ndarray = float(af_concentration.item())
    else:
        af_concentration_value = af_concentration.astype(np.float64, copy=False)

    sex_prior = tuple(
        np.asarray(
            map_estimates.get("model_sex_prior", np.asarray([0.5, 0.5]))
        ).astype(np.float32).tolist()
    )
    return CNVModel(
        n_states=int(np.asarray(map_estimates.get("model_n_states", 6)).item()),
        autosome_prior_mode=str(
            np.asarray(
                map_estimates.get("model_autosome_prior_mode", "dirichlet")
            ).item()
        ),
        alpha_ref=float(np.asarray(map_estimates.get("model_alpha_ref", 50.0)).item()),
        alpha_non_ref=float(
            np.asarray(map_estimates.get("model_alpha_non_ref", 1.0)).item()
        ),
        autosome_nonref_mean_alpha=float(
            np.asarray(
                map_estimates.get("model_autosome_nonref_mean_alpha", 1.0)
            ).item()
        ),
        autosome_nonref_mean_beta=float(
            np.asarray(
                map_estimates.get("model_autosome_nonref_mean_beta", 19.0)
            ).item()
        ),
        autosome_nonref_concentration=float(
            np.asarray(
                map_estimates.get("model_autosome_nonref_concentration", 20.0)
            ).item()
        ),
        var_sample=float(
            np.asarray(map_estimates.get("model_var_sample", 0.001)).item()
        ),
        raw_variance_power=float(
            np.asarray(map_estimates.get("model_raw_variance_power", 2.0)).item()
        ),
        epsilon_mean=float(np.asarray(map_estimates.get("epsilon_mean", 0.0)).item()),
        epsilon_concentration=float(
            np.asarray(map_estimates.get("model_epsilon_concentration", 1.0)).item()
        ),
        device=device,
        dtype=torch.float64,
        af_concentration=af_concentration_value,
        af_weight=float(np.asarray(map_estimates.get("model_af_weight", 0.0)).item()),
        af_outlier_weight=float(
            np.asarray(map_estimates.get("model_af_outlier_weight", 0.05)).item()
        ),
        af_background_concentration=float(
            np.asarray(
                map_estimates.get("model_af_background_concentration", 0.5)
            ).item()
        ),
        learn_af_temperature=bool(
            np.asarray(map_estimates.get("model_learn_af_temperature", False)).item()
        ),
        af_temperature_prior_scale=float(
            np.asarray(
                map_estimates.get("model_af_temperature_prior_scale", 0.5)
            ).item()
        ),
        alpha_sex_ref=float(
            np.asarray(map_estimates.get("model_alpha_sex_ref", 1.0)).item()
        ),
        alpha_sex_non_ref=float(
            np.asarray(map_estimates.get("model_alpha_sex_non_ref", 1.0)).item()
        ),
        sex_prior=sex_prior,
        sex_cn_weight=float(
            np.asarray(map_estimates.get("model_sex_cn_weight", 3.0)).item()
        ),
        sample_depth_max=float(
            np.asarray(map_estimates.get("sample_depth_max", 10000.0)).item()
        ),
        autosomal_baseline_cn=np.asarray(
            map_estimates.get("autosomal_baseline_cn", None)
        ) if "autosomal_baseline_cn" in map_estimates else None,
    )


def _sample_ppd_observation(
    data: DepthData,
    rng: np.random.RandomState,
    discrete_posterior: Dict[str, np.ndarray],
    continuous_maps: Dict[str, np.ndarray],
) -> np.ndarray:
    """Sample one posterior predictive observation matrix."""
    cn_probs = discrete_posterior["cn_posterior"]
    n_bins, n_samples, n_states = cn_probs.shape
    bin_bias = _compose_effective_bin_bias_numpy(
        n_bins,
        n_samples,
        fixed_bias=np.asarray(
            continuous_maps.get("bin_bias_matrix", continuous_maps["bin_bias"])
        ),
    )
    sample_var = np.atleast_1d(np.asarray(continuous_maps["sample_var"]).squeeze())
    if "bin_var" in continuous_maps:
        bin_var = np.atleast_1d(np.asarray(continuous_maps["bin_var"]).squeeze())
    else:
        bin_var = np.zeros(n_bins, dtype=np.float64)
    additive_background = compose_additive_background_matrix(
        continuous_maps.get("bin_epsilon"),
        n_bins,
        n_samples,
        dtype=np.float64,
    )

    flat_probs = cn_probs.reshape(-1, n_states)
    cum = np.cumsum(flat_probs, axis=1)
    u = rng.uniform(size=(flat_probs.shape[0], 1))
    cn_sample = (u < cum).argmax(axis=1).reshape(n_bins, n_samples)

    sample_depth = np.atleast_1d(
        np.asarray(continuous_maps["sample_depth"]).squeeze()
    )
    overdispersion = bin_var[:, np.newaxis] + sample_var[np.newaxis, :]
    mean = data.bin_length_kb.detach().cpu().numpy()[:, np.newaxis]
    mean = mean * sample_depth[np.newaxis, :]
    autosomal_baseline_cn = continuous_maps.get("autosomal_baseline_cn", 2.0)
    mean = mean * np.maximum(
        _raw_expected_depth_units(
            cn_sample.astype(np.float32),
            bin_bias,
            additive_background,
            autosomal_baseline_copy_number=autosomal_baseline_cn,
        ),
        0.0,
    )
    raw_variance_power = float(
        np.asarray(
            continuous_maps.get("model_raw_variance_power", 2.0)
        ).item()
    )
    effective_overdispersion = _effective_negative_binomial_overdispersion_numpy(
        overdispersion,
        mean,
        raw_variance_power,
    )
    concentration = 1.0 / np.maximum(effective_overdispersion, 1e-8)
    latent_rate = rng.gamma(
        shape=concentration,
        scale=np.where(mean > 0.0, mean / concentration, 0.0),
    )
    draw = rng.poisson(latent_rate)

    return draw.astype(np.float32)


# ── posterior predictive sampling ───────────────────────────────────────────


def generate_ppd_depth(
    data: DepthData,
    map_estimates: Dict[str, np.ndarray],
    cn_posterior: Dict[str, np.ndarray],
    n_draws: int = 100,
    seed: int = 42,
) -> np.ndarray:
    """Generate posterior predictive draws of observed depth.

    For each bin *i* and sample *j* with CN posterior *p(c|x)* and
    fitted continuous parameters, draw from the mixture:

    1. Sample a CN state *c* from the posterior.
    2. Draw depth from the configured observation family.

    For the negative-binomial observation family on raw counts, the
    predictive mean is
    ``bin_length_kb_i × sample_depth_j × (c × bin_bias_{ij} + I(c=0) × additive_background_{ij}) / 2`` and
    the predictive variance is
    ``μ + (bin_var_i + sample_var_j) × μ^raw_variance_power``.

    Args:
        data: :class:`DepthData` instance.
        map_estimates: MAP estimates dict (from ``inference_artifacts.npz``).
        cn_posterior: Discrete posterior dict.
        n_draws: Number of posterior predictive draws per bin/sample.
        seed: Random seed.
    Returns:
        Array of shape ``(n_draws, n_bins, n_samples)`` with simulated depths.
    """
    rng = np.random.RandomState(seed)
    cn_probs = np.asarray(cn_posterior["cn_posterior"])
    if cn_probs.ndim != 3:
        raise ValueError("cn_posterior must have shape (n_bins, n_samples, n_states).")
    expected_shape = (int(data.n_bins), int(data.n_samples))
    observed_shape = (int(cn_probs.shape[0]), int(cn_probs.shape[1]))
    if observed_shape != expected_shape:
        raise ValueError(
            "cn_posterior shape does not match the supplied depth matrix. "
            f"Expected bins×samples {expected_shape}, got {observed_shape}. "
            "Align the ppd input depth/site data to the fitted infer sample subset."
        )
    draws = np.empty((n_draws, data.n_bins, data.n_samples), dtype=np.float32)
    for d in range(n_draws):
        draws[d] = _sample_ppd_observation(
            data,
            rng,
            cn_posterior,
            map_estimates,
        )
    return draws


def _compute_randomized_pit(
    ppd_draws: np.ndarray,
    observed: np.ndarray,
    seed: int = 0,
) -> np.ndarray:
    """Compute randomized PIT values from posterior predictive draws.

    For discrete predictive distributions, ties at the observed value create
    non-uniform one-sided tail histograms even when the model is calibrated.
    This helper randomizes within the predictive mass at the observed value so
    the resulting PIT is the appropriate uniform-target diagnostic.
    """
    observed = np.asarray(observed)
    n_draws = int(ppd_draws.shape[0])
    n_less = (ppd_draws < observed[np.newaxis, :, :]).sum(axis=0)
    n_equal = (ppd_draws == observed[np.newaxis, :, :]).sum(axis=0)
    rng = np.random.RandomState(seed)
    tie_break = rng.uniform(size=observed.shape)
    pit = (n_less + tie_break * n_equal) / max(n_draws, 1)
    return np.clip(pit, 0.0, 1.0)


# ── summary statistics ──────────────────────────────────────────────────────


def compute_ppd_bin_summary(
    data: DepthData,
    ppd_draws: np.ndarray,
    map_estimates: Dict[str, np.ndarray],
    cn_posterior: Dict[str, np.ndarray],
) -> pd.DataFrame:
    """Compute per-bin, per-sample PPD summary statistics.

    For each bin/sample, computes:
    - ``ppd_mean``, ``ppd_std``: mean and std of predictive draws.
    - ``residual``: observed − ppd_mean.
    - ``z_score``: residual / ppd_std.
    - ``randomized_pit``: randomized probability integral transform value.
    - ``tail_prob``: fraction of PPD draws ≥ observed (one-sided).
    - ``two_tail_prob``: P(|draw − ppd_mean| ≥ |obs − ppd_mean|).

    Args:
        data: :class:`DepthData` instance.
        ppd_draws: Array of shape ``(n_draws, n_bins, n_samples)``.
        map_estimates: MAP estimates dict.
        cn_posterior: CN posterior dict.

    Returns:
        DataFrame with one row per (bin, sample).
    """
    obs = data.depth.detach().cpu().numpy()  # (n_bins, n_samples)
    cn_probs = cn_posterior["cn_posterior"]

    ppd_mean = ppd_draws.mean(axis=0)
    ppd_std = ppd_draws.std(axis=0)
    ppd_std_safe = np.maximum(ppd_std, 1e-10)

    residual = obs - ppd_mean
    z_score = residual / ppd_std_safe

    randomized_pit = _compute_randomized_pit(ppd_draws, obs, seed=0)

    # Tail probability: fraction of draws ≥ observed
    n_draws = ppd_draws.shape[0]
    tail_prob = (ppd_draws >= obs[np.newaxis, :, :]).sum(axis=0) / n_draws

    # Two-sided tail probability
    abs_dev_obs = np.abs(obs - ppd_mean)
    abs_dev_draws = np.abs(ppd_draws - ppd_mean[np.newaxis, :, :])
    two_tail_prob = (abs_dev_draws >= abs_dev_obs[np.newaxis, :, :]).sum(axis=0) / n_draws

    # Quantiles of PPD
    ppd_q05 = np.percentile(ppd_draws, 5, axis=0)
    ppd_q25 = np.percentile(ppd_draws, 25, axis=0)
    ppd_q50 = np.percentile(ppd_draws, 50, axis=0)
    ppd_q75 = np.percentile(ppd_draws, 75, axis=0)
    ppd_q95 = np.percentile(ppd_draws, 95, axis=0)
    outside_90pct_interval = np.logical_or(obs < ppd_q05, obs > ppd_q95)
    outside_50pct_interval = np.logical_or(obs < ppd_q25, obs > ppd_q75)

    rows: list[dict] = []
    for i in range(data.n_bins):
        for j in range(data.n_samples):
            cn_map = int(np.argmax(cn_probs[i, j, :]))
            rows.append({
                "chr": data.chr[i],
                "start": int(data.start[i]),
                "end": int(data.end[i]),
                "sample": data.sample_ids[j],
                "observed_depth": float(obs[i, j]),
                "cn_map": cn_map,
                "ppd_mean": float(ppd_mean[i, j]),
                "ppd_std": float(ppd_std[i, j]),
                "ppd_q05": float(ppd_q05[i, j]),
                "ppd_q25": float(ppd_q25[i, j]),
                "ppd_q50": float(ppd_q50[i, j]),
                "ppd_q75": float(ppd_q75[i, j]),
                "ppd_q95": float(ppd_q95[i, j]),
                "residual": float(residual[i, j]),
                "z_score": float(z_score[i, j]),
                "randomized_pit": float(randomized_pit[i, j]),
                "tail_prob": float(tail_prob[i, j]),
                "two_tail_prob": float(two_tail_prob[i, j]),
                "outside_90pct_interval": bool(outside_90pct_interval[i, j]),
                "outside_50pct_interval": bool(outside_50pct_interval[i, j]),
            })

    return pd.DataFrame(rows)


def compute_ppd_bin_quality_summary(
    ppd_bin_df: pd.DataFrame,
    outside_90_thresholds: tuple[float, float] = (0.15, 0.20),
) -> pd.DataFrame:
    """Aggregate bin-level PPD fit into per-bin quality metrics.

    For each bin, estimate the posterior probability that the true fraction of
    samples outside the central 90% posterior predictive interval exceeds a
    user-chosen tolerance threshold. BINQ scores are Phred-scaled versions of
    those bad-bin probabilities, so higher is better.
    """
    thresholds = tuple(float(t) for t in outside_90_thresholds)
    if any(t <= 0.0 or t >= 1.0 for t in thresholds):
        raise ValueError("outside_90_thresholds must lie strictly between 0 and 1.")

    rows: list[dict] = []
    group_cols = ["chr", "start", "end"]
    for (chrom, start, end), sdf in ppd_bin_df.groupby(group_cols, sort=False):
        n_samples = int(len(sdf))
        n_outside_90 = int(sdf["outside_90pct_interval"].sum())
        n_outside_50 = int(sdf["outside_50pct_interval"].sum())
        frac_outside_90 = n_outside_90 / max(n_samples, 1)
        frac_outside_50 = n_outside_50 / max(n_samples, 1)
        alpha = 1.0 + n_outside_90
        beta = 1.0 + (n_samples - n_outside_90)

        row = {
            "chr": chrom,
            "start": int(start),
            "end": int(end),
            "n_samples": n_samples,
            "n_outside_90pct_interval": n_outside_90,
            "frac_outside_90pct_interval": float(frac_outside_90),
            "n_outside_50pct_interval": n_outside_50,
            "frac_outside_50pct_interval": float(frac_outside_50),
            "median_two_tail_prob": float(sdf["two_tail_prob"].median()),
            "mean_two_tail_prob": float(sdf["two_tail_prob"].mean()),
        }

        for threshold in thresholds:
            label = int(round(100 * threshold))
            bad_prob = float(1.0 - sp_stats.beta.cdf(threshold, alpha, beta))
            row[f"binq{label}_bad_prob"] = bad_prob
            row[f"BINQ{label}"] = _phred_scale_error_probability(bad_prob)

        rows.append(row)

    return pd.DataFrame(rows).sort_values(group_cols).reset_index(drop=True)


def compute_call_stability_quality_summary(
    data: DepthData,
    cn_posterior: Dict[str, np.ndarray],
    instability_thresholds: tuple[float, float] = (0.15, 0.20),
) -> pd.DataFrame:
    """Aggregate per-bin call quality from posterior confidence and stability.

    ``BINQ`` measures posterior-predictive model fit. That is useful for model
    diagnostics, but it is not the right default for deciding whether a copy-
    number call is reliable. ``CALLQ`` instead estimates the probability that a
    bin has too large a per-sample call-error burden, where per-sample error is
    the larger of posterior CN uncertainty and multi-draw MAP instability.
    """
    thresholds = tuple(float(t) for t in instability_thresholds)
    if any(t <= 0.0 or t >= 1.0 for t in thresholds):
        raise ValueError("instability_thresholds must lie strictly between 0 and 1.")

    cn_probs = np.asarray(cn_posterior["cn_posterior"], dtype=np.float64)
    if cn_probs.shape[:2] != (data.n_bins, data.n_samples):
        raise ValueError("cn_posterior must have shape (n_bins, n_samples, n_states).")
    posterior_confidence = np.clip(cn_probs.max(axis=2), 0.0, 1.0)
    posterior_error = 1.0 - posterior_confidence

    stability = cn_posterior.get("cn_map_stability")
    if stability is None:
        stability = np.ones((data.n_bins, data.n_samples), dtype=np.float64)
    else:
        stability = np.asarray(stability, dtype=np.float64)
    if stability.shape != (data.n_bins, data.n_samples):
        raise ValueError("cn_map_stability must have shape (n_bins, n_samples).")

    rows: list[dict] = []
    for i in range(data.n_bins):
        bin_stability = np.clip(stability[i, :], 0.0, 1.0)
        bin_instability = 1.0 - bin_stability
        bin_posterior_error = np.clip(posterior_error[i, :], 0.0, 1.0)
        bin_call_error = np.maximum(bin_instability, bin_posterior_error)
        alpha = 1.0 + float(bin_call_error.sum())
        beta = 1.0 + float((1.0 - bin_call_error).sum())
        row = {
            "chr": data.chr[i],
            "start": int(data.start[i]),
            "end": int(data.end[i]),
            "mean_call_stability": float(bin_stability.mean()),
            "median_call_stability": float(np.median(bin_stability)),
            "mean_call_instability": float(bin_instability.mean()),
            "mean_posterior_call_error": float(bin_posterior_error.mean()),
            "mean_call_error": float(bin_call_error.mean()),
        }
        for threshold in thresholds:
            label = int(round(100 * threshold))
            bad_prob = float(1.0 - sp_stats.beta.cdf(threshold, alpha, beta))
            row[f"callq{label}_bad_prob"] = bad_prob
            row[f"CALLQ{label}"] = _phred_scale_error_probability(bad_prob)
        rows.append(row)

    return pd.DataFrame(rows).sort_values(["chr", "start", "end"]).reset_index(drop=True)


def compute_ppd_chromosome_summary(
    data: DepthData,
    ppd_draws: np.ndarray,
    cn_posterior: Dict[str, np.ndarray],
) -> pd.DataFrame:
    """Compute per-chromosome, per-sample PPD summary statistics.

    Aggregates bin-level PPD stats to the chromosome level: mean residual,
    mean absolute residual, RMSE, mean z-score, mean tail probability,
    and Bayesian p-value (fraction of draws where sum-of-squared-residuals
    for the replicated data exceeds that for the observed data).

    Args:
        data: :class:`DepthData` instance.
        ppd_draws: Array of shape ``(n_draws, n_bins, n_samples)``.
        cn_posterior: CN posterior dict.

    Returns:
        DataFrame with one row per (sample, chromosome).
    """
    obs = data.depth.detach().cpu().numpy()
    cn_probs = cn_posterior["cn_posterior"]

    ppd_mean = ppd_draws.mean(axis=0)
    ppd_std_safe = np.maximum(ppd_draws.std(axis=0), 1e-10)

    residual = obs - ppd_mean
    z_score = residual / ppd_std_safe
    n_draws = ppd_draws.shape[0]
    tail_prob = (ppd_draws >= obs[np.newaxis, :, :]).sum(axis=0) / n_draws

    unique_chrs = np.unique(data.chr)
    rows: list[dict] = []

    for si in range(data.n_samples):
        for chr_name in unique_chrs:
            mask = data.chr == chr_name
            n_bins = int(mask.sum())
            if n_bins == 0:
                continue

            chr_obs = obs[mask, si]
            chr_res = residual[mask, si]
            chr_z = z_score[mask, si]
            chr_tail = tail_prob[mask, si]
            chr_ppd = ppd_draws[:, mask, si]  # (n_draws, n_chr_bins)

            # Bayesian p-value: compare chi-squared-like discrepancy
            # D(y, θ) = Σ (y_i − E[y_i|θ])² / Var[y_i|θ]
            chr_ppd_mean = ppd_mean[mask, si]
            chr_ppd_std = ppd_std_safe[mask, si]

            disc_obs = np.sum(((chr_obs - chr_ppd_mean) / chr_ppd_std) ** 2)
            disc_rep = np.sum(
                ((chr_ppd - chr_ppd_mean[np.newaxis, :]) / chr_ppd_std[np.newaxis, :]) ** 2,
                axis=1,
            )
            bayesian_pval = float(np.mean(disc_rep >= disc_obs))

            # Dominant CN
            chr_cn_probs = cn_probs[mask, si, :]
            cn_map = np.argmax(chr_cn_probs, axis=-1)
            counts = np.bincount(cn_map, minlength=6)
            dominant_cn = int(np.argmax(counts))

            rows.append({
                "sample": data.sample_ids[si],
                "chromosome": chr_name,
                "dominant_cn": dominant_cn,
                "n_bins": n_bins,
                "mean_residual": float(np.mean(chr_res)),
                "mean_abs_residual": float(np.mean(np.abs(chr_res))),
                "rmse": float(np.sqrt(np.mean(chr_res ** 2))),
                "mean_z_score": float(np.mean(chr_z)),
                "std_z_score": float(np.std(chr_z)),
                "mean_tail_prob": float(np.mean(chr_tail)),
                "bayesian_pvalue": bayesian_pval,
                "mean_observed_depth": float(np.mean(chr_obs)),
                "mean_predicted_depth": float(np.mean(chr_ppd_mean)),
            })

    df = pd.DataFrame(rows).sort_values(["sample", "chromosome"])
    return df


def compute_ppd_global_summary(
    ppd_bin_df: pd.DataFrame,
) -> pd.DataFrame:
    """Compute global (whole-genome) PPD calibration statistics.

    Args:
        ppd_bin_df: Per-bin PPD summary DataFrame.

    Returns:
        Single-row DataFrame with global calibration metrics.
    """
    tail_probs = ppd_bin_df["tail_prob"].values
    z_scores = ppd_bin_df["z_score"].values

    # KS test: tail_probs should be Uniform(0, 1) if well-calibrated
    ks_stat, ks_pval = sp_stats.kstest(tail_probs, "uniform")

    # Z-scores should be approximately Normal(0, 1)
    z_mean = float(np.mean(z_scores))
    z_std = float(np.std(z_scores))

    # Fraction of bins outside 90% predictive interval
    outside_90 = float(
        np.mean(
            np.logical_or(
                ppd_bin_df["observed_depth"] < ppd_bin_df["ppd_q05"],
                ppd_bin_df["observed_depth"] > ppd_bin_df["ppd_q95"],
            )
        )
    )
    # Fraction outside 50% predictive interval
    outside_50 = float(
        np.mean(
            np.logical_or(
                ppd_bin_df["observed_depth"] < ppd_bin_df["ppd_q25"],
                ppd_bin_df["observed_depth"] > ppd_bin_df["ppd_q75"],
            )
        )
    )

    return pd.DataFrame([{
        "n_bins_x_samples": len(ppd_bin_df),
        "mean_residual": float(ppd_bin_df["residual"].mean()),
        "mean_abs_residual": float(ppd_bin_df["residual"].abs().mean()),
        "rmse": float(np.sqrt((ppd_bin_df["residual"] ** 2).mean())),
        "z_score_mean": z_mean,
        "z_score_std": z_std,
        "tail_prob_ks_stat": float(ks_stat),
        "tail_prob_ks_pval": float(ks_pval),
        "frac_outside_90pct_interval": outside_90,
        "frac_outside_50pct_interval": outside_50,
        "expected_frac_outside_90pct": 0.10,
        "expected_frac_outside_50pct": 0.50,
    }])


# ── CLI ─────────────────────────────────────────────────────────────────────


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the ppd subcommand."""
    p = argparse.ArgumentParser(
        description="Posterior predictive check: generate replicated data "
                    "and compare to observed",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "-i", "--input", required=True,
        help="Preprocessed depth TSV (output of 'preprocess')",
    )
    p.add_argument(
        "-a", "--artifacts", required=True,
        help="inference_artifacts.npz (output of 'infer')",
    )
    p.add_argument(
        "-o", "--output-dir", required=True,
        help="Output directory",
    )
    p.add_argument(
        "--draws", type=int, default=100,
        help="Number of posterior predictive draws per bin/sample",
    )
    p.add_argument(
        "--site-data", default=None,
        help="Per-site allele data .npz (output of 'preprocess')",
    )
    p.add_argument(
        "--autosomal-baseline-cn-tsv",
        default=None,
        help=(
            "Optional TSV used to exclude samples omitted by infer before PPD. "
            "If omitted, ppd also looks for sample_autosomal_baseline_cn.tsv "
            "next to --artifacts."
        ),
    )
    p.add_argument(
        "--seed", type=int, default=42,
        help="Random seed for PPD sampling",
    )
    p.add_argument(
        "--device", choices=["cpu", "cuda"], default="cpu",
    )
    return p.parse_args()


def _run_ppd(args: argparse.Namespace, logger) -> None:
    """Run posterior predictive diagnostics after logging is configured."""

    # ── load inference artifacts ────────────────────────────────────────
    logger.info("Loading inference artifacts")
    map_est, cn_post = load_inference_artifacts(args.artifacts)
    depth_space = "raw"

    # ── load data ───────────────────────────────────────────────────────
    df = pd.read_csv(args.input, sep="\t", index_col=0)
    logger.info(
        "Loaded depth matrix for PPD: bins=%d samples=%d",
        len(df),
        len(get_sample_columns(df)),
    )

    sd = None
    if args.site_data:
        sd = load_site_data(args.site_data)
        logger.info("Loaded per-site allele tensors for PPD")

    df, sd = _prepare_ppd_inputs(
        df,
        sd,
        map_est,
        cn_post,
        args.artifacts,
        args.autosomal_baseline_cn_tsv,
    )

    if sd is not None:
        sd = apply_effective_site_pop_af(sd, map_est)

    data = DepthData(
        df, device=args.device, dtype=torch.float32,
        clamp_threshold=None,
        depth_space=depth_space,
        site_data=sd,
    )

    # ── generate posterior predictive draws ──────────────────────────────
    logger.info("Generating posterior predictive draws: draws=%d", args.draws)
    ppd_draws = generate_ppd_depth(
        data,
        map_est,
        cn_post,
        n_draws=args.draws,
        seed=args.seed,
    )

    # ── save raw PPD draws (compressed) ─────────────────────────────────
    draws_path = os.path.join(args.output_dir, "ppd_draws.npz")
    np.savez_compressed(draws_path, ppd_draws=ppd_draws)
    output_artifacts = [draws_path]

    # ── per-bin summary ─────────────────────────────────────────────────
    ppd_bin_df = compute_ppd_bin_summary(data, ppd_draws, map_est, cn_post)
    ppd_bin_path = os.path.join(args.output_dir, "ppd_bin_summary.tsv.gz")
    ppd_bin_df.to_csv(ppd_bin_path, sep="\t", index=False, compression="gzip")
    output_artifacts.append(ppd_bin_path)

    ppd_quality_df = compute_ppd_bin_quality_summary(ppd_bin_df)
    call_quality_df = compute_call_stability_quality_summary(data, cn_post)
    if not call_quality_df.empty:
        ppd_quality_df = ppd_quality_df.merge(
            call_quality_df,
            on=["chr", "start", "end"],
            how="left",
        )
    ppd_quality_path = os.path.join(args.output_dir, "ppd_bin_quality.tsv")
    ppd_quality_df.to_csv(ppd_quality_path, sep="\t", index=False)
    output_artifacts.append(ppd_quality_path)

    # ── per-chromosome summary ──────────────────────────────────────────
    ppd_chr_df = compute_ppd_chromosome_summary(data, ppd_draws, cn_post)
    ppd_chr_path = os.path.join(args.output_dir, "ppd_chromosome_summary.tsv")
    ppd_chr_df.to_csv(ppd_chr_path, sep="\t", index=False)
    output_artifacts.append(ppd_chr_path)

    # ── global calibration summary ──────────────────────────────────────
    ppd_global_df = compute_ppd_global_summary(ppd_bin_df)
    ppd_global_path = os.path.join(args.output_dir, "ppd_global_summary.tsv")
    ppd_global_df.to_csv(ppd_global_path, sep="\t", index=False)
    output_artifacts.append(ppd_global_path)
    logger.info(
        "PPD summaries written: bin_rows=%d chromosome_rows=%d",
        len(ppd_bin_df),
        len(ppd_chr_df),
    )
    log_output_artifacts(logger, output_artifacts)


def main() -> None:
    """Entry point for ``gatk-sv-ploidy ppd``."""
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    with tool_logging_context(
        tool_name="ppd",
        output_dir=args.output_dir,
        args=args,
        random_seeds={"ppd": args.seed},
    ) as logger:
        _run_ppd(args, logger)


if __name__ == "__main__":
    main()
