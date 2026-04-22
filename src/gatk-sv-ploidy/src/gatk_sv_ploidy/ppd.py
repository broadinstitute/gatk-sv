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
        \cdot \mathcal{N}(\tilde{d} \mid c \cdot b_i + \epsilon_i,
                           \sigma^2_{ij})

where *b_i* is the MAP ``bin_bias``, *\epsilon_i* is the MAP additive
background depth, and *p(c | x)* is the analytical CN posterior.

For ``negative_binomial`` fits on raw counts, the predictive mean becomes

.. math::

    \mu_{ij}(c) = \ell_i \cdot s_j \cdot (c \cdot b_i + \epsilon_i) / 2

with *\ell_i* the bin length in kilobases and *s_j* the MAP sample-depth
latent, and the predictive variance follows

.. math::

    \mu_{ij}(c) + (\alpha_i + \beta_j) \mu_{ij}(c)^2
"""

from __future__ import annotations

import argparse
import logging
import os
from typing import Dict

import numpy as np
import pandas as pd
import torch
from scipy import stats as sp_stats

from gatk_sv_ploidy._util import (
    DEPTH_SPACES,
    read_observation_type,
    validate_depth_space,
)
from gatk_sv_ploidy.data import DepthData, load_site_data
from gatk_sv_ploidy.infer import load_inference_artifacts
from gatk_sv_ploidy.models import (
    CNVModel,
    _matched_residual_scale,
    _normalize_obs_likelihood_name,
    _precompute_af_table,
    _raw_expected_depth_units,
)

logger = logging.getLogger(__name__)


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
    logger.info(
        "Using effective site AFs from inference artifacts for PPD "
        "(applied=%s).",
        bool(np.asarray(map_estimates.get("site_af_estimation_applied", False)).item()),
    )
    return updated


def _phred_scale_error_probability(
    error_prob: float,
    max_quality: float = 99.0,
) -> float:
    """Convert an error probability into a Phred-scaled quality score."""
    return float(min(max_quality, -10.0 * np.log10(max(error_prob, 1e-300))))


def _clamp_threshold_for_depth_space(depth_space: str) -> float | None:
    """Use clamping only for normalized-depth posterior predictive checks."""
    return 5.0 if depth_space == "normalized" else None


def _resolve_input_depth_space(
    requested_depth_space: str,
    obs_likelihood: str,
    input_path: str,
    map_estimates: Dict[str, np.ndarray],
) -> str:
    """Resolve PPD input depth space from preprocess marker or artifacts."""
    marker_depth_space = read_observation_type(input_path)
    if marker_depth_space is not None:
        return validate_depth_space(marker_depth_space, obs_likelihood)

    artifact_depth_space = requested_depth_space
    if artifact_depth_space == "auto" and "depth_space" in map_estimates:
        artifact_depth_space = np.asarray(map_estimates["depth_space"]).item()
    return validate_depth_space(artifact_depth_space, obs_likelihood)


def _extract_saved_posterior_draws(
    map_estimates: Dict[str, np.ndarray],
) -> Dict[str, np.ndarray]:
    """Return saved continuous latent draws from inference artifacts."""
    return {
        key[len("posterior_draws_"):]: value
        for key, value in map_estimates.items()
        if key.startswith("posterior_draws_")
    }


def _build_model_from_artifacts(
    map_estimates: Dict[str, np.ndarray],
    device: str,
    obs_likelihood: str,
    obs_df: float,
) -> CNVModel:
    """Reconstruct a model instance for analytical inference in PPD."""
    sex_prior = tuple(
        np.asarray(
            map_estimates.get("model_sex_prior", np.asarray([0.5, 0.5]))
        ).astype(np.float32).tolist()
    )
    dtype = torch.float64 if obs_likelihood == "negative_binomial" else torch.float32
    return CNVModel(
        n_states=int(np.asarray(map_estimates.get("model_n_states", 6)).item()),
        autosome_prior_mode=str(
            np.asarray(
                map_estimates.get("model_autosome_prior_mode", "shrinkage")
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
        var_bias_bin=float(
            np.asarray(map_estimates.get("model_var_bias_bin", 0.05)).item()
        ),
        var_sample=float(
            np.asarray(map_estimates.get("model_var_sample", 0.001)).item()
        ),
        var_bin=float(np.asarray(map_estimates.get("model_var_bin", 0.001)).item()),
        epsilon_mean=float(np.asarray(map_estimates.get("epsilon_mean", 0.0)).item()),
        device=device,
        dtype=dtype,
        guide_type=str(np.asarray(map_estimates.get("model_guide_type", "delta")).item()),
        af_concentration=float(
            np.asarray(map_estimates.get("model_af_concentration", 50.0)).item()
        ),
        af_weight=float(np.asarray(map_estimates.get("model_af_weight", 0.0)).item()),
        af_evidence_mode=str(
            np.asarray(map_estimates.get("model_af_evidence_mode", "absolute")).item()
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
        obs_likelihood=obs_likelihood,
        obs_df=obs_df,
        sample_depth_max=float(
            np.asarray(map_estimates.get("sample_depth_max", 10000.0)).item()
        ),
    )


def _sample_ppd_observation(
    data: DepthData,
    rng: np.random.RandomState,
    cn_probs: np.ndarray,
    continuous_maps: Dict[str, np.ndarray],
    obs_likelihood: str,
    obs_df: float,
) -> np.ndarray:
    """Sample one posterior predictive observation matrix."""
    bin_bias = np.atleast_1d(np.asarray(continuous_maps["bin_bias"]).squeeze())
    sample_var = np.atleast_1d(np.asarray(continuous_maps["sample_var"]).squeeze())
    bin_var = np.atleast_1d(np.asarray(continuous_maps["bin_var"]).squeeze())
    bin_epsilon = np.atleast_1d(
        np.asarray(continuous_maps.get("bin_epsilon", np.zeros_like(bin_bias))).squeeze()
    )

    n_bins, n_samples, n_states = cn_probs.shape
    flat_probs = cn_probs.reshape(-1, n_states)
    cum = np.cumsum(flat_probs, axis=1)
    u = rng.uniform(size=(flat_probs.shape[0], 1))
    cn_sample = (u < cum).argmax(axis=1).reshape(n_bins, n_samples)

    expected = cn_sample.astype(np.float32) * bin_bias[:, np.newaxis]
    expected = expected + bin_epsilon[:, np.newaxis]

    if obs_likelihood == "negative_binomial":
        sample_depth = np.atleast_1d(
            np.asarray(continuous_maps["sample_depth"]).squeeze()
        )
        overdispersion = bin_var[:, np.newaxis] + sample_var[np.newaxis, :]
        mean = data.bin_length_kb.detach().cpu().numpy()[:, np.newaxis]
        mean = mean * sample_depth[np.newaxis, :]
        mean = mean * np.maximum(
            _raw_expected_depth_units(
                cn_sample.astype(np.float32),
                bin_bias[:, np.newaxis],
                bin_epsilon[:, np.newaxis],
            ),
            0.0,
        )
        concentration = 1.0 / np.maximum(overdispersion, 1e-8)
        latent_rate = rng.gamma(
            shape=concentration,
            scale=np.where(mean > 0.0, mean / concentration, 0.0),
        )
        draw = rng.poisson(latent_rate)
    else:
        variance = bin_var[:, np.newaxis] + sample_var[np.newaxis, :]
        variance = variance * np.maximum(expected, 1e-3)
        scale = _matched_residual_scale(variance, obs_likelihood, obs_df)

        if obs_likelihood == "normal":
            draw = rng.normal(expected, scale)
        elif obs_likelihood == "laplace":
            draw = rng.laplace(expected, scale)
        else:
            draw = expected + rng.standard_t(obs_df, size=expected.shape) * scale

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
    MAP continuous parameters, draw from the mixture:

    1. Sample a CN state *c* from the posterior.
    2. Draw depth from the configured observation family.

    For continuous observation families, the predictive variance is
    ``(bin_var_i + sample_var_j) × max(c × bin_bias_i + epsilon_i, 1e-3)``.

    For the negative-binomial observation family on raw counts, the
    predictive mean is
    ``bin_length_kb_i × sample_depth_j × (c × bin_bias_i + epsilon_i) / 2`` and
    the predictive variance is
    ``μ + (bin_var_i + sample_var_j) × μ²``.

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

    obs_likelihood = _normalize_obs_likelihood_name(
        np.asarray(map_estimates.get("obs_likelihood", "normal")).item()
    )
    obs_df = float(np.asarray(map_estimates.get("obs_df", 6.0)).item())
    saved_draws = _extract_saved_posterior_draws(map_estimates)
    required_sites = {"bin_bias", "sample_var", "bin_var"}
    if obs_likelihood == "negative_binomial":
        required_sites.add("sample_depth")

    if required_sites.issubset(saved_draws):
        draw_count = int(next(iter(saved_draws.values())).shape[0])
        if all(int(value.shape[0]) == draw_count for value in saved_draws.values()):
            logger.info(
                "Using %d saved continuous posterior draws for PPD integration.",
                draw_count,
            )
            model = _build_model_from_artifacts(
                map_estimates,
                data.depth.device.type,
                obs_likelihood,
                obs_df,
            )
            af_table_np = None
            if data.site_alt is not None and model.af_weight > 0:
                with torch.no_grad():
                    af_table_np = _precompute_af_table(
                        data.site_alt,
                        data.site_total,
                        data.site_pop_af,
                        data.site_mask,
                        n_states=model.n_states,
                        concentration=model.af_concentration,
                    ).cpu().numpy()

            draw_indices = rng.randint(draw_count, size=n_draws)
            draws = np.empty((n_draws, data.n_bins, data.n_samples), dtype=np.float32)
            for out_idx, draw_idx in enumerate(draw_indices):
                draw_maps = {
                    site: np.asarray(value[draw_idx])
                    for site, value in saved_draws.items()
                }
                draw_post = model._run_discrete_inference_fixed_latents(
                    data,
                    draw_maps,
                    af_table=af_table_np,
                )
                draws[out_idx] = _sample_ppd_observation(
                    data,
                    rng,
                    draw_post["cn_posterior"],
                    draw_maps,
                    obs_likelihood,
                    obs_df,
                )
            return draws

        logger.warning(
            "Saved posterior-draw arrays have inconsistent leading dimensions; "
            "falling back to plug-in PPD.",
        )

    cn_probs = cn_posterior["cn_posterior"]
    draws = np.empty((n_draws, data.n_bins, data.n_samples), dtype=np.float32)
    for d in range(n_draws):
        draws[d] = _sample_ppd_observation(
            data,
            rng,
            cn_probs,
            map_estimates,
            obs_likelihood,
            obs_df,
        )
    return draws


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
    """Aggregate multi-draw call stability into per-bin quality metrics."""
    thresholds = tuple(float(t) for t in instability_thresholds)
    if any(t <= 0.0 or t >= 1.0 for t in thresholds):
        raise ValueError("instability_thresholds must lie strictly between 0 and 1.")

    stability = cn_posterior.get("cn_map_stability")
    if stability is None:
        return pd.DataFrame(columns=["chr", "start", "end"])

    stability = np.asarray(stability, dtype=np.float64)
    if stability.shape != (data.n_bins, data.n_samples):
        raise ValueError("cn_map_stability must have shape (n_bins, n_samples).")

    rows: list[dict] = []
    for i in range(data.n_bins):
        bin_stability = np.clip(stability[i, :], 0.0, 1.0)
        bin_instability = 1.0 - bin_stability
        alpha = 1.0 + float(bin_instability.sum())
        beta = 1.0 + float(bin_stability.sum())
        row = {
            "chr": data.chr[i],
            "start": int(data.start[i]),
            "end": int(data.end[i]),
            "mean_call_stability": float(bin_stability.mean()),
            "median_call_stability": float(np.median(bin_stability)),
            "mean_call_instability": float(bin_instability.mean()),
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
        "--depth-space", choices=["auto", *DEPTH_SPACES], default="auto",
        help="Interpret the input matrix as normalized depth or raw counts. 'auto' first consults preprocess observation_type.txt, then saved inference metadata, and finally the observation likelihood.",
    )
    p.add_argument(
        "--seed", type=int, default=42,
        help="Random seed for PPD sampling",
    )
    p.add_argument(
        "--device", choices=["cpu", "cuda"], default="cpu",
    )
    return p.parse_args()


def main() -> None:
    """Entry point for ``gatk-sv-ploidy ppd``."""
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    # ── load inference artifacts ────────────────────────────────────────
    logger.info("Loading inference artifacts: %s", args.artifacts)
    map_est, cn_post = load_inference_artifacts(args.artifacts)
    obs_likelihood = _normalize_obs_likelihood_name(
        np.asarray(map_est.get("obs_likelihood", "normal")).item()
    )
    depth_space = _resolve_input_depth_space(
        args.depth_space,
        obs_likelihood,
        args.input,
        map_est,
    )

    # ── load data ───────────────────────────────────────────────────────
    logger.info("Loading preprocessed depth: %s", args.input)
    df = pd.read_csv(args.input, sep="\t", index_col=0)

    sd = None
    if args.site_data:
        sd = load_site_data(args.site_data)
        sd = apply_effective_site_pop_af(sd, map_est)

    data = DepthData(
        df, device=args.device, dtype=torch.float32,
        clamp_threshold=_clamp_threshold_for_depth_space(depth_space),
        depth_space=depth_space,
        site_data=sd,
    )
    logger.info("Using %s input with %s observation model.", depth_space, obs_likelihood)

    # ── generate posterior predictive draws ──────────────────────────────
    logger.info("Generating %d posterior predictive draws …", args.draws)
    ppd_draws = generate_ppd_depth(
        data, map_est, cn_post, n_draws=args.draws, seed=args.seed,
    )
    logger.info("PPD draws shape: %s", ppd_draws.shape)

    # ── save raw PPD draws (compressed) ─────────────────────────────────
    draws_path = os.path.join(args.output_dir, "ppd_draws.npz")
    np.savez_compressed(draws_path, ppd_draws=ppd_draws)
    logger.info("PPD draws saved to %s", draws_path)

    # ── per-bin summary ─────────────────────────────────────────────────
    logger.info("Computing per-bin PPD summary …")
    ppd_bin_df = compute_ppd_bin_summary(data, ppd_draws, map_est, cn_post)
    ppd_bin_path = os.path.join(args.output_dir, "ppd_bin_summary.tsv.gz")
    ppd_bin_df.to_csv(ppd_bin_path, sep="\t", index=False, compression="gzip")
    logger.info("PPD bin summary saved to %s", ppd_bin_path)

    logger.info("Computing per-bin PPD quality summary …")
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
    logger.info("PPD bin quality summary saved to %s", ppd_quality_path)

    # ── per-chromosome summary ──────────────────────────────────────────
    logger.info("Computing per-chromosome PPD summary …")
    ppd_chr_df = compute_ppd_chromosome_summary(data, ppd_draws, cn_post)
    ppd_chr_path = os.path.join(args.output_dir, "ppd_chromosome_summary.tsv")
    ppd_chr_df.to_csv(ppd_chr_path, sep="\t", index=False)
    logger.info("PPD chromosome summary saved to %s", ppd_chr_path)

    # ── global calibration summary ──────────────────────────────────────
    logger.info("Computing global calibration summary …")
    ppd_global_df = compute_ppd_global_summary(ppd_bin_df)
    ppd_global_path = os.path.join(args.output_dir, "ppd_global_summary.tsv")
    ppd_global_df.to_csv(ppd_global_path, sep="\t", index=False)
    logger.info("PPD global summary saved to %s", ppd_global_path)

    # ── print key calibration metrics ───────────────────────────────────
    row = ppd_global_df.iloc[0]
    logger.info("")
    logger.info("=== Posterior Predictive Check Summary ===")
    logger.info("  RMSE:                    %.4f", row["rmse"])
    logger.info("  Mean z-score:            %.4f (ideal: 0)", row["z_score_mean"])
    logger.info("  Std z-score:             %.4f (ideal: 1)", row["z_score_std"])
    logger.info("  Tail prob KS p-value:    %.4f (>0.05 = well-calibrated)",
                row["tail_prob_ks_pval"])
    logger.info("  Frac outside 90%% PI:    %.3f (ideal: 0.10)",
                row["frac_outside_90pct_interval"])
    logger.info("  Frac outside 50%% PI:    %.3f (ideal: 0.50)",
                row["frac_outside_50pct_interval"])

    logger.info("\nDone.")


if __name__ == "__main__":
    main()
