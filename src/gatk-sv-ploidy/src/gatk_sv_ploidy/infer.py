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
    DEFAULT_AF_WEIGHT,
    DEPTH_SPACES,
    compute_cnq_from_probabilities,
    default_obs_likelihood_for_depth_space,
    read_observation_type,
    summarize_contig_ploidy_from_bin_calls,
    validate_depth_space,
)
from gatk_sv_ploidy.data import DepthData, load_site_data
from gatk_sv_ploidy.models import (
    AUTOSOME_PRIOR_MODES,
    OBS_LIKELIHOODS,
    CNVModel,
)

logger = logging.getLogger(__name__)
SITE_AF_ESTIMATORS = ("off", "auto", "naive-bayes")


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
            if "bin_epsilon" in map_estimates:
                row["bin_epsilon"] = float(map_estimates["bin_epsilon"].flatten()[i])
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
    g.add_argument("--alpha-ref", type=float, default=50.0,
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
    g.add_argument("--var-bias-bin", type=float, default=0.02,
                   help="LogNormal scale for per-bin bias")
    g.add_argument("--var-sample", type=float, default=0.00025,
                   help="Exponential mean for per-sample variance")
    g.add_argument("--var-bin", type=float, default=0.001,
                   help="Exponential mean for per-bin variance")
    g.add_argument("--epsilon-mean", type=float, default=1e-2,
                   help="Mean of the Exponential prior on per-bin additive background depth; set to 0 to disable")
    g.add_argument("--guide-type", choices=["delta", "diagonal", "lowrank"], default="delta",
                   help="Variational guide type")
    g.add_argument("--obs-likelihood", choices=["auto", *list(OBS_LIKELIHOODS)], default="auto",
                   help="Observation likelihood for the depth/count matrix. 'auto' resolves from preprocess observation_type.txt when present.")
    g.add_argument("--obs-df", type=float, default=3.5,
                   help="Student-t degrees of freedom when --obs-likelihood=studentt")
    g.add_argument("--sample-depth-max", type=float, default=10000.0,
                   help="Upper bound of the Uniform prior on per-sample depth scale when --obs-likelihood=negative_binomial")

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
                   help="Global scale/temperature applied to the summed per-bin allele fraction log-likelihood (0 to disable)")
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
    g.add_argument("--learn-af-temperature", action="store_true", default=False,
                   help="Learn a single global AF temperature instead of keeping --af-weight fixed")
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
    g.add_argument("--log-freq", type=int, default=50)
    g.add_argument("--jit", action="store_true", default=False)
    g.add_argument("--early-stopping", action="store_true", default=True)
    g.add_argument("--no-early-stopping", dest="early_stopping",
                   action="store_false")
    g.add_argument("--patience", type=int, default=50)
    g.add_argument("--min-delta", type=float, default=1000.0)

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
    logger.info("Loading preprocessed depth: %s", args.input)
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

    learn_af_temperature = args.learn_af_temperature and data.site_alt is not None
    learn_site_af = args.learn_site_af and data.site_alt is not None

    if data.site_alt is not None:
        logger.info(
            "Allele-fraction evidence enabled (af_weight=%.2f, af_concentration=%.1f, mode=%s, summed over observed sites/bin, learn_af_temperature=%s)",
            args.af_weight,
            args.af_concentration,
            args.af_evidence_mode,
            learn_af_temperature,
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
    )

    model.train(
        data,
        max_iter=args.max_iter,
        lr_init=args.lr_init,
        lr_min=args.lr_min,
        lr_decay=args.lr_decay,
        log_freq=args.log_freq,
        jit=args.jit,
        early_stopping=args.early_stopping,
        patience=args.patience,
        min_delta=args.min_delta,
    )

    # ── save training loss ──────────────────────────────────────────────
    loss_df = pd.DataFrame(model.loss_history)
    loss_path = os.path.join(args.output_dir, "training_loss.tsv")
    loss_df.to_csv(loss_path, sep="\t", index=False)
    logger.info("Training loss saved to %s", loss_path)

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

    # ── log sex karyotype results ───────────────────────────────────────
    if "sex_posterior" in cn_post:
        sex_post = cn_post["sex_posterior"]  # (n_samples, 2)
        n_xx = int((sex_post[:, 0] > 0.5).sum())
        n_xy = int((sex_post[:, 1] > 0.5).sum())
        logger.info(
            "Sex karyotype assignment: %d XX, %d XY (of %d samples)",
            n_xx, n_xy, data.n_samples,
        )
    if "sex_karyotype" in map_est:
        sk = map_est["sex_karyotype"]
        logger.info(
            "Point-estimate sex karyotype: %d XX, %d XY",
            int((sk == 0).sum()), int((sk == 1).sum()),
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
        logger.info("Site AF estimates saved to %s", site_af_path)

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
    logger.info("Inference artifacts saved to %s", artifacts_path)

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
    logger.info("Bin statistics saved to %s", bin_path)

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
    logger.info("Chromosome statistics saved to %s", chr_path)


if __name__ == "__main__":
    main()
