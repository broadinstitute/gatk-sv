"""
Hierarchical Bayesian model for copy-number inference (Pyro).

Provides :class:`CNVModel` which wraps the probabilistic model, the
variational guide, training via SVI, MAP estimation, and exact discrete
posterior inference over copy-number states.

The allele-fraction likelihood analytically marginalises over SNP genotype
states at every site in each bin, avoiding the circular-reasoning problem
of pre-filtering sites by observed allele fraction.
"""

from __future__ import annotations

import copy
import math
import logging
import random
from typing import Dict, List, Optional, Sequence

import numpy as np
import torch
import pyro
import pyro.distributions as dist
from pyro import poutine
from pyro.ops.indexing import Vindex
from pyro.infer import (
    SVI,
    JitTraceEnum_ELBO,
    TraceEnum_ELBO,
    config_enumerate,
    infer_discrete,
)
from pyro.infer.autoguide import (
    AutoDelta,
    AutoDiagonalNormal,
    AutoGuideList,
    AutoLowRankMultivariateNormal,
)
from pyro.infer.autoguide.initialization import init_to_median, init_to_sample
from tqdm import tqdm

from gatk_sv_ploidy._util import (
    compose_additive_background_matrix,
    DEFAULT_AF_CONCENTRATION,
    DEFAULT_AF_WEIGHT,
    NEGATIVE_BINOMIAL_OBS_LIKELIHOOD,
)
from gatk_sv_ploidy.data import DepthData

logger = logging.getLogger(__name__)

_SAMPLE_DEPTH_PRIOR_BOOTSTRAP_DRAWS = 256
_SAMPLE_DEPTH_PRIOR_BOOTSTRAP_SEED = 0
_SAMPLE_DEPTH_PRIOR_SD_FLOOR_FRAC = 0.05


OBS_LIKELIHOODS = (
    "normal",
    "studentt",
    "laplace",
    NEGATIVE_BINOMIAL_OBS_LIKELIHOOD,
)

AF_EVIDENCE_MODES = (
    "absolute",
    "relative",
)

AUTOSOME_PRIOR_MODES = (
    "dirichlet",
    "shrinkage",
)

DEFAULT_EPSILON_MEAN = 1e-4
DEFAULT_EPSILON_CONCENTRATION = 0.5
DEFAULT_BACKGROUND_FACTORS = 2
DEFAULT_BACKGROUND_SAMPLE_SCALE = 0.02
DEFAULT_BACKGROUND_BIN_SCALE = 0.5


def _windowed_relative_elbo_change(
    loss_history: Sequence[float],
    window: int,
) -> Optional[float]:
    """Return the relative ELBO shift between the two latest rolling windows."""
    if window < 1:
        raise ValueError("convergence_window must be at least 1.")
    if len(loss_history) < 2 * window:
        return None

    previous_window = np.asarray(loss_history[-2 * window:-window], dtype=np.float64)
    current_window = np.asarray(loss_history[-window:], dtype=np.float64)
    previous_mean = float(previous_window.mean())
    current_mean = float(current_window.mean())
    baseline = max(abs(previous_mean), np.finfo(np.float64).eps)
    return abs(current_mean - previous_mean) / baseline


# ── marginalized allele-fraction likelihood ─────────────────────────────────


def _marginalized_af_log_lik(
    site_alt: torch.Tensor,
    site_total: torch.Tensor,
    site_pop_af: torch.Tensor,
    site_mask: torch.Tensor,
    cn: torch.Tensor,
    n_states: int,
    concentration: float,
) -> torch.Tensor:
    """Compute per-bin AF log-likelihood marginalised over genotype states.

    For each SNP site *s* inside a bin with copy number *c*, the genotype
    (number of alt alleles) *k* ranges from 0 to *c*.  We marginalise:

    .. math::

        P(a_s | c, p_s, n_s) = \\sum_{k=0}^{c}
            \\mathrm{Binom}(k | c, p_s) \\cdot
            \\mathrm{BetaBinom}(a_s | \\alpha_k, \\beta_k, n_s)

    where *p_s* is the population alt-allele frequency, *a_s* / *n_s* are the
    observed alt / total counts, and ``α_k = conc * (k/c) + ε``,
    ``β_k = conc * (1 − k/c) + ε``.

    This is injected into the model via :func:`pyro.factor`.

    Args:
        site_alt: ``(n_bins, max_sites, n_samples)`` observed alt counts.
        site_total: ``(n_bins, max_sites, n_samples)`` observed total counts.
        site_pop_af: Population alt-allele frequencies with shape
            ``(n_bins, max_sites)`` or sample-specific leave-one-out values with
            shape ``(n_bins, max_sites, n_samples)``.
        site_mask: ``(n_bins, max_sites, n_samples)`` boolean validity mask.
        cn: ``(n_bins, n_samples)`` integer copy-number state per bin/sample
            (may have extra leading enumeration dimensions from Pyro).
        n_states: Maximum CN state (6 for CN 0–5).
        concentration: Beta-Binomial concentration parameter.

    Returns:
        Per-bin, per-sample log-likelihood averaged over valid sites.
    """
    eps = 1e-6

    # cn may have enumeration dims from config_enumerate: (..., n_bins, n_samples)
    # Insert a sites dim so it broadcasts: (..., n_bins, 1, n_samples)
    cn_e = cn.unsqueeze(-2)
    cn_f = cn_e.float()

    site_total_safe = torch.clamp(site_total, min=1)
    if site_pop_af.dim() == 2:
        pop_af = site_pop_af.unsqueeze(-1)  # (n_bins, max_sites, 1)
    elif site_pop_af.dim() == 3:
        pop_af = site_pop_af
    else:
        raise ValueError(
            "site_pop_af must have shape (n_bins, max_sites) or "
            "(n_bins, max_sites, n_samples)."
        )
    pop_af_safe = torch.clamp(pop_af, min=eps, max=1.0 - eps)

    # Iterate over genotype states k = 0 .. n_states-1
    log_terms = []
    for k in range(n_states):
        k_f = float(k)

        # Genotype prior: Binom(k | c, p_s)
        valid_geno = (cn_e >= k) & (cn_f > 0)

        # Manual log Binom(k | c, p) to avoid distribution validation issues
        # and explicit expand calls — rely on broadcasting throughout.
        cn_safe = cn_f.clamp(min=1)
        log_comb = torch.lgamma(cn_safe + 1)
        log_comb = log_comb - torch.lgamma(
            torch.tensor(k_f + 1, device=cn_f.device),
        )
        log_comb = log_comb - torch.lgamma(cn_safe - k_f + 1)
        log_p = k_f * torch.log(pop_af_safe) + (cn_safe - k_f) * torch.log(
            1.0 - pop_af_safe,
        )
        log_binom = log_comb + log_p

        log_binom_prior = torch.where(
            valid_geno,
            log_binom,
            torch.where(
                (cn_e == 0) & (k == 0),
                torch.zeros_like(cn_f),  # log(1) = 0 for k=0, c=0
                torch.full_like(cn_f, -1e10),
            ),
        )

        # Allele-fraction likelihood: BetaBinom(a_s | alpha, beta, n_s)
        # Expected allele fraction for genotype k at copy number c: k/c
        # For c=0: uninformative (flat); clamp to [0,1] for k > cn safety
        eaf = torch.where(
            cn_f > 0,
            (torch.tensor(k_f, device=cn_f.device) / cn_safe).clamp(0.0, 1.0),
            torch.tensor(0.5, device=cn_f.device),
        )
        alpha = concentration * eaf + eps
        beta = concentration * (1.0 - eaf) + eps

        log_bb = dist.BetaBinomial(
            alpha, beta, site_total_safe, validate_args=False,
        ).log_prob(site_alt)

        log_terms.append(log_binom_prior + log_bb)

    # logsumexp over genotype states
    log_stack = torch.stack(log_terms, dim=0)  # (n_states, ..., n_bins, max_sites, n_samples)
    log_marginal = torch.logsumexp(log_stack, dim=0)  # (..., n_bins, max_sites, n_samples)

    # Mask out invalid sites and sum over observed sites so bins with more
    # informative loci contribute proportionally more AF evidence.
    log_marginal = torch.where(site_mask, log_marginal, torch.zeros_like(log_marginal))
    return log_marginal.sum(dim=-2)


def _marginalized_af_log_lik_numpy(
    site_alt: np.ndarray,
    site_total: np.ndarray,
    site_pop_af: np.ndarray,
    site_mask: np.ndarray,
    cn_state: int,
    n_states: int,
    concentration: float,
) -> np.ndarray:
    """NumPy version for analytical discrete inference.

    Computes the per-bin AF log-likelihood for a single CN state, marginalised
    over genotypes.

    Args:
        site_alt: ``(n_bins, max_sites, n_samples)``
        site_total: ``(n_bins, max_sites, n_samples)``
        site_pop_af: ``(n_bins, max_sites)`` or sample-specific values with
            shape ``(n_bins, max_sites, n_samples)``.
        site_mask: ``(n_bins, max_sites, n_samples)``
        cn_state: The CN value (integer).
        n_states: Max CN.
        concentration: Beta-Binomial concentration.

    Returns:
        ``(n_bins, n_samples)`` log-likelihood summed over observed sites
        within each bin.
    """
    return _site_level_marginalized_af_log_lik_numpy(
        site_alt,
        site_total,
        site_pop_af,
        site_mask,
        cn_state=cn_state,
        n_states=n_states,
        concentration=concentration,
    ).sum(axis=1)


def _site_level_marginalized_af_log_lik_numpy(
    site_alt: np.ndarray,
    site_total: np.ndarray,
    site_pop_af: np.ndarray,
    site_mask: np.ndarray,
    cn_state: int,
    n_states: int,
    concentration: float,
) -> np.ndarray:
    """Return the per-site marginalized AF log-likelihood for one CN state.

    This matches :func:`_marginalized_af_log_lik_numpy` before the final sum over
    sites, so diagnostics can inspect site-level AF contributions without
    changing the inference math.
    """
    from scipy.stats import betabinom as betabinom_scipy
    from scipy.stats import binom as binom_scipy

    eps = 1e-6
    site_total_safe = np.maximum(site_total, 1)
    if site_pop_af.ndim == 2:
        pop_af = site_pop_af[:, :, np.newaxis]   # (n_bins, max_sites, 1)
    elif site_pop_af.ndim == 3:
        pop_af = site_pop_af
    else:
        raise ValueError(
            "site_pop_af must have shape (n_bins, max_sites) or "
            "(n_bins, max_sites, n_samples)."
        )
    pop_af_safe = np.clip(pop_af, eps, 1.0 - eps)

    log_terms = []
    for k in range(cn_state + 1):
        # Genotype prior: Binom(k | cn_state, pop_af)
        if cn_state > 0:
            log_bp = binom_scipy.logpmf(k, cn_state, pop_af_safe)
        else:
            # cn=0: only k=0 is valid
            log_bp = np.where(k == 0, 0.0, -1e10) * np.ones_like(site_alt, dtype=np.float64)

        # BetaBinomial likelihood
        eaf = k / max(cn_state, 1)
        alpha = concentration * eaf + eps
        beta = concentration * (1.0 - eaf) + eps
        log_bb = betabinom_scipy.logpmf(site_alt, site_total_safe, alpha, beta)
        log_bb = np.nan_to_num(log_bb, nan=0.0, posinf=0.0, neginf=-100.0)

        log_terms.append(log_bp + log_bb)

    # logsumexp over genotype states
    log_stack = np.stack(log_terms, axis=0)  # (n_geno, n_bins, max_sites, n_samples)
    max_val = np.max(log_stack, axis=0, keepdims=True)
    log_marginal = max_val.squeeze(0) + np.log(
        np.sum(np.exp(log_stack - max_val), axis=0) + 1e-30
    )

    return np.where(site_mask, log_marginal, 0.0)


def _precompute_af_table(
    site_alt: torch.Tensor,
    site_total: torch.Tensor,
    site_pop_af: torch.Tensor,
    site_mask: torch.Tensor,
    n_states: int,
    concentration: float,
) -> torch.Tensor:
    """Precompute AF log-likelihood for every CN state.

    The result is constant across SVI iterations (it depends only on observed
    data and the discrete CN state, not on any continuous latent variables),
    so computing it once and reusing it via a simple table lookup eliminates
    the dominant per-step cost.

    Args:
        site_alt: ``(n_bins, max_sites, n_samples)``
        site_total: ``(n_bins, max_sites, n_samples)``
        site_pop_af: ``(n_bins, max_sites)``
        site_mask: ``(n_bins, max_sites, n_samples)``
        n_states: Number of CN states.
        concentration: Beta-Binomial concentration.

    Returns:
        Tensor of shape ``(n_states, n_bins, n_samples)`` where entry
        ``[c, b, s]`` is the marginalized AF log-likelihood for CN state *c*
        at bin *b* and sample *s*, summed over observed sites.
    """
    n_bins = site_alt.shape[0]
    n_samples = site_alt.shape[2]
    table = torch.empty(n_states, n_bins, n_samples,
                        device=site_alt.device, dtype=site_alt.dtype)
    for c in range(n_states):
        cn_c = torch.full((n_bins, n_samples), c, dtype=torch.long,
                          device=site_alt.device)
        table[c] = _marginalized_af_log_lik(
            site_alt, site_total, site_pop_af, site_mask,
            cn=cn_c, n_states=n_states, concentration=concentration,
        )
    return table


def _center_af_table_torch(
    af_table: torch.Tensor,
    reference_cn_probs: torch.Tensor,
) -> torch.Tensor:
    """Convert AF log-likelihoods into log-likelihood ratios vs. a fixed reference mixture."""
    ref_log_probs = torch.log(
        torch.clamp(reference_cn_probs, min=1e-10),
    ).transpose(0, 1).unsqueeze(-1)
    baseline = torch.logsumexp(ref_log_probs + af_table, dim=0, keepdim=True)
    return af_table - baseline


def _center_af_table_numpy(
    af_table: np.ndarray,
    reference_cn_probs: np.ndarray,
) -> np.ndarray:
    """NumPy counterpart of :func:`_center_af_table_torch`."""
    ref_log_probs = np.log(np.maximum(reference_cn_probs, 1e-10)).T[:, :, np.newaxis]
    raw = ref_log_probs + af_table
    max_val = np.max(raw, axis=0, keepdims=True)
    baseline = max_val + np.log(
        np.sum(np.exp(raw - max_val), axis=0, keepdims=True) + 1e-30,
    )
    return af_table - baseline


def _cohort_site_pop_af_from_counts_torch(
    site_alt: torch.Tensor,
    site_total: torch.Tensor,
    site_mask: torch.Tensor,
) -> torch.Tensor:
    """Estimate cohort AFs from the observed counts using Jeffreys smoothing."""
    valid_counts = site_mask & (site_total > 0)
    valid_alt = torch.where(valid_counts, site_alt, torch.zeros_like(site_alt))
    valid_total = torch.where(valid_counts, site_total, torch.zeros_like(site_total))
    pooled_alt = valid_alt.sum(dim=-1)
    pooled_total = valid_total.sum(dim=-1)
    return torch.where(
        pooled_total > 0,
        (pooled_alt + 0.5) / (pooled_total + 1.0),
        torch.full_like(pooled_total, 0.5),
    )


def _cohort_site_pop_af_from_counts_numpy(
    site_alt: np.ndarray,
    site_total: np.ndarray,
    site_mask: np.ndarray,
) -> np.ndarray:
    """NumPy counterpart of :func:`_cohort_site_pop_af_from_counts_torch`."""
    valid_counts = np.asarray(site_mask, dtype=bool) & (site_total > 0)
    valid_alt = np.where(valid_counts, site_alt, 0.0)
    valid_total = np.where(valid_counts, site_total, 0.0)
    pooled_alt = valid_alt.sum(axis=-1, dtype=np.float64)
    pooled_total = valid_total.sum(axis=-1, dtype=np.float64)
    return np.where(
        pooled_total > 0,
        (pooled_alt + 0.5) / (pooled_total + 1.0),
        0.5,
    )


def _resolve_fixed_site_pop_af_torch(
    site_alt: torch.Tensor,
    site_total: torch.Tensor,
    site_pop_af: torch.Tensor,
    site_mask: torch.Tensor,
    *,
    atol: float = 1e-3,
) -> tuple[torch.Tensor, bool]:
    """Expand self-pooled AF inputs into sample-excluded leave-one-out AFs."""
    if site_pop_af.dim() != 2:
        return site_pop_af, False

    observed_sites = site_mask.any(dim=-1)
    if not bool(observed_sites.any()):
        return site_pop_af, False

    pooled_site_pop_af = _cohort_site_pop_af_from_counts_torch(
        site_alt,
        site_total,
        site_mask,
    )
    if not torch.allclose(
        site_pop_af[observed_sites],
        pooled_site_pop_af[observed_sites],
        atol=atol,
        rtol=0.0,
    ):
        return site_pop_af, False

    valid_counts = site_mask & (site_total > 0)
    valid_alt = torch.where(valid_counts, site_alt, torch.zeros_like(site_alt))
    valid_total = torch.where(valid_counts, site_total, torch.zeros_like(site_total))
    pooled_alt = valid_alt.sum(dim=-1, keepdim=True)
    pooled_total = valid_total.sum(dim=-1, keepdim=True)
    leave_one_out_alt = pooled_alt - valid_alt
    leave_one_out_total = pooled_total - valid_total
    fallback_af = site_pop_af.unsqueeze(-1).expand_as(site_total)
    leave_one_out_af = torch.where(
        leave_one_out_total > 0,
        (leave_one_out_alt + 0.5) / (leave_one_out_total + 1.0),
        fallback_af,
    )
    return leave_one_out_af, True


def _resolve_fixed_site_pop_af_numpy(
    site_alt: np.ndarray,
    site_total: np.ndarray,
    site_pop_af: np.ndarray,
    site_mask: np.ndarray,
    *,
    atol: float = 1e-3,
) -> tuple[np.ndarray, bool]:
    """NumPy counterpart of :func:`_resolve_fixed_site_pop_af_torch`."""
    if site_pop_af.ndim != 2:
        return site_pop_af, False

    observed_sites = np.any(site_mask, axis=-1)
    if not np.any(observed_sites):
        return site_pop_af, False

    pooled_site_pop_af = _cohort_site_pop_af_from_counts_numpy(
        site_alt,
        site_total,
        site_mask,
    )
    if not np.allclose(
        site_pop_af[observed_sites],
        pooled_site_pop_af[observed_sites],
        atol=atol,
        rtol=0.0,
    ):
        return site_pop_af, False

    valid_counts = np.asarray(site_mask, dtype=bool) & (site_total > 0)
    valid_alt = np.where(valid_counts, site_alt, 0.0)
    valid_total = np.where(valid_counts, site_total, 0.0)
    pooled_alt = valid_alt.sum(axis=-1, keepdims=True, dtype=np.float64)
    pooled_total = valid_total.sum(axis=-1, keepdims=True, dtype=np.float64)
    leave_one_out_alt = pooled_alt - valid_alt
    leave_one_out_total = pooled_total - valid_total
    fallback_af = np.broadcast_to(site_pop_af[:, :, np.newaxis], site_total.shape)
    leave_one_out_af = np.where(
        leave_one_out_total > 0,
        (leave_one_out_alt + 0.5) / (leave_one_out_total + 1.0),
        fallback_af,
    )
    return leave_one_out_af, True


def _normalize_obs_likelihood_name(name: str) -> str:
    """Validate and normalize observation likelihood names."""
    normalized = str(name).strip().lower()
    if normalized not in OBS_LIKELIHOODS:
        raise ValueError(
            f"Unknown obs_likelihood: {name!r}. "
            f"Choose one of {OBS_LIKELIHOODS}."
        )
    return normalized


def _normalize_af_evidence_mode(name: str) -> str:
    """Validate and normalize AF evidence mode names."""
    normalized = str(name).strip().lower()
    if normalized not in AF_EVIDENCE_MODES:
        raise ValueError(
            f"Unknown af_evidence_mode: {name!r}. "
            f"Choose one of {AF_EVIDENCE_MODES}."
        )
    return normalized


def _negative_binomial_total_count(
    overdispersion: torch.Tensor,
) -> torch.Tensor:
    """Convert NB overdispersion to a Gamma-Poisson concentration."""
    return 1.0 / torch.clamp(overdispersion, min=1e-8)


def _make_negative_binomial_distribution(
    mean: torch.Tensor,
    overdispersion: torch.Tensor,
) -> dist.Distribution:
    """Build a Gamma-Poisson distribution with ``Var = μ + α μ²``."""
    concentration = _negative_binomial_total_count(overdispersion)
    mean_safe = torch.clamp(mean, min=1e-8)
    rate = concentration / mean_safe
    return dist.GammaPoisson(concentration, rate)


def _negative_binomial_log_lik_numpy(
    obs: np.ndarray,
    mean: np.ndarray,
    overdispersion: np.ndarray,
) -> np.ndarray:
    """Analytical Gamma-Poisson log-likelihood with ``Var = μ + α μ²``."""
    from scipy.special import gammaln

    obs_rounded = np.rint(obs)
    if not np.allclose(obs, obs_rounded, atol=1e-6):
        raise ValueError(
            "negative_binomial observation likelihood requires integer-valued raw counts."
        )

    obs_safe = obs_rounded.astype(np.float64, copy=False)
    mean_safe = np.maximum(mean, 1e-10)
    concentration = 1.0 / np.maximum(overdispersion, 1e-8)
    log_prob = gammaln(obs_safe + concentration)
    log_prob = log_prob - gammaln(concentration)
    log_prob = log_prob - gammaln(obs_safe + 1.0)
    log_prob = log_prob + concentration * np.log(
        concentration / (concentration + mean_safe)
    )
    log_prob = log_prob + obs_safe * np.log(
        mean_safe / (concentration + mean_safe)
    )

    zero_mean = mean <= 0.0
    log_prob = np.where(zero_mean & (obs_safe == 0.0), 0.0, log_prob)
    log_prob = np.where(zero_mean & (obs_safe > 0.0), -1e10, log_prob)
    return log_prob


def _raw_expected_depth_units(
    copy_number,
    bin_bias,
    additive_background,
):
    """Map normalized-depth expectations onto a fixed diploid raw-count baseline.

    The continuous-depth model uses ``CN * bin_bias + additive_background``
    with autosomal bins centered near diploid copy number. The additive
    background may include the contig-shared epsilon floor and any fitted
    low-rank background factors after expansion onto bins. Raw counts are
    parameterized by each sample's diploid baseline depth per kb, so we
    convert the normalized-depth expectation back to that baseline scale by
    dividing by 2 before multiplying by ``sample_depth`` and bin length.
    """
    return (copy_number * bin_bias + additive_background) / 2.0


def _normalize_background_bin_factors_torch(values: torch.Tensor) -> torch.Tensor:
    """Normalize per-factor bin loadings to unit mean across bins."""
    if values.dim() != 2:
        raise ValueError("background_bin_factors must be 2D.")
    if values.shape[1] == 0:
        return values
    eps = torch.finfo(values.dtype).eps
    return values / torch.clamp(values.mean(dim=0, keepdim=True), min=eps)


def _compose_additive_background_torch(
    bin_epsilon: torch.Tensor,
    background_bin_factors: torch.Tensor | None = None,
    background_sample_factors: torch.Tensor | None = None,
) -> torch.Tensor:
    """Compose the effective additive background matrix for the Pyro model."""
    additive = bin_epsilon
    if background_bin_factors is None and background_sample_factors is None:
        return additive
    if background_bin_factors is None or background_sample_factors is None:
        raise ValueError(
            "background_bin_factors and background_sample_factors must both be provided."
        )
    return additive + torch.matmul(
        _normalize_background_bin_factors_torch(background_bin_factors),
        background_sample_factors,
    )


def _matched_residual_scale(
    target_variance: torch.Tensor | np.ndarray | float,
    obs_likelihood: str,
    obs_df: float,
) -> torch.Tensor | np.ndarray | float:
    """Convert target variance to the scale parameter of the residual family."""
    obs_likelihood = _normalize_obs_likelihood_name(obs_likelihood)
    if obs_likelihood == NEGATIVE_BINOMIAL_OBS_LIKELIHOOD:
        raise ValueError("negative_binomial does not use a residual scale parameter.")
    if obs_likelihood == "normal":
        return target_variance ** 0.5
    if obs_likelihood == "laplace":
        return (target_variance / 2.0) ** 0.5
    if obs_df <= 2.0:
        raise ValueError("Student-t observation likelihood requires obs_df > 2.")
    return (target_variance * (obs_df - 2.0) / obs_df) ** 0.5


def _make_depth_distribution(
    expected: torch.Tensor,
    variance: torch.Tensor,
    obs_likelihood: str,
    obs_df: float,
) -> dist.Distribution:
    """Build the selected observation distribution with matched variance."""
    scale = _matched_residual_scale(variance, obs_likelihood, obs_df)
    if obs_likelihood == "normal":
        return dist.Normal(expected, scale)
    if obs_likelihood == "laplace":
        return dist.Laplace(expected, scale)
    return dist.StudentT(obs_df, loc=expected, scale=scale)


def _depth_log_lik_numpy(
    obs: np.ndarray,
    expected: np.ndarray,
    variance: np.ndarray,
    obs_likelihood: str,
    obs_df: float,
) -> np.ndarray:
    """Analytical residual log-likelihood matching the selected family."""
    from scipy import stats as sp_stats

    obs_likelihood = _normalize_obs_likelihood_name(obs_likelihood)
    if obs_likelihood == NEGATIVE_BINOMIAL_OBS_LIKELIHOOD:
        raise ValueError("negative_binomial uses the count likelihood helper.")
    safe_variance = np.maximum(variance, 1e-10)
    scale = _matched_residual_scale(safe_variance, obs_likelihood, obs_df)
    if obs_likelihood == "normal":
        return sp_stats.norm.logpdf(obs, loc=expected, scale=scale)
    if obs_likelihood == "laplace":
        return sp_stats.laplace.logpdf(obs, loc=expected, scale=scale)
    return sp_stats.t.logpdf((obs - expected) / scale, df=obs_df) - np.log(scale)


class CNVModel:
    """Hierarchical Bayesian model for whole-genome CN detection.

    The generative model uses per-bin CN-state priors, per-bin bias and
    variance, and per-sample variance. Autosomes can either use the historical
    Dirichlet-Categorical prior or a shrinkage prior that separates total
    non-reference mass from the distribution over alternative CN states. For
    the historical continuous observation families, observed normalised depth
    is drawn from a configurable residual model whose mean is
    ``CN × bin_bias + additive_background`` and whose variance scales with
    expected depth as
    ``(bin_var + sample_var) × max(CN × bin_bias + additive_background, 1e-3)``.

    When ``obs_likelihood='negative_binomial'``, the model instead consumes
    raw integer counts. It introduces a per-sample diploid-baseline
    depth scale ``sample_depth`` with an empirical-Bayes LogNormal prior
    centered on that sample's autosomal median counts per kb and scaled by
    a bootstrap estimate of the median's uncertainty, and uses the bin-length
    exposure in kilobases to model the count mean as

    ``bin_length_kb × sample_depth × (CN × bin_bias + additive_background) / 2``

    with negative-binomial variance

    ``μ + (bin_var + sample_var) × μ²``.

    When per-site allele data is provided, the model additionally includes
    a marginalized allele-fraction likelihood that jointly considers all
    possible genotype states at each SNP site, weighted by population
    allele frequencies and summed over observed sites per bin.

    Parameters
    ----------
    n_states : int
        Number of copy-number states (default 6 for CN 0–5).
    autosome_prior_mode : str
        ``'shrinkage'`` to use a hierarchical autosomal non-reference mass
        prior, or ``'dirichlet'`` for the historical autosomal
        Dirichlet-Categorical prior.
    alpha_ref : float
        Dirichlet concentration for the reference state (CN = 2) in autosomal
        ``dirichlet`` mode.
    alpha_non_ref : float
        Dirichlet concentration for non-reference states in autosomal
        ``dirichlet`` mode and for the alternative-state simplex in autosomal
        ``shrinkage`` mode.
    autosome_nonref_mean_alpha : float
        Beta prior alpha for the cohort-wide autosomal non-reference mean in
        ``shrinkage`` mode.
    autosome_nonref_mean_beta : float
        Beta prior beta for the cohort-wide autosomal non-reference mean in
        ``shrinkage`` mode.
    autosome_nonref_concentration : float
        Fixed concentration controlling how tightly each autosomal bin's
        non-reference prevalence is shrunk toward the learned cohort-wide mean
        in ``shrinkage`` mode.
    var_bias_bin : float
        Scale of the LogNormal prior on per-bin bias.
    var_sample : float
        Mean of the Exponential prior on per-sample variance.
    var_bin : float
        Scale of the Exponential prior on per-bin variance.
    epsilon_mean : float
        Prior mean of the global scale controlling the additive background
        depth term for each contig / sample pair. Defaults to ``1e-4``. Set
        to 0 to disable the epsilon floor.
    epsilon_concentration : float
        Gamma concentration for the additive background-depth prior. Values
        below 1 keep most bin/sample pairs near zero while retaining a heavy
        positive tail; ``1.0`` recovers the original Exponential prior.
    background_factors : int
        Number of low-rank additive background factors. Set to 0 to disable
        the structured background term entirely.
    background_sample_scale : float
        HalfNormal scale for the per-sample amplitudes of each background
        factor.
    background_bin_scale : float
        LogNormal scale for the raw per-bin loadings of each background
        factor before unit-mean normalization across bins.
    device, dtype : str, torch.dtype
        Torch device and floating-point type.
    guide_type : str
        ``'delta'`` for :class:`AutoDelta` or ``'diagonal'`` for
        :class:`AutoDiagonalNormal`, or ``'lowrank'`` for
        :class:`AutoLowRankMultivariateNormal`.
    af_concentration : float
        Beta-Binomial concentration for the allele-fraction model.
    af_weight : float
        Fixed global scale applied to the summed per-bin allele-fraction
        log-likelihood when ``learn_af_temperature`` is disabled; otherwise
        the prior median for the learned global AF temperature (0 to disable).
    learn_af_temperature : bool
        Whether to learn a single global AF temperature instead of keeping
        ``af_weight`` fixed.
    learn_site_pop_af : bool
        Whether to infer per-site population AFs inside the model using a
        Delta guide while keeping the requested guide family for the remaining
        latent variables.
    site_af_prior_strength : float
        Strength of the Beta prior on latent site AFs. The current
        ``site_pop_af`` values are used as prior means.
    af_temperature_prior_scale : float
        LogNormal prior scale for the learned AF temperature.
    alpha_sex_ref : float
        Dirichlet concentration for CN = 2 on sex chromosomes (default 1.0,
        flat).  The sex–CN coupling factor provides the main sex-dependent
        prior, so this can remain uninformative.
    alpha_sex_non_ref : float
        Dirichlet concentration for non-reference states on sex chromosomes.
    sex_prior : tuple
        Prior probabilities ``(P(XX), P(XY))`` for the per-sample sex
        karyotype latent variable.
    sex_cn_weight : float
        Weight of the per-bin sex–CN coupling factor.  Normalised by
        chromosome-type bin count so the total coupling per chromosome equals
        this value.  Set to 0 to disable the sex karyotype latent entirely.
    obs_likelihood : str
        Observation family for depth residuals: ``normal``, ``studentt``, or
        ``laplace`` for normalized depth, or ``negative_binomial`` for raw
        counts.
    obs_df : float
        Degrees of freedom for the Student-t observation likelihood.
    sample_depth_max : float
        Upper clip applied when deriving the empirical-Bayes prior center and
        guide initialisation for per-sample depth scale in the
        negative-binomial observation model.
    """

    _latent_sites = ["bin_bias", "sample_var", "bin_var"]

    def __init__(
        self,
        n_states: int = 6,
        autosome_prior_mode: str = "dirichlet",
        alpha_ref: float = 50.0,
        alpha_non_ref: float = 1.0,
        autosome_nonref_mean_alpha: float = 1.0,
        autosome_nonref_mean_beta: float = 19.0,
        autosome_nonref_concentration: float = 20.0,
        var_bias_bin: float = 0.02,
        var_sample: float = 0.00025,
        var_bin: float = 0.001,
        epsilon_mean: float = DEFAULT_EPSILON_MEAN,
        epsilon_concentration: float = DEFAULT_EPSILON_CONCENTRATION,
        background_factors: int = DEFAULT_BACKGROUND_FACTORS,
        background_sample_scale: float = DEFAULT_BACKGROUND_SAMPLE_SCALE,
        background_bin_scale: float = DEFAULT_BACKGROUND_BIN_SCALE,
        device: str = "cpu",
        dtype: torch.dtype = torch.float32,
        guide_type: str = "delta",
        af_concentration: float = DEFAULT_AF_CONCENTRATION,
        af_weight: float = DEFAULT_AF_WEIGHT,
        af_evidence_mode: str = "relative",
        learn_af_temperature: bool = True,
        learn_site_pop_af: bool = False,
        site_af_prior_strength: float = 20.0,
        af_temperature_prior_scale: float = 0.5,
        alpha_sex_ref: float = 1.0,
        alpha_sex_non_ref: float = 1.0,
        sex_prior: tuple = (0.5, 0.5),
        sex_cn_weight: float = 3.0,
        obs_likelihood: str = "normal",
        obs_df: float = 3.5,
        sample_depth_max: float = 10000.0,
        freeze_bin_bias: bool = False,
        freeze_cn_prior: bool = False,
        freeze_sample_var: bool = False,
        freeze_sample_depth: bool = False,
    ) -> None:
        self.n_states = n_states
        self.autosome_prior_mode = autosome_prior_mode
        self.alpha_ref = alpha_ref
        self.alpha_non_ref = alpha_non_ref
        self.autosome_nonref_mean_alpha = autosome_nonref_mean_alpha
        self.autosome_nonref_mean_beta = autosome_nonref_mean_beta
        self.autosome_nonref_concentration = autosome_nonref_concentration
        self.var_bias_bin = var_bias_bin
        self.var_sample = var_sample
        self.var_bin = var_bin
        self.epsilon_mean = epsilon_mean
        self.epsilon_concentration = epsilon_concentration
        self.background_factors = background_factors
        self.background_sample_scale = background_sample_scale
        self.background_bin_scale = background_bin_scale
        self.device = device
        self.dtype = dtype
        self.guide_type = guide_type
        self.af_concentration = af_concentration
        self.af_weight = af_weight
        self.af_evidence_mode = _normalize_af_evidence_mode(af_evidence_mode)
        self.learn_af_temperature = learn_af_temperature
        self.learn_site_pop_af = learn_site_pop_af
        self.site_af_prior_strength = site_af_prior_strength
        self.af_temperature_prior_scale = af_temperature_prior_scale
        self.alpha_sex_ref = alpha_sex_ref
        self.alpha_sex_non_ref = alpha_sex_non_ref
        self.sex_prior = sex_prior
        self.sex_cn_weight = sex_cn_weight
        self.obs_likelihood = _normalize_obs_likelihood_name(obs_likelihood)
        self.obs_df = obs_df
        self.sample_depth_max = sample_depth_max
        self.freeze_bin_bias = bool(freeze_bin_bias)
        self.freeze_cn_prior = bool(freeze_cn_prior)
        self.freeze_sample_var = bool(freeze_sample_var)
        self.freeze_sample_depth = bool(freeze_sample_depth)

        if self.autosome_prior_mode not in AUTOSOME_PRIOR_MODES:
            raise ValueError(
                f"Unknown autosome_prior_mode: {self.autosome_prior_mode!r}. "
                f"Choose one of {AUTOSOME_PRIOR_MODES}."
            )
        if self.freeze_cn_prior and self.autosome_prior_mode == "shrinkage":
            raise ValueError(
                "freeze_cn_prior is only supported with autosome_prior_mode='dirichlet'."
            )
        if self.autosome_nonref_mean_alpha <= 0 or self.autosome_nonref_mean_beta <= 0:
            raise ValueError(
                "autosome_nonref_mean_alpha and autosome_nonref_mean_beta must be positive."
            )
        if self.autosome_nonref_concentration <= 0:
            raise ValueError("autosome_nonref_concentration must be positive.")
        if self.af_temperature_prior_scale <= 0:
            raise ValueError("af_temperature_prior_scale must be positive.")
        if self.learn_af_temperature and self.af_weight <= 0:
            raise ValueError("learn_af_temperature requires af_weight > 0.")
        if self.learn_site_pop_af and self.af_weight <= 0:
            raise ValueError("learn_site_pop_af requires af_weight > 0.")
        if self.site_af_prior_strength < 0:
            raise ValueError("site_af_prior_strength must be non-negative.")

        if self.epsilon_mean < 0:
            raise ValueError("epsilon_mean must be non-negative.")
        if self.epsilon_concentration <= 0:
            raise ValueError("epsilon_concentration must be positive.")
        if int(self.background_factors) != self.background_factors:
            raise ValueError("background_factors must be an integer.")
        self.background_factors = int(self.background_factors)
        if self.background_factors < 0:
            raise ValueError("background_factors must be non-negative.")
        if self.background_sample_scale <= 0:
            raise ValueError("background_sample_scale must be positive.")
        if self.background_bin_scale <= 0:
            raise ValueError("background_bin_scale must be positive.")

        if self.obs_likelihood == "studentt" and self.obs_df <= 2.0:
            raise ValueError("Student-t observation likelihood requires obs_df > 2.")
        if self.obs_likelihood == NEGATIVE_BINOMIAL_OBS_LIKELIHOOD and self.sample_depth_max <= 0:
            raise ValueError("sample_depth_max must be positive.")

        self._latent_sites = list(self._latent_sites)
        if self.freeze_bin_bias:
            self._latent_sites = [
                site for site in self._latent_sites if site != "bin_bias"
            ]
        if self.freeze_sample_var:
            self._latent_sites = [
                site for site in self._latent_sites if site != "sample_var"
            ]
        if self.autosome_prior_mode == "shrinkage":
            self._latent_sites.extend(
                [
                    "autosome_nonref_mean",
                    "autosome_nonref_prob",
                    "autosome_alt_cn_probs",
                    "sex_cn_probs",
                ]
            )
        else:
            if not self.freeze_cn_prior:
                self._latent_sites.append("cn_probs")
        if all(
            (
                self.obs_likelihood == NEGATIVE_BINOMIAL_OBS_LIKELIHOOD,
                not self.freeze_sample_depth,
            )
        ):
            self._latent_sites.append("sample_depth")
        if self.epsilon_mean > 0:
            self._latent_sites.append("bin_epsilon")
        if self.background_factors > 0:
            self._latent_sites.extend(
                ["background_bin_factors", "background_sample_factors"]
            )
        if self.learn_af_temperature:
            self._latent_sites.append("af_temperature")
        if self.learn_site_pop_af:
            self._latent_sites.append("site_pop_af_latent")

        self.ref_state = 2
        self.nonref_state_indices = tuple(
            idx for idx in range(self.n_states) if idx != self.ref_state
        )
        self.nonref_state_indices_torch = torch.tensor(
            self.nonref_state_indices,
            dtype=torch.long,
            device=device,
        )
        self.nonref_state_indices_numpy = np.asarray(
            self.nonref_state_indices,
            dtype=np.int64,
        )

        # Per-chromosome-type Dirichlet alpha table: (3, n_states)
        # Row 0: autosomes — strong CN=2 prior
        # Row 1: chrX — flat prior (sex_cn coupling handles sex-dependent CN)
        # Row 2: chrY — flat prior
        alpha_auto = [alpha_non_ref] * n_states
        alpha_auto[2] = alpha_ref
        alpha_sex = [alpha_sex_non_ref] * n_states
        alpha_sex[2] = alpha_sex_ref
        self.alpha_table = torch.tensor(
            [alpha_auto, alpha_sex, alpha_sex], dtype=dtype, device=device,
        )

        # Sex-CN coupling score table: (2, 3, n_states)
        # Indexed by [sex_state, chr_type, cn_state].
        # 0 for the expected CN state; -1 penalty for unexpected states.
        # Autosome (chr_type=0): all zeros — no sex effect.
        sex_cn_score = torch.zeros(2, 3, n_states, dtype=dtype, device=device)
        # chrX (type 1): XX→CN 2, XY→CN 1
        sex_cn_score[0, 1, :] = -1.0
        sex_cn_score[0, 1, 2] = 0.0
        sex_cn_score[1, 1, :] = -1.0
        sex_cn_score[1, 1, 1] = 0.0
        # chrY (type 2): XX→CN 0, XY→CN 1
        sex_cn_score[0, 2, :] = -1.0
        sex_cn_score[0, 2, 0] = 0.0
        sex_cn_score[1, 2, :] = -1.0
        sex_cn_score[1, 2, 1] = 0.0
        self.sex_cn_score = sex_cn_score

        self.loss_history: Dict[str, List[float]] = {"epoch": [], "elbo": []}

        self.guide = self._build_guide()

    def _build_base_guide(self, blocked_model, guide_kwargs: dict):
        """Construct the requested autoguide for a blocked model."""
        if self.guide_type == "delta":
            return AutoDelta(blocked_model, **guide_kwargs)
        if self.guide_type == "diagonal":
            return AutoDiagonalNormal(blocked_model, **guide_kwargs)
        if self.guide_type == "lowrank":
            return AutoLowRankMultivariateNormal(blocked_model, **guide_kwargs)
        raise ValueError(
            f"Unknown guide_type: {self.guide_type!r}. "
            "Choose 'delta', 'diagonal', or 'lowrank'."
        )

    def _build_guide(self, init_loc_fn=None):
        """Construct the variational guide, optionally with custom init."""
        guide_kwargs = {}
        if init_loc_fn is not None:
            guide_kwargs["init_loc_fn"] = init_loc_fn

        if self.learn_site_pop_af:
            guide = AutoGuideList(self.model)
            guide.append(
                AutoDelta(
                    poutine.block(self.model, expose=["site_pop_af_latent"]),
                    **guide_kwargs,
                )
            )
            remaining_latents = [
                site for site in self._latent_sites if site != "site_pop_af_latent"
            ]
            if remaining_latents:
                guide.append(
                    self._build_base_guide(
                        poutine.block(self.model, expose=remaining_latents),
                        guide_kwargs,
                    )
                )
            return guide

        blocked = poutine.block(self.model, expose=self._latent_sites)
        return self._build_base_guide(blocked, guide_kwargs)

    def _compose_autosome_cn_probs_torch(
        self,
        autosome_nonref_prob: torch.Tensor,
        autosome_alt_cn_probs: torch.Tensor,
    ) -> torch.Tensor:
        """Expand autosomal non-reference mass and alternative simplex."""
        cn_probs = torch.zeros(
            autosome_alt_cn_probs.shape[:-1] + (self.n_states,),
            dtype=autosome_alt_cn_probs.dtype,
            device=autosome_alt_cn_probs.device,
        )
        autosome_nonref_prob = torch.clamp(
            autosome_nonref_prob,
            min=1e-6,
            max=1.0 - 1e-6,
        )
        cn_probs[..., self.ref_state] = 1.0 - autosome_nonref_prob
        cn_probs[..., self.nonref_state_indices_torch] = (
            autosome_nonref_prob.unsqueeze(-1) * autosome_alt_cn_probs
        )
        return cn_probs

    def _compose_autosome_cn_probs_numpy(
        self,
        autosome_nonref_prob: np.ndarray,
        autosome_alt_cn_probs: np.ndarray,
    ) -> np.ndarray:
        """NumPy version of autosomal non-reference/simplex expansion."""
        autosome_nonref_prob = np.asarray(autosome_nonref_prob, dtype=np.float64)
        autosome_alt_cn_probs = np.asarray(autosome_alt_cn_probs, dtype=np.float64)
        autosome_nonref_prob = np.clip(autosome_nonref_prob, 1e-6, 1.0 - 1e-6)
        cn_probs = np.zeros(
            autosome_alt_cn_probs.shape[:-1] + (self.n_states,),
            dtype=autosome_alt_cn_probs.dtype,
        )
        cn_probs[..., self.ref_state] = 1.0 - autosome_nonref_prob
        cn_probs[..., self.nonref_state_indices_numpy] = (
            autosome_nonref_prob[..., np.newaxis] * autosome_alt_cn_probs
        )
        return cn_probs

    def _default_cn_probs_torch(
        self,
        chr_type: torch.Tensor | None,
        n_bins: int,
    ) -> torch.Tensor:
        """Return the historical Dirichlet mean simplex for each bin."""
        if chr_type is not None:
            alpha = self.alpha_table[chr_type]
        else:
            alpha = torch.full(
                (n_bins, self.n_states),
                self.alpha_non_ref,
                dtype=self.dtype,
                device=self.device,
            )
            alpha[:, self.ref_state] = self.alpha_ref
        return alpha / alpha.sum(dim=-1, keepdim=True)

    def _af_reference_cn_probs_torch(
        self,
        chr_type: torch.Tensor | None,
        n_bins: int,
    ) -> torch.Tensor:
        """Build the fixed CN reference mixture used for relative AF evidence."""
        if self.autosome_prior_mode != "shrinkage":
            return self._default_cn_probs_torch(chr_type, n_bins)

        autosome_nonref_mean = self.autosome_nonref_mean_alpha / (
            self.autosome_nonref_mean_alpha + self.autosome_nonref_mean_beta
        )
        autosome_nonref_prob = torch.full(
            (n_bins,),
            autosome_nonref_mean,
            dtype=self.dtype,
            device=self.device,
        )
        autosome_alt_cn_probs = torch.full(
            (n_bins, self.n_states - 1),
            1.0 / (self.n_states - 1),
            dtype=self.dtype,
            device=self.device,
        )
        autosome_probs = self._compose_autosome_cn_probs_torch(
            autosome_nonref_prob,
            autosome_alt_cn_probs,
        )

        if chr_type is None:
            return autosome_probs

        cn_probs = self._default_cn_probs_torch(chr_type, n_bins)
        autosome_mask = chr_type == 0
        if autosome_mask.any():
            cn_probs[autosome_mask] = autosome_probs[autosome_mask]
        return cn_probs

    def _build_cn_probs_from_estimates_torch(
        self,
        source: Dict[str, torch.Tensor],
        chr_type: torch.Tensor | None,
    ) -> torch.Tensor:
        """Build per-bin CN simplexes from sampled or derived latent values."""
        if "cn_probs" in source:
            return source["cn_probs"]

        n_bins = int(chr_type.shape[0]) if chr_type is not None else int(source["bin_bias"].shape[0])
        if self.autosome_prior_mode != "shrinkage":
            return self._default_cn_probs_torch(chr_type, n_bins)

        if chr_type is None:
            return self._compose_autosome_cn_probs_torch(
                source["autosome_nonref_prob"],
                source["autosome_alt_cn_probs"],
            )

        cn_probs = torch.empty(
            (n_bins, self.n_states),
            dtype=self.dtype,
            device=self.device,
        )
        autosome_mask = chr_type == 0
        if autosome_mask.any():
            cn_probs[autosome_mask] = self._compose_autosome_cn_probs_torch(
                source["autosome_nonref_prob"],
                source["autosome_alt_cn_probs"],
            )
        if (~autosome_mask).any():
            if "sex_cn_probs" in source and source["sex_cn_probs"].numel() > 0:
                cn_probs[~autosome_mask] = source["sex_cn_probs"]
            else:
                cn_probs[~autosome_mask] = self._default_cn_probs_torch(
                    chr_type[~autosome_mask],
                    int((~autosome_mask).sum().item()),
                )
        return cn_probs

    def _build_cn_probs_from_estimates_numpy(
        self,
        maps: Dict[str, np.ndarray],
        chr_type: np.ndarray | None,
        n_bins: int,
    ) -> np.ndarray:
        """NumPy mirror of per-bin CN simplex reconstruction."""
        if "cn_probs" in maps:
            cn_probs = np.asarray(maps["cn_probs"]).squeeze()
            if cn_probs.ndim == 1:
                cn_probs = cn_probs.reshape(1, -1)
            return cn_probs

        if self.autosome_prior_mode != "shrinkage" or (
            "autosome_nonref_prob" not in maps or "autosome_alt_cn_probs" not in maps
        ):
            if chr_type is not None:
                alpha = self.alpha_table.detach().cpu().numpy()[chr_type]
            else:
                alpha = np.full(
                    (n_bins, self.n_states),
                    self.alpha_non_ref,
                    dtype=np.float64,
                )
                alpha[:, self.ref_state] = self.alpha_ref
            return alpha / alpha.sum(axis=1, keepdims=True)

        autosome_nonref_prob = np.atleast_1d(
            np.asarray(maps["autosome_nonref_prob"]).squeeze()
        )
        autosome_alt_cn_probs = np.asarray(maps["autosome_alt_cn_probs"]).squeeze()
        if autosome_alt_cn_probs.ndim == 1:
            autosome_alt_cn_probs = autosome_alt_cn_probs.reshape(1, -1)

        if chr_type is None:
            return self._compose_autosome_cn_probs_numpy(
                autosome_nonref_prob,
                autosome_alt_cn_probs,
            )

        cn_probs = np.empty((n_bins, self.n_states), dtype=np.float64)
        autosome_mask = chr_type == 0
        if np.any(autosome_mask):
            cn_probs[autosome_mask] = self._compose_autosome_cn_probs_numpy(
                autosome_nonref_prob,
                autosome_alt_cn_probs,
            )
        if np.any(~autosome_mask):
            if "sex_cn_probs" in maps:
                sex_cn_probs = np.asarray(maps["sex_cn_probs"]).squeeze()
                if sex_cn_probs.ndim == 1:
                    sex_cn_probs = sex_cn_probs.reshape(1, -1)
                cn_probs[~autosome_mask] = sex_cn_probs
            else:
                sex_alpha = self.alpha_table.detach().cpu().numpy()[chr_type[~autosome_mask]]
                cn_probs[~autosome_mask] = sex_alpha / sex_alpha.sum(axis=1, keepdims=True)
        return cn_probs

    def _af_reference_cn_probs_numpy(
        self,
        chr_type: np.ndarray | None,
        n_bins: int,
    ) -> np.ndarray:
        """NumPy counterpart of :meth:`_af_reference_cn_probs_torch`."""
        if self.autosome_prior_mode != "shrinkage":
            alpha = self.alpha_table.detach().cpu().numpy()
            if chr_type is not None:
                selected = alpha[chr_type]
            else:
                selected = np.full(
                    (n_bins, self.n_states),
                    self.alpha_non_ref,
                    dtype=np.float64,
                )
                selected[:, self.ref_state] = self.alpha_ref
            return selected / selected.sum(axis=1, keepdims=True)

        autosome_nonref_mean = self.autosome_nonref_mean_alpha / (
            self.autosome_nonref_mean_alpha + self.autosome_nonref_mean_beta
        )
        autosome_nonref_prob = np.full((n_bins,), autosome_nonref_mean, dtype=np.float64)
        autosome_alt_cn_probs = np.full(
            (n_bins, self.n_states - 1),
            1.0 / (self.n_states - 1),
            dtype=np.float64,
        )
        autosome_probs = self._compose_autosome_cn_probs_numpy(
            autosome_nonref_prob,
            autosome_alt_cn_probs,
        )

        if chr_type is None:
            return autosome_probs

        alpha = self.alpha_table.detach().cpu().numpy()[chr_type]
        cn_probs = alpha / alpha.sum(axis=1, keepdims=True)
        autosome_mask = chr_type == 0
        if np.any(autosome_mask):
            cn_probs[autosome_mask] = autosome_probs[autosome_mask]
        return cn_probs

    def _apply_af_evidence_mode_torch(
        self,
        af_table: torch.Tensor,
        chr_type: torch.Tensor | None,
    ) -> torch.Tensor:
        """Apply the configured AF evidence transformation to a precomputed table."""
        if self.af_evidence_mode != "relative":
            return af_table
        reference_cn_probs = self._af_reference_cn_probs_torch(
            chr_type,
            af_table.shape[1],
        )
        return _center_af_table_torch(af_table, reference_cn_probs)

    def _apply_af_evidence_mode_numpy(
        self,
        af_table: np.ndarray,
        chr_type: np.ndarray | None,
    ) -> np.ndarray:
        """NumPy counterpart of :meth:`_apply_af_evidence_mode_torch`."""
        if self.af_evidence_mode != "relative":
            return af_table
        reference_cn_probs = self._af_reference_cn_probs_numpy(
            chr_type,
            af_table.shape[1],
        )
        return _center_af_table_numpy(af_table, reference_cn_probs)

    def _prepare_af_table_torch(
        self,
        site_alt: torch.Tensor,
        site_total: torch.Tensor,
        site_pop_af: torch.Tensor,
        site_mask: torch.Tensor,
        chr_type: torch.Tensor | None,
        auto_leave_one_out: bool = True,
    ) -> torch.Tensor:
        """Precompute and transform the AF table according to the configured evidence mode."""
        if auto_leave_one_out:
            site_pop_af, used_leave_one_out = _resolve_fixed_site_pop_af_torch(
                site_alt,
                site_total,
                site_pop_af,
                site_mask,
            )
            if used_leave_one_out:
                logger.info(
                    "AF table: detected self-pooled site_pop_af; using leave-one-out population AFs per sample."
                )
        af_table = _precompute_af_table(
            site_alt,
            site_total,
            site_pop_af,
            site_mask,
            n_states=self.n_states,
            concentration=self.af_concentration,
        )
        return self._apply_af_evidence_mode_torch(af_table, chr_type)

    def _af_scale_torch(
        self,
        source: Dict[str, torch.Tensor] | None = None,
    ) -> torch.Tensor:
        """Return the fixed or learned AF scale used to temper summed AF evidence."""
        if source is not None and "af_temperature" in source:
            return source["af_temperature"]
        return torch.tensor(self.af_weight, dtype=self.dtype, device=self.device)

    def _af_scale_numpy(
        self,
        maps: Dict[str, np.ndarray] | None = None,
    ) -> float:
        """NumPy counterpart of :meth:`_af_scale_torch`."""
        if maps is not None and "af_temperature" in maps:
            return float(np.asarray(maps["af_temperature"]).item())
        return float(self.af_weight)

    def _get_continuous_estimates(
        self,
        data: DepthData,
        estimate_method: str = "current",
        model_kw: dict | None = None,
    ) -> Dict[str, torch.Tensor]:
        """Return continuous latent estimates from the fitted guide.

        ``estimate_method='current'`` preserves the historical behavior of
        taking one guide draw. ``estimate_method='median'`` plugs in the guide
        medians, which is deterministic for non-delta guides.
        """
        model_kw = self._model_kwargs(data) if model_kw is None else model_kw

        if estimate_method == "current":
            guide_trace = poutine.trace(self.guide).get_trace(**model_kw)
            source = {
                site: guide_trace.nodes[site]["value"] for site in self._latent_sites
            }
        elif estimate_method == "median":
            guide_median = self.guide.median(**model_kw)
            source = {site: guide_median[site] for site in self._latent_sites}
        else:
            raise ValueError(
                f"Unknown estimate_method: {estimate_method!r}. "
                "Choose 'current' or 'median'."
            )

        estimates: Dict[str, torch.Tensor] = {}
        for site, value in source.items():
            if isinstance(value, torch.Tensor):
                estimates[site] = value.detach().clone()
            else:
                estimates[site] = torch.as_tensor(value, device=self.device)
        if all(
            (
                self.freeze_bin_bias,
                "bin_bias" not in estimates,
                model_kw.get("bin_bias_fixed") is not None,
            )
        ):
            estimates["bin_bias"] = model_kw["bin_bias_fixed"].detach().clone()
        if all(
            (
                self.freeze_sample_var,
                "sample_var" not in estimates,
                model_kw.get("sample_var_fixed") is not None,
            )
        ):
            estimates["sample_var"] = (
                model_kw["sample_var_fixed"].detach().clone()
            )
        if "cn_probs" not in estimates:
            estimates["cn_probs"] = self._build_cn_probs_from_estimates_torch(
                estimates,
                model_kw.get("chr_type"),
            )
        if all(
            (
                self.obs_likelihood == NEGATIVE_BINOMIAL_OBS_LIKELIHOOD,
                self.freeze_sample_depth,
                "sample_depth" not in estimates,
                model_kw.get("sample_depth_fixed") is not None,
            )
        ):
            estimates["sample_depth"] = (
                model_kw["sample_depth_fixed"].detach().clone()
            )
        return estimates

    def _infer_discrete_assignments(
        self,
        data: DepthData,
        continuous_estimates: Dict[str, torch.Tensor],
        model_kw: dict | None = None,
    ) -> Dict[str, np.ndarray]:
        """Infer discrete latents conditioned on fixed continuous estimates."""
        model_kw = self._model_kwargs(data) if model_kw is None else model_kw
        has_sex = (
            self.sex_cn_weight > 0 and
            "chr_type" in model_kw and
            model_kw["chr_type"] is not None
        )
        first_dim = -3 - int(has_sex)

        conditioned = poutine.condition(
            self.model,
            data={
                site: continuous_estimates[site]
                for site in self._latent_sites
                if site in continuous_estimates
            },
        )
        inferred = infer_discrete(
            conditioned, temperature=0, first_available_dim=first_dim,
        )
        trace = poutine.trace(inferred).get_trace(**model_kw)

        assignments = {
            "cn": trace.nodes["cn"]["value"].detach().cpu().numpy(),
        }
        if has_sex and "sex_karyotype" in trace.nodes:
            assignments["sex_karyotype"] = (
                trace.nodes["sex_karyotype"]["value"].detach().cpu().numpy()
            )
        return assignments

    def _estimate_sample_depth_init(self, data: DepthData) -> torch.Tensor:
        """Estimate neutral-baseline counts-per-kb from autosomal bins."""
        autosome_median_per_kb = self._estimate_sample_depth_center(data)
        return torch.tensor(
            autosome_median_per_kb,
            dtype=self.dtype,
            device=self.device,
        )

    def _get_autosomal_depth_per_kb(self, data: DepthData) -> np.ndarray:
        """Return per-sample autosomal counts-per-kb matrix."""
        autosome_mask = data.chr_type.detach().cpu().numpy() == 0
        if not np.any(autosome_mask):
            autosome_mask = np.ones(data.n_bins, dtype=bool)

        depth_np = data.depth.detach().cpu().numpy()[autosome_mask, :]
        bin_length_kb_np = data.bin_length_kb.detach().cpu().numpy()[autosome_mask]
        return np.divide(
            depth_np,
            bin_length_kb_np[:, np.newaxis],
            out=np.zeros_like(depth_np),
            where=bin_length_kb_np[:, np.newaxis] > 0,
        )

    def _estimate_sample_depth_center(self, data: DepthData) -> np.ndarray:
        """Estimate per-sample autosomal median counts-per-kb."""
        per_kb = self._get_autosomal_depth_per_kb(data)
        center_np = np.median(per_kb, axis=0)
        return np.clip(center_np, 1e-3, self.sample_depth_max * 0.95)

    def _estimate_sample_depth_bootstrap_sd(
        self,
        per_kb: np.ndarray,
        center_np: np.ndarray,
    ) -> np.ndarray:
        """Estimate per-sample uncertainty of the autosomal median via bootstrap."""
        floor_np = np.maximum(_SAMPLE_DEPTH_PRIOR_SD_FLOOR_FRAC * center_np, 1e-3)
        n_bins, n_samples = per_kb.shape
        if n_bins <= 1:
            return floor_np

        rng = np.random.RandomState(_SAMPLE_DEPTH_PRIOR_BOOTSTRAP_SEED)
        bootstrap_indices = rng.randint(
            0,
            n_bins,
            size=(_SAMPLE_DEPTH_PRIOR_BOOTSTRAP_DRAWS, n_bins),
        )
        raw_sd_np = np.empty(n_samples, dtype=np.float64)
        for sample_idx in range(n_samples):
            sample_values = per_kb[:, sample_idx]
            bootstrap_medians = np.median(sample_values[bootstrap_indices], axis=1)
            raw_sd_np[sample_idx] = np.std(bootstrap_medians, ddof=1)
        return np.maximum(raw_sd_np, floor_np)

    def _estimate_sample_depth_prior(self, data: DepthData) -> Dict[str, torch.Tensor]:
        """Estimate empirical-Bayes LogNormal prior parameters for sample depth."""
        per_kb = self._get_autosomal_depth_per_kb(data)
        center_np = np.clip(
            np.median(per_kb, axis=0),
            1e-3,
            self.sample_depth_max * 0.95,
        )
        raw_sd_np = self._estimate_sample_depth_bootstrap_sd(per_kb, center_np)

        sigma2_np = np.log1p((raw_sd_np / np.maximum(center_np, 1e-6)) ** 2)
        scale_np = np.sqrt(np.maximum(sigma2_np, 1e-12))
        loc_np = np.log(np.maximum(center_np, 1e-6)) - 0.5 * sigma2_np
        return {
            "center": torch.tensor(center_np, dtype=self.dtype, device=self.device),
            "raw_sd": torch.tensor(raw_sd_np, dtype=self.dtype, device=self.device),
            "loc": torch.tensor(loc_np, dtype=self.dtype, device=self.device),
            "scale": torch.tensor(scale_np, dtype=self.dtype, device=self.device),
        }

    def _make_init_loc_fn(self, data: DepthData, *, randomize: bool = False):
        """Build a guide init function anchored to raw autosomal depth."""
        if not randomize and self.obs_likelihood != NEGATIVE_BINOMIAL_OBS_LIKELIHOOD:
            return None

        base_init = init_to_sample if randomize else init_to_median(num_samples=50)

        if self.obs_likelihood != NEGATIVE_BINOMIAL_OBS_LIKELIHOOD:
            return base_init

        sample_depth_prior = self._estimate_sample_depth_prior(data)
        sample_depth_init = sample_depth_prior["center"]
        if not randomize:
            logger.info(
                "Initial sample_depth from autosomal counts/kb: median=%.3f range=[%.3f, %.3f]",
                float(sample_depth_init.median().item()),
                float(sample_depth_init.min().item()),
                float(sample_depth_init.max().item()),
            )
            logger.info(
                "Sample_depth LogNormal prior from autosomal medians: bootstrap raw_sd median=%.3f range=[%.3f, %.3f], log_scale median=%.4f range=[%.4f, %.4f]",
                float(sample_depth_prior["raw_sd"].median().item()),
                float(sample_depth_prior["raw_sd"].min().item()),
                float(sample_depth_prior["raw_sd"].max().item()),
                float(sample_depth_prior["scale"].median().item()),
                float(sample_depth_prior["scale"].min().item()),
                float(sample_depth_prior["scale"].max().item()),
            )

        def init_loc_fn(site):
            if site["name"] == "sample_depth":
                if randomize:
                    sampled_depth = dist.LogNormal(
                        sample_depth_prior["loc"],
                        sample_depth_prior["scale"],
                    ).sample()
                    return torch.clamp(
                        sampled_depth,
                        min=1e-3,
                        max=self.sample_depth_max * 0.95,
                    )
                return sample_depth_init
            if site["name"] == "site_pop_af_latent" and data.site_pop_af is not None:
                return torch.clamp(data.site_pop_af, min=1e-6, max=1.0 - 1e-6)
            return base_init(site)

        return init_loc_fn

    @staticmethod
    def _capture_rng_state() -> dict[str, object]:
        """Capture Python, NumPy, and Torch RNG state for restart search."""
        state: dict[str, object] = {
            "python": random.getstate(),
            "numpy": np.random.get_state(),
            "torch": torch.random.get_rng_state(),
        }
        if torch.cuda.is_available():
            state["torch_cuda"] = torch.cuda.get_rng_state_all()
        return state

    @staticmethod
    def _restore_rng_state(state: dict[str, object]) -> None:
        """Restore RNG state after Monte Carlo guide initialization."""
        random.setstate(state["python"])
        np.random.set_state(state["numpy"])
        torch.random.set_rng_state(state["torch"])
        torch_cuda_state = state.get("torch_cuda")
        if torch_cuda_state is not None and torch.cuda.is_available():
            torch.cuda.set_rng_state_all(torch_cuda_state)

    def _select_svi_initialization(
        self,
        data: DepthData,
        *,
        model_kw: dict,
        elbo,
        init_restarts: int,
    ) -> None:
        """Choose the best initial guide state before optimization."""
        if init_restarts < 1:
            raise ValueError("init_restarts must be at least 1.")

        default_init_loc_fn = self._make_init_loc_fn(data)
        if init_restarts == 1:
            pyro.clear_param_store()
            self.guide = self._build_guide(default_init_loc_fn)
            return

        logger.info(
            "Evaluating %d SVI initializations before gradient descent (%d anchored, %d random).",
            init_restarts,
            1,
            init_restarts - 1,
        )

        rng_state = self._capture_rng_state()
        seed_modulus = np.iinfo(np.uint32).max + 1
        base_seed = int(torch.initial_seed()) % seed_modulus
        best_loss = float("inf")
        best_restart_idx: Optional[int] = None
        best_state = None
        best_guide = None
        candidate_losses: list[float] = []

        try:
            for restart_idx in range(init_restarts):
                pyro.clear_param_store()
                if restart_idx == 0:
                    init_loc_fn = default_init_loc_fn
                else:
                    pyro.set_rng_seed((base_seed + restart_idx) % seed_modulus)
                    init_loc_fn = self._make_init_loc_fn(data, randomize=True)

                guide = self._build_guide(init_loc_fn)
                try:
                    loss = float(elbo.loss(self.model, guide, **model_kw))
                except Exception as exc:
                    logger.warning(
                        "SVI init candidate %d/%d failed during loss evaluation: %s",
                        restart_idx + 1,
                        init_restarts,
                        exc,
                    )
                    continue

                if not np.isfinite(loss):
                    logger.warning(
                        "SVI init candidate %d/%d produced non-finite loss %s and will be skipped.",
                        restart_idx + 1,
                        init_restarts,
                        loss,
                    )
                    continue

                candidate_losses.append(loss)
                if loss < best_loss:
                    best_loss = loss
                    best_restart_idx = restart_idx
                    best_state = copy.deepcopy(pyro.get_param_store().get_state())
                    best_guide = guide
        finally:
            self._restore_rng_state(rng_state)

        pyro.clear_param_store()
        if best_state is None or best_guide is None or best_restart_idx is None:
            logger.warning(
                "All SVI initialization candidates failed; falling back to the anchored initialization."
            )
            self.guide = self._build_guide(default_init_loc_fn)
            return

        pyro.get_param_store().set_state(best_state)
        self.guide = best_guide

        loss_arr = np.asarray(candidate_losses, dtype=np.float64)
        logger.info(
            "Selected SVI init candidate %d/%d with initial loss %.4f (median=%.4f, p10=%.4f, p90=%.4f).",
            best_restart_idx + 1,
            init_restarts,
            best_loss,
            float(np.median(loss_arr)),
            float(np.quantile(loss_arr, 0.10)),
            float(np.quantile(loss_arr, 0.90)),
        )

    # ── probabilistic model ─────────────────────────────────────────────

    @config_enumerate(default="parallel")
    def model(
        self,
        depth: torch.Tensor,
        n_bins: int = None,
        n_samples: int = None,
        bin_length_kb: Optional[torch.Tensor] = None,
        bin_bias_fixed: Optional[torch.Tensor] = None,
        sample_var_fixed: Optional[torch.Tensor] = None,
        sample_depth_prior_loc: Optional[torch.Tensor] = None,
        sample_depth_prior_scale: Optional[torch.Tensor] = None,
        sample_depth_fixed: Optional[torch.Tensor] = None,
        n_contigs: Optional[int] = None,
        contig_index: Optional[torch.Tensor] = None,
        chr_type: Optional[torch.Tensor] = None,
        sex_cn_weight_per_bin: Optional[torch.Tensor] = None,
        site_alt: Optional[torch.Tensor] = None,
        site_total: Optional[torch.Tensor] = None,
        site_pop_af: Optional[torch.Tensor] = None,
        site_mask: Optional[torch.Tensor] = None,
        af_table: Optional[torch.Tensor] = None,
    ) -> None:
        """Pyro generative model.

        Args:
            depth: Observed normalized depth or raw-count tensor
                (n_bins × n_samples).
            n_bins: Number of genomic bins.
            n_samples: Number of samples.
            bin_length_kb: Per-bin exposure in kilobases. Required when
                ``obs_likelihood='negative_binomial'``.
            bin_bias_fixed: Fixed per-bin bias values used when
                ``freeze_bin_bias=True``.
            sample_var_fixed: Fixed per-sample overdispersion values used when
                ``freeze_sample_var=True``.
            sample_depth_prior_loc: LogNormal location parameter for the
                per-sample diploid baseline depth scale.
            sample_depth_prior_scale: LogNormal scale parameter for the
                per-sample diploid baseline depth scale.
            n_contigs: Number of unique contigs represented by the bins.
            contig_index: Per-bin index into the contig axis.
            chr_type: Per-bin chromosome type (0=autosome, 1=chrX, 2=chrY).
            sex_cn_weight_per_bin: Per-bin normalised sex–CN coupling weight.
            site_alt: Per-site alt counts (n_bins, max_sites, n_samples).
                Only used when *af_table* is ``None`` (falls back to
                on-the-fly computation).
            site_total: Per-site total counts (n_bins, max_sites, n_samples).
            site_pop_af: Per-site population AFs (n_bins, max_sites).
            site_mask: Per-site validity mask (n_bins, max_sites, n_samples).
            af_table: Precomputed AF log-likelihood table
                ``(n_states, n_bins, n_samples)`` from
                :func:`_precompute_af_table`.  When provided, the model
                performs a cheap table lookup instead of recomputing
                BetaBinomial log-probs at every SVI step.
        """
        zero = torch.zeros(1, device=self.device, dtype=self.dtype)
        one = torch.ones(1, device=self.device, dtype=self.dtype)
        cn_states = torch.arange(
            0, self.n_states, device=self.device, dtype=self.dtype
        )
        sample_var_mean = torch.tensor(
            self.var_sample,
            device=self.device,
            dtype=self.dtype,
        )
        bin_var_rate = torch.tensor(
            1.0 / self.var_bin,
            device=self.device,
            dtype=self.dtype,
        )

        plate_bins = pyro.plate("bins", n_bins, dim=-2, device=self.device)
        plate_samples = pyro.plate("samples", n_samples, dim=-1, device=self.device)
        if self.learn_af_temperature:
            af_temperature = pyro.sample(
                "af_temperature",
                dist.LogNormal(
                    torch.tensor(
                        math.log(max(self.af_weight, 1e-6)),
                        device=self.device,
                        dtype=self.dtype,
                    ),
                    torch.tensor(
                        self.af_temperature_prior_scale,
                        device=self.device,
                        dtype=self.dtype,
                    ),
                ),
            )
        else:
            af_temperature = self._af_scale_torch()

        if self.learn_site_pop_af and site_pop_af is not None and self.af_weight > 0:
            prior_mean = torch.clamp(site_pop_af, min=1e-6, max=1.0 - 1e-6)
            prior_strength = torch.tensor(
                self.site_af_prior_strength,
                device=self.device,
                dtype=self.dtype,
            )
            site_pop_af = pyro.sample(
                "site_pop_af_latent",
                dist.Beta(
                    1.0 + prior_strength * prior_mean,
                    1.0 + prior_strength * (1.0 - prior_mean),
                ).to_event(2),
            )

        # Per-sample variance and sex karyotype
        sample_depth = None
        with plate_samples:
            if self.freeze_sample_var:
                if sample_var_fixed is None:
                    raise ValueError(
                        "sample_var_fixed is required when freeze_sample_var=True."
                    )
                sample_var = pyro.sample(
                    "sample_var",
                    dist.Delta(
                        sample_var_fixed,
                        log_density=torch.zeros_like(sample_var_fixed),
                    ),
                    obs=sample_var_fixed,
                )
            else:
                sample_var = pyro.sample(
                    "sample_var",
                    dist.Exponential(1.0 / sample_var_mean),
                )
            if self.obs_likelihood == NEGATIVE_BINOMIAL_OBS_LIKELIHOOD:
                if self.freeze_sample_depth:
                    if sample_depth_fixed is None:
                        raise ValueError(
                            "sample_depth_fixed is required when freeze_sample_depth=True."
                        )
                    sample_depth = pyro.sample(
                        "sample_depth",
                        dist.Delta(
                            sample_depth_fixed,
                            log_density=torch.zeros_like(sample_depth_fixed),
                        ),
                        obs=sample_depth_fixed,
                    )
                else:
                    if sample_depth_prior_loc is None or sample_depth_prior_scale is None:
                        raise ValueError(
                            "sample_depth_prior_loc and sample_depth_prior_scale are required "
                            "for negative_binomial observation likelihood."
                        )
                    sample_depth = pyro.sample(
                        "sample_depth",
                        dist.LogNormal(
                            sample_depth_prior_loc,
                            sample_depth_prior_scale,
                        ),
                    )
            if self.sex_cn_weight > 0 and chr_type is not None:
                sex_prior = torch.tensor(
                    self.sex_prior, device=self.device, dtype=self.dtype,
                )
                sex_karyotype = pyro.sample(
                    "sex_karyotype", dist.Categorical(sex_prior),
                )

        # Per-bin latent variables
        with plate_bins:
            if self.freeze_bin_bias:
                if bin_bias_fixed is None:
                    raise ValueError(
                        "bin_bias_fixed is required when freeze_bin_bias=True."
                    )
                bin_bias = pyro.sample(
                    "bin_bias",
                    dist.Delta(
                        bin_bias_fixed,
                        log_density=torch.zeros_like(bin_bias_fixed),
                    ),
                    obs=bin_bias_fixed,
                )
            else:
                bin_bias = pyro.sample(
                    "bin_bias", dist.LogNormal(zero, self.var_bias_bin)
                )
            bin_var = pyro.sample("bin_var", dist.Exponential(bin_var_rate))

        if self.epsilon_mean > 0:
            if n_contigs is None or contig_index is None:
                raise ValueError(
                    "n_contigs and contig_index are required when epsilon_mean > 0."
                )
            epsilon_concentration = torch.tensor(
                self.epsilon_concentration,
                device=self.device,
                dtype=self.dtype,
            )
            contig_epsilon = pyro.sample(
                "bin_epsilon",
                dist.Gamma(
                    epsilon_concentration,
                    epsilon_concentration / torch.tensor(
                        self.epsilon_mean,
                        device=self.device,
                        dtype=self.dtype,
                    ),
                ).expand([n_contigs, n_samples]).to_event(2),
            )
            bin_epsilon = contig_epsilon[contig_index, :]
        else:
            bin_epsilon = torch.zeros(
                (n_bins, n_samples),
                dtype=self.dtype,
                device=self.device,
            )

        if self.background_factors > 0:
            background_bin_factors = pyro.sample(
                "background_bin_factors",
                dist.LogNormal(
                    zero,
                    torch.tensor(
                        self.background_bin_scale,
                        device=self.device,
                        dtype=self.dtype,
                    ),
                ).expand([n_bins, self.background_factors]).to_event(2),
            )
            background_sample_factors = pyro.sample(
                "background_sample_factors",
                dist.HalfNormal(
                    torch.tensor(
                        self.background_sample_scale,
                        device=self.device,
                        dtype=self.dtype,
                    )
                ).expand([self.background_factors, n_samples]).to_event(2),
            )
            additive_background = _compose_additive_background_torch(
                bin_epsilon,
                background_bin_factors,
                background_sample_factors,
            )
        else:
            additive_background = bin_epsilon

        if self.autosome_prior_mode == "shrinkage":
            autosome_nonref_mean = pyro.sample(
                "autosome_nonref_mean",
                dist.Beta(
                    torch.tensor(
                        self.autosome_nonref_mean_alpha,
                        device=self.device,
                        dtype=self.dtype,
                    ),
                    torch.tensor(
                        self.autosome_nonref_mean_beta,
                        device=self.device,
                        dtype=self.dtype,
                    ),
                ),
            )
            autosome_mask = (
                chr_type == 0 if chr_type is not None else torch.ones(
                    n_bins,
                    dtype=torch.bool,
                    device=self.device,
                )
            )
            autosome_idx = torch.nonzero(autosome_mask, as_tuple=False).squeeze(-1)
            sex_idx = torch.nonzero(~autosome_mask, as_tuple=False).squeeze(-1)

            autosome_nonref_prob = torch.empty(0, dtype=self.dtype, device=self.device)
            autosome_alt_cn_probs = torch.empty(
                (0, self.n_states - 1),
                dtype=self.dtype,
                device=self.device,
            )
            if autosome_idx.numel() > 0:
                nonref_alpha = torch.clamp(
                    autosome_nonref_mean * self.autosome_nonref_concentration,
                    min=1e-4,
                )
                ref_alpha = torch.clamp(
                    (1.0 - autosome_nonref_mean) * self.autosome_nonref_concentration,
                    min=1e-4,
                )
                with pyro.plate("autosome_bins_prior", autosome_idx.numel()):
                    autosome_nonref_prob = pyro.sample(
                        "autosome_nonref_prob",
                        dist.Beta(nonref_alpha, ref_alpha),
                    )
                    autosome_alt_cn_probs = pyro.sample(
                        "autosome_alt_cn_probs",
                        dist.Dirichlet(
                            torch.full(
                                (autosome_idx.numel(), self.n_states - 1),
                                self.alpha_non_ref,
                                dtype=self.dtype,
                                device=self.device,
                            )
                        ),
                    )

            sex_cn_probs = torch.empty(
                (0, self.n_states),
                dtype=self.dtype,
                device=self.device,
            )
            if sex_idx.numel() > 0:
                sex_alpha = self.alpha_table[chr_type[sex_idx]]
                with pyro.plate("sex_bins_prior", sex_idx.numel()):
                    sex_cn_probs = pyro.sample(
                        "sex_cn_probs",
                        dist.Dirichlet(sex_alpha),
                    )

            cn_probs = torch.empty(
                (n_bins, self.n_states),
                dtype=self.dtype,
                device=self.device,
            )
            if autosome_idx.numel() > 0:
                cn_probs[autosome_idx] = self._compose_autosome_cn_probs_torch(
                    autosome_nonref_prob,
                    autosome_alt_cn_probs,
                )
            if sex_idx.numel() > 0:
                cn_probs[sex_idx] = sex_cn_probs
            cn_probs = pyro.deterministic("cn_probs", cn_probs)
        else:
            if self.freeze_cn_prior:
                cn_probs = pyro.deterministic(
                    "cn_probs",
                    self._default_cn_probs_torch(chr_type, n_bins),
                )
            else:
                with plate_bins:
                    # Dirichlet-Categorical CN prior — per-chromosome-type alphas.
                    if chr_type is not None:
                        alpha = self.alpha_table[chr_type].unsqueeze(-2)
                    else:
                        alpha = self.alpha_non_ref * one.expand(self.n_states)
                        alpha[2] = self.alpha_ref
                    cn_probs = pyro.sample("cn_probs", dist.Dirichlet(alpha))

        # Per-bin, per-sample observations
        with plate_bins, plate_samples:
            # The legacy Dirichlet path produces an explicit singleton sample
            # dimension. Restore that axis for derived shrinkage priors so the
            # categorical probabilities broadcast across the samples plate.
            cn_probs_for_samples = (
                cn_probs.unsqueeze(-2) if cn_probs.dim() == 2 else cn_probs
            )
            cn = pyro.sample("cn", dist.Categorical(cn_probs_for_samples))
            expected_units = Vindex(cn_states)[cn] * bin_bias + additive_background

            if self.obs_likelihood == NEGATIVE_BINOMIAL_OBS_LIKELIHOOD:
                if bin_length_kb is None:
                    raise ValueError(
                        "bin_length_kb is required for negative_binomial observation likelihood."
                    )
                mean = _raw_expected_depth_units(
                    Vindex(cn_states)[cn],
                    bin_bias,
                    additive_background,
                )
                mean = torch.clamp(mean, min=0.0)
                mean = mean * sample_depth
                mean = mean * bin_length_kb.unsqueeze(-1)
                overdispersion = sample_var + bin_var
                pyro.sample(
                    "obs",
                    _make_negative_binomial_distribution(mean, overdispersion),
                    obs=depth,
                )
            else:
                base_variance = sample_var + bin_var
                variance = base_variance * torch.clamp(expected_units, min=1e-3)
                pyro.sample(
                    "obs",
                    _make_depth_distribution(
                        expected_units,
                        variance,
                        self.obs_likelihood,
                        self.obs_df,
                    ),
                    obs=depth,
                )

            # ── sex–CN coupling factor ────────────────────────────────────────
            if self.sex_cn_weight > 0 and chr_type is not None:
                chr_type_exp = chr_type.unsqueeze(-1)       # (n_bins, 1)
                score = Vindex(self.sex_cn_score)[
                    sex_karyotype, chr_type_exp, cn,
                ]
                pyro.factor(
                    "sex_cn_prior", sex_cn_weight_per_bin * score,
                )

            # ── marginalized allele-fraction likelihood ───────────────────
            if af_table is not None and self.af_weight > 0:
                # Fast path: precomputed table lookup.
                # af_table shape: (n_states, n_bins, n_samples).
                # cn has enumeration dims: (..., 1, 1) with values 0..n_states-1.
                # Select the row matching each enumerated cn state via torch.where.
                af_log_lik = torch.tensor(0.0, device=self.device)
                for c in range(self.n_states):
                    af_log_lik = af_log_lik + torch.where(
                        cn == c, af_table[c], zero,
                    )
                pyro.factor("af_lik", af_temperature * af_log_lik)
            elif site_alt is not None and site_total is not None and self.af_weight > 0:
                # Slow fallback: precompute once inside the model call.
                af_table = self._prepare_af_table_torch(
                    site_alt,
                    site_total,
                    site_pop_af,
                    site_mask,
                    chr_type,
                    auto_leave_one_out=not self.learn_site_pop_af,
                )
                af_log_lik = torch.tensor(0.0, device=self.device)
                for c in range(self.n_states):
                    af_log_lik = af_log_lik + torch.where(
                        cn == c, af_table[c], zero,
                    )
                pyro.factor("af_lik", af_temperature * af_log_lik)

    # ── training ────────────────────────────────────────────────────────

    def _model_kwargs(self, data: DepthData) -> dict:
        """Build keyword arguments dict for model/guide calls.

        When per-site allele data is present and ``af_weight > 0``, the
        expensive BetaBinomial marginalisation is done **once** here and
        passed to the model as a precomputed lookup table (``af_table``).
        The stored AF term is summed over observed sites within each
        bin/sample so bins with more informative loci contribute
        proportionally more evidence.
        """
        kw: dict = {
            "depth": data.depth,
            "n_bins": data.n_bins,
            "n_samples": data.n_samples,
        }
        if self.epsilon_mean > 0:
            kw["n_contigs"] = data.n_contigs
            kw["contig_index"] = data.contig_index
        if self.freeze_sample_var:
            kw["sample_var_fixed"] = torch.zeros(
                (data.n_samples,),
                device=self.device,
                dtype=self.dtype,
            )
        if self.freeze_bin_bias:
            kw["bin_bias_fixed"] = torch.ones(
                (data.n_bins, 1),
                device=self.device,
                dtype=self.dtype,
            )
        if self.obs_likelihood == NEGATIVE_BINOMIAL_OBS_LIKELIHOOD:
            kw["bin_length_kb"] = data.bin_length_kb
            sample_depth_prior = self._estimate_sample_depth_prior(data)
            kw["sample_depth_prior_loc"] = sample_depth_prior["loc"]
            kw["sample_depth_prior_scale"] = sample_depth_prior["scale"]
            if self.freeze_sample_depth:
                kw["sample_depth_fixed"] = torch.tensor(
                    self._estimate_sample_depth_center(data),
                    device=self.device,
                    dtype=self.dtype,
                )
        # Chromosome-type tensor and per-bin sex-CN coupling weight
        if hasattr(data, "chr_type") and data.chr_type is not None:
            kw["chr_type"] = data.chr_type
            if self.sex_cn_weight > 0:
                ct = data.chr_type
                counts = torch.zeros(3, device=ct.device, dtype=self.dtype)
                for t in range(3):
                    counts[t] = (ct == t).sum().float().clamp(min=1)
                kw["sex_cn_weight_per_bin"] = (
                    self.sex_cn_weight / counts[ct]
                ).unsqueeze(-1)  # (n_bins, 1)
        if data.site_alt is not None and self.af_weight > 0:
            if self.learn_site_pop_af:
                kw["site_alt"] = data.site_alt
                kw["site_total"] = data.site_total
                kw["site_pop_af"] = data.site_pop_af
                kw["site_mask"] = data.site_mask
            else:
                with torch.no_grad():
                    kw["af_table"] = self._prepare_af_table_torch(
                        data.site_alt, data.site_total,
                        data.site_pop_af, data.site_mask,
                        kw.get("chr_type"),
                    )
        return kw

    def train(
        self,
        data: DepthData,
        *,
        max_iter: int = 1000,
        lr_init: float = 0.01,
        lr_min: float = 0.001,
        lr_decay: float = 500,
        adam_beta1: float = 0.9,
        adam_beta2: float = 0.999,
        grad_clip_norm: float | None = 10.0,
        init_restarts: int = 1,
        log_freq: int = 50,
        jit: bool = False,
        early_stopping: bool = True,
        patience: int = 50,
        convergence_window: int = 50,
        convergence_rtol: float = 1e-4,
    ) -> List[float]:
        """Train the model with stochastic variational inference.

        Args:
            data: :class:`DepthData` instance.
            max_iter: Maximum SVI iterations.
            lr_init: Initial learning rate.
            lr_min: Minimum learning rate.
            lr_decay: Exponential decay constant.
            adam_beta1, adam_beta2: Adam beta parameters.
            grad_clip_norm: Global gradient-norm clip applied by the
                Pyro optimizer. Set to ``None`` or a non-positive value to
                disable clipping.
            init_restarts: Number of candidate guide initializations to
                evaluate before gradient descent. The first candidate uses
                the existing anchored initialization and the remainder sample
                random starts from the prior.
            log_freq: Logging frequency (iterations).
            jit: Whether to JIT-compile the ELBO.
            early_stopping: Enable early stopping.
            patience: Consecutive window comparisons below the relative
                tolerance before stopping.
            convergence_window: Number of iterations per rolling ELBO window.
            convergence_rtol: Relative tolerance between successive rolling
                ELBO-window means.

        Returns:
            List of per-epoch ELBO losses.
        """
        logger.info("Initialising training ...")
        if init_restarts < 1:
            raise ValueError("init_restarts must be at least 1.")
        if early_stopping:
            if patience < 1:
                raise ValueError("patience must be at least 1.")
            if convergence_window < 1:
                raise ValueError("convergence_window must be at least 1.")
            if convergence_rtol < 0:
                raise ValueError("convergence_rtol must be non-negative.")
        if grad_clip_norm is not None:
            if not np.isfinite(grad_clip_norm):
                raise ValueError("grad_clip_norm must be finite when provided.")
            if grad_clip_norm <= 0:
                grad_clip_norm = None

        self.loss_history = {"epoch": [], "elbo": []}

        pyro.clear_param_store()
        model_kw = self._model_kwargs(data)

        clip_args = None
        if grad_clip_norm is not None:
            clip_args = {"clip_norm": float(grad_clip_norm)}

        scheduler = pyro.optim.LambdaLR(
            {
                "optimizer": torch.optim.Adam,
                "optim_args": {"lr": 1.0, "betas": (adam_beta1, adam_beta2)},
                "lr_lambda": lambda k: (
                    lr_min + (lr_init - lr_min) * np.exp(-k / lr_decay)
                ),
            },
            clip_args=clip_args,
        )

        elbo = JitTraceEnum_ELBO() if jit else TraceEnum_ELBO()
        self._select_svi_initialization(
            data,
            model_kw=model_kw,
            elbo=elbo,
            init_restarts=init_restarts,
        )
        svi = SVI(self.model, self.guide, optim=scheduler, loss=elbo)

        if grad_clip_norm is None:
            logger.info("Training guardrails: gradient clipping disabled.")
        else:
            logger.info(
                "Training guardrails: clip gradient norm at %.2f.",
                grad_clip_norm,
            )
        logger.info("Training for up to %d iterations ...", max_iter)

        patience_ctr = 0

        with tqdm(range(max_iter), desc="Training", unit="epoch") as pbar:
            for epoch in pbar:
                loss = svi.step(**model_kw)
                scheduler.step()

                self.loss_history["epoch"].append(epoch)
                self.loss_history["elbo"].append(loss)

                relative_change = None
                if early_stopping:
                    relative_change = _windowed_relative_elbo_change(
                        self.loss_history["elbo"],
                        convergence_window,
                    )
                if relative_change is None:
                    pbar.set_postfix(loss=f"{loss:.4f}")
                else:
                    pbar.set_postfix(loss=f"{loss:.4f}", rel=f"{relative_change:.2e}")

                if (epoch + 1) % log_freq == 0:
                    if relative_change is None:
                        tqdm.write(f"[epoch {epoch + 1:04d}]  loss: {loss:.4f}")
                    else:
                        tqdm.write(
                            f"[epoch {epoch + 1:04d}]  loss: {loss:.4f}  "
                            f"rel_change: {relative_change:.2e}"
                        )

                if early_stopping and relative_change is not None:
                    if relative_change < convergence_rtol:
                        patience_ctr += 1
                    else:
                        patience_ctr = 0
                    if patience_ctr >= patience:
                        tqdm.write(
                            f"\nEarly stopping at epoch {epoch + 1} "
                            f"(window={convergence_window}, "
                            f"rel_change={relative_change:.2e})"
                        )
                        break

        logger.info("Training complete.  Final loss: %.4f", loss)
        return self.loss_history["elbo"]

    # ── posterior estimation ─────────────────────────────────────────────

    def get_map_estimates(
        self,
        data: DepthData,
        estimate_method: str = "current",
        model_kw: dict | None = None,
    ) -> Dict[str, np.ndarray]:
        """Compute point estimates for all latent variables.

        Args:
            data: :class:`DepthData` used during training.
            estimate_method: ``'current'`` for the historical single guide
                draw or ``'median'`` for deterministic guide medians.

        Returns:
            Dictionary mapping site names to numpy arrays.
        """
        logger.info("Computing point estimates (method=%s) ...", estimate_method)
        model_kw = self._model_kwargs(data) if model_kw is None else model_kw
        continuous_estimates = self._get_continuous_estimates(
            data,
            estimate_method=estimate_method,
            model_kw=model_kw,
        )

        estimates: Dict[str, np.ndarray] = {}
        for site, value in continuous_estimates.items():
            estimates[site] = value.detach().cpu().numpy()
        if all(
            (
                self.freeze_bin_bias,
                "bin_bias" not in estimates,
                "bin_bias_fixed" in model_kw,
            )
        ):
            estimates["bin_bias"] = model_kw["bin_bias_fixed"].detach().cpu().numpy()
        if all(
            (
                self.freeze_sample_var,
                "sample_var" not in estimates,
                "sample_var_fixed" in model_kw,
            )
        ):
            estimates["sample_var"] = (
                model_kw["sample_var_fixed"].detach().cpu().numpy()
            )
        if all(
            (
                self.obs_likelihood == NEGATIVE_BINOMIAL_OBS_LIKELIHOOD,
                self.freeze_sample_depth,
                "sample_depth" not in estimates,
                "sample_depth_fixed" in model_kw,
            )
        ):
            estimates["sample_depth"] = (
                model_kw["sample_depth_fixed"].detach().cpu().numpy()
            )
        estimates.update(
            self._infer_discrete_assignments(
                data,
                continuous_estimates,
                model_kw=model_kw,
            )
        )
        return estimates

    def _run_discrete_inference_fixed_latents(
        self,
        data: DepthData,
        maps: Dict[str, np.ndarray],
        af_table: np.ndarray | None = None,
    ) -> Dict[str, np.ndarray]:
        """Compute exact discrete CN posteriors with fixed continuous latents."""
        bin_bias = np.atleast_1d(np.asarray(maps["bin_bias"]).squeeze())
        sample_var = np.atleast_1d(np.asarray(maps["sample_var"]).squeeze())
        bin_var = np.atleast_1d(np.asarray(maps["bin_var"]).squeeze())
        additive_background = compose_additive_background_matrix(
            maps.get("bin_epsilon"),
            data.n_bins,
            data.n_samples,
            dtype=np.float64,
            contig_index=data.contig_index.detach().cpu().numpy(),
            n_contigs=data.n_contigs,
            background_bin_factors=maps.get("background_bin_factors"),
            background_sample_factors=maps.get("background_sample_factors"),
        )
        chr_type_np = None
        if hasattr(data, "chr_type") and data.chr_type is not None:
            chr_type_np = data.chr_type.detach().cpu().numpy()
        cn_probs = self._build_cn_probs_from_estimates_numpy(
            maps,
            chr_type_np,
            data.n_bins,
        )
        obs = data.depth.detach().cpu().numpy()
        cn_states = np.arange(self.n_states, dtype=obs.dtype).reshape(-1, 1, 1)
        obs_b = obs[np.newaxis, :, :]

        if self.obs_likelihood == NEGATIVE_BINOMIAL_OBS_LIKELIHOOD:
            sample_depth = np.atleast_1d(
                np.asarray(maps["sample_depth"]).squeeze()
            )
            overdispersion = (
                sample_var[np.newaxis, :] + bin_var[:, np.newaxis]
            )
            expected_units = _raw_expected_depth_units(
                cn_states,
                bin_bias[np.newaxis, :, np.newaxis],
                additive_background[np.newaxis, :, :],
            )
            mean = data.bin_length_kb.detach().cpu().numpy()[np.newaxis, :, np.newaxis]
            mean = mean * sample_depth[np.newaxis, np.newaxis, :]
            mean = mean * np.maximum(expected_units, 0.0)
            log_lik = _negative_binomial_log_lik_numpy(
                obs_b,
                mean,
                overdispersion[np.newaxis, :, :],
            )
        else:
            base_variance = sample_var[np.newaxis, :] + bin_var[:, np.newaxis]
            expected_depth = cn_states * bin_bias[np.newaxis, :, np.newaxis]
            expected_depth = expected_depth + additive_background[np.newaxis, :, :]
            variance = base_variance[np.newaxis, :, :]
            variance = variance * np.maximum(expected_depth, 1e-3)
            log_lik = _depth_log_lik_numpy(
                obs_b,
                expected_depth,
                variance,
                self.obs_likelihood,
                self.obs_df,
            )

        log_prior = np.log(np.maximum(cn_probs.T[:, :, np.newaxis], 1e-10))
        base_log_unnorm = log_lik + log_prior  # (n_states, n_bins, n_samples)

        af_scale = self._af_scale_numpy(maps)
        if af_table is not None and af_scale > 0:
            base_log_unnorm += af_scale * af_table
        elif data.site_alt is not None and af_scale > 0:
            site_alt_np = data.site_alt.detach().cpu().numpy()
            site_total_np = data.site_total.detach().cpu().numpy()
            site_pop_af_np = np.asarray(
                maps.get(
                    "site_pop_af_latent",
                    data.site_pop_af.detach().cpu().numpy(),
                ),
                dtype=np.float64,
            )
            site_mask_np = data.site_mask.detach().cpu().numpy()
            if "site_pop_af_latent" not in maps:
                site_pop_af_np, used_leave_one_out = _resolve_fixed_site_pop_af_numpy(
                    site_alt_np,
                    site_total_np,
                    site_pop_af_np,
                    site_mask_np,
                )
                if used_leave_one_out:
                    logger.info(
                        "AF inference: detected self-pooled site_pop_af; using leave-one-out population AFs per sample."
                    )
            raw_af_table = np.empty(
                (self.n_states, data.n_bins, data.n_samples),
                dtype=np.float64,
            )

            for c in range(self.n_states):
                raw_af_table[c] = _marginalized_af_log_lik_numpy(
                    site_alt_np, site_total_np, site_pop_af_np, site_mask_np,
                    cn_state=c, n_states=self.n_states,
                    concentration=self.af_concentration,
                )

            af_ll_table = self._apply_af_evidence_mode_numpy(
                raw_af_table,
                chr_type_np,
            )
            base_log_unnorm += af_scale * af_ll_table

            n_sites = int(site_mask_np.any(axis=2).sum())
            logger.info(
                "AF marginalised likelihood: %d total sites, summed per "
                "bin/sample (scale=%.2f, mode=%s)",
                n_sites, af_scale,
                self.af_evidence_mode,
            )

        has_sex = (
            self.sex_cn_weight > 0 and
            hasattr(data, "chr_type") and
            data.chr_type is not None
        )
        if has_sex:
            chr_type_np = data.chr_type.cpu().numpy()
            sex_cn_np = self.sex_cn_score.cpu().numpy()  # (2, 3, n_states)
            sex_prior_np = np.array(self.sex_prior, dtype=np.float64)

            ct_counts = np.array(
                [max(int((chr_type_np == t).sum()), 1) for t in range(3)],
                dtype=np.float64,
            )
            w_per_bin = self.sex_cn_weight / ct_counts[chr_type_np]

            n_bins, n_samp = data.n_bins, data.n_samples
            log_evidence = np.zeros((2, n_samp))
            cn_post_given_sex = np.zeros(
                (2, self.n_states, n_bins, n_samp),
            )

            for s in range(2):
                penalty = np.zeros((self.n_states, n_bins, 1))
                for c in range(self.n_states):
                    penalty[c, :, 0] = (
                        w_per_bin * sex_cn_np[s, chr_type_np, c]
                    )
                log_unnorm_s = base_log_unnorm + penalty

                mx = np.max(log_unnorm_s, axis=0, keepdims=True)
                ev = np.exp(log_unnorm_s - mx)
                log_Z_bin = mx.squeeze(0) + np.log(ev.sum(axis=0) + 1e-30)

                log_evidence[s] = log_Z_bin.sum(axis=0)
                cn_post_given_sex[s] = ev / ev.sum(axis=0, keepdims=True)

            log_sex = (
                np.log(np.maximum(sex_prior_np, 1e-10))[:, np.newaxis] +
                log_evidence
            )
            mx_s = np.max(log_sex, axis=0, keepdims=True)
            sex_post = np.exp(log_sex - mx_s)
            sex_post /= sex_post.sum(axis=0, keepdims=True)

            cn_marginal = np.zeros((self.n_states, n_bins, n_samp))
            for s in range(2):
                cn_marginal += (
                    sex_post[s][np.newaxis, np.newaxis, :] *
                    cn_post_given_sex[s]
                )

            cn_posterior = np.transpose(cn_marginal, (1, 2, 0)).astype(
                np.float32, copy=False,
            )
            return {
                "cn_posterior": cn_posterior,
                "sex_posterior": sex_post.T.astype(np.float32, copy=False),
            }

        max_log = np.max(base_log_unnorm, axis=0, keepdims=True)
        exp_vals = np.exp(base_log_unnorm - max_log)
        posterior = exp_vals / np.sum(exp_vals, axis=0, keepdims=True)
        cn_posterior = np.transpose(posterior, (1, 2, 0)).astype(
            np.float32, copy=False,
        )
        return {"cn_posterior": cn_posterior}

    def run_discrete_inference(
        self,
        data: DepthData,
        map_estimates: Dict[str, np.ndarray] | None = None,
        af_table: np.ndarray | None = None,
    ) -> Dict[str, np.ndarray]:
        """Compute exact discrete CN posteriors analytically.

        Once the continuous latents are fixed at point estimates, the
        hidden copy-number state for each bin / sample pair has a small,
        finite state space and a closed-form likelihood. That means the
        posterior can be computed directly with Bayes' rule instead of via
        repeated Monte Carlo calls to :func:`infer_discrete`.

        When per-site allele data is available, the marginalized genotype
        likelihood for each CN state is added to the log-posterior.

        Args:
            data: :class:`DepthData` instance.
            map_estimates: Optional precomputed output from
                :meth:`get_map_estimates` to avoid recomputation.
            af_table: Optional precomputed AF log-likelihood table of shape
                ``(n_states, n_bins, n_samples)``.

        Returns:
            Dictionary with ``'cn_posterior'`` of shape
            ``(n_bins, n_samples, n_states)`` and, when sex–CN coupling is
            active, ``'sex_posterior'`` of shape ``(n_samples, 2)`` giving
            ``[P(XX), P(XY)]`` per sample.
        """
        logger.info("Running discrete inference ...")

        maps = map_estimates if map_estimates is not None else self.get_map_estimates(data)
        posterior = self._run_discrete_inference_fixed_latents(
            data,
            maps,
            af_table=af_table,
        )
        logger.info("Exact analytical discrete inference complete.")
        return posterior

    def run_discrete_inference_multi_draw(
        self,
        data: DepthData,
        n_draws: int = 30,
        af_table: np.ndarray | None = None,
        draw_estimate_collector: Dict[str, List[np.ndarray]] | None = None,
    ) -> Dict[str, np.ndarray]:
        """Average discrete CN posteriors over repeated guide draws.

        When ``draw_estimate_collector`` is provided, each guide draw's
        continuous latent values are appended to that dictionary so downstream
        consumers can form a true posterior predictive distribution instead of
        relying only on plug-in point estimates.
        """
        if n_draws < 1:
            raise ValueError("n_draws must be at least 1.")

        if self.guide_type == "delta":
            logger.info(
                "Guide type 'delta' is deterministic; multi-draw inference "
                "collapses to a single point-estimate pass.",
            )
            point_estimates = self.get_map_estimates(data, estimate_method="current")
            return self.run_discrete_inference(
                data,
                map_estimates=point_estimates,
                af_table=af_table,
            )

        logger.info(
            "Running multi-draw discrete inference (%d guide draws) ...",
            n_draws,
        )
        model_kw = self._model_kwargs(data)
        posterior_sum: Dict[str, np.ndarray] = {}
        cn_map_counts: np.ndarray | None = None
        sex_map_counts: np.ndarray | None = None
        initialized = False
        for _ in range(n_draws):
            draw_estimates = self._get_continuous_estimates(
                data,
                estimate_method="current",
                model_kw=model_kw,
            )
            draw_maps = {
                site: value.detach().cpu().numpy()
                for site, value in draw_estimates.items()
            }
            if draw_estimate_collector is not None:
                for site in draw_maps:
                    if site == "site_pop_af_latent":
                        continue
                    draw_estimate_collector.setdefault(site, [])
                for site, value in draw_maps.items():
                    if site == "site_pop_af_latent":
                        continue
                    draw_estimate_collector[site].append(np.array(value, copy=True))
            draw_post = self._run_discrete_inference_fixed_latents(
                data,
                draw_maps,
                af_table=af_table,
            )
            if not initialized:
                posterior_sum = {
                    key: np.asarray(value, dtype=np.float64)
                    for key, value in draw_post.items()
                }
                cn_map_counts = np.zeros_like(
                    draw_post["cn_posterior"],
                    dtype=np.float64,
                )
                if "sex_posterior" in draw_post:
                    sex_map_counts = np.zeros_like(
                        draw_post["sex_posterior"],
                        dtype=np.float64,
                    )
                initialized = True
            else:
                for key, value in draw_post.items():
                    if key not in posterior_sum:
                        continue
                    posterior_sum[key] += np.asarray(value, dtype=np.float64)
            cn_map = np.argmax(draw_post["cn_posterior"], axis=2)
            np.add.at(
                cn_map_counts,
                (
                    np.arange(data.n_bins)[:, np.newaxis],
                    np.arange(data.n_samples)[np.newaxis, :],
                    cn_map,
                ),
                1.0,
            )
            if sex_map_counts is not None and "sex_posterior" in draw_post:
                sex_map = np.argmax(draw_post["sex_posterior"], axis=1)
                np.add.at(
                    sex_map_counts,
                    (np.arange(data.n_samples), sex_map),
                    1.0,
                )

        averaged = {
            key: (value / float(n_draws)).astype(np.float32, copy=False)
            for key, value in posterior_sum.items()
        }
        if cn_map_counts is not None:
            cn_map_counts = cn_map_counts / float(n_draws)
            final_cn_map = np.argmax(averaged["cn_posterior"], axis=2)
            averaged["cn_map_stability"] = np.take_along_axis(
                cn_map_counts,
                final_cn_map[..., np.newaxis],
                axis=2,
            ).squeeze(axis=2).astype(np.float32, copy=False)
        if sex_map_counts is not None and "sex_posterior" in averaged:
            sex_map_counts = sex_map_counts / float(n_draws)
            final_sex_map = np.argmax(averaged["sex_posterior"], axis=1)
            averaged["sex_map_stability"] = np.take_along_axis(
                sex_map_counts,
                final_sex_map[:, np.newaxis],
                axis=1,
            ).squeeze(axis=1).astype(np.float32, copy=False)
        logger.info("Exact analytical multi-draw discrete inference complete.")
        return averaged
