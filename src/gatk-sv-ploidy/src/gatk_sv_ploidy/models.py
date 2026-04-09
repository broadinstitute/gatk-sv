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

import logging
from typing import Dict, List, Optional

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
from pyro.infer.autoguide import AutoDelta, AutoDiagonalNormal
from tqdm import tqdm

from gatk_sv_ploidy._util import DEFAULT_AF_CONCENTRATION, DEFAULT_AF_WEIGHT
from gatk_sv_ploidy.data import DepthData

logger = logging.getLogger(__name__)


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
        site_pop_af: ``(n_bins, max_sites)`` population alt-allele frequencies.
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
    pop_af = site_pop_af.unsqueeze(-1)  # (n_bins, max_sites, 1)
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

    # Mask out invalid sites and average over observed sites so bins with more
    # retained loci do not automatically outweigh bins with fewer loci.
    log_marginal = torch.where(site_mask, log_marginal, torch.zeros_like(log_marginal))
    site_counts = site_mask.to(log_marginal.dtype).sum(dim=-2)
    site_counts_safe = torch.clamp(site_counts, min=1.0)
    return log_marginal.sum(dim=-2) / site_counts_safe


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
        site_pop_af: ``(n_bins, max_sites)``
        site_mask: ``(n_bins, max_sites, n_samples)``
        cn_state: The CN value (integer).
        n_states: Max CN.
        concentration: Beta-Binomial concentration.

    Returns:
        ``(n_bins, n_samples)`` log-likelihood averaged over observed sites
        within each bin.
    """
    from scipy.stats import betabinom as betabinom_scipy
    from scipy.stats import binom as binom_scipy

    eps = 1e-6
    site_total_safe = np.maximum(site_total, 1)
    pop_af = site_pop_af[:, :, np.newaxis]   # (n_bins, max_sites, 1)
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

    # Mask and average over observed sites so AF contribution is not a direct
    # function of how many loci were retained in each bin.
    log_marginal = np.where(site_mask, log_marginal, 0.0)
    site_counts = site_mask.sum(axis=1, dtype=np.float64)
    site_counts_safe = np.maximum(site_counts, 1.0)
    return log_marginal.sum(axis=1) / site_counts_safe


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
        at bin *b* and sample *s*, normalized by observed site count.
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
    logger.info("Precomputed AF table: shape %s", list(table.shape))
    return table


class CNVModel:
    """Hierarchical Bayesian model for whole-genome CN detection.

    The generative model uses a Dirichlet-Categorical prior over CN states
    (0–5), per-bin bias and variance, and per-sample variance.  Observed
    normalised depth is drawn from a Normal likelihood whose mean is
    ``CN × bin_bias`` and whose variance is ``bin_var + sample_var``.

    When per-site allele data is provided, the model additionally includes
    a marginalized allele-fraction likelihood that jointly considers all
    possible genotype states at each SNP site, weighted by population
    allele frequencies and normalized by observed site count per bin.

    Parameters
    ----------
    n_states : int
        Number of copy-number states (default 6 for CN 0–5).
    alpha_ref : float
        Dirichlet concentration for the reference state (CN = 2).
    alpha_non_ref : float
        Dirichlet concentration for non-reference states.
    var_bias_bin : float
        Scale of the LogNormal prior on per-bin bias.
    var_sample : float
        Scale of the Exponential prior on per-sample variance.
    var_bin : float
        Scale of the Exponential prior on per-bin variance.
    device, dtype : str, torch.dtype
        Torch device and floating-point type.
    guide_type : str
        ``'delta'`` for :class:`AutoDelta` or ``'diagonal'`` for
        :class:`AutoDiagonalNormal`.
    af_concentration : float
        Beta-Binomial concentration for the allele-fraction model.
    af_weight : float
        Relative weight of the per-bin, site-count-normalized allele-fraction
        likelihood (0 to disable).
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
    """

    # The list of continuous latent sites exposed to the guide.
    _latent_sites = ["bin_bias", "sample_var", "bin_var", "cn_probs"]

    def __init__(
        self,
        n_states: int = 6,
        alpha_ref: float = 50.0,
        alpha_non_ref: float = 1.0,
        var_bias_bin: float = 0.1,
        var_sample: float = 0.2,
        var_bin: float = 0.2,
        device: str = "cpu",
        dtype: torch.dtype = torch.float32,
        guide_type: str = "diagonal",
        af_concentration: float = DEFAULT_AF_CONCENTRATION,
        af_weight: float = DEFAULT_AF_WEIGHT,
        alpha_sex_ref: float = 1.0,
        alpha_sex_non_ref: float = 1.0,
        sex_prior: tuple = (0.5, 0.5),
        sex_cn_weight: float = 5.0,
    ) -> None:
        self.n_states = n_states
        self.alpha_ref = alpha_ref
        self.alpha_non_ref = alpha_non_ref
        self.var_bias_bin = var_bias_bin
        self.var_sample = var_sample
        self.var_bin = var_bin
        self.device = device
        self.dtype = dtype
        self.guide_type = guide_type
        self.af_concentration = af_concentration
        self.af_weight = af_weight
        self.alpha_sex_ref = alpha_sex_ref
        self.alpha_sex_non_ref = alpha_sex_non_ref
        self.sex_prior = sex_prior
        self.sex_cn_weight = sex_cn_weight

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

        # Build variational guide
        blocked = poutine.block(self.model, expose=self._latent_sites)
        if guide_type == "delta":
            self.guide = AutoDelta(blocked)
        elif guide_type == "diagonal":
            self.guide = AutoDiagonalNormal(blocked)
        else:
            raise ValueError(
                f"Unknown guide_type: {guide_type!r}. Choose 'diagonal' or 'delta'."
            )

    # ── probabilistic model ─────────────────────────────────────────────

    @config_enumerate(default="parallel")
    def model(
        self,
        depth: torch.Tensor,
        n_bins: int = None,
        n_samples: int = None,
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
            depth: Observed normalised depth tensor (n_bins × n_samples).
            n_bins: Number of genomic bins.
            n_samples: Number of samples.
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

        plate_bins = pyro.plate("bins", n_bins, dim=-2, device=self.device)
        plate_samples = pyro.plate("samples", n_samples, dim=-1, device=self.device)

        # Per-sample variance and sex karyotype
        with plate_samples:
            sample_var = pyro.sample(
                "sample_var", dist.Exponential(1.0 / self.var_sample)
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
            bin_bias = pyro.sample(
                "bin_bias", dist.LogNormal(zero, self.var_bias_bin)
            )
            bin_var = pyro.sample(
                "bin_var", dist.Exponential(1.0 / self.var_bin)
            )

            # Dirichlet-Categorical CN prior — per-chromosome-type alphas
            if chr_type is not None:
                # (n_bins, n_states) → (n_bins, 1, n_states) so the bins
                # dimension sits at batch dim -2, matching plate_bins.
                alpha = self.alpha_table[chr_type].unsqueeze(-2)
            else:
                alpha = self.alpha_non_ref * one.expand(self.n_states)
                alpha[2] = self.alpha_ref
            cn_probs = pyro.sample("cn_probs", dist.Dirichlet(alpha))

        # Per-bin, per-sample observations
        with plate_bins, plate_samples:
            cn = pyro.sample("cn", dist.Categorical(cn_probs))
            expected = Vindex(cn_states)[cn] * bin_bias
            variance = sample_var + bin_var
            pyro.sample("obs", dist.Normal(expected, variance.sqrt()), obs=depth)

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
                pyro.factor("af_lik", self.af_weight * af_log_lik)
            elif site_alt is not None and site_total is not None and self.af_weight > 0:
                # Slow fallback: compute BetaBinomial on the fly.
                af_log_lik = _marginalized_af_log_lik(
                    site_alt, site_total, site_pop_af, site_mask,
                    cn=cn, n_states=self.n_states,
                    concentration=self.af_concentration,
                )
                pyro.factor("af_lik", self.af_weight * af_log_lik)

    # ── training ────────────────────────────────────────────────────────

    def _model_kwargs(self, data: DepthData) -> dict:
        """Build keyword arguments dict for model/guide calls.

        When per-site allele data is present and ``af_weight > 0``, the
        expensive BetaBinomial marginalisation is done **once** here and
        passed to the model as a precomputed lookup table (``af_table``).
        The stored AF term is normalized by observed site count within each
        bin/sample so changing ``max_sites_per_bin`` does not require a
        compensating change to ``af_weight``.
        """
        kw: dict = {
            "depth": data.depth,
            "n_bins": data.n_bins,
            "n_samples": data.n_samples,
        }
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
            with torch.no_grad():
                kw["af_table"] = _precompute_af_table(
                    data.site_alt, data.site_total,
                    data.site_pop_af, data.site_mask,
                    n_states=self.n_states,
                    concentration=self.af_concentration,
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
        log_freq: int = 50,
        jit: bool = False,
        early_stopping: bool = True,
        patience: int = 50,
        min_delta: float = 1000.0,
    ) -> List[float]:
        """Train the model with stochastic variational inference.

        Args:
            data: :class:`DepthData` instance.
            max_iter: Maximum SVI iterations.
            lr_init: Initial learning rate.
            lr_min: Minimum learning rate.
            lr_decay: Exponential decay constant.
            adam_beta1, adam_beta2: Adam beta parameters.
            log_freq: Logging frequency (iterations).
            jit: Whether to JIT-compile the ELBO.
            early_stopping: Enable early stopping.
            patience: Epochs without improvement before stopping.
            min_delta: Minimum loss improvement.

        Returns:
            List of per-epoch ELBO losses.
        """
        logger.info("Initialising training …")
        pyro.clear_param_store()

        scheduler = pyro.optim.LambdaLR(
            {
                "optimizer": torch.optim.Adam,
                "optim_args": {"lr": 1.0, "betas": (adam_beta1, adam_beta2)},
                "lr_lambda": lambda k: (
                    lr_min + (lr_init - lr_min) * np.exp(-k / lr_decay)
                ),
            }
        )

        elbo = JitTraceEnum_ELBO() if jit else TraceEnum_ELBO()
        svi = SVI(self.model, self.guide, optim=scheduler, loss=elbo)

        logger.info("Training for up to %d iterations …", max_iter)

        best_loss = float("inf")
        patience_ctr = 0

        model_kw = self._model_kwargs(data)

        with tqdm(range(max_iter), desc="Training", unit="epoch") as pbar:
            for epoch in pbar:
                loss = svi.step(**model_kw)
                scheduler.step()

                self.loss_history["epoch"].append(epoch)
                self.loss_history["elbo"].append(loss)
                pbar.set_postfix(loss=f"{loss:.4f}")

                if (epoch + 1) % log_freq == 0:
                    tqdm.write(f"[epoch {epoch + 1:04d}]  loss: {loss:.4f}")

                if early_stopping:
                    if loss < best_loss - min_delta:
                        best_loss = loss
                        patience_ctr = 0
                    else:
                        patience_ctr += 1
                    if patience_ctr >= patience:
                        tqdm.write(
                            f"\nEarly stopping at epoch {epoch + 1} "
                            f"(best loss: {best_loss:.4f})"
                        )
                        break

        logger.info("Training complete.  Final loss: %.4f", loss)
        return self.loss_history["elbo"]

    # ── posterior estimation ─────────────────────────────────────────────

    def get_map_estimates(self, data: DepthData) -> Dict[str, np.ndarray]:
        """Compute MAP estimates for all latent variables.

        Args:
            data: :class:`DepthData` used during training.

        Returns:
            Dictionary mapping site names to numpy arrays.
        """
        logger.info("Computing MAP estimates …")
        model_kw = self._model_kwargs(data)

        has_sex = (
            self.sex_cn_weight > 0
            and "chr_type" in model_kw
            and model_kw["chr_type"] is not None
        )
        first_dim = -4 if has_sex else -3

        guide_trace = poutine.trace(self.guide).get_trace(**model_kw)
        replayed = poutine.replay(self.model, trace=guide_trace)
        inferred = infer_discrete(
            replayed, temperature=0, first_available_dim=first_dim,
        )
        trace = poutine.trace(inferred).get_trace(**model_kw)

        estimates: Dict[str, np.ndarray] = {}
        for site in ("bin_bias", "sample_var", "bin_var", "cn_probs"):
            estimates[site] = (
                guide_trace.nodes[site]["value"].detach().cpu().numpy()
            )
        estimates["cn"] = trace.nodes["cn"]["value"].detach().cpu().numpy()
        if has_sex and "sex_karyotype" in trace.nodes:
            estimates["sex_karyotype"] = (
                trace.nodes["sex_karyotype"]["value"].detach().cpu().numpy()
            )
        return estimates

    def run_discrete_inference(
        self,
        data: DepthData,
        map_estimates: Dict[str, np.ndarray] | None = None,
    ) -> Dict[str, np.ndarray]:
        """Compute exact discrete CN posteriors analytically.

        Once the continuous latents are fixed at their MAP estimates, the
        hidden copy-number state for each bin / sample pair has a small,
        finite state space and a Gaussian likelihood.  That means the
        posterior can be computed directly with Bayes' rule instead of via
        repeated Monte Carlo calls to :func:`infer_discrete`.

        When per-site allele data is available, the marginalized genotype
        likelihood for each CN state is added to the log-posterior.

        Args:
            data: :class:`DepthData` instance.
            map_estimates: Optional precomputed output from
                :meth:`get_map_estimates` to avoid recomputation.

        Returns:
            Dictionary with ``'cn_posterior'`` of shape
            ``(n_bins, n_samples, n_states)`` and, when sex–CN coupling is
            active, ``'sex_posterior'`` of shape ``(n_samples, 2)`` giving
            ``[P(XX), P(XY)]`` per sample.
        """
        logger.info("Running exact analytical discrete inference …")

        maps = map_estimates if map_estimates is not None else self.get_map_estimates(data)

        bin_bias = np.asarray(maps["bin_bias"]).squeeze()
        sample_var = np.asarray(maps["sample_var"]).squeeze()
        bin_var = np.asarray(maps["bin_var"]).squeeze()
        cn_probs = np.asarray(maps["cn_probs"]).squeeze()

        if cn_probs.ndim == 1:
            cn_probs = cn_probs.reshape(1, -1)

        obs = data.depth.detach().cpu().numpy()

        variance = sample_var[np.newaxis, :] + bin_var[:, np.newaxis]
        std = np.sqrt(np.maximum(variance, 1e-10))

        cn_states = np.arange(self.n_states, dtype=obs.dtype).reshape(-1, 1, 1)
        expected_depth = cn_states * bin_bias[np.newaxis, :, np.newaxis]

        obs_b = obs[np.newaxis, :, :]
        std_b = std[np.newaxis, :, :]

        log_lik = -0.5 * np.log(2 * np.pi * std_b ** 2) - (
            (obs_b - expected_depth) ** 2
        ) / (2 * std_b ** 2)

        log_prior = np.log(np.maximum(cn_probs.T[:, :, np.newaxis], 1e-10))
        base_log_unnorm = log_lik + log_prior  # (n_states, n_bins, n_samples)

        # ── marginalized allele-fraction likelihood ──────────────────────
        if data.site_alt is not None and self.af_weight > 0:
            site_alt_np = data.site_alt.detach().cpu().numpy()
            site_total_np = data.site_total.detach().cpu().numpy()
            site_pop_af_np = data.site_pop_af.detach().cpu().numpy()
            site_mask_np = data.site_mask.detach().cpu().numpy()

            for c in range(self.n_states):
                af_ll = _marginalized_af_log_lik_numpy(
                    site_alt_np, site_total_np, site_pop_af_np, site_mask_np,
                    cn_state=c, n_states=self.n_states,
                    concentration=self.af_concentration,
                )
                base_log_unnorm[c] += self.af_weight * af_ll

            n_sites = int(site_mask_np.any(axis=2).sum())
            logger.info(
                "AF marginalised likelihood: %d total sites, normalized "
                "per bin/sample (weight=%.2f)",
                n_sites, self.af_weight,
            )

        # ── sex-aware marginalisation ────────────────────────────────
        has_sex = (
            self.sex_cn_weight > 0
            and hasattr(data, "chr_type")
            and data.chr_type is not None
        )
        if has_sex:
            chr_type_np = data.chr_type.cpu().numpy()
            sex_cn_np = self.sex_cn_score.cpu().numpy()  # (2, 3, n_states)
            sex_prior_np = np.array(self.sex_prior, dtype=np.float64)

            # Per-bin weight normalised by chr-type bin count
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

            # Posterior over sex per sample
            log_sex = (
                np.log(np.maximum(sex_prior_np, 1e-10))[:, np.newaxis]
                + log_evidence
            )
            mx_s = np.max(log_sex, axis=0, keepdims=True)
            sex_post = np.exp(log_sex - mx_s)
            sex_post /= sex_post.sum(axis=0, keepdims=True)

            # Marginal CN posterior: weighted sum over sex states
            cn_marginal = np.zeros((self.n_states, n_bins, n_samp))
            for s in range(2):
                cn_marginal += (
                    sex_post[s][np.newaxis, np.newaxis, :]
                    * cn_post_given_sex[s]
                )

            cn_posterior = np.transpose(cn_marginal, (1, 2, 0)).astype(
                np.float32, copy=False,
            )
            logger.info(
                "Sex posterior (mean): P(XX)=%.3f  P(XY)=%.3f",
                sex_post[0].mean(), sex_post[1].mean(),
            )
            logger.info("Exact analytical discrete inference complete.")
            return {
                "cn_posterior": cn_posterior,
                "sex_posterior": sex_post.T.astype(np.float32, copy=False),
            }

        # ── fallback: no sex coupling ────────────────────────────────
        max_log = np.max(base_log_unnorm, axis=0, keepdims=True)
        exp_vals = np.exp(base_log_unnorm - max_log)
        posterior = exp_vals / np.sum(exp_vals, axis=0, keepdims=True)
        cn_posterior = np.transpose(posterior, (1, 2, 0)).astype(
            np.float32, copy=False,
        )

        logger.info("Exact analytical discrete inference complete.")
        return {"cn_posterior": cn_posterior}
