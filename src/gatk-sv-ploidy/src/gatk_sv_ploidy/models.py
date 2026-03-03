"""
Hierarchical Bayesian model for copy-number inference (Pyro).

Provides :class:`CNVModel` which wraps the probabilistic model, the
variational guide, training via SVI, MAP estimation, and discrete posterior
inference over copy-number states.
"""

from __future__ import annotations

import logging
from typing import Dict, List

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

from gatk_sv_ploidy.data import DepthData

logger = logging.getLogger(__name__)


class CNVModel:
    """Hierarchical Bayesian model for whole-genome CN detection.

    The generative model uses a Dirichlet-Categorical prior over CN states
    (0–5), per-bin bias and variance, and per-sample variance.  Observed
    normalised depth is drawn from a Normal likelihood whose mean is
    ``CN × bin_bias`` and whose variance is ``bin_var + sample_var``.

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
    ) -> None:
        """Pyro generative model.

        Args:
            depth: Observed normalised depth tensor (n_bins × n_samples).
            n_bins: Number of genomic bins.
            n_samples: Number of samples.
        """
        zero = torch.zeros(1, device=self.device, dtype=self.dtype)
        one = torch.ones(1, device=self.device, dtype=self.dtype)
        cn_states = torch.arange(
            0, self.n_states, device=self.device, dtype=self.dtype
        )

        plate_bins = pyro.plate("bins", n_bins, dim=-2, device=self.device)
        plate_samples = pyro.plate("samples", n_samples, dim=-1, device=self.device)

        # Per-sample variance
        with plate_samples:
            sample_var = pyro.sample(
                "sample_var", dist.Exponential(1.0 / self.var_sample)
            )

        # Per-bin latent variables
        with plate_bins:
            bin_bias = pyro.sample(
                "bin_bias", dist.LogNormal(zero, self.var_bias_bin)
            )
            bin_var = pyro.sample(
                "bin_var", dist.Exponential(1.0 / self.var_bin)
            )

            # Dirichlet-Categorical CN prior (heavily weights CN = 2)
            alpha = self.alpha_non_ref * one.expand(self.n_states)
            alpha[2] = self.alpha_ref
            cn_probs = pyro.sample("cn_probs", dist.Dirichlet(alpha))

        # Per-bin, per-sample observations
        with plate_bins, plate_samples:
            cn = pyro.sample("cn", dist.Categorical(cn_probs))
            expected = Vindex(cn_states)[cn] * bin_bias
            variance = sample_var + bin_var
            pyro.sample("obs", dist.Normal(expected, variance.sqrt()), obs=depth)

    # ── training ────────────────────────────────────────────────────────

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
                "lr_lambda": lambda k: lr_min
                + (lr_init - lr_min) * np.exp(-k / lr_decay),
            }
        )

        elbo = JitTraceEnum_ELBO() if jit else TraceEnum_ELBO()
        svi = SVI(self.model, self.guide, optim=scheduler, loss=elbo)

        logger.info("Training for up to %d iterations …", max_iter)

        best_loss = float("inf")
        patience_ctr = 0

        with tqdm(range(max_iter), desc="Training", unit="epoch") as pbar:
            for epoch in pbar:
                loss = svi.step(
                    depth=data.depth, n_bins=data.n_bins, n_samples=data.n_samples
                )
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

        guide_trace = poutine.trace(self.guide).get_trace(
            depth=data.depth, n_bins=data.n_bins, n_samples=data.n_samples
        )
        replayed = poutine.replay(self.model, trace=guide_trace)
        inferred = infer_discrete(replayed, temperature=0, first_available_dim=-3)
        trace = poutine.trace(inferred).get_trace(
            depth=data.depth, n_bins=data.n_bins, n_samples=data.n_samples
        )

        estimates: Dict[str, np.ndarray] = {}
        for site in ("bin_bias", "sample_var", "bin_var", "cn_probs"):
            estimates[site] = (
                guide_trace.nodes[site]["value"].detach().cpu().numpy()
            )
        estimates["cn"] = trace.nodes["cn"]["value"].detach().cpu().numpy()
        return estimates

    def run_discrete_inference(
        self,
        data: DepthData,
        n_samples: int = 1000,
    ) -> Dict[str, np.ndarray]:
        """Sample discrete CN posterior by repeated forward passes.

        Args:
            data: :class:`DepthData` instance.
            n_samples: Number of posterior samples.

        Returns:
            Dictionary with key ``'cn_posterior'`` mapping to an array of shape
            ``(n_bins, n_data_samples, n_states)`` containing per-state
            frequencies.
        """
        logger.info("Running discrete inference (%d samples) …", n_samples)

        guide_trace = poutine.trace(self.guide).get_trace(
            depth=data.depth, n_bins=data.n_bins, n_samples=data.n_samples
        )
        replayed = poutine.replay(self.model, trace=guide_trace)

        cn_draws: list[np.ndarray] = []
        with torch.no_grad():
            for _ in tqdm(range(n_samples), desc="Discrete inference", unit="sample"):
                inferred = infer_discrete(
                    replayed, temperature=1, first_available_dim=-3
                )
                trace = poutine.trace(inferred).get_trace(
                    depth=data.depth,
                    n_bins=data.n_bins,
                    n_samples=data.n_samples,
                )
                cn_draws.append(trace.nodes["cn"]["value"].detach().cpu().numpy())

        cn_stack = np.array(cn_draws)  # (n_samples, n_bins, n_data_samples)
        cn_freq = np.zeros((data.n_bins, data.n_samples, self.n_states))
        for state in range(self.n_states):
            cn_freq[..., state] = (cn_stack == state).mean(axis=0)

        logger.info("Discrete inference complete.")
        return {"cn_posterior": cn_freq}
