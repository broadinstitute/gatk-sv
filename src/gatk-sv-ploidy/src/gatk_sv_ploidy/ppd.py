"""
Posterior predictive check (ppd) subcommand.

Generates posterior predictive draws from the trained model and computes
summary statistics comparing the predictive distribution to observed data.
This is a model-checking diagnostic: if the model fits well, the observed
data should look like a typical sample from the posterior predictive
distribution.

For each bin × sample pair the predictive depth is a mixture of Gaussians:

.. math::

    p(\\tilde{d} | x) = \\sum_{c=0}^{5} p(c | x)
        \\cdot \\mathcal{N}(\\tilde{d} \\mid c \\cdot b_i,
                           \\sigma^2_{ij})

where *b_i* is the MAP ``bin_bias``, *σ²_{ij}* = ``bin_var_i`` +
``sample_var_j``, and *p(c | x)* is the analytical CN posterior.
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

from gatk_sv_ploidy.data import DepthData, load_site_data
from gatk_sv_ploidy.infer import load_inference_artifacts

logger = logging.getLogger(__name__)


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
    2. Draw depth from Normal(c * bin_bias_i, sqrt(bin_var_i + sample_var_j)).

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

    cn_probs = cn_posterior["cn_posterior"]  # (n_bins, n_samples, n_states)
    bin_bias = np.asarray(map_estimates["bin_bias"]).squeeze()
    sample_var = np.asarray(map_estimates["sample_var"]).squeeze()
    bin_var = np.asarray(map_estimates["bin_var"]).squeeze()

    n_bins, n_samples, n_states = cn_probs.shape

    # Precompute per-bin/sample variance and std
    variance = bin_var[:, np.newaxis] + sample_var[np.newaxis, :]  # (n_bins, n_samples)
    std = np.sqrt(np.maximum(variance, 1e-10))

    draws = np.empty((n_draws, n_bins, n_samples), dtype=np.float32)

    for d in range(n_draws):
        # Vectorised: sample CN states from the categorical posterior
        # Flatten probs to (n_bins * n_samples, n_states), draw, reshape
        flat_probs = cn_probs.reshape(-1, n_states)
        # Cumulative sum trick for vectorised categorical sampling
        cum = np.cumsum(flat_probs, axis=1)
        u = rng.uniform(size=(flat_probs.shape[0], 1))
        cn_sample = (u < cum).argmax(axis=1).reshape(n_bins, n_samples)

        # Expected depth = cn * bin_bias
        expected = cn_sample.astype(np.float32) * bin_bias[:, np.newaxis]

        # Draw from Normal
        draws[d] = rng.normal(expected, std).astype(np.float32)

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
            })

    return pd.DataFrame(rows)


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
            (ppd_bin_df["observed_depth"] < ppd_bin_df["ppd_q05"])
            | (ppd_bin_df["observed_depth"] > ppd_bin_df["ppd_q95"])
        )
    )
    # Fraction outside 50% predictive interval
    outside_50 = float(
        np.mean(
            (ppd_bin_df["observed_depth"] < ppd_bin_df["ppd_q25"])
            | (ppd_bin_df["observed_depth"] > ppd_bin_df["ppd_q75"])
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

    # ── load data ───────────────────────────────────────────────────────
    logger.info("Loading preprocessed depth: %s", args.input)
    df = pd.read_csv(args.input, sep="\t", index_col=0)

    sd = None
    if args.site_data:
        sd = load_site_data(args.site_data)

    data = DepthData(
        df, device=args.device, dtype=torch.float32,
        site_data=sd,
    )

    # ── load inference artifacts ────────────────────────────────────────
    logger.info("Loading inference artifacts: %s", args.artifacts)
    map_est, cn_post = load_inference_artifacts(args.artifacts)

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
