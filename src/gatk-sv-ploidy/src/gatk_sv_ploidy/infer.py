"""
Infer subcommand — train model and run copy-number inference.

Loads preprocessed depth data, trains the Pyro CNV model, obtains MAP
estimates and discrete posterior samples, then writes per-bin and
per-chromosome summary statistics.
"""

from __future__ import annotations

import argparse
import logging
import os
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import pyro
import torch
from scipy import stats

from gatk_sv_ploidy.data import DepthData
from gatk_sv_ploidy.models import CNVModel

logger = logging.getLogger(__name__)


# ── summary statistics ──────────────────────────────────────────────────────


def _get_chr_info(
    data: DepthData,
    cn: np.ndarray,
    cn_probs: np.ndarray,
    chr_name: str,
    sample_idx: int,
    n_states: int = 6,
) -> Tuple[int, float, int]:
    """Return (most_common_cn, mean_prob, n_bins) for one chromosome/sample."""
    mask = data.chr == chr_name
    n_bins = int(mask.sum())
    if n_bins == 0:
        return 0, 0.0, 0
    chr_cn = cn[mask, sample_idx]
    counts = np.bincount(chr_cn, minlength=n_states)
    best_cn = int(np.argmax(counts))
    mean_prob = float(cn_probs[mask, sample_idx, best_cn].mean())
    return best_cn, mean_prob, n_bins


def detect_aneuploidies(
    data: DepthData,
    map_estimates: Dict[str, np.ndarray],
    cn_posterior: Dict[str, np.ndarray],
    prob_threshold: float = 0.5,
) -> Dict[int, List[Tuple[str, int, float]]]:
    """Detect per-chromosome aneuploidies for every sample.

    Args:
        data: :class:`DepthData` used for inference.
        map_estimates: MAP estimates from :meth:`CNVModel.get_map_estimates`.
        cn_posterior: Posterior from :meth:`CNVModel.run_discrete_inference`.
        prob_threshold: Minimum mean CN probability to call an aneuploidy.

    Returns:
        Dictionary mapping sample index → list of
        ``(chr_name, cn_state, mean_prob)`` tuples for aneuploid chromosomes.
    """
    cn = map_estimates["cn"]
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
            best_cn, mean_prob, _ = _get_chr_info(data, cn, cn_probs, chr_name, si)
            if best_cn != 2 and mean_prob > prob_threshold:
                aneuploid[si].append((chr_name, best_cn, mean_prob))

    # Sex chromosomes: aneuploidy when karyotype is not XX or XY
    for si in range(data.n_samples):
        x_cn, x_prob, x_bins = _get_chr_info(data, cn, cn_probs, "chrX", si)
        y_cn, y_prob, y_bins = _get_chr_info(data, cn, cn_probs, "chrY", si)

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
) -> pd.DataFrame:
    """Build a per-bin, per-sample results DataFrame.

    Args:
        data: :class:`DepthData` instance.
        map_estimates: MAP estimates dict.
        cn_posterior: Discrete posterior dict.

    Returns:
        DataFrame with columns for chr, start, end, sample, observed depth,
        MAP CN, per-state probabilities, and model parameters.
    """
    cn_probs = cn_posterior["cn_posterior"]
    rows: list[dict] = []

    for i in range(data.n_bins):
        for j in range(data.n_samples):
            prob = cn_probs[i, j, :]
            rows.append(
                {
                    "chr": data.chr[i],
                    "start": int(data.start[i]),
                    "end": int(data.end[i]),
                    "sample": data.sample_ids[j],
                    "observed_depth": float(data.depth[i, j].cpu().numpy()),
                    "cn_map": int(map_estimates["cn"][i, j]),
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
            )

    return pd.DataFrame(rows)


def build_chromosome_stats(
    data: DepthData,
    map_estimates: Dict[str, np.ndarray],
    cn_posterior: Dict[str, np.ndarray],
    aneuploid_map: Dict[int, List[Tuple[str, int, float]]],
) -> pd.DataFrame:
    """Build per-chromosome, per-sample summary statistics.

    Args:
        data: :class:`DepthData` instance.
        map_estimates: MAP estimates dict.
        cn_posterior: Discrete posterior dict.
        aneuploid_map: Output of :func:`detect_aneuploidies`.

    Returns:
        DataFrame with one row per (sample, chromosome).
    """
    cn = map_estimates["cn"]
    cn_probs = cn_posterior["cn_posterior"]
    unique_chrs = np.unique(data.chr)
    rows: list[dict] = []

    for si in range(data.n_samples):
        aneu_set = {c for c, _, _ in aneuploid_map.get(si, [])}
        for chr_name in unique_chrs:
            mask = data.chr == chr_name
            n_bins = int(mask.sum())
            if n_bins == 0:
                continue

            chr_cn = cn[mask, si]
            counts = np.bincount(chr_cn, minlength=6)
            best_cn = int(np.argmax(counts))
            mean_prob = float(cn_probs[mask, si, best_cn].mean())
            depths = data.depth[mask, si].cpu().numpy()

            rows.append(
                {
                    "sample": data.sample_ids[si],
                    "chromosome": chr_name,
                    "copy_number": best_cn,
                    "mean_cn_probability": mean_prob,
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
            )

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

    # Model priors
    g = p.add_argument_group("model priors")
    g.add_argument("--alpha-ref", type=float, default=1.0,
                   help="Dirichlet concentration for CN=2")
    g.add_argument("--alpha-non-ref", type=float, default=1.0,
                   help="Dirichlet concentration for other CN states")
    g.add_argument("--var-bias-bin", type=float, default=0.01,
                   help="LogNormal scale for per-bin bias")
    g.add_argument("--var-sample", type=float, default=0.001,
                   help="Exponential rate for per-sample variance")
    g.add_argument("--var-bin", type=float, default=0.001,
                   help="Exponential rate for per-bin variance")
    g.add_argument("--guide-type", choices=["delta", "diagonal"], default="delta",
                   help="Variational guide type")

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
    g.add_argument("--n-discrete-samples", type=int, default=1000)
    g.add_argument("--prob-threshold", type=float, default=0.5,
                   help="Min mean CN probability for aneuploidy call")

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

    device = args.device
    data = DepthData(df, device=device, dtype=torch.float32)

    # ── build & train model ─────────────────────────────────────────────
    model = CNVModel(
        n_states=6,
        alpha_ref=args.alpha_ref,
        alpha_non_ref=args.alpha_non_ref,
        var_bias_bin=args.var_bias_bin,
        var_sample=args.var_sample,
        var_bin=args.var_bin,
        device=device,
        dtype=torch.float32,
        guide_type=args.guide_type,
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

    # ── MAP + discrete inference ────────────────────────────────────────
    map_est = model.get_map_estimates(data)
    cn_post = model.run_discrete_inference(
        data, n_samples=args.n_discrete_samples
    )

    # ── detect aneuploidies ─────────────────────────────────────────────
    aneuploid_map = detect_aneuploidies(
        data, map_est, cn_post, prob_threshold=args.prob_threshold
    )

    # ── write bin stats ─────────────────────────────────────────────────
    bin_df = build_bin_stats(data, map_est, cn_post)
    bin_path = os.path.join(args.output_dir, "bin_stats.tsv.gz")
    bin_df.to_csv(bin_path, sep="\t", index=False, compression="gzip")
    logger.info("Bin statistics saved to %s", bin_path)

    # ── write chromosome stats ──────────────────────────────────────────
    chr_df = build_chromosome_stats(data, map_est, cn_post, aneuploid_map)
    chr_path = os.path.join(args.output_dir, "chromosome_stats.tsv")
    chr_df.to_csv(chr_path, sep="\t", index=False)
    logger.info("Chromosome statistics saved to %s", chr_path)


if __name__ == "__main__":
    main()
