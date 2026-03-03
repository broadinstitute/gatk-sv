"""
Preprocess subcommand — read, normalise, and filter depth data.

Reads a raw depth TSV, normalises each sample by its autosomal median so that
diploid depth ≈ 2.0, filters low-quality bins by median/MAD thresholds, and
writes the preprocessed matrix to disk.
"""

from __future__ import annotations

import argparse
import logging
import os

import numpy as np
import pandas as pd
from scipy import stats

from gatk_sv_ploidy._util import get_sample_columns
from gatk_sv_ploidy.data import read_depth_tsv

logger = logging.getLogger(__name__)


# ── normalisation ───────────────────────────────────────────────────────────


def normalise_depth(df: pd.DataFrame) -> pd.DataFrame:
    """Normalise per-sample depth so that autosomal CN = 2 corresponds to 2.0.

    For each sample the median depth across autosomal bins is computed and all
    bins (including sex chromosomes) are scaled by ``2 / median``.

    Args:
        df: DataFrame with ``Chr`` metadata column and per-sample depth columns.

    Returns:
        Copy of *df* with normalised depth values.
    """
    sample_cols = get_sample_columns(df)
    autosome_mask = ~df["Chr"].isin(["chrX", "chrY"])
    medians = np.median(df.loc[autosome_mask, sample_cols].values, axis=0)

    logger.info(
        "Autosomal medians: min=%.3f, max=%.3f, mean=%.3f",
        medians.min(),
        medians.max(),
        medians.mean(),
    )

    df = df.copy()
    df[sample_cols] = 2.0 * df[sample_cols].values / medians[np.newaxis, :]
    return df


# ── bin quality filtering ───────────────────────────────────────────────────


def filter_low_quality_bins(
    df: pd.DataFrame,
    *,
    autosome_median_min: float = 1.0,
    autosome_median_max: float = 3.0,
    autosome_mad_max: float = 2.0,
    chrX_median_min: float = 0.0,
    chrX_median_max: float = 3.0,
    chrX_mad_max: float = 2.0,
    chrY_median_min: float = 0.0,
    chrY_median_max: float = 3.0,
    chrY_mad_max: float = 2.0,
    min_bins_per_chr: int = 10,
) -> pd.DataFrame:
    """Filter bins whose cross-sample statistics fall outside thresholds.

    Thresholds are applied separately for autosomes, chrX, and chrY because
    expected depth distributions differ across chromosome types.

    Args:
        df: Normalised depth DataFrame.
        autosome_median_min/max: Median depth range for autosomal bins.
        autosome_mad_max: Maximum MAD for autosomal bins.
        chrX_median_min/max: Median depth range for chrX bins.
        chrX_mad_max: Maximum MAD for chrX bins.
        chrY_median_min/max: Median depth range for chrY bins.
        chrY_mad_max: Maximum MAD for chrY bins.
        min_bins_per_chr: Raise if any chromosome has fewer bins after
            filtering.

    Returns:
        Filtered DataFrame.

    Raises:
        ValueError: If a chromosome has fewer than *min_bins_per_chr* bins.
    """
    sample_cols = get_sample_columns(df)
    depths = df[sample_cols].values
    medians = np.median(depths, axis=1)
    mads = stats.median_abs_deviation(depths, axis=1)

    logger.info("Starting bins: %d", len(df))

    keep = np.ones(len(df), dtype=bool)

    # -- per-chromosome-type thresholds --
    _thresholds = {
        "autosome": (
            ~df["Chr"].isin(["chrX", "chrY"]),
            autosome_median_min,
            autosome_median_max,
            autosome_mad_max,
        ),
        "chrX": (
            df["Chr"] == "chrX",
            chrX_median_min,
            chrX_median_max,
            chrX_mad_max,
        ),
        "chrY": (
            df["Chr"] == "chrY",
            chrY_median_min,
            chrY_median_max,
            chrY_mad_max,
        ),
    }

    for label, (mask, med_min, med_max, mad_max) in _thresholds.items():
        if not mask.any():
            continue
        ok = (medians >= med_min) & (medians <= med_max) & (mads <= mad_max)
        n_before = int(mask.sum())
        keep[mask.values] &= ok[mask.values]
        n_after = int((mask & keep).sum())
        logger.info(
            "%s: median [%.1f, %.1f], MAD ≤ %.1f → %d / %d bins kept",
            label,
            med_min,
            med_max,
            mad_max,
            n_after,
            n_before,
        )

    df_out = df[keep].copy()
    logger.info("Bins after filtering: %d (removed %d)", len(df_out), len(df) - len(df_out))

    # Verify sufficient bins per chromosome
    counts = df_out.groupby("Chr").size()
    bad = counts[counts < min_bins_per_chr]
    if len(bad) > 0:
        msg = "; ".join(f"{c}={n}" for c, n in bad.items())
        raise ValueError(
            f"Insufficient bins after filtering ({msg}). "
            "Consider --skip-bin-filter or relaxing thresholds."
        )

    return df_out


# ── CLI ─────────────────────────────────────────────────────────────────────


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the preprocess subcommand."""
    p = argparse.ArgumentParser(
        description="Preprocess depth data for aneuploidy inference",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "-i", "--input", required=True,
        help="Input TSV file with raw read-depth (bins × samples)",
    )
    p.add_argument(
        "-o", "--output-dir", required=True,
        help="Output directory",
    )
    p.add_argument(
        "--viable-only", action="store_true", default=False,
        help="Subset to chr13, chr18, chr21, chrX, chrY only",
    )
    p.add_argument(
        "--skip-bin-filter", action="store_true", default=False,
        help="Skip bin quality filtering",
    )

    # Filter thresholds
    g = p.add_argument_group("bin-filter thresholds")
    g.add_argument("--autosome-median-min", type=float, default=1.0)
    g.add_argument("--autosome-median-max", type=float, default=3.0)
    g.add_argument("--autosome-mad-max", type=float, default=2.0)
    g.add_argument("--chrX-median-min", type=float, default=0.0)
    g.add_argument("--chrX-median-max", type=float, default=3.0)
    g.add_argument("--chrX-mad-max", type=float, default=2.0)
    g.add_argument("--chrY-median-min", type=float, default=0.0)
    g.add_argument("--chrY-median-max", type=float, default=3.0)
    g.add_argument("--chrY-mad-max", type=float, default=2.0)

    return p.parse_args()


def main() -> None:
    """Entry point for ``gatk-sv-ploidy preprocess``."""
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    # 1. Read
    df = read_depth_tsv(args.input)

    # 2. Optionally subset to viable trisomy chromosomes
    if args.viable_only:
        viable = {"chr13", "chr18", "chr21", "chrX", "chrY"}
        n_before = len(df)
        df = df[df["Chr"].isin(viable)]
        logger.info(
            "Viable-only filter: %d → %d bins (%s)",
            n_before,
            len(df),
            sorted(df["Chr"].unique()),
        )

    # 3. Normalise
    df = normalise_depth(df)

    # 4. Filter
    if args.skip_bin_filter:
        logger.info("Skipping bin quality filtering (--skip-bin-filter)")
    else:
        df = filter_low_quality_bins(
            df,
            autosome_median_min=args.autosome_median_min,
            autosome_median_max=args.autosome_median_max,
            autosome_mad_max=args.autosome_mad_max,
            chrX_median_min=args.chrX_median_min,
            chrX_median_max=args.chrX_median_max,
            chrX_mad_max=args.chrX_mad_max,
            chrY_median_min=args.chrY_median_min,
            chrY_median_max=args.chrY_median_max,
            chrY_mad_max=args.chrY_mad_max,
        )

    # 5. Write preprocessed depth
    out_path = os.path.join(args.output_dir, "preprocessed_depth.tsv")
    df.to_csv(out_path, sep="\t")
    logger.info("Preprocessed depth written to %s", out_path)


if __name__ == "__main__":
    main()
