"""
Shared utilities for gatk-sv-ploidy.

Small helper functions used by multiple subcommands.
"""

from __future__ import annotations

import logging
import os
from typing import List, Optional, Sequence

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from gatk_sv_ploidy._plot_style import DEFAULT_RASTER_DPI, save_publication_figure

logger = logging.getLogger(__name__)

# ── column / chromosome constants ───────────────────────────────────────────

METADATA_COLS = frozenset({
    "Chr",
    "Start",
    "End",
    "BinLengthBp",
    "source_file",
    "Bin",
})
"""Columns that are *not* per-sample depth values."""

CHR_ORDER = {f"chr{i}": i for i in range(1, 23)}
CHR_ORDER["chrX"] = 23
CHR_ORDER["chrY"] = 24

AUTOSOME_NAMES = [f"chr{i}" for i in range(1, 23)]
SEX_CHR_NAMES = ["chrX", "chrY"]

MAX_GENOTYPE_STATES = 6
"""Maximum number of CN states modelled (CN 0 … 5)."""

DEFAULT_AF_CONCENTRATION = 50.0
"""Default Beta-Binomial concentration for the per-site allele model."""

DEFAULT_AF_WEIGHT = 0.25
"""Default relative weight of the allele-fraction likelihood term."""
DEPTH_SPACES = ("normalized", "raw")
"""Supported depth/count spaces for the depth matrix input."""
BINQ_FIELD_OPTIONS = ("auto", "BINQ15", "BINQ20", "CALLQ15", "CALLQ20")
"""Supported BINQ/CALLQ selector values for filtering and plotting."""

NEGATIVE_BINOMIAL_OBS_LIKELIHOOD = "negative_binomial"
"""Observation likelihood name used for raw-count negative-binomial fits."""
MAX_CNQ = 99
"""Maximum phred-scaled CN quality score reported by the ploidy model."""

_BASELINE_PLOIDY_LABELS = {
    1: "HAPLOID",
    2: "DIPLOID",
    3: "TRIPLOID",
    4: "TETRAPLOID",
}


def expected_allosome_copy_number_pairs(
    autosomal_baseline_cn: int,
) -> tuple[tuple[int, int], tuple[int, int]]:
    """Return female-like and male-like chrX/chrY CN pairs for a baseline CN."""
    baseline = int(autosomal_baseline_cn)
    if baseline < 1:
        raise ValueError("autosomal_baseline_cn must be at least 1.")

    female_like = (baseline, 0)
    male_y_cn = max(1, baseline // 2)
    male_y_cn = min(male_y_cn, baseline)
    male_like = (baseline - male_y_cn, male_y_cn)
    return female_like, male_like


def is_expected_allosome_copy_number_pair(
    x_cn: int,
    y_cn: int,
    autosomal_baseline_cn: int,
) -> bool:
    """Return whether chrX/chrY CN matches a baseline-aware allosome pair."""
    return (
        int(x_cn),
        int(y_cn),
    ) in set(expected_allosome_copy_number_pairs(autosomal_baseline_cn))


def baseline_ploidy_label(autosomal_baseline_cn: int) -> str:
    """Return the standard label for a sample's autosomal baseline CN."""
    baseline = int(autosomal_baseline_cn)
    return _BASELINE_PLOIDY_LABELS.get(baseline, f"BASELINE_CN{baseline}")


def resolve_binq_field(
    bin_quality_df: pd.DataFrame,
    requested_field: str,
) -> str:
    """Resolve a quality field, defaulting auto to BINQ20 when available."""
    if requested_field != "auto":
        return requested_field
    if "BINQ20" in bin_quality_df.columns:
        return "BINQ20"
    if "CALLQ20" in bin_quality_df.columns:
        return "CALLQ20"
    return "BINQ20"


# ── dataframe helpers ───────────────────────────────────────────────────────

def get_sample_columns(df: pd.DataFrame) -> List[str]:
    """Return column names that represent per-sample depth values.

    Every column whose name is not in :data:`METADATA_COLS` is considered a
    sample column.

    Args:
        df: DataFrame whose columns include both metadata and sample IDs.

    Returns:
        Ordered list of sample column names.
    """
    return [c for c in df.columns if c not in METADATA_COLS]


def get_chromosome_type(chromosome: str) -> str:
    """Classify *chromosome* as ``'Autosomal'``, ``'chrX'``, or ``'chrY'``.

    Args:
        chromosome: Chromosome name (e.g. ``'chr1'``, ``'chrX'``).

    Returns:
        One of ``'chrX'``, ``'chrY'``, or ``'Autosomal'``.
    """
    if chromosome in ("chrX", "X"):
        return "chrX"
    if chromosome in ("chrY", "Y"):
        return "chrY"
    return "Autosomal"


def format_column_name(col: str) -> str:
    """Format a column name for display (``'mean_depth'`` → ``'Mean Depth'``).

    Args:
        col: Raw column name.

    Returns:
        Title-cased string with underscores replaced by spaces.
    """
    return col.replace("_", " ").title()


def load_exclusion_ids(path: str) -> List[str]:
    """Load sample IDs to exclude (one per line).

    Args:
        path: Path to a text file.

    Returns:
        List of sample-ID strings.
    """
    with open(path) as fh:
        return [line.strip() for line in fh if line.strip()]


def safe_path_label(path: str) -> str:
    """Return a display-safe path label for logs."""
    path_str = str(path)
    stripped = path_str.rstrip(os.sep)
    if not stripped:
        return path_str
    return os.path.basename(stripped) or stripped


def summarize_numeric_array(
    values: np.ndarray | Sequence[float],
) -> dict[str, float | int]:
    """Summarize numeric values with privacy-safe quantiles."""
    arr = np.asarray(values, dtype=np.float64).reshape(-1)
    finite = arr[np.isfinite(arr)]
    summary: dict[str, float | int] = {
        "n_total": int(arr.size),
        "n_finite": int(finite.size),
    }
    if finite.size == 0:
        return summary

    quantiles = np.quantile(finite, [0.05, 0.25, 0.50, 0.75, 0.95])
    summary.update(
        {
            "mean": float(finite.mean()),
            "sd": float(finite.std(ddof=0)),
            "p05": float(quantiles[0]),
            "p25": float(quantiles[1]),
            "median": float(quantiles[2]),
            "p75": float(quantiles[3]),
            "p95": float(quantiles[4]),
        }
    )
    return summary


def coerce_bin_epsilon_matrix(
    values: np.ndarray | Sequence[float] | float | None,
    n_bins: int,
    n_samples: int,
    dtype: np.dtype | type = np.float64,
) -> np.ndarray:
    """Return ``bin_epsilon`` as a ``(n_bins, n_samples)`` matrix.

    Current runs store either a scalar epsilon floor, one epsilon value per
    bin, or a full ``(n_bins, n_samples)`` matrix. This helper expands those
    forms to the effective per-bin matrix.
    """
    if values is None:
        return np.zeros((n_bins, n_samples), dtype=dtype)

    arr = np.asarray(values, dtype=dtype)
    if arr.ndim > 2:
        arr = np.squeeze(arr)
    if arr.ndim == 0:
        return np.full((n_bins, n_samples), float(arr), dtype=dtype)

    if arr.ndim == 1:
        if n_bins == 1 and arr.shape[0] == n_samples:
            return arr.reshape(1, n_samples).astype(dtype, copy=False)
        if arr.shape[0] == n_bins:
            return np.broadcast_to(arr[:, np.newaxis], (n_bins, n_samples)).astype(
                dtype,
                copy=False,
            )
        raise ValueError(
            "bin_epsilon must have length n_bins or shape "
            "(n_bins, n_samples)."
        )

    if arr.ndim == 2:
        if arr.shape == (n_bins, n_samples):
            return arr.astype(dtype, copy=False)
        if arr.shape == (n_bins, 1):
            return np.broadcast_to(arr, (n_bins, n_samples)).astype(
                dtype,
                copy=False,
            )
        if arr.shape == (1, n_samples):
            return np.broadcast_to(arr, (n_bins, n_samples)).astype(
                dtype,
                copy=False,
            )
        raise ValueError(
            "bin_epsilon must have shape (n_bins, n_samples)."
        )

    raise ValueError("bin_epsilon must be scalar, 1D, or 2D.")


def _coerce_background_bin_factor_matrix(
    values: np.ndarray | Sequence[float],
    n_bins: int,
    dtype: np.dtype | type,
) -> np.ndarray:
    """Return normalized per-bin background factors as an ``(n_bins, k)`` matrix."""
    arr = np.asarray(values, dtype=dtype)
    if arr.ndim > 2:
        arr = np.squeeze(arr)
    if arr.ndim == 1:
        if arr.shape[0] != n_bins:
            raise ValueError("background_bin_factors must have length n_bins.")
        arr = arr[:, np.newaxis]
    if arr.ndim != 2 or arr.shape[0] != n_bins:
        raise ValueError("background_bin_factors must have shape (n_bins, k).")
    if arr.shape[1] == 0:
        return arr.astype(dtype, copy=False)
    column_means = np.maximum(arr.mean(axis=0, keepdims=True), 1e-8)
    return (arr / column_means).astype(dtype, copy=False)


def _coerce_background_sample_factor_matrix(
    values: np.ndarray | Sequence[float],
    n_samples: int,
    n_factors: int,
    dtype: np.dtype | type,
) -> np.ndarray:
    """Return per-sample background factors as a ``(k, n_samples)`` matrix."""
    arr = np.asarray(values, dtype=dtype)
    if arr.ndim > 2:
        arr = np.squeeze(arr)
    if arr.ndim == 1:
        if n_factors != 1 or arr.shape[0] != n_samples:
            raise ValueError(
                "background_sample_factors must have shape (k, n_samples)."
            )
        arr = arr[np.newaxis, :]
    if arr.ndim != 2:
        raise ValueError("background_sample_factors must be 1D or 2D.")
    if arr.shape == (n_factors, n_samples):
        return arr.astype(dtype, copy=False)
    if arr.shape == (n_samples, n_factors) and n_samples != n_factors:
        return arr.T.astype(dtype, copy=False)
    raise ValueError("background_sample_factors must have shape (k, n_samples).")


def compose_additive_background_matrix(
    bin_epsilon: np.ndarray | Sequence[float] | float | None,
    n_bins: int,
    n_samples: int,
    dtype: np.dtype | type = np.float64,
    background_bin_factors: np.ndarray | Sequence[float] | None = None,
    background_sample_factors: np.ndarray | Sequence[float] | None = None,
) -> np.ndarray:
    """Compose the effective additive background matrix.

    This includes the ``bin_epsilon`` floor, plus any low-rank additive
    background factors stored as ``background_bin_factors`` and
    ``background_sample_factors``.
    """
    additive = coerce_bin_epsilon_matrix(
        bin_epsilon,
        n_bins,
        n_samples,
        dtype=dtype,
    )
    if background_bin_factors is None and background_sample_factors is None:
        return additive
    if background_bin_factors is None or background_sample_factors is None:
        raise ValueError(
            "background_bin_factors and background_sample_factors must both be provided."
        )
    bin_factor_matrix = _coerce_background_bin_factor_matrix(
        background_bin_factors,
        n_bins,
        dtype,
    )
    sample_factor_matrix = _coerce_background_sample_factor_matrix(
        background_sample_factors,
        n_samples,
        bin_factor_matrix.shape[1],
        dtype,
    )
    if bin_factor_matrix.shape[1] == 0:
        return additive
    # ``np.matmul`` emits spurious floating-point RuntimeWarnings on some
    # macOS NumPy/BLAS builds here even when both inputs and the result are
    # finite. ``np.dot`` is equivalent for these 2D arrays and avoids that
    # false-positive warning path.
    return additive + np.dot(bin_factor_matrix, sample_factor_matrix)


def format_numeric_summary(
    label: str,
    values: np.ndarray | Sequence[float],
    precision: int = 3,
) -> str:
    """Format a privacy-safe numeric summary string for logs."""
    summary = summarize_numeric_array(values)
    n_total = int(summary["n_total"])
    n_finite = int(summary["n_finite"])
    if n_total == 0:
        return f"{label}: n=0"
    if n_finite == 0:
        return f"{label}: n={n_total}, finite=0"

    fmt = f"{{:.{precision}f}}"
    return (
        f"{label}: n={n_total}, finite={n_finite}, "
        f"mean={fmt.format(float(summary['mean']))}, "
        f"sd={fmt.format(float(summary['sd']))}, "
        f"p05={fmt.format(float(summary['p05']))}, "
        f"p25={fmt.format(float(summary['p25']))}, "
        f"median={fmt.format(float(summary['median']))}, "
        f"p75={fmt.format(float(summary['p75']))}, "
        f"p95={fmt.format(float(summary['p95']))}"
    )


def format_count_summary(
    label: str,
    counts: np.ndarray | Sequence[int],
    names: Optional[Sequence[str]] = None,
) -> str:
    """Format count and fraction summaries for categorical diagnostics."""
    count_arr = np.asarray(counts, dtype=np.int64).reshape(-1)
    total = int(count_arr.sum())
    if names is None:
        names = [str(i) for i in range(count_arr.size)]
    if len(names) != count_arr.size:
        raise ValueError("names must match the number of counts")
    if total == 0:
        return f"{label}: total=0"

    parts = [
        f"{name}={int(count)} ({count / total:.1%})"
        for name, count in zip(names, count_arr)
    ]
    return f"{label}: total={total}, " + ", ".join(parts)


def compute_cnq_from_probabilities(
    probabilities: np.ndarray,
    max_score: int = MAX_CNQ,
) -> np.ndarray:
    """Compute per-bin CN quality from the top-two posterior probabilities.

    CNQ is defined as ``round(-10 * log10(1 - p_diff))``, where
    ``p_diff`` is the difference between the largest and second-largest CN
    posterior probabilities. Scores are clipped to ``[0, max_score]``.

    Args:
        probabilities: Array whose last dimension enumerates CN states.
        max_score: Maximum reported CNQ.

    Returns:
        Integer NumPy array with the same leading dimensions as
        *probabilities* and one score per posterior vector.
    """
    probs = np.asarray(probabilities, dtype=np.float64)
    if probs.shape[-1] == 0:
        raise ValueError("CN posterior array must include at least one state")

    if probs.shape[-1] == 1:
        best = probs[..., 0]
        second = np.zeros_like(best)
    else:
        kth = probs.shape[-1] - 2
        top2 = np.partition(probs, kth=kth, axis=-1)[..., -2:]
        best = top2[..., 1]
        second = top2[..., 0]

    p_diff = np.clip(best - second, 0.0, 1.0)
    error_prob = np.clip(1.0 - p_diff, np.finfo(np.float64).tiny, 1.0)
    cnq = np.rint(-10.0 * np.log10(error_prob))
    return np.clip(cnq, 0, max_score).astype(np.int16, copy=False)


def compute_plq_from_probabilities(
    probabilities: np.ndarray,
    max_score: int = MAX_CNQ,
) -> np.ndarray:
    """Compute ploidy quality from the top-two posterior-probability ratio.

    PLQ is the phred-scaled log-likelihood ratio between the most likely and
    second most likely ploidy states:

    ``round(10 * log10(p_best / p_second))``

    Scores are clipped to ``[0, max_score]``.

    Args:
        probabilities: Array whose last dimension enumerates ploidy states.
        max_score: Maximum reported PLQ.

    Returns:
        Integer NumPy array with the same leading dimensions as
        *probabilities* and one score per posterior vector.
    """
    probs = np.asarray(probabilities, dtype=np.float64)
    if probs.shape[-1] == 0:
        raise ValueError(
            "PLQ posterior array must include at least one state"
        )

    if probs.shape[-1] == 1:
        return np.full(probs.shape[:-1], max_score, dtype=np.int16)

    kth = probs.shape[-1] - 2
    top2 = np.partition(probs, kth=kth, axis=-1)[..., -2:]
    best = top2[..., 1]
    second = top2[..., 0]

    tiny = np.finfo(np.float64).tiny
    ratio = np.maximum(best, tiny) / np.maximum(second, tiny)
    plq = np.rint(10.0 * np.log10(ratio))
    plq = np.where(best > tiny, plq, 0.0)
    return np.clip(plq, 0, max_score).astype(np.int16, copy=False)


def summarize_contig_ploidy_from_bin_calls(
    cn_map: np.ndarray,
    cnq: np.ndarray,
    n_states: int = MAX_GENOTYPE_STATES,
    max_score: int = MAX_CNQ,
) -> tuple[int, float, np.ndarray, int]:
    """Summarize contig ploidy from per-bin MAP calls and CNQ values.

    The contig ploidy is the most frequent bin-level MAP copy number. The
    returned support value is the fraction of bins assigned to that winning CN.
    PLQ is the rounded mean CNQ across all bins on the contig.
    """
    cn_map_arr = np.asarray(cn_map, dtype=np.int64).reshape(-1)
    cnq_arr = np.asarray(cnq, dtype=np.float64).reshape(-1)

    if cn_map_arr.size == 0:
        raise ValueError("cn_map must include at least one bin")
    if cnq_arr.size != cn_map_arr.size:
        raise ValueError("cnq must have one value per bin")
    if np.any(cn_map_arr < 0) or np.any(cn_map_arr >= n_states):
        raise ValueError("cn_map contains copy-number states outside the model range")

    counts = np.bincount(cn_map_arr, minlength=n_states).astype(np.float64, copy=False)
    total_bins = float(cn_map_arr.size)
    ploidy_fractions = counts / total_bins
    best_cn = int(np.argmax(counts))
    best_fraction = float(ploidy_fractions[best_cn])

    finite_cnq = cnq_arr[np.isfinite(cnq_arr)]
    if finite_cnq.size == 0:
        plq = 0
    else:
        plq = int(np.clip(np.rint(finite_cnq.mean()), 0, max_score))

    return best_cn, best_fraction, ploidy_fractions.astype(np.float32, copy=False), plq


def compute_contig_posterior_from_bin_posteriors(
    bin_posteriors: np.ndarray,
) -> np.ndarray:
    """Aggregate per-bin CN posteriors into a contig-wide constant-state posterior.

    For a chromosome with ``n`` bins, each ploidy hypothesis fixes all bins to
    the same copy-number state. The hypothesis score is therefore the product
    of the corresponding per-bin posterior probabilities. This helper performs
    the calculation in log space and normalizes the resulting state scores.

    Args:
        bin_posteriors: Array of shape ``(n_bins, n_states)`` containing
            per-bin posterior probabilities.

    Returns:
        Normalized posterior probabilities of shape ``(n_states,)``.
    """
    probs = np.asarray(bin_posteriors, dtype=np.float64)
    if probs.ndim != 2:
        raise ValueError(
            "bin_posteriors must have shape (n_bins, n_states)"
        )
    if probs.shape[0] == 0 or probs.shape[1] == 0:
        raise ValueError("bin_posteriors must include at least one bin and state")

    log_joint = np.log(np.clip(probs, np.finfo(np.float64).tiny, 1.0)).sum(
        axis=0,
    )
    log_joint -= np.max(log_joint)
    posterior = np.exp(log_joint)
    posterior /= np.sum(posterior)
    return posterior.astype(np.float32, copy=False)


# ── plotting helpers ────────────────────────────────────────────────────────

def save_and_close_plot(
    output_dir: str,
    filename: str,
    *,
    subdir: str = "diagnostics",
    dpi: int = DEFAULT_RASTER_DPI,
) -> None:
    """Save the current matplotlib figure and close it.

    Args:
        output_dir: Top-level output directory.
        filename: File name (e.g. ``'plot.png'``).
        subdir: Sub-directory under *output_dir* (default ``'diagnostics'``).
        dpi: Resolution for raster formats.
    """
    plt.tight_layout()
    dest = os.path.join(output_dir, subdir)
    os.makedirs(dest, exist_ok=True)
    out = os.path.join(dest, filename)
    actual_out = save_publication_figure(plt.gcf(), out, dpi=dpi)
    plt.close()
    logger.debug("Saved plot: %s", os.path.relpath(actual_out, output_dir))


def add_chromosome_labels(
    ax: plt.Axes,
    chr_array: np.ndarray,
    x_transformed: Optional[np.ndarray] = None,
) -> None:
    """Add chromosome boundary lines and labels to an axis.

    Draws vertical dashed lines at chromosome boundaries and labels each
    chromosome at its midpoint on the x-axis.

    Args:
        ax: Matplotlib axes instance.
        chr_array: 1-D array of chromosome names (one per bin, in order).
        x_transformed: Optional transformed x-coordinates (for equal-width
            chromosome layout).  When *None*, raw bin indices are used.
    """
    chr_changes = np.where(chr_array[:-1] != chr_array[1:])[0] + 1
    boundaries = np.concatenate([[0], chr_changes, [len(chr_array)]])

    labels: list[str] = []
    positions: list[float] = []
    boundary_positions: list[float] = []

    for i in range(len(boundaries) - 1):
        start = boundaries[i]
        end = boundaries[i + 1]

        if x_transformed is not None:
            positions.append(i + 0.5)
            if i > 0:
                boundary_positions.append(float(i))
        else:
            positions.append((start + end) / 2)
            if start > 0:
                boundary_positions.append(float(start))

        labels.append(str(chr_array[start]).replace("chr", ""))

    for b in boundary_positions:
        ax.axvline(
            b,
            color="#424242",
            linestyle="--",
            alpha=0.9,
            linewidth=1.1,
            zorder=8,
            clip_on=False,
        )

    ax.set_xticks(positions)
    ax.set_xticklabels(labels, rotation=0, ha="center")
    ax.set_xlabel("Chromosome")
