"""
Shared utilities for gatk-sv-ploidy.

Small helper functions used by multiple subcommands.
"""

from __future__ import annotations

import logging
import os
from typing import List, Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

logger = logging.getLogger(__name__)

# ── column / chromosome constants ───────────────────────────────────────────

METADATA_COLS = frozenset({"Chr", "Start", "End", "source_file", "Bin"})
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
MAX_CNQ = 99
"""Maximum phred-scaled CN quality score reported by the ploidy model."""


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


# ── plotting helpers ────────────────────────────────────────────────────────

def save_and_close_plot(
    output_dir: str,
    filename: str,
    *,
    subdir: str = "diagnostics",
    dpi: int = 300,
) -> None:
    """Save the current matplotlib figure and close it.

    Args:
        output_dir: Top-level output directory.
        filename: File name (e.g. ``'training_loss.png'``).
        subdir: Sub-directory under *output_dir* (default ``'diagnostics'``).
        dpi: Resolution for raster formats.
    """
    plt.tight_layout()
    dest = os.path.join(output_dir, subdir)
    os.makedirs(dest, exist_ok=True)
    out = os.path.join(dest, filename)
    plt.savefig(out, dpi=dpi, bbox_inches="tight")
    plt.close()
    logger.info("Saved plot: %s/%s", subdir, filename)


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
        ax.axvline(b, color="gray", linestyle="--", alpha=1, linewidth=1)

    ax.set_xticks(positions)
    ax.set_xticklabels(labels, rotation=0, ha="center")
    ax.set_xlabel("Chromosome")
