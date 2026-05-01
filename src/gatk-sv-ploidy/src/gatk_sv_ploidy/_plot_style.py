"""Shared Nature-oriented visual style helpers for ploidy plots."""

from __future__ import annotations

from contextlib import contextmanager
import logging
from pathlib import Path
from typing import Iterable

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns


MM_PER_INCH = 25.4
NATURE_SINGLE_COLUMN_MM = 89
NATURE_DOUBLE_COLUMN_MM = 183
NATURE_MAX_HEIGHT_MM = 170
DEFAULT_RASTER_DPI = 450
AXIS_LABEL_SIZE_PT = 8
AXIS_TITLE_SIZE_PT = 8
TICK_LABEL_SIZE_PT = 8
LEGEND_SIZE_PT = 8
PANEL_LABEL_SIZE_PT = 8
ANNOTATION_SIZE_PT = 8

NATURE_PALETTE = [
    "#000000",
    "#E69F00",
    "#56B4E9",
    "#009E73",
    "#F0E442",
    "#0072B2",
    "#D55E00",
    "#CC79A7",
]

BASE_PALETTE = {
    "black": "#000000",
    "orange": "#E69F00",
    "sky_blue": "#56B4E9",
    "teal": "#009E73",
    "yellow": "#F0E442",
    "blue": "#0072B2",
    "vermillion": "#D55E00",
    "purple": "#CC79A7",
    "gray": "#666666",
    "light_gray": "#D9D9D9",
    "ink": "#000000",
}

CHR_TYPE_PALETTE = {
    "Autosomal": BASE_PALETTE["blue"],
    "chrX": BASE_PALETTE["orange"],
    "chrY": BASE_PALETTE["teal"],
}

CN_STATE_PALETTE = {
    0: BASE_PALETTE["black"],
    1: BASE_PALETTE["orange"],
    2: BASE_PALETTE["blue"],
    3: BASE_PALETTE["teal"],
    4: BASE_PALETTE["purple"],
    5: BASE_PALETTE["vermillion"],
}

SEX_PALETTE = {
    "MALE": BASE_PALETTE["blue"],
    "FEMALE": BASE_PALETTE["purple"],
    "UNKNOWN": BASE_PALETTE["gray"],
}

HIGHLIGHT_COLOR = BASE_PALETTE["vermillion"]
REFERENCE_COLOR = BASE_PALETTE["black"]
MUTED_GRID_COLOR = BASE_PALETTE["light_gray"]
_PLOT_OUTPUT_FORMATS = {"png", "pdf"}
_ACTIVE_PLOT_OUTPUT_FORMAT = "png"


def _suppress_chatty_plot_dependency_loggers() -> None:
    """Keep verbose plotting dependencies from inheriting CLI INFO logging.

    Matplotlib can emit very noisy fontTools subsetting logs when figures are
    exported as PDF. The plot CLI intentionally uses INFO-level output for
    user-facing progress, so we raise the fontTools logger threshold here to
    keep those internal messages out of normal runs.
    """
    for logger_name in ("fontTools", "fontTools.subset"):
        dep_logger = logging.getLogger(logger_name)
        if dep_logger.level in (logging.NOTSET, logging.DEBUG, logging.INFO):
            dep_logger.setLevel(logging.WARNING)


def set_plot_output_format(output_format: str) -> None:
    """Set the active figure format for plot-saving helpers."""
    normalized = str(output_format).lower()
    if normalized not in _PLOT_OUTPUT_FORMATS:
        allowed = ", ".join(sorted(_PLOT_OUTPUT_FORMATS))
        raise ValueError(
            f"Unsupported plot output format: {output_format!r}. "
            f"Choose one of: {allowed}."
        )

    global _ACTIVE_PLOT_OUTPUT_FORMAT
    _ACTIVE_PLOT_OUTPUT_FORMAT = normalized


def get_plot_output_format() -> str:
    """Return the active figure format for plot-saving helpers."""
    return _ACTIVE_PLOT_OUTPUT_FORMAT


@contextmanager
def plot_output_format(output_format: str):
    """Temporarily set the active plot output format."""
    previous = get_plot_output_format()
    set_plot_output_format(output_format)
    try:
        yield
    finally:
        set_plot_output_format(previous)


def resolve_plot_output_path(output_path: str | Path) -> Path:
    """Resolve a requested figure path to the active output format."""
    path = Path(output_path)
    suffix = f".{get_plot_output_format()}"
    if path.suffix.lower() != suffix:
        path = path.with_suffix(suffix)
    return path


def figure_size(width_mm: float, height_mm: float) -> tuple[float, float]:
    """Return a matplotlib ``figsize`` tuple from physical millimetres."""
    return width_mm / MM_PER_INCH, height_mm / MM_PER_INCH


def bounded_figure_size(width_mm: float, height_mm: float) -> tuple[float, float]:
    """Return a figure size constrained to Nature one/two-column bounds."""
    bounded_width_mm = min(float(width_mm), float(NATURE_DOUBLE_COLUMN_MM))
    bounded_height_mm = min(float(height_mm), float(NATURE_MAX_HEIGHT_MM))
    return figure_size(bounded_width_mm, bounded_height_mm)


def single_column_size(height_mm: float) -> tuple[float, float]:
    """Return a Nature single-column figure size with bounded height."""
    return bounded_figure_size(NATURE_SINGLE_COLUMN_MM, height_mm)


def double_column_size(height_mm: float) -> tuple[float, float]:
    """Return a Nature double-column figure size with bounded height."""
    return bounded_figure_size(NATURE_DOUBLE_COLUMN_MM, height_mm)


def apply_theme() -> None:
    """Apply Nature-style matplotlib/seaborn defaults globally."""
    _suppress_chatty_plot_dependency_loggers()
    sns.set_theme(
        context="paper",
        style="ticks",
        palette=NATURE_PALETTE,
        rc={
            "axes.edgecolor": REFERENCE_COLOR,
            "axes.labelcolor": REFERENCE_COLOR,
            "axes.linewidth": 0.5,
            "axes.titlesize": AXIS_TITLE_SIZE_PT,
            "axes.titleweight": "normal",
            "axes.labelsize": AXIS_LABEL_SIZE_PT,
            "figure.facecolor": "white",
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": AXIS_LABEL_SIZE_PT,
            "grid.color": MUTED_GRID_COLOR,
            "grid.linewidth": 0.35,
            "legend.edgecolor": REFERENCE_COLOR,
            "legend.fontsize": LEGEND_SIZE_PT,
            "legend.frameon": False,
            "savefig.facecolor": "white",
            "xtick.color": REFERENCE_COLOR,
            "xtick.labelsize": TICK_LABEL_SIZE_PT,
            "xtick.major.width": 0.5,
            "ytick.color": REFERENCE_COLOR,
            "ytick.labelsize": TICK_LABEL_SIZE_PT,
            "ytick.major.width": 0.5,
        },
    )
    mpl.rcParams.update(
        {
            "axes.grid": False,
            "axes.axisbelow": True,
            "figure.dpi": 150,
            "lines.linewidth": 1.0,
            "patch.linewidth": 0.5,
            "savefig.dpi": DEFAULT_RASTER_DPI,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def save_publication_figure(
    fig: plt.Figure,
    output_path: str | Path,
    *,
    dpi: int = DEFAULT_RASTER_DPI,
) -> Path:
    """Save a figure in the currently selected plot output format."""
    _suppress_chatty_plot_dependency_loggers()
    path = resolve_plot_output_path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    return path


def polish_axes(axes: plt.Axes | Iterable[plt.Axes]) -> None:
    """Apply small, consistent finishing touches to one or more axes."""
    if isinstance(axes, plt.Axes):
        axes_iter = [axes]
    else:
        axes_iter = axes
    for ax in axes_iter:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(axis="both", width=0.5, length=2.5)


def compact_int(value: int | float) -> str:
    """Format large integer-like counts for report labels."""
    value = float(value)
    if abs(value) >= 1_000_000:
        return f"{value / 1_000_000:.1f}M"
    if abs(value) >= 1_000:
        return f"{value / 1_000:.1f}K"
    return f"{value:.0f}"
