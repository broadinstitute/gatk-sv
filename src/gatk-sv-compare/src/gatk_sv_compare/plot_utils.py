"""Reusable plotting utilities for analysis modules."""

from __future__ import annotations

from math import sqrt
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

plt.style.use("seaborn-v0_8-white")

SVTYPE_COLORS = {
    "DEL": "#D43925",
    "DUP": "#2376B2",
    "INS": "#D474E0",
    "INV": "#FA931E",
    "BND": "#397246",
    "CTX": "#397246",
    "CPX": "#71E38C",
    "CNV": "#7459B2",
    "INS:MEI": "#D474E0",
    "OTH": "#397246",
}

CALLSET_COLORS = {"a": "#4C78A8", "b": "#F58518"}
STATE_COLORS = {
    "matched": "#54A24B",
    "unmatched": "#9D9D9D",
    "pass": "#54A24B",
    "nominal": "#EECA3B",
    "bonferroni": "#B279A2",
}
SUMMARY_COLORS = {
    "primary": "#4C78A8",
    "secondary": "#F58518",
    "tertiary": "#72B7B2",
    "quaternary": "#E45756",
    "neutral": "#7A7A7A",
    "edge": "#D9D9D9",
}
MEI_PEAK_COLORS = {"alu": "#F3C969", "sva": "#E8A24D", "l1": "#D66D75"}
OVERLAP_COLORS = {"matched": STATE_COLORS["matched"], "unmatched": STATE_COLORS["unmatched"]}
HWE_COLORS = {"pass": STATE_COLORS["pass"], "nominal": STATE_COLORS["nominal"], "bonferroni": STATE_COLORS["bonferroni"]}
FIGURE_DPI = 450
NATURE_SINGLE_COLUMN_FIGSIZE = (3.5, 2.4)
NATURE_DOUBLE_COLUMN_FIGSIZE = (7.2, 4.5)
DEFAULT_FIGSIZE = NATURE_DOUBLE_COLUMN_FIGSIZE
TITLE_FONTSIZE = 7
LABEL_FONTSIZE = 6
TICK_FONTSIZE = 5.5
LEGEND_FONTSIZE = 5.5

plt.rcParams.update(
    {
        "figure.figsize": DEFAULT_FIGSIZE,
        "figure.dpi": FIGURE_DPI,
        "savefig.dpi": FIGURE_DPI,
        "savefig.facecolor": "white",
        "savefig.edgecolor": "white",
        "savefig.transparent": False,
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": LABEL_FONTSIZE,
        "axes.titlesize": TITLE_FONTSIZE,
        "axes.titleweight": "bold",
        "axes.labelsize": LABEL_FONTSIZE,
        "axes.labelcolor": "black",
        "axes.edgecolor": "black",
        "axes.linewidth": 0.5,
        "xtick.labelsize": TICK_FONTSIZE,
        "ytick.labelsize": TICK_FONTSIZE,
        "xtick.color": "black",
        "ytick.color": "black",
        "xtick.major.width": 0.5,
        "ytick.major.width": 0.5,
        "xtick.major.size": 3.0,
        "ytick.major.size": 3.0,
        "legend.fontsize": LEGEND_FONTSIZE,
        "legend.frameon": False,
        "text.color": "black",
        "lines.linewidth": 0.8,
        "patch.linewidth": 0.5,
        "axes.grid": False,
        "grid.alpha": 0.0,
        "grid.linewidth": 0.0,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "image.cmap": "viridis",
    }
)


def _style_axis(ax) -> None:
    if not ax.axison:
        return
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(0.5)
    ax.spines["bottom"].set_linewidth(0.5)
    ax.set_axisbelow(True)
    ax.grid(False)
    title = ax.get_title()
    if title:
        ax.set_title(title, fontsize=TITLE_FONTSIZE, fontweight="bold")
    xlabel = ax.get_xlabel()
    if xlabel:
        ax.set_xlabel(xlabel, fontsize=LABEL_FONTSIZE)
    ylabel = ax.get_ylabel()
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=LABEL_FONTSIZE)
    ax.tick_params(axis="both", labelsize=TICK_FONTSIZE)
    legend = ax.get_legend()
    if legend is not None:
        legend.set_frame_on(False)
        for text in legend.get_texts():
            text.set_fontsize(LEGEND_FONTSIZE)


def apply_publication_style(fig) -> None:
    for ax in fig.axes:
        _style_axis(ax)


def single_column_figsize(height: float | None = None) -> tuple[float, float]:
    width, default_height = NATURE_SINGLE_COLUMN_FIGSIZE
    return (width, default_height if height is None else height)


def double_column_figsize(height: float | None = None) -> tuple[float, float]:
    width, default_height = NATURE_DOUBLE_COLUMN_FIGSIZE
    return (width, default_height if height is None else height)


def _sliding_window_median_trend(x_values: np.ndarray, y_values: np.ndarray, bandwidth: float = 0.05):
    valid = np.isfinite(x_values) & np.isfinite(y_values)
    x = np.asarray(x_values[valid], dtype=float)
    y = np.asarray(y_values[valid], dtype=float)
    if x.size < 5:
        return None, None
    x_min = float(np.min(x))
    x_max = float(np.max(x))
    if x_min == x_max:
        return None, None
    centers = np.arange(max(0.0, x_min), min(1.0, x_max) + 0.5 * bandwidth, bandwidth / 2.0)
    medians_x = []
    medians = []
    half_bandwidth = bandwidth / 2.0
    for center in centers:
        start = center - half_bandwidth
        end = center + half_bandwidth
        mask = (x >= start) & (x <= end)
        if not np.any(mask):
            continue
        medians_x.append(float(np.median(x[mask])))
        medians.append(float(np.median(y[mask])))
    if len(medians_x) < 2:
        return None, None
    return np.asarray(medians_x, dtype=float), np.asarray(medians, dtype=float)


def plot_stacked_bars(ax, matrix, colors, labels, scaled=True, annotate_counts=True):
    values = np.asarray(matrix, dtype=float)
    if values.ndim != 2:
        raise ValueError("matrix must be 2-dimensional")
    totals = values.sum(axis=1, keepdims=True)
    heights = values / np.where(totals == 0, 1, totals) if scaled else values
    bottoms = np.zeros(values.shape[0], dtype=float)
    x_positions = np.arange(values.shape[0])
    for index in range(values.shape[1]):
        ax.bar(x_positions, heights[:, index], bottom=bottoms, color=colors[index], label=labels[index])
        if annotate_counts:
            for row_idx, value in enumerate(values[:, index]):
                if value:
                    ax.text(x_positions[row_idx], bottoms[row_idx] + heights[row_idx, index] / 2, str(int(value)), ha="center", va="center", fontsize=8)
        bottoms += heights[:, index]
    return ax


def plot_scatter_af(ax, x, y, label_x, label_y, show_lm=True, show_rolling=True):
    x_values = np.asarray(x, dtype=float)
    y_values = np.asarray(y, dtype=float)
    ax.scatter(x_values, y_values, alpha=0.5, s=6, color=SUMMARY_COLORS["neutral"], edgecolors="none")
    min_val = float(np.nanmin([x_values.min(initial=0), y_values.min(initial=0)])) if x_values.size and y_values.size else 0.0
    max_val = float(np.nanmax([x_values.max(initial=1), y_values.max(initial=1)])) if x_values.size and y_values.size else 1.0
    ax.plot([min_val, max_val], [min_val, max_val], color="gray", linestyle="--", linewidth=1)
    if show_rolling and x_values.size > 4:
        trend_x, trend_y = _sliding_window_median_trend(x_values, y_values, bandwidth=0.05)
        if trend_x is not None and trend_y is not None:
            ax.plot(trend_x, trend_y, color="red", linewidth=1.2, alpha=0.9)
    elif show_lm and x_values.size > 1 and y_values.size > 1:
        slope, intercept = np.polyfit(x_values, y_values, 1)
        fit_x = np.linspace(min_val, max_val, num=100)
        ax.plot(fit_x, slope * fit_x + intercept, color="red", linewidth=1.0, alpha=0.8)
    ax.set_xlabel(label_x)
    ax.set_ylabel(label_y)
    return ax


def plot_ternary(ax, aa, ab, bb, colors, draw_hwe_curve=True, alpha=0.05):
    del alpha
    ab_values = np.asarray(ab, dtype=float)
    bb_values = np.asarray(bb, dtype=float)
    x_values = 0.5 * (2 * bb_values + ab_values)
    y_values = (sqrt(3) / 2.0) * ab_values
    ax.scatter(x_values, y_values, c=colors, s=5, alpha=0.6)
    triangle = np.array([[0, 0], [1, 0], [0.5, sqrt(3) / 2.0], [0, 0]])
    ax.plot(triangle[:, 0], triangle[:, 1], color="black", linewidth=1)
    if draw_hwe_curve:
        p = np.linspace(0, 1, 200)
        ab_curve = 2 * p * (1 - p)
        bb_curve = p**2
        curve_x = 0.5 * (2 * bb_curve + ab_curve)
        curve_y = (sqrt(3) / 2.0) * ab_curve
        ax.plot(curve_x, curve_y, color="gray", linestyle="--", linewidth=1)
    ax.text(0.0, -0.045, "REF", ha="left", va="top")
    ax.text(0.5, (sqrt(3) / 2.0) + 0.035, "HET", ha="center", va="bottom")
    ax.text(1.0, -0.045, "HOMVAR", ha="right", va="top")
    ax.set_axis_off()
    return ax


def plot_heatmap_annotated(ax, matrix, row_labels, col_labels, fmt="{value}"):
    values = np.asarray(matrix)
    image = ax.imshow(values, aspect="auto", cmap="viridis")
    ax.set_xticks(np.arange(len(col_labels)))
    ax.set_xticklabels(col_labels)
    ax.set_yticks(np.arange(len(row_labels)))
    ax.set_yticklabels(row_labels)
    for row_idx in range(values.shape[0]):
        for col_idx in range(values.shape[1]):
            value = values[row_idx, col_idx]
            ax.text(
                col_idx,
                row_idx,
                fmt.format(value=value),
                ha="center",
                va="center",
                fontsize=TICK_FONTSIZE,
                color="black",
            )
    return image


def plot_histogram(ax, values, bins, color, label):
    ax.hist(values, bins=bins, color=color, label=label, alpha=0.7, edgecolor=SUMMARY_COLORS["edge"], linewidth=0.6)
    if label:
        ax.legend()
    return ax


def plot_boxplot_grouped(ax, data_dict, colors):
    labels = list(data_dict.keys())
    values = [data_dict[label] for label in labels]
    boxplot = ax.boxplot(values, labels=labels, patch_artist=True)
    for patch, label in zip(boxplot["boxes"], labels):
        patch.set_facecolor(colors.get(label, SUMMARY_COLORS["neutral"]))
    return ax


def plot_beeswarm_horizontal(ax, values_list, colors, labels, show_mean=True):
    for index, values in enumerate(values_list):
        y_coords = np.full(len(values), index, dtype=float)
        ax.scatter(values, y_coords, color=colors[index], alpha=0.6, s=10)
        if show_mean and len(values):
            ax.axvline(float(np.mean(values)), color=colors[index], linewidth=1, linestyle="--")
    ax.set_yticks(np.arange(len(labels)))
    ax.set_yticklabels(labels)
    return ax


def plot_dnr_vs_continuous(ax, bins, dnr_matrix, svtype_colors, k_smooth=4, log_x=True):
    x_values = np.asarray(bins, dtype=float)
    for label, values in dnr_matrix.items():
        y_values = np.asarray(values, dtype=float)
        if k_smooth > 1 and y_values.size >= k_smooth:
            kernel = np.ones(k_smooth) / float(k_smooth)
            y_values = np.convolve(y_values, kernel, mode="same")
        ax.plot(x_values, y_values, label=label, color=svtype_colors.get(label, SUMMARY_COLORS["neutral"]))
    if log_x:
        ax.set_xscale("log")
    ax.legend()
    return ax


def plot_peak_histogram(ax, values, bins, peak_regions=None, log_x=True):
    ax.hist(values, bins=bins, color=SUMMARY_COLORS["primary"], alpha=0.7, edgecolor=SUMMARY_COLORS["edge"], linewidth=0.6)
    if peak_regions:
        for start, end, color in peak_regions:
            ax.axvspan(start, end, color=color, alpha=0.15)
    if log_x:
        ax.set_xscale("log")
    return ax


def save_figure(fig, path: Path, dpi=FIGURE_DPI):
    path.parent.mkdir(parents=True, exist_ok=True)
    apply_publication_style(fig)
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
