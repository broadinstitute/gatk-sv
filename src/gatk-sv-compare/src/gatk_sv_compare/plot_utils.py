"""Reusable plotting utilities for analysis modules."""

from __future__ import annotations

from math import sqrt
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

SVTYPE_COLORS = {
    "DEL": "#D43925",
    "DUP": "#2376B2",
    "INS": "#7B2D8E",
    "INV": "#F57E20",
    "BND": "#4DAF4A",
    "CTX": "#4DAF4A",
    "CPX": "#E31A8B",
    "CNV": "#984EA3",
    "INS:MEI": "#AF5FA0",
    "OTH": "#999999",
}

OVERLAP_COLORS = {"matched": "#4DAC26", "unmatched": "#696969"}
HWE_COLORS = {"pass": "#4DAC26", "nominal": "#81F850", "bonferroni": "#AC26A1"}
FIGURE_DPI = 300
DEFAULT_FIGSIZE = (12, 8)


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
    ax.scatter(x_values, y_values, alpha=0.5, s=12)
    min_val = float(np.nanmin([x_values.min(initial=0), y_values.min(initial=0)])) if x_values.size and y_values.size else 0.0
    max_val = float(np.nanmax([x_values.max(initial=1), y_values.max(initial=1)])) if x_values.size and y_values.size else 1.0
    ax.plot([min_val, max_val], [min_val, max_val], color="gray", linestyle="--", linewidth=1)
    if show_lm and x_values.size > 1 and y_values.size > 1:
        slope, intercept = np.polyfit(x_values, y_values, 1)
        fit_x = np.linspace(min_val, max_val, num=100)
        ax.plot(fit_x, slope * fit_x + intercept, color="red", linewidth=1.5)
    if show_rolling and x_values.size > 4:
        order = np.argsort(x_values)
        sorted_x = x_values[order]
        sorted_y = y_values[order]
        kernel = np.ones(5) / 5.0
        smooth_y = np.convolve(sorted_y, kernel, mode="same")
        ax.plot(sorted_x, smooth_y, color="darkred", linewidth=1.2)
    ax.set_xlabel(label_x)
    ax.set_ylabel(label_y)
    return ax


def plot_ternary(ax, aa, ab, bb, colors, draw_hwe_curve=True, alpha=0.05):
    del alpha
    ab_values = np.asarray(ab, dtype=float)
    bb_values = np.asarray(bb, dtype=float)
    x_values = 0.5 * (2 * bb_values + ab_values)
    y_values = (sqrt(3) / 2.0) * ab_values
    ax.scatter(x_values, y_values, c=colors, s=10, alpha=0.6)
    triangle = np.array([[0, 0], [1, 0], [0.5, sqrt(3) / 2.0], [0, 0]])
    ax.plot(triangle[:, 0], triangle[:, 1], color="black", linewidth=1)
    if draw_hwe_curve:
        p = np.linspace(0, 1, 200)
        ab_curve = 2 * p * (1 - p)
        bb_curve = p**2
        curve_x = 0.5 * (2 * bb_curve + ab_curve)
        curve_y = (sqrt(3) / 2.0) * ab_curve
        ax.plot(curve_x, curve_y, color="gray", linestyle="--", linewidth=1)
    ax.set_axis_off()
    return ax


def plot_heatmap_annotated(ax, matrix, row_labels, col_labels, fmt="{value}"):
    values = np.asarray(matrix)
    image = ax.imshow(values, aspect="auto")
    ax.set_xticks(np.arange(len(col_labels)))
    ax.set_xticklabels(col_labels)
    ax.set_yticks(np.arange(len(row_labels)))
    ax.set_yticklabels(row_labels)
    for row_idx in range(values.shape[0]):
        for col_idx in range(values.shape[1]):
            ax.text(col_idx, row_idx, fmt.format(value=values[row_idx, col_idx]), ha="center", va="center", fontsize=8)
    return image


def plot_histogram(ax, values, bins, color, label):
    ax.hist(values, bins=bins, color=color, label=label, alpha=0.7)
    if label:
        ax.legend()
    return ax


def plot_boxplot_grouped(ax, data_dict, colors):
    labels = list(data_dict.keys())
    values = [data_dict[label] for label in labels]
    boxplot = ax.boxplot(values, labels=labels, patch_artist=True)
    for patch, label in zip(boxplot["boxes"], labels):
        patch.set_facecolor(colors.get(label, "#999999"))
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
        ax.plot(x_values, y_values, label=label, color=svtype_colors.get(label, "#999999"))
    if log_x:
        ax.set_xscale("log")
    ax.legend()
    return ax


def plot_peak_histogram(ax, values, bins, peak_regions=None, log_x=True):
    ax.hist(values, bins=bins, color="#4C72B0", alpha=0.7)
    if peak_regions:
        for start, end, color in peak_regions:
            ax.axvspan(start, end, color=color, alpha=0.15)
    if log_x:
        ax.set_xscale("log")
    return ax


def save_figure(fig, path: Path, dpi=FIGURE_DPI):
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
