"""Overall count plots for single-dimension site distributions."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ..aggregate import AggregatedData
from ..config import AnalysisConfig
from ..dimensions import ordered_contexts, ordered_svtypes
from ..plot_utils import MEI_PEAK_COLORS, SUMMARY_COLORS, SVTYPE_COLORS, double_column_figsize, save_figure
from .base import AnalysisModule

_SIZE_BIN_EDGES = np.logspace(np.log10(50.0), 7, 101)
_MEI_PEAK_REGIONS = [(200, 400, MEI_PEAK_COLORS["alu"]), (1500, 3000, MEI_PEAK_COLORS["sva"]), (5000, 7000, MEI_PEAK_COLORS["l1"])]
_MEI_PEAK_LABELS = [("Alu", 200, 400), ("SVA", 1500, 3000), ("LINE1", 5000, 7000)]


def _filtered_sites(sites: pd.DataFrame, pass_only: bool) -> pd.DataFrame:
    if not pass_only:
        return sites.copy()
    return sites.loc[sites["in_filtered_pass_view"]].copy()


def _af_bin_edges(sites: pd.DataFrame, low_ac_max: int = 20, high_bin_count: int = 80) -> np.ndarray:
    valid_sites = sites.loc[sites["af"].notna()].copy()
    positive_af = valid_sites.loc[valid_sites["af"].astype(float) > 0, "af"].astype(float)
    if positive_af.empty:
        return np.logspace(-4, 0, 101)

    min_af = float(positive_af.min())
    if min_af >= 1.0:
        min_af = max(1.0e-4, min_af / 10.0)

    low_site_mask = (
        valid_sites["ac"].notna()
        & valid_sites["af"].notna()
        & (valid_sites["ac"].astype(float) > 0)
        & (valid_sites["ac"].astype(float) <= float(low_ac_max))
        & (valid_sites["af"].astype(float) > 0)
    )
    low_sites = valid_sites.loc[low_site_mask].copy()
    low_values = np.unique(low_sites["af"].astype(float).to_numpy())
    low_values = low_values[np.isfinite(low_values)]
    low_values.sort()

    if low_values.size == 0:
        return np.logspace(np.log10(min_af), 0, 101)

    if low_values.size == 1:
        first_edge = max(min_af * 0.8, low_values[0] * 0.8, np.nextafter(0.0, 1.0))
        low_edges = np.asarray([first_edge, low_values[0] * 1.2], dtype=float)
    else:
        midpoints = np.sqrt(low_values[:-1] * low_values[1:])
        first_edge = max(low_values[0] / np.sqrt(low_values[1] / low_values[0]), np.nextafter(0.0, 1.0))
        low_edges = np.concatenate(([first_edge], midpoints))

    low_max = float(low_values[-1])
    high_values = positive_af.loc[positive_af > low_max].to_numpy(dtype=float)
    if high_values.size:
        transition_edge = float(np.sqrt(low_max * np.min(high_values)))
    else:
        transition_edge = min(1.0, low_max * 1.25)
    transition_edge = max(transition_edge, low_max * 1.0001)
    low_edges = np.concatenate((low_edges, [transition_edge]))

    if transition_edge >= 1.0:
        edges = np.concatenate((low_edges[:-1], [1.0]))
        return np.unique(edges)

    high_edges = np.logspace(np.log10(transition_edge), 0, high_bin_count + 1)
    edges = np.concatenate((low_edges[:-1], high_edges))
    return np.unique(edges)


def _plot_type_counts(sites: pd.DataFrame, output_path: Path, title: str) -> None:
    counts = sites["svtype"].value_counts()
    counts = counts.reindex(ordered_svtypes(counts.index), fill_value=0)
    fig, ax = plt.subplots(figsize=double_column_figsize(3.0))
    ax.bar(counts.index, counts.values, color=[SVTYPE_COLORS.get(label, SUMMARY_COLORS["neutral"]) for label in counts.index], edgecolor=SUMMARY_COLORS["edge"], linewidth=0.8)
    ax.set_ylabel("Variant count")
    ax.set_xlabel("SV type")
    ax.set_title(title)
    save_figure(fig, output_path)


def _plot_size_distribution(sites: pd.DataFrame, output_path: Path, title: str) -> None:
    fig, ax = plt.subplots(figsize=double_column_figsize(3.0))
    for start, end, color in _MEI_PEAK_REGIONS:
        ax.axvspan(start, end, color=color, alpha=0.12, zorder=0)
    valid_sites = sites.loc[sites["svlen"].notna() & (sites["svlen"] > 0)]
    for svtype in ordered_svtypes(valid_sites["svtype"].unique()):
        frame = valid_sites.loc[valid_sites["svtype"] == svtype]
        values = np.clip(frame["svlen"].astype(float).to_numpy(), 50.0, 1.0e7)
        counts, edges = np.histogram(values, bins=_SIZE_BIN_EDGES)
        centers = np.sqrt(edges[:-1] * edges[1:])
        ax.plot(centers, counts, linewidth=1.5, label=svtype, color=SVTYPE_COLORS.get(svtype, SUMMARY_COLORS["neutral"]))
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(50.0, 1.0e7)
    ax.set_xlabel("SVLEN")
    ax.set_ylabel("Count")
    ax.legend(fontsize=8)
    ax.set_title(title)
    y_max = ax.get_ylim()[1]
    for label, start, end in _MEI_PEAK_LABELS:
        ax.text(np.sqrt(start * end), y_max * 0.92, label, ha="center", va="top", fontsize=10)
    save_figure(fig, output_path)


def _plot_af_distribution(sites: pd.DataFrame, output_path: Path, title: str) -> None:
    fig, ax = plt.subplots(figsize=double_column_figsize(3.0))
    valid_sites = sites.loc[sites["af"].notna()]
    af_bin_edges = _af_bin_edges(valid_sites)
    min_af = float(af_bin_edges[0])
    for svtype in ordered_svtypes(valid_sites["svtype"].unique()):
        frame = valid_sites.loc[valid_sites["svtype"] == svtype]
        values = np.clip(frame["af"].astype(float).to_numpy(), min_af, 1.0)
        counts, edges = _normalized_histogram(values, af_bin_edges)
        centers = np.sqrt(edges[:-1] * edges[1:])
        ax.plot(centers, counts, linewidth=1.5, label=svtype, color=SVTYPE_COLORS.get(svtype, SUMMARY_COLORS["neutral"]))
    ax.set_xscale("log")
    ax.set_xlim(min_af, 1.0)
    ax.set_xlabel("Allele frequency")
    ax.set_ylabel("Proportion of SV type")
    ax.legend(fontsize=8)
    ax.set_title(title)
    save_figure(fig, output_path)


def _normalized_histogram(values: np.ndarray, bins: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    counts, edges = np.histogram(values, bins=bins)
    total = int(counts.sum())
    if total == 0:
        return counts.astype(float), edges
    return counts.astype(float) / float(total), edges


def _plot_context_by_type(sites: pd.DataFrame, output_path: Path, title: str) -> None:
    context_table = pd.crosstab(sites["svtype"], sites["genomic_context"], normalize="index")
    context_table = context_table.reindex(index=ordered_svtypes(context_table.index), columns=ordered_contexts(context_table.columns), fill_value=0.0)
    fig, ax = plt.subplots(figsize=double_column_figsize(3.0))
    bottom = None
    for context in context_table.columns:
        values = context_table[context].values
        ax.bar(context_table.index, values, bottom=bottom, label=context)
        for patch in ax.patches[-len(values):]:
            patch.set_edgecolor(SUMMARY_COLORS["edge"])
            patch.set_linewidth(0.8)
        bottom = values if bottom is None else bottom + values
    ax.set_ylabel("Fraction")
    ax.set_xlabel("SV type")
    ax.set_title(title)
    ax.legend(fontsize=8, ncol=max(1, len(context_table.columns)), loc="lower center", bbox_to_anchor=(0.5, 1.02), frameon=False)
    save_figure(fig, output_path)


class OverallCountsModule(AnalysisModule):
    @property
    def name(self) -> str:
        return "overall_counts"

    def run(self, data: AggregatedData, config: AnalysisConfig) -> None:
        output_dir = self.output_dir(config)
        for label, sites in ((data.label_a, data.sites_a), (data.label_b, data.sites_b)):
            filtered = _filtered_sites(sites, config.pass_only)
            _plot_type_counts(filtered, output_dir / f"sv_count_by_type.{label}.png", label)
            _plot_size_distribution(filtered, output_dir / f"size_distribution.{label}.png", label)
            _plot_af_distribution(filtered, output_dir / f"af_distribution.{label}.png", label)
            _plot_context_by_type(filtered, output_dir / f"genomic_context_by_type.{label}.png", label)
