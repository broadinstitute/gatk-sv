"""Overall count plots for single-dimension site distributions."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ..aggregate import AggregatedData
from ..config import AnalysisConfig
from ..dimensions import explode_algorithm_buckets, ordered_algorithms, ordered_contexts, ordered_evidence_buckets, ordered_svtypes
from ..plot_utils import MEI_PEAK_COLORS, SUMMARY_COLORS, SVTYPE_COLORS, double_column_figsize, save_figure
from .base import AnalysisModule

_SIZE_BIN_EDGES = np.logspace(np.log10(50.0), 7, 101)
_MEI_PEAK_REGIONS = [(200, 400, MEI_PEAK_COLORS["alu"]), (1500, 3000, MEI_PEAK_COLORS["sva"]), (5000, 7000, MEI_PEAK_COLORS["l1"])]
_MEI_PEAK_LABELS = [("Alu", 200, 400), ("SVA", 1500, 3000), ("LINE1", 5000, 7000)]


def _ordered_categories(values, field: str) -> list[str]:
    if field == "algorithm":
        return ordered_algorithms(values)
    if field == "evidence_bucket":
        return ordered_evidence_buckets(values)
    if field == "genomic_context":
        return ordered_contexts(values)
    return ordered_svtypes(values)


def _category_colors(labels: list[str], field: str) -> dict[str, tuple[float, float, float, float] | str]:
    if field == "svtype":
        return {label: SVTYPE_COLORS.get(label, SUMMARY_COLORS["neutral"]) for label in labels}
    cmap = plt.get_cmap("tab20")
    return {label: cmap(index % cmap.N) for index, label in enumerate(labels)}


def _filtered_sites(sites: pd.DataFrame, pass_only: bool) -> pd.DataFrame:
    if not pass_only:
        return sites.copy()
    return sites.loc[sites["in_filtered_pass_view"]].copy()


def _infer_reference_an(sites: pd.DataFrame, low_ac_max: int) -> int | None:
    if "an" in sites.columns:
        an_values = pd.to_numeric(sites["an"], errors="coerce")
        an_values = an_values.loc[np.isfinite(an_values) & (an_values > 0)]
        if not an_values.empty:
            rounded = np.rint(an_values.to_numpy(dtype=float)).astype(int)
            values, counts = np.unique(rounded, return_counts=True)
            return int(values[np.argmax(counts)])

    if "ac" in sites.columns:
        ac_values = pd.to_numeric(sites["ac"], errors="coerce")
        af_values = pd.to_numeric(sites["af"], errors="coerce")
        ratio_mask = np.isfinite(ac_values) & np.isfinite(af_values)
        ratio_mask = ratio_mask & (ac_values > 0) & (af_values > 0)
        ratio_mask = ratio_mask & (ac_values <= float(low_ac_max))
        if ratio_mask.any():
            ratios = np.rint(ac_values.loc[ratio_mask].to_numpy(dtype=float) / af_values.loc[ratio_mask].to_numpy(dtype=float)).astype(int)
            ratios = ratios[ratios > 0]
            if ratios.size:
                values, counts = np.unique(ratios, return_counts=True)
                return int(values[np.argmax(counts)])

    return None


def _af_bin_layout(sites: pd.DataFrame, low_ac_max: int = 20, high_bin_count: int = 80) -> tuple[np.ndarray, np.ndarray]:
    valid_sites = sites.loc[sites["af"].notna()].copy()
    positive_af = valid_sites.loc[valid_sites["af"].astype(float) > 0, "af"].astype(float)
    if positive_af.empty:
        edges = np.logspace(-4, 0, 101)
        centers = np.sqrt(edges[:-1] * edges[1:])
        return edges, centers

    reference_an = _infer_reference_an(valid_sites, low_ac_max)
    if reference_an is None:
        min_af = float(positive_af.min())
        if min_af >= 1.0:
            min_af = max(1.0e-4, min_af / 10.0)
        edges = np.logspace(np.log10(min_af), 0, 101)
        centers = np.sqrt(edges[:-1] * edges[1:])
        return edges, centers

    low_ac_limit = min(low_ac_max, reference_an)
    if low_ac_limit < 1:
        edges = np.logspace(-4, 0, 101)
        centers = np.sqrt(edges[:-1] * edges[1:])
        return edges, centers

    low_centers = np.arange(1, low_ac_limit + 1, dtype=float) / float(reference_an)
    low_edges = (np.arange(0, low_ac_limit + 1, dtype=float) + 0.5) / float(reference_an)
    low_edges[0] = max(low_edges[0], np.nextafter(0.0, 1.0))
    transition_edge = float(min(1.0, (low_ac_limit + 0.5) / float(reference_an)))

    if transition_edge >= 1.0 or high_bin_count <= 0:
        edges = np.concatenate((low_edges[:-1], [1.0]))
        centers = low_centers
        return np.unique(edges), centers

    high_edges = np.linspace(transition_edge, 1.0, high_bin_count + 1, dtype=float)
    high_centers = 0.5 * (high_edges[:-1] + high_edges[1:])
    edges = np.concatenate((low_edges[:-1], high_edges))
    centers = np.concatenate((low_centers, high_centers))
    return edges, centers


def _af_bin_edges(sites: pd.DataFrame, low_ac_max: int = 20, high_bin_count: int = 80) -> np.ndarray:
    return _af_bin_layout(sites, low_ac_max=low_ac_max, high_bin_count=high_bin_count)[0]


def _plot_type_counts(sites: pd.DataFrame, output_path: Path, title: str, field: str = "svtype", xlabel: str = "SV type") -> None:
    if field == "algorithm":
        sites = explode_algorithm_buckets(sites)
    counts = sites[field].value_counts()
    counts = counts.reindex(_ordered_categories(counts.index, field), fill_value=0)
    colors = _category_colors(list(counts.index.astype(str)), field)
    fig, ax = plt.subplots(figsize=double_column_figsize(3.0))
    ax.bar(counts.index, counts.values, color=[colors[str(label)] for label in counts.index], edgecolor=SUMMARY_COLORS["edge"], linewidth=0.8)
    ax.set_ylabel("Variant count")
    ax.set_xlabel(xlabel)
    ax.set_title(title)
    save_figure(fig, output_path)


def _plot_size_distribution(sites: pd.DataFrame, output_path: Path, title: str, group_field: str = "svtype") -> None:
    fig, ax = plt.subplots(figsize=double_column_figsize(3.0))
    for start, end, color in _MEI_PEAK_REGIONS:
        ax.axvspan(start, end, color=color, alpha=0.12, zorder=0)
    if group_field == "algorithm":
        sites = explode_algorithm_buckets(sites)
    valid_sites = sites.loc[sites["svlen"].notna() & (sites["svlen"] > 0)]
    group_values = _ordered_categories(valid_sites[group_field].unique(), group_field)
    colors = _category_colors(group_values, group_field)
    for group_value in group_values:
        frame = valid_sites.loc[valid_sites[group_field] == group_value]
        if frame.empty:
            continue
        values = np.clip(frame["svlen"].astype(float).to_numpy(), 50.0, 1.0e7)
        counts, edges = np.histogram(values, bins=_SIZE_BIN_EDGES)
        centers = np.sqrt(edges[:-1] * edges[1:])
        ax.plot(centers, counts, linewidth=1.5, label=str(group_value), color=colors[str(group_value)])
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


def _plot_af_distribution(sites: pd.DataFrame, output_path: Path, title: str, group_field: str = "svtype") -> None:
    fig, ax = plt.subplots(figsize=double_column_figsize(3.0))
    if group_field == "algorithm":
        sites = explode_algorithm_buckets(sites)
    valid_sites = sites.loc[sites["af"].notna()]
    af_bin_edges, af_bin_centers = _af_bin_layout(valid_sites)
    min_af = float(af_bin_edges[0])
    group_values = _ordered_categories(valid_sites[group_field].unique(), group_field)
    colors = _category_colors(group_values, group_field)
    for group_value in group_values:
        frame = valid_sites.loc[valid_sites[group_field] == group_value]
        if frame.empty:
            continue
        values = np.clip(frame["af"].astype(float).to_numpy(), min_af, 1.0)
        counts, edges = _normalized_histogram(values, af_bin_edges)
        centers = af_bin_centers if counts.shape == af_bin_centers.shape else np.sqrt(edges[:-1] * edges[1:])
        ax.plot(centers, counts, linewidth=1.5, label=str(group_value), color=colors[str(group_value)])
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


def _plot_context_by_type(sites: pd.DataFrame, output_path: Path, title: str, group_field: str = "svtype", xlabel: str = "SV type") -> None:
    if group_field == "algorithm":
        sites = explode_algorithm_buckets(sites)
    context_table = pd.crosstab(sites[group_field], sites["genomic_context"], normalize="index")
    context_table = context_table.reindex(index=_ordered_categories(context_table.index, group_field), columns=ordered_contexts(context_table.columns), fill_value=0.0)
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
    ax.set_xlabel(xlabel)
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
            _plot_type_counts(filtered, output_dir / f"sv_count_by_algorithm.{label}.png", label, field="algorithm", xlabel="Algorithm")
            _plot_type_counts(filtered, output_dir / f"sv_count_by_evidence.{label}.png", label, field="evidence_bucket", xlabel="Evidence")
            _plot_type_counts(filtered, output_dir / f"sv_count_by_context.{label}.png", label, field="genomic_context", xlabel="Genomic context")
            _plot_size_distribution(filtered, output_dir / f"size_distribution.{label}.png", label)
            _plot_size_distribution(filtered, output_dir / f"size_distribution.by_algorithm.{label}.png", label, group_field="algorithm")
            _plot_size_distribution(filtered, output_dir / f"size_distribution.by_evidence.{label}.png", label, group_field="evidence_bucket")
            _plot_size_distribution(filtered, output_dir / f"size_distribution.by_context.{label}.png", label, group_field="genomic_context")
            _plot_af_distribution(filtered, output_dir / f"af_distribution.{label}.png", label)
            _plot_af_distribution(filtered, output_dir / f"af_distribution.by_algorithm.{label}.png", label, group_field="algorithm")
            _plot_af_distribution(filtered, output_dir / f"af_distribution.by_evidence.{label}.png", label, group_field="evidence_bucket")
            _plot_af_distribution(filtered, output_dir / f"af_distribution.by_context.{label}.png", label, group_field="genomic_context")
            _plot_context_by_type(filtered, output_dir / f"genomic_context_by_type.{label}.png", label)
            _plot_context_by_type(filtered, output_dir / f"genomic_context_by_algorithm.{label}.png", label, group_field="algorithm", xlabel="Algorithm")
            _plot_context_by_type(filtered, output_dir / f"genomic_context_by_evidence.{label}.png", label, group_field="evidence_bucket", xlabel="Evidence")
