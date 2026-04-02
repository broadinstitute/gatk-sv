"""Overall count plots for single-dimension site distributions."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ..aggregate import AggregatedData
from ..config import AnalysisConfig
from ..plot_utils import SVTYPE_COLORS, save_figure
from .base import AnalysisModule

_SIZE_BIN_EDGES = np.logspace(np.log10(50.0), 7, 101)
_AF_BIN_EDGES = np.linspace(0.0, 1.0, 101)


def _filtered_sites(sites: pd.DataFrame, pass_only: bool) -> pd.DataFrame:
    if not pass_only:
        return sites.copy()
    return sites.loc[sites["in_filtered_pass_view"]].copy()


def _plot_type_counts(sites: pd.DataFrame, output_path: Path) -> None:
    counts = sites["svtype"].value_counts().sort_index()
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.bar(counts.index, counts.values, color=[SVTYPE_COLORS.get(label, "#999999") for label in counts.index], edgecolor="#E5E5E5", linewidth=0.8)
    ax.set_ylabel("Variant count")
    ax.set_xlabel("SV type")
    ax.set_title("SV count by type")
    save_figure(fig, output_path)


def _plot_size_distribution(sites: pd.DataFrame, output_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 4))
    for svtype, frame in sites.loc[sites["svlen"].notna() & (sites["svlen"] > 0)].groupby("svtype"):
        values = np.clip(frame["svlen"].astype(float).to_numpy(), 50.0, 1.0e7)
        counts, edges = np.histogram(values, bins=_SIZE_BIN_EDGES)
        centers = np.sqrt(edges[:-1] * edges[1:])
        ax.plot(centers, counts, linewidth=1.5, label=svtype, color=SVTYPE_COLORS.get(svtype, "#999999"))
    ax.set_xscale("log")
    ax.set_xlim(50.0, 1.0e7)
    ax.set_xlabel("SVLEN")
    ax.set_ylabel("Count")
    ax.legend(fontsize=8)
    ax.set_title("Size distribution by type")
    save_figure(fig, output_path)


def _plot_af_distribution(sites: pd.DataFrame, output_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 4))
    for svtype, frame in sites.loc[sites["af"].notna()].groupby("svtype"):
        values = np.clip(frame["af"].astype(float).to_numpy(), 0.0, 1.0)
        counts, edges = np.histogram(values, bins=_AF_BIN_EDGES)
        centers = (edges[:-1] + edges[1:]) / 2.0
        ax.plot(centers, counts, linewidth=1.5, label=svtype, color=SVTYPE_COLORS.get(svtype, "#999999"))
    ax.set_xlim(0.0, 1.0)
    ax.set_xlabel("Allele frequency")
    ax.set_ylabel("Count")
    ax.legend(fontsize=8)
    ax.set_title("AF distribution by type")
    save_figure(fig, output_path)


def _plot_context_by_type(sites: pd.DataFrame, output_path: Path) -> None:
    context_table = pd.crosstab(sites["svtype"], sites["genomic_context"], normalize="index")
    fig, ax = plt.subplots(figsize=(8, 4))
    bottom = None
    for context in context_table.columns:
        values = context_table[context].values
        ax.bar(context_table.index, values, bottom=bottom, label=context)
        for patch in ax.patches[-len(values):]:
            patch.set_edgecolor("#E5E5E5")
            patch.set_linewidth(0.8)
        bottom = values if bottom is None else bottom + values
    ax.set_ylabel("Fraction")
    ax.set_xlabel("SV type")
    ax.legend(fontsize=8)
    ax.set_title("Genomic context by type")
    save_figure(fig, output_path)


class OverallCountsModule(AnalysisModule):
    @property
    def name(self) -> str:
        return "overall_counts"

    def run(self, data: AggregatedData, config: AnalysisConfig) -> None:
        output_dir = self.output_dir(config)
        for label, sites in ((data.label_a, data.sites_a), (data.label_b, data.sites_b)):
            filtered = _filtered_sites(sites, config.pass_only)
            _plot_type_counts(filtered, output_dir / f"sv_count_by_type.{label}.png")
            _plot_size_distribution(filtered, output_dir / f"size_distribution.{label}.png")
            _plot_af_distribution(filtered, output_dir / f"af_distribution.{label}.png")
            _plot_context_by_type(filtered, output_dir / f"genomic_context_by_type.{label}.png")
