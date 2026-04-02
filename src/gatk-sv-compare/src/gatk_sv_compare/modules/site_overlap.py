"""Overlap summary plots and tables based on STATUS annotations."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ..aggregate import AggregatedData
from ..config import AnalysisConfig
from ..dimensions import STATUS_MATCHED
from ..plot_utils import OVERLAP_COLORS, plot_heatmap_annotated, save_figure
from .base import AnalysisModule


def _filtered_sites(sites: pd.DataFrame, pass_only: bool) -> pd.DataFrame:
    if not pass_only:
        return sites.copy()
    return sites.loc[sites["in_filtered_pass_view"]].copy()


def _overlap_metrics(sites: pd.DataFrame, suffix: str) -> pd.DataFrame:
    grouped = sites.groupby(["svtype", "size_bucket", "af_bucket", "genomic_context"], dropna=False)
    metrics = grouped.agg(
        **{f"n_total_{suffix}": ("variant_id", "count")},
        **{f"n_matched_{suffix}": ("status", lambda values: int((values == STATUS_MATCHED).sum()))},
    ).reset_index()
    metrics[f"pct_matched_{suffix}"] = np.where(
        metrics[f"n_total_{suffix}"] > 0,
        metrics[f"n_matched_{suffix}"] / metrics[f"n_total_{suffix}"],
        np.nan,
    )
    return metrics


def build_overlap_metrics(sites_a: pd.DataFrame, sites_b: pd.DataFrame, pass_only: bool = False) -> pd.DataFrame:
    metrics_a = _overlap_metrics(_filtered_sites(sites_a, pass_only), "a")
    metrics_b = _overlap_metrics(_filtered_sites(sites_b, pass_only), "b")
    return metrics_a.merge(metrics_b, on=["svtype", "size_bucket", "af_bucket", "genomic_context"], how="outer").fillna(0)


def _plot_overlap_bar(sites: pd.DataFrame, field: str, output_path: Path) -> None:
    grouped = sites.groupby(field, dropna=False).agg(
        matched=("status", lambda values: int((values == STATUS_MATCHED).sum())),
        total=("variant_id", "count"),
    )
    unmatched = grouped["total"] - grouped["matched"]
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.bar(grouped.index.astype(str), grouped["matched"], color=OVERLAP_COLORS["matched"], label="matched")
    ax.bar(grouped.index.astype(str), unmatched, bottom=grouped["matched"], color=OVERLAP_COLORS["unmatched"], label="unmatched")
    ax.set_ylabel("Variant count")
    ax.legend(fontsize=8)
    ax.set_title(f"Overlap by {field}")
    save_figure(fig, output_path)


def _plot_heatmap(sites: pd.DataFrame, row_field: str, col_field: str, output_path: Path) -> None:
    grouped = sites.groupby([row_field, col_field], dropna=False).agg(
        matched=("status", lambda values: int((values == STATUS_MATCHED).sum())),
        total=("variant_id", "count"),
    ).reset_index()
    grouped["pct"] = np.where(grouped["total"] > 0, grouped["matched"] / grouped["total"], 0.0)
    matrix = grouped.pivot(index=row_field, columns=col_field, values="pct").fillna(0.0)
    fig, ax = plt.subplots(figsize=(8, 5))
    image = plot_heatmap_annotated(ax, matrix.values, list(matrix.index), list(matrix.columns), fmt="{value:.2f}")
    fig.colorbar(image, ax=ax)
    ax.set_title(f"Overlap rate: {row_field} × {col_field}")
    save_figure(fig, output_path)


class SiteOverlapModule(AnalysisModule):
    @property
    def name(self) -> str:
        return "site_overlap"

    @property
    def requires_concordance(self) -> bool:
        return True

    def run(self, data: AggregatedData, config: AnalysisConfig) -> None:
        output_dir = self.output_dir(config)
        tables_dir = output_dir / "tables"
        tables_dir.mkdir(parents=True, exist_ok=True)
        metrics = build_overlap_metrics(data.sites_a, data.sites_b, pass_only=config.pass_only)
        metrics.to_csv(tables_dir / "overlap_metrics.tsv", sep="\t", index=False)
        metrics.to_parquet(tables_dir / "overlap_metrics.parquet", index=False)

        for label, sites in ((data.label_a, data.sites_a), (data.label_b, data.sites_b)):
            filtered = _filtered_sites(sites, config.pass_only)
            _plot_overlap_bar(filtered, "svtype", output_dir / f"overlap.by_class.{label}.png")
            _plot_overlap_bar(filtered, "size_bucket", output_dir / f"overlap.by_size.{label}.png")
            _plot_heatmap(filtered, "size_bucket", "svtype", output_dir / f"heatmap.size_x_class.{label}.png")
            _plot_heatmap(filtered, "af_bucket", "svtype", output_dir / f"heatmap.freq_x_class.{label}.png")
