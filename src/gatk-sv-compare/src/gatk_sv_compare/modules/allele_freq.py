"""Allele-frequency correlation analysis for matched sites."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

from ..aggregate import AggregatedData
from ..config import AnalysisConfig
from ..plot_utils import plot_scatter_af, save_figure
from .base import AnalysisModule


def build_af_correlation_table(data: AggregatedData) -> pd.DataFrame:
    matched = data.matched_pairs.copy()
    if matched.empty:
        return pd.DataFrame(columns=["group", "n_matched", "pearson_r", "pearson_p", "spearman_rho", "spearman_p", "mean_abs_diff"])
    rows = []
    for group_name, frame in [("overall", matched)] + [(svtype, group) for svtype, group in matched.groupby("svtype_a")]:
        if frame.empty:
            continue
        x_values = frame["af_a"].astype(float)
        y_values = frame["af_b"].astype(float)
        if len(frame) >= 2:
            pearson_r, pearson_p = pearsonr(x_values, y_values)
            spearman_rho, spearman_p = spearmanr(x_values, y_values)
        else:
            pearson_r = pearson_p = spearman_rho = spearman_p = np.nan
        rows.append(
            {
                "group": group_name,
                "n_matched": len(frame),
                "pearson_r": pearson_r,
                "pearson_p": pearson_p,
                "spearman_rho": spearman_rho,
                "spearman_p": spearman_p,
                "mean_abs_diff": float((x_values - y_values).abs().mean()),
            }
        )
    return pd.DataFrame(rows)


class AlleleFreqModule(AnalysisModule):
    @property
    def name(self) -> str:
        return "allele_freq"

    def run(self, data: AggregatedData, config: AnalysisConfig) -> None:
        output_dir = self.output_dir(config)
        tables_dir = output_dir / "tables"
        tables_dir.mkdir(parents=True, exist_ok=True)
        stats = build_af_correlation_table(data)
        stats.to_csv(tables_dir / "af_correlation_stats.tsv", sep="\t", index=False)

        overall = data.matched_pairs.copy()
        fig, ax = plt.subplots(figsize=(6, 6))
        if not overall.empty:
            plot_scatter_af(ax, overall["af_a"], overall["af_b"], data.label_a, data.label_b)
        else:
            ax.text(0.5, 0.5, "No matched sites", ha="center", va="center")
        ax.set_title("AF correlation")
        save_figure(fig, output_dir / "af_correlation.overall.png")

        if not overall.empty:
            svtypes = list(overall["svtype_a"].dropna().unique())
            fig, axes = plt.subplots(1, len(svtypes), figsize=(max(6, 4 * len(svtypes)), 4), squeeze=False)
            for axis, svtype in zip(axes[0], svtypes):
                frame = overall.loc[overall["svtype_a"] == svtype]
                plot_scatter_af(axis, frame["af_a"], frame["af_b"], data.label_a, data.label_b)
                axis.set_title(str(svtype))
            save_figure(fig, output_dir / "af_correlation.by_type.png")
