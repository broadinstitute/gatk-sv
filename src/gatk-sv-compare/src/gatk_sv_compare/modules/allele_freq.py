"""Allele-frequency correlation analysis for matched sites."""

from __future__ import annotations

import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

from ..aggregate import AggregatedData
from ..config import AnalysisConfig
from ..dimensions import ordered_contexts, ordered_plot_af_buckets, ordered_plot_size_buckets, ordered_svtypes
from ..plot_utils import double_column_figsize, plot_scatter_af, save_figure, single_column_figsize
from .base import AnalysisModule, write_tsv_gz


def _compute_correlation_stats(x_values: pd.Series, y_values: pd.Series) -> tuple[float, float, float, float]:
    valid = pd.DataFrame({"x": x_values, "y": y_values}).dropna()
    if len(valid) < 2:
        return np.nan, np.nan, np.nan, np.nan
    x = valid["x"].astype(float).to_numpy()
    y = valid["y"].astype(float).to_numpy()
    if np.nanstd(x) == 0.0 or np.nanstd(y) == 0.0:
        if np.allclose(x, y):
            return 1.0, 0.0, 1.0, 0.0
        return np.nan, np.nan, np.nan, np.nan
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pearson_r, pearson_p = pearsonr(x, y)
        spearman_rho, spearman_p = spearmanr(x, y)
    return float(pearson_r), float(pearson_p), float(spearman_rho), float(spearman_p)


def _annotate_r_squared(ax: plt.Axes, pearson_r: float) -> None:
    r_squared = pearson_r * pearson_r if np.isfinite(pearson_r) else np.nan
    text = "$r^2$=NA" if not np.isfinite(r_squared) else f"$r^2$={r_squared:.3f}"
    ax.text(
        0.03,
        0.97,
        text,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9,
        bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "none"},
    )


def build_af_correlation_table(data: AggregatedData) -> pd.DataFrame:
    matched = data.matched_pairs.copy()
    if matched.empty:
        return pd.DataFrame(columns=["group", "n_matched", "pearson_r", "pearson_p", "spearman_rho", "spearman_p", "mean_abs_diff"])
    rows = []
    for group_name, frame in [("overall", matched)] + [(svtype, group) for svtype, group in matched.groupby("svtype_a")]:
        if frame.empty:
            continue
        valid = frame[["af_a", "af_b"]].dropna()
        x_values = valid["af_a"].astype(float)
        y_values = valid["af_b"].astype(float)
        pearson_r, pearson_p, spearman_rho, spearman_p = _compute_correlation_stats(x_values, y_values)
        rows.append(
            {
                "group": group_name,
                "n_matched": len(valid),
                "pearson_r": pearson_r,
                "pearson_p": pearson_p,
                "spearman_rho": spearman_rho,
                "spearman_p": spearman_p,
                "mean_abs_diff": float((x_values - y_values).abs().mean()) if len(valid) else np.nan,
            }
        )
    return pd.DataFrame(rows)


def _matched_pairs_with_categories(data: AggregatedData) -> pd.DataFrame:
    matched = data.matched_pairs.copy()
    if matched.empty:
        return matched
    site_meta = data.sites_a[["variant_id", "size_bucket", "af_bucket", "genomic_context"]].rename(
        columns={
            "variant_id": "variant_id_a",
            "size_bucket": "size_bucket_a",
            "af_bucket": "af_bucket_a",
            "genomic_context": "genomic_context_a",
        }
    )
    return matched.merge(site_meta, on="variant_id_a", how="left")


def _ordered_group_values(frame: pd.DataFrame, group_field: str) -> list[str]:
    values = frame[group_field].dropna().astype(str).unique()
    if group_field == "svtype_a":
        return [value for value in ordered_svtypes(values) if value != "CTX"]
    if group_field == "size_bucket_a":
        return ordered_plot_size_buckets(values)
    if group_field == "af_bucket_a":
        return ordered_plot_af_buckets(values)
    if group_field == "genomic_context_a":
        return ordered_contexts(values)
    return sorted(values)


def _plot_grouped_af_correlation(
    matched_pairs: pd.DataFrame,
    group_field: str,
    output_path,
    label_a: str,
    label_b: str,
) -> None:
    group_values = _ordered_group_values(matched_pairs, group_field)
    panel_count = max(len(group_values), 1)
    ncols = 3
    nrows = int(np.ceil(panel_count / ncols))
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=double_column_figsize(max(3.6, 2.1 * nrows)),
        squeeze=False,
        sharex=True,
        sharey=True,
        constrained_layout=True,
    )
    flat_axes = axes.flatten()
    for index, (axis, group_value) in enumerate(zip(flat_axes, group_values)):
        frame = matched_pairs.loc[matched_pairs[group_field].astype(str) == str(group_value), ["af_a", "af_b"]].dropna()
        if not frame.empty:
            plot_scatter_af(axis, frame["af_a"], frame["af_b"], label_a, label_b)
            axis.set_xlim(0.0, 1.0)
            axis.set_ylim(0.0, 1.0)
            pearson_r, _, _, _ = _compute_correlation_stats(frame["af_a"], frame["af_b"])
            _annotate_r_squared(axis, pearson_r)
        else:
            axis.text(0.5, 0.5, "No AF data", ha="center", va="center")
        axis.set_title(str(group_value))
        if index // ncols < nrows - 1:
            axis.set_xlabel("")
        if index % ncols != 0:
            axis.set_ylabel("")
    for axis in flat_axes[len(group_values):]:
        axis.set_axis_off()
    save_figure(fig, output_path)


class AlleleFreqModule(AnalysisModule):
    @property
    def name(self) -> str:
        return "allele_freq"

    def run(self, data: AggregatedData, config: AnalysisConfig) -> None:
        output_dir = self.output_dir(config)
        tables_dir = output_dir / "tables"
        tables_dir.mkdir(parents=True, exist_ok=True)
        stats = build_af_correlation_table(data)
        write_tsv_gz(stats, tables_dir / "af_correlation_stats.tsv")

        overall = _matched_pairs_with_categories(data)
        fig, ax = plt.subplots(figsize=single_column_figsize(3.2))
        overall_valid = overall[["af_a", "af_b"]].dropna() if not overall.empty else pd.DataFrame()
        if not overall_valid.empty:
            plot_scatter_af(ax, overall_valid["af_a"], overall_valid["af_b"], data.label_a, data.label_b)
            ax.set_xlim(0.0, 1.0)
            ax.set_ylim(0.0, 1.0)
            overall_row = stats.loc[stats["group"] == "overall"]
            pearson_r = float(overall_row.iloc[0]["pearson_r"]) if not overall_row.empty else np.nan
            _annotate_r_squared(ax, pearson_r)
        else:
            ax.text(0.5, 0.5, "No matched sites", ha="center", va="center")
        ax.set_title("AF correlation")
        save_figure(fig, output_dir / "af_correlation.overall.png")

        if not overall.empty:
            _plot_grouped_af_correlation(overall, "svtype_a", output_dir / "af_correlation.by_type.png", data.label_a, data.label_b)
            _plot_grouped_af_correlation(overall, "size_bucket_a", output_dir / "af_correlation.by_size.png", data.label_a, data.label_b)
            _plot_grouped_af_correlation(overall, "af_bucket_a", output_dir / "af_correlation.by_af.png", data.label_a, data.label_b)
            _plot_grouped_af_correlation(overall, "genomic_context_a", output_dir / "af_correlation.by_context.png", data.label_a, data.label_b)
