"""Per-sample site and allele counts from the source VCFs."""

from __future__ import annotations

from pathlib import Path
from typing import Dict

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam
from matplotlib.gridspec import GridSpec
from matplotlib.transforms import offset_copy

from ..config import AnalysisConfig
from ..dimensions import normalize_svtype, ordered_svtypes
from ..plot_utils import CALLSET_COLORS, double_column_figsize, save_figure
from ..vcf_format import filter_values, safe_info_get
from .base import AnalysisModule, write_tsv_gz


def collect_per_sample_counts(vcf_path: Path, pass_only: bool = False) -> pd.DataFrame:
    with pysam.VariantFile(str(vcf_path)) as vcf:
        sample_names = list(vcf.header.samples)
        per_sample: Dict[tuple, Dict[str, int]] = {}
        for record in vcf:
            filters = filter_values(record)
            if pass_only and not ({"PASS", "MULTIALLELIC"} & filters):
                continue
            svtype = normalize_svtype(str(safe_info_get(record, "SVTYPE", "UNKNOWN")), ",".join(record.alts or ()))
            for sample_name in sample_names:
                sample = record.samples[sample_name]
                key = (sample_name, svtype)
                accumulator = per_sample.setdefault(key, {"sites": 0, "alleles": 0})
                gt = sample.get("GT")
                if gt and not any(allele is None for allele in gt):
                    alt_count = sum(1 for allele in gt if allele and allele > 0)
                    if alt_count > 0:
                        accumulator["sites"] += 1
                        accumulator["alleles"] += alt_count
                        continue
                cn = sample.get("CN")
                ecn = sample.get("ECN")
                if cn not in (None, ".") and ecn not in (None, ".") and int(cn) != int(ecn):
                    accumulator["sites"] += 1
                    accumulator["alleles"] += abs(int(cn) - int(ecn))
        rows = [
            {
                "sample": sample,
                "svtype": svtype,
                "sites": values["sites"],
                "alleles": values["alleles"],
            }
            for (sample, svtype), values in per_sample.items()
        ]
        return pd.DataFrame(rows)


def _sample_distribution(counts: pd.DataFrame, value_column: str, svtype: str | None = None) -> np.ndarray:
    if counts.empty:
        return np.asarray([], dtype=float)
    if svtype is None:
        grouped = counts.groupby("sample", dropna=False)[value_column].sum()
    else:
        grouped = counts.loc[counts["svtype"] == svtype].groupby("sample", dropna=False)[value_column].sum()
    return grouped.to_numpy(dtype=float)


def _draw_swarm_with_summary(
    ax: plt.Axes,
    values: np.ndarray,
    x_center: float,
    color: str,
    rng: np.random.Generator,
    median_label_dx: float = 0.05,
    median_label_dy_fraction: float = 0.035,
    median_label_dx_pixels: float = 10.0,
) -> None:
    if values.size == 0:
        return
    jitter = rng.uniform(-0.045, 0.045, size=values.size)
    ax.scatter(np.full(values.size, x_center) + jitter, values, s=12, alpha=0.8, color=color, edgecolors="none", zorder=2)
    q1, median, q3 = np.quantile(values, [0.25, 0.5, 0.75])
    y_span = max(float(np.max(values) - np.min(values)), 1.0)
    ax.vlines(x_center, q1, q3, colors="black", linewidth=1.2, zorder=3)
    ax.hlines([q1, q3], x_center - 0.045, x_center + 0.045, colors="black", linewidth=1.0, zorder=3)
    ax.hlines(median, x_center - 0.07, x_center + 0.07, colors="black", linewidth=1.8, zorder=4)
    label_transform = offset_copy(ax.transData, fig=ax.figure, x=median_label_dx_pixels, y=0.0, units="dots")
    ax.text(
        x_center + median_label_dx,
        median - y_span * median_label_dy_fraction,
        f"{median:,.0f}",
        transform=label_transform,
        color="black",
        fontsize=8,
        va="center",
        ha="left",
        fontweight="bold",
    )


def _plot_side_by_side_swarms(
    counts_a: pd.DataFrame,
    counts_b: pd.DataFrame,
    value_column: str,
    output_path: Path,
    y_label: str,
    label_a: str,
    label_b: str,
) -> None:
    svtypes = ordered_svtypes(set(counts_a["svtype"]).union(set(counts_b["svtype"])))
    fig = plt.figure(figsize=double_column_figsize(3.8))
    grid = GridSpec(1, 2, figure=fig, width_ratios=[1.35, max(4.8, len(svtypes))], wspace=0.24)
    ax_all = fig.add_subplot(grid[0, 0])
    ax_types = fig.add_subplot(grid[0, 1])
    if not svtypes:
        ax_types.text(0.5, 0.5, "No sample counts", ha="center", va="center")
        ax_all.set_axis_off()
        ax_types.set_axis_off()
        save_figure(fig, output_path)
        return
    rng = np.random.default_rng(0)
    width = 0.18

    overall_a = _sample_distribution(counts_a, value_column)
    overall_b = _sample_distribution(counts_b, value_column)
    _draw_swarm_with_summary(ax_all, overall_a, -width, CALLSET_COLORS["a"], rng, median_label_dx=0.07)
    _draw_swarm_with_summary(ax_all, overall_b, width, CALLSET_COLORS["b"], rng, median_label_dx=0.07)
    ax_all.set_xlim(-0.45, 0.45)
    ax_all.set_xticks([0.0])
    ax_all.set_xticklabels(["ALL"], rotation=90)
    ax_all.set_ylabel(y_label)
    ax_all.grid(axis="y", alpha=0.25)

    base_positions = np.arange(len(svtypes), dtype=float)
    for index, svtype in enumerate(svtypes):
        values_a = _sample_distribution(counts_a, value_column, svtype)
        values_b = _sample_distribution(counts_b, value_column, svtype)
        _draw_swarm_with_summary(ax_types, values_a, base_positions[index] - width, CALLSET_COLORS["a"], rng, median_label_dx=0.045)
        _draw_swarm_with_summary(ax_types, values_b, base_positions[index] + width, CALLSET_COLORS["b"], rng, median_label_dx=0.045)
    ax_types.set_xticks(base_positions)
    ax_types.set_xticklabels(svtypes, rotation=90)
    ax_types.set_xlabel("SV type")
    ax_types.grid(axis="y", alpha=0.25)
    ax_types.set_ylabel(y_label)
    ax_types.legend(
        handles=[
            plt.Line2D([0], [0], marker="o", linestyle="None", color=CALLSET_COLORS["a"], label=label_a, markersize=6),
            plt.Line2D([0], [0], marker="o", linestyle="None", color=CALLSET_COLORS["b"], label=label_b, markersize=6),
        ],
        loc="upper right",
        fontsize=8,
    )
    title_prefix = "SV Sites per Genome" if value_column == "sites" else "SV Alleles per Genome"
    sample_count = max(len(set(counts_a.get("sample", pd.Series(dtype=object)))), len(set(counts_b.get("sample", pd.Series(dtype=object)))))
    fig.suptitle(title_prefix, fontsize=16, fontweight="bold", y=0.98)
    ax_types.set_title(f"N={sample_count} Samples", fontsize=12, pad=6)
    save_figure(fig, output_path)


class CountsPerGenomeModule(AnalysisModule):
    @property
    def name(self) -> str:
        return "counts_per_genome"

    @property
    def requires_genotype_pass(self) -> bool:
        return True

    def run(self, data, config: AnalysisConfig) -> None:
        del data
        output_dir = self.output_dir(config)
        tables_dir = output_dir / "tables"
        tables_dir.mkdir(parents=True, exist_ok=True)
        counts_by_label: Dict[str, pd.DataFrame] = {}
        for label, vcf_path in ((config.vcf_a_label, config.vcf_a_path), (config.vcf_b_label, config.vcf_b_path)):
            if vcf_path is None:
                continue
            counts = collect_per_sample_counts(vcf_path, pass_only=config.pass_only)
            counts_by_label[label] = counts
            if config.per_sample_counts_table:
                write_tsv_gz(counts, tables_dir / f"per_sample_counts.{label}.tsv")

        counts_a = counts_by_label.get(config.vcf_a_label, pd.DataFrame(columns=["sample", "svtype", "sites", "alleles"]))
        counts_b = counts_by_label.get(config.vcf_b_label, pd.DataFrame(columns=["sample", "svtype", "sites", "alleles"]))
        _plot_side_by_side_swarms(
            counts_a,
            counts_b,
            "sites",
            output_dir / "sites_per_genome.by_type.png",
            "Sites per sample",
            config.vcf_a_label,
            config.vcf_b_label,
        )
        _plot_side_by_side_swarms(
            counts_a,
            counts_b,
            "alleles",
            output_dir / "alleles_per_genome.by_type.png",
            "Alleles per sample",
            config.vcf_a_label,
            config.vcf_b_label,
        )
