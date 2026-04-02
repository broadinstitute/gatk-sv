"""Hardy-Weinberg ternary and carrier-frequency diagnostics."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import chisquare

from ..aggregate import AggregatedData
from ..config import AnalysisConfig
from ..dimensions import af_bucket_sort_key, genomic_context_sort_key, size_bucket_sort_key, svtype_sort_key
from ..plot_utils import HWE_COLORS, SUMMARY_COLORS, plot_scatter_af, plot_ternary, save_figure, single_column_figsize
from .base import AnalysisModule, column_safe_label, write_tsv_gz


def _filtered_sites(sites: pd.DataFrame, pass_only: bool) -> pd.DataFrame:
    filtered = sites.loc[sites["svtype"].isin(["DEL", "DUP", "INS", "INS:MEI", "INV"])].copy()
    filtered = filtered.loc[filtered["n_bi_genos"] > 0]
    if pass_only:
        filtered = filtered.loc[filtered["in_filtered_pass_view"]]
    return filtered


def _hwe_p_value(row: pd.Series) -> float:
    total = int(row["n_bi_genos"])
    if total <= 0:
        return float("nan")
    hom_ref = float(row["n_hom_ref"])
    het = float(row["n_het"])
    hom_alt = float(row["n_hom_alt"])
    p_alt = (2.0 * hom_alt + het) / (2.0 * total)
    expected = np.asarray([
        total * (1.0 - p_alt) ** 2,
        total * 2.0 * p_alt * (1.0 - p_alt),
        total * p_alt**2,
    ])
    if np.any(expected == 0):
        return float("nan")
    observed = np.asarray([hom_ref, het, hom_alt])
    return float(chisquare(observed, expected).pvalue)


def build_hwe_table(sites: pd.DataFrame, pass_only: bool = False) -> pd.DataFrame:
    filtered = _filtered_sites(sites, pass_only)
    if filtered.empty:
        return pd.DataFrame(columns=["variant_id", "svtype", "size_bucket", "af_bucket", "genomic_context", "aa", "ab", "bb", "carrier_freq", "af", "hwe_p", "hwe_class"])
    bonferroni_threshold = 0.05 / max(len(filtered), 1)
    table = filtered[["variant_id", "svtype", "size_bucket", "af_bucket", "genomic_context", "af", "n_bi_genos", "n_hom_ref", "n_het", "n_hom_alt"]].copy()
    table["aa"] = table["n_hom_ref"] / table["n_bi_genos"]
    table["ab"] = table["n_het"] / table["n_bi_genos"]
    table["bb"] = table["n_hom_alt"] / table["n_bi_genos"]
    table["carrier_freq"] = (table["n_het"] + table["n_hom_alt"]) / table["n_bi_genos"]
    table["hwe_p"] = table.apply(_hwe_p_value, axis=1)
    table["hwe_class"] = np.where(
        table["hwe_p"] < bonferroni_threshold,
        "bonferroni",
        np.where(table["hwe_p"] < 0.05, "nominal", "pass"),
    )
    return table[["variant_id", "svtype", "size_bucket", "af_bucket", "genomic_context", "aa", "ab", "bb", "carrier_freq", "af", "hwe_p", "hwe_class"]]


def summarize_hwe_by_bucket(table: pd.DataFrame, label: str) -> pd.DataFrame:
    token = column_safe_label(label)
    columns = [
        "svtype",
        "size_bucket",
        "af_bucket",
        "genomic_context",
        f"n_variants_{token}",
        f"frac_pass_{token}",
        f"frac_nominal_{token}",
        f"frac_bonferroni_{token}",
        f"mean_carrier_freq_{token}",
        f"mean_af_{token}",
    ]
    if table.empty:
        return pd.DataFrame(columns=columns)
    grouped = table.groupby(["svtype", "size_bucket", "af_bucket", "genomic_context"], dropna=False).agg(
        n_variants=("variant_id", "count"),
        frac_pass=("hwe_class", lambda values: float((values == "pass").mean())),
        frac_nominal=("hwe_class", lambda values: float((values == "nominal").mean())),
        frac_bonferroni=("hwe_class", lambda values: float((values == "bonferroni").mean())),
        mean_carrier_freq=("carrier_freq", "mean"),
        mean_af=("af", "mean"),
    ).reset_index()
    return grouped.rename(columns={
        "n_variants": f"n_variants_{token}",
        "frac_pass": f"frac_pass_{token}",
        "frac_nominal": f"frac_nominal_{token}",
        "frac_bonferroni": f"frac_bonferroni_{token}",
        "mean_carrier_freq": f"mean_carrier_freq_{token}",
        "mean_af": f"mean_af_{token}",
    })[columns]
    return grouped.sort_values(
        by=["svtype", "size_bucket", "af_bucket", "genomic_context"],
        key=lambda series: series.map(
            lambda value: (
                svtype_sort_key(value) if series.name == "svtype" else
                size_bucket_sort_key(value) if series.name == "size_bucket" else
                af_bucket_sort_key(value) if series.name == "af_bucket" else
                genomic_context_sort_key(value)
            )
        ),
    ).reset_index(drop=True)


def _annotate_hwe_summary(ax: plt.Axes, table: pd.DataFrame) -> None:
    if table.empty:
        return
    nominal_in_hwe_fraction = float((table["hwe_p"] >= 0.05).mean())
    bonf_in_hwe_fraction = float((table["hwe_class"] != "bonferroni").mean())
    ax.text(
        0.02,
        1.04,
        f"Nom. HWE={nominal_in_hwe_fraction:.1%}\nBonf. HWE={bonf_in_hwe_fraction:.1%}",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=6,
        clip_on=False,
    )


class GenotypeDistModule(AnalysisModule):
    @property
    def name(self) -> str:
        return "genotype_dist"

    def _run_one(self, sites: pd.DataFrame, label: str, config: AnalysisConfig, output_dir: Path) -> pd.DataFrame:
        table = build_hwe_table(sites, pass_only=config.pass_only)

        fig, ax = plt.subplots(figsize=single_column_figsize(3.0))
        if not table.empty:
            colors = [HWE_COLORS.get(value, SUMMARY_COLORS["neutral"]) for value in table["hwe_class"]]
            plot_ternary(ax, table["aa"], table["ab"], table["bb"], colors)
            _annotate_hwe_summary(ax, table)
        else:
            ax.text(0.5, 0.5, "No eligible variants", ha="center", va="center")
        ax.set_title(label)
        save_figure(fig, output_dir / f"ternary.all.{label}.png")

        fig, ax = plt.subplots(figsize=single_column_figsize(3.0))
        if not table.empty:
            plot_scatter_af(ax, table["carrier_freq"], table["af"].fillna(0.0), "Carrier frequency", "Allele frequency")
        else:
            ax.text(0.5, 0.5, "No eligible variants", ha="center", va="center")
        ax.set_title(label)
        save_figure(fig, output_dir / f"carrier_freq_vs_af.{label}.png")
        return table.assign(label=label)

    def run(self, data: AggregatedData, config: AnalysisConfig) -> None:
        output_dir = self.output_dir(config)
        tables_dir = output_dir / "tables"
        tables_dir.mkdir(parents=True, exist_ok=True)
        table_a = self._run_one(data.sites_a, data.label_a, config, output_dir)
        table_b = self._run_one(data.sites_b, data.label_b, config, output_dir)
        summary_a = summarize_hwe_by_bucket(table_a, data.label_a)
        summary_b = summarize_hwe_by_bucket(table_b, data.label_b)
        combined = summary_a.merge(summary_b, on=["svtype", "size_bucket", "af_bucket", "genomic_context"], how="outer")
        write_tsv_gz(combined, tables_dir / "hwe_stats.tsv")
