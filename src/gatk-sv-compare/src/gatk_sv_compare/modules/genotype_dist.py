"""Hardy-Weinberg ternary and carrier-frequency diagnostics."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import chisquare

from ..aggregate import AggregatedData
from ..config import AnalysisConfig
from ..plot_utils import HWE_COLORS, plot_scatter_af, plot_ternary, save_figure
from .base import AnalysisModule


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
        return pd.DataFrame(columns=["variant_id", "svtype", "size_bucket", "aa", "ab", "bb", "carrier_freq", "af", "hwe_p", "hwe_class"])
    bonferroni_threshold = 0.05 / max(len(filtered), 1)
    table = filtered[["variant_id", "svtype", "size_bucket", "af", "n_bi_genos", "n_hom_ref", "n_het", "n_hom_alt"]].copy()
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
    return table[["variant_id", "svtype", "size_bucket", "aa", "ab", "bb", "carrier_freq", "af", "hwe_p", "hwe_class"]]


class GenotypeDistModule(AnalysisModule):
    @property
    def name(self) -> str:
        return "genotype_dist"

    def _run_one(self, sites: pd.DataFrame, label: str, config: AnalysisConfig, output_dir: Path) -> None:
        table = build_hwe_table(sites, pass_only=config.pass_only)
        tables_dir = output_dir / "tables"
        tables_dir.mkdir(parents=True, exist_ok=True)
        out_name = "hwe_stats.tsv" if label == "combined" else f"hwe_stats.{label}.tsv"
        table.to_csv(tables_dir / out_name, sep="\t", index=False)

        fig, ax = plt.subplots(figsize=(6, 5))
        if not table.empty:
            colors = [HWE_COLORS.get(value, "#999999") for value in table["hwe_class"]]
            plot_ternary(ax, table["aa"], table["ab"], table["bb"], colors)
        else:
            ax.text(0.5, 0.5, "No eligible variants", ha="center", va="center")
        ax.set_title(f"Genotype ternary: {label}")
        save_figure(fig, output_dir / f"ternary.all.{label}.png")

        fig, ax = plt.subplots(figsize=(6, 5))
        if not table.empty:
            plot_scatter_af(ax, table["carrier_freq"], table["af"].fillna(0.0), "Carrier frequency", "Allele frequency")
        else:
            ax.text(0.5, 0.5, "No eligible variants", ha="center", va="center")
        ax.set_title(f"Carrier frequency vs AF: {label}")
        save_figure(fig, output_dir / f"carrier_freq_vs_af.{label}.png")

    def run(self, data: AggregatedData, config: AnalysisConfig) -> None:
        output_dir = self.output_dir(config)
        self._run_one(data.sites_a, data.label_a, config, output_dir)
        self._run_one(data.sites_b, data.label_b, config, output_dir)
