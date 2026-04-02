"""Size-signature diagnostics for insertion peaks and implausible variants."""

from __future__ import annotations

from typing import Dict, List

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.signal import find_peaks
from scipy.stats import ks_2samp

from ..aggregate import AggregatedData
from ..config import AnalysisConfig
from ..plot_utils import plot_peak_histogram, save_figure
from .base import AnalysisModule

_PEAK_WINDOWS = {
    "alu": (200.0, 400.0),
    "sva": (1500.0, 3000.0),
    "l1": (5000.0, 7000.0),
}


def quantify_retrotransposon_peaks(ins_svlens: np.ndarray) -> Dict[str, object]:
    cleaned = np.asarray([value for value in ins_svlens if value and value > 0], dtype=float)
    if cleaned.size == 0:
        return {
            "alu_peak_detected": False,
            "alu_peak_center": np.nan,
            "alu_peak_prominence": np.nan,
            "sva_peak_detected": False,
            "sva_peak_center": np.nan,
            "sva_peak_prominence": np.nan,
            "l1_peak_detected": False,
            "l1_peak_center": np.nan,
            "l1_peak_prominence": np.nan,
            "n_insertions": 0,
        }
    bins = np.logspace(np.log10(max(50.0, cleaned.min())), np.log10(max(100000.0, cleaned.max())), num=200)
    counts, edges = np.histogram(cleaned, bins=bins)
    peaks, properties = find_peaks(counts, prominence=1)
    metrics: Dict[str, object] = {"n_insertions": int(cleaned.size)}
    centers = (edges[:-1] + edges[1:]) / 2.0
    for prefix, (lower, upper) in _PEAK_WINDOWS.items():
        selected = [index for index in peaks if lower <= centers[index] <= upper]
        if not selected:
            metrics[f"{prefix}_peak_detected"] = False
            metrics[f"{prefix}_peak_center"] = np.nan
            metrics[f"{prefix}_peak_prominence"] = np.nan
            continue
        best = max(selected, key=lambda idx: properties["prominences"][list(peaks).index(idx)])
        prominence = float(properties["prominences"][list(peaks).index(best)])
        metrics[f"{prefix}_peak_detected"] = True
        metrics[f"{prefix}_peak_center"] = float(centers[best])
        metrics[f"{prefix}_peak_prominence"] = prominence
    return metrics


def flag_implausible_variants(sites: pd.DataFrame, contig_lengths: Dict[str, int], max_chrom_fraction: float = 0.20) -> pd.DataFrame:
    rows: List[dict] = []
    for _, row in sites.iterrows():
        contig_length = contig_lengths.get(str(row["contig"]))
        svlen = row.get("svlen")
        if contig_length is None or pd.isna(svlen) or float(svlen) <= 0:
            continue
        chrom_fraction = float(svlen) / float(contig_length)
        if chrom_fraction > max_chrom_fraction:
            rows.append(
                {
                    "variant_id": row["variant_id"],
                    "contig": row["contig"],
                    "svtype": row["svtype"],
                    "svlen": row["svlen"],
                    "chrom_fraction": chrom_fraction,
                }
            )
    return pd.DataFrame(rows)


def summarize_mei_subtypes(sites: pd.DataFrame) -> pd.DataFrame:
    mei_sites = sites.loc[sites["svtype"] == "INS:MEI"].copy()
    if mei_sites.empty:
        return pd.DataFrame(columns=["subtype", "n", "median_svlen", "iqr_svlen"])
    mei_sites["subtype"] = mei_sites["alt_allele"].astype(str)
    grouped = mei_sites.groupby("subtype")
    return grouped.agg(
        n=("variant_id", "count"),
        median_svlen=("svlen", "median"),
        iqr_svlen=("svlen", lambda values: float(values.quantile(0.75) - values.quantile(0.25))),
    ).reset_index()


class SizeSignaturesModule(AnalysisModule):
    @property
    def name(self) -> str:
        return "size_signatures"

    def run(self, data: AggregatedData, config: AnalysisConfig) -> None:
        output_dir = self.output_dir(config)
        tables_dir = output_dir / "tables"
        tables_dir.mkdir(parents=True, exist_ok=True)

        peak_rows = []
        for label, sites in ((data.label_a, data.sites_a), (data.label_b, data.sites_b)):
            filtered = sites.loc[sites["in_filtered_pass_view"]] if config.pass_only else sites
            ins_svlens = filtered.loc[filtered["svtype"].isin(["INS", "INS:MEI"]), "svlen"].dropna().astype(float).values
            peak_metrics = quantify_retrotransposon_peaks(ins_svlens)
            peak_metrics["label"] = label
            peak_rows.append(peak_metrics)

            fig, ax = plt.subplots(figsize=(7, 4))
            if ins_svlens.size:
                plot_peak_histogram(
                    ax,
                    ins_svlens,
                    np.logspace(np.log10(max(50.0, ins_svlens.min())), np.log10(max(100000.0, ins_svlens.max())), num=40),
                    peak_regions=[(200, 400, "#FFD54F"), (1500, 3000, "#FFB74D"), (5000, 7000, "#EF5350")],
                )
            else:
                ax.text(0.5, 0.5, "No insertions", ha="center", va="center")
            ax.set_title(f"Insertion size distribution: {label}")
            save_figure(fig, output_dir / f"retrotransposon_peaks.{label}.png")

            implausible = flag_implausible_variants(filtered, config.contig_lengths)
            implausible.to_csv(tables_dir / f"implausible_variants.{label}.tsv", sep="\t", index=False)
            fig, ax = plt.subplots(figsize=(7, 4))
            if not implausible.empty:
                ax.scatter(implausible["variant_id"], implausible["svlen"], color="red")
                ax.tick_params(axis="x", rotation=45)
            else:
                ax.text(0.5, 0.5, "No implausible variants", ha="center", va="center")
            ax.set_ylabel("SVLEN")
            ax.set_title(f"Implausible variants: {label}")
            save_figure(fig, output_dir / f"implausible_variants.{label}.png")

            mei_summary = summarize_mei_subtypes(filtered)
            mei_summary.to_csv(tables_dir / f"mei_subtype_summary.{label}.tsv", sep="\t", index=False)

        peak_table = pd.DataFrame(peak_rows)
        peak_table.to_csv(tables_dir / "retrotransposon_peaks.tsv", sep="\t", index=False)

        ins_a = data.sites_a.loc[data.sites_a["svtype"].isin(["INS", "INS:MEI"]), "svlen"].dropna().astype(float)
        ins_b = data.sites_b.loc[data.sites_b["svtype"].isin(["INS", "INS:MEI"]), "svlen"].dropna().astype(float)
        if len(ins_a) and len(ins_b):
            ks_statistic, ks_pvalue = ks_2samp(ins_a, ins_b)
        else:
            ks_statistic, ks_pvalue = np.nan, np.nan
        pd.DataFrame([
            {"svtype": "INS+INS:MEI", "ks_statistic": ks_statistic, "p_value": ks_pvalue}
        ]).to_csv(tables_dir / "size_distribution_comparison.tsv", sep="\t", index=False)

        fig, ax = plt.subplots(figsize=(7, 4))
        if len(ins_a):
            ax.hist(ins_a, bins=20, alpha=0.5, label=data.label_a)
        if len(ins_b):
            ax.hist(ins_b, bins=20, alpha=0.5, label=data.label_b)
        ax.set_xscale("log")
        ax.legend()
        ax.set_title("Insertion size overlay")
        save_figure(fig, output_dir / "ins_size_overlay.png")
