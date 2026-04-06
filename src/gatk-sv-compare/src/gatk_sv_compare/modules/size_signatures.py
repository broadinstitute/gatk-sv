"""Size-signature diagnostics for implausible variants and MEI subtype summaries."""

from __future__ import annotations

from typing import Dict, List

import numpy as np
import pandas as pd
from scipy.signal import find_peaks

from ..aggregate import AggregatedData
from ..config import AnalysisConfig
from .base import AnalysisModule, write_tsv_gz

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

        for label, sites in ((data.label_a, data.sites_a), (data.label_b, data.sites_b)):
            filtered = sites.loc[sites["in_filtered_pass_view"]] if config.pass_only else sites
            implausible = flag_implausible_variants(filtered, config.contig_lengths)
            write_tsv_gz(implausible, output_dir / f"implausible_variants.{label}.tsv")

            mei_summary = summarize_mei_subtypes(filtered)
            write_tsv_gz(mei_summary, output_dir / f"mei_subtype_summary.{label}.tsv")

