"""GQ distribution summaries derived directly from the source VCFs."""

from __future__ import annotations

from pathlib import Path
from typing import List

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam

from ..config import AnalysisConfig
from ..dimensions import normalize_svtype, ordered_svtypes
from ..plot_utils import SUMMARY_COLORS, SVTYPE_COLORS, save_figure, single_column_figsize
from ..vcf_format import filter_values, safe_info_get
from .base import AnalysisModule, write_tsv_gz


def _iter_alt_gq_rows(vcf_path: Path, pass_only: bool) -> List[dict]:
    rows: List[dict] = []
    with pysam.VariantFile(str(vcf_path)) as vcf:
        for record in vcf:
            filters = filter_values(record)
            if pass_only and not ({"PASS", "MULTIALLELIC"} & filters):
                continue
            svtype = normalize_svtype(str(safe_info_get(record, "SVTYPE", "UNKNOWN")), ",".join(record.alts or ()))
            for sample in record.samples.values():
                gt = sample.get("GT")
                gq = sample.get("GQ")
                if gq in (None, ".") or not gt or any(allele is None for allele in gt):
                    continue
                alt_count = sum(1 for allele in gt if allele and allele > 0)
                if alt_count > 0:
                    rows.append({"svtype": svtype, "gq": float(gq)})
    return rows


def summarize_gq(vcf_path: Path, pass_only: bool = False) -> pd.DataFrame:
    rows = _iter_alt_gq_rows(vcf_path, pass_only)
    if not rows:
        return pd.DataFrame(columns=["group", "n", "mean_gq", "median_gq", "q25_gq", "q75_gq"])
    frame = pd.DataFrame(rows)
    summaries = []
    for group_name, group in [("overall", frame)] + [(svtype, frame.loc[frame["svtype"] == svtype]) for svtype in ordered_svtypes(frame["svtype"].unique())]:
        summaries.append(
            {
                "group": group_name,
                "n": len(group),
                "mean_gq": float(group["gq"].mean()),
                "median_gq": float(group["gq"].median()),
                "q25_gq": float(group["gq"].quantile(0.25)),
                "q75_gq": float(group["gq"].quantile(0.75)),
            }
        )
    return pd.DataFrame(summaries)


def _nice_histogram_step(max_value: float, target_bins: int = 20) -> float:
    if not np.isfinite(max_value) or max_value <= 0:
        return 5.0
    rough_step = max_value / max(target_bins, 1)
    magnitude = 10 ** np.floor(np.log10(max(rough_step, 1.0)))
    normalized = rough_step / magnitude
    if normalized <= 1.0:
        multiplier = 1.0
    elif normalized <= 2.0:
        multiplier = 2.0
    elif normalized <= 5.0:
        multiplier = 5.0
    else:
        multiplier = 10.0
    return float(multiplier * magnitude)


def _gq_histogram_bins(values: np.ndarray) -> np.ndarray:
    finite_values = np.asarray(values, dtype=float)
    finite_values = finite_values[np.isfinite(finite_values)]
    if finite_values.size == 0:
        return np.asarray([0.0, 5.0], dtype=float)
    max_value = float(finite_values.max())
    step = _nice_histogram_step(max_value)
    upper = max(step, np.ceil(max_value / step) * step)
    return np.arange(0.0, upper + step, step, dtype=float)


class GenotypeQualityModule(AnalysisModule):
    @property
    def name(self) -> str:
        return "genotype_quality"

    @property
    def requires_gq(self) -> bool:
        return True

    def run(self, data, config: AnalysisConfig) -> None:
        del data
        output_dir = self.output_dir(config)
        tables_dir = output_dir / "tables"
        tables_dir.mkdir(parents=True, exist_ok=True)
        for label, vcf_path in ((config.vcf_a_label, config.vcf_a_path), (config.vcf_b_label, config.vcf_b_path)):
            if vcf_path is None:
                continue
            summary = summarize_gq(vcf_path, pass_only=config.pass_only)
            write_tsv_gz(summary, tables_dir / f"gq_summary.{label}.tsv")

            rows = _iter_alt_gq_rows(vcf_path, config.pass_only)
            fig, ax = plt.subplots(figsize=single_column_figsize(2.8))
            if rows:
                frame = pd.DataFrame(rows)
                svtypes = ordered_svtypes(frame["svtype"].dropna().unique())
                values = [frame.loc[frame["svtype"] == svtype, "gq"].to_numpy(dtype=float) for svtype in svtypes]
                histogram_bins = _gq_histogram_bins(frame["gq"].to_numpy(dtype=float))
                colors = [SVTYPE_COLORS.get(str(svtype), SUMMARY_COLORS["neutral"]) for svtype in svtypes]
                ax.hist(
                    values,
                    bins=histogram_bins,
                    stacked=True,
                    color=colors,
                    alpha=0.85,
                    edgecolor=SUMMARY_COLORS["edge"],
                    linewidth=0.6,
                    label=[str(svtype) for svtype in svtypes],
                )
                ax.set_xlim(float(histogram_bins[0]), float(histogram_bins[-1]))
                if svtypes:
                    ax.legend(fontsize=8)
            else:
                ax.text(0.5, 0.5, "No alt genotypes", ha="center", va="center")
            ax.set_xlabel("GQ")
            ax.set_ylabel("Count")
            ax.set_title(label)
            save_figure(fig, output_dir / f"gq_histogram.overall.{label}.png")
