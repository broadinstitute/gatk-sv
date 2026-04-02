"""GQ distribution summaries derived directly from the source VCFs."""

from __future__ import annotations

from pathlib import Path
from typing import List

import matplotlib.pyplot as plt
import pandas as pd
import pysam

from ..config import AnalysisConfig
from ..dimensions import normalize_svtype
from ..vcf_format import filter_values
from .base import AnalysisModule, write_tsv_gz


def _iter_alt_gq_rows(vcf_path: Path, pass_only: bool) -> List[dict]:
    rows: List[dict] = []
    with pysam.VariantFile(str(vcf_path)) as vcf:
        for record in vcf:
            filters = filter_values(record)
            if pass_only and not ({"PASS", "MULTIALLELIC"} & filters):
                continue
            svtype = normalize_svtype(str(record.info.get("SVTYPE", "UNKNOWN")), ",".join(record.alts or ()))
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
    for group_name, group in [("overall", frame)] + [(svtype, group) for svtype, group in frame.groupby("svtype")]:
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
            fig, ax = plt.subplots(figsize=(6, 4))
            if rows:
                frame = pd.DataFrame(rows)
                ax.hist(frame["gq"], bins=list(range(0, 105, 5)), color="#4C72B0", alpha=0.8)
            else:
                ax.text(0.5, 0.5, "No alt genotypes", ha="center", va="center")
            ax.set_xlabel("GQ")
            ax.set_ylabel("Count")
            ax.set_title(f"GQ histogram: {label}")
            fig.savefig(output_dir / f"gq_histogram.overall.{label}.png", dpi=300, bbox_inches="tight")
            plt.close(fig)
