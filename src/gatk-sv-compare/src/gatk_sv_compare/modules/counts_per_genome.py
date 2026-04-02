"""Per-sample site and allele counts from the source VCFs."""

from __future__ import annotations

from pathlib import Path
from typing import Dict

import matplotlib.pyplot as plt
import pandas as pd
import pysam

from ..config import AnalysisConfig
from ..dimensions import normalize_svtype
from ..vcf_format import filter_values
from .base import AnalysisModule


def collect_per_sample_counts(vcf_path: Path, pass_only: bool = False) -> pd.DataFrame:
    with pysam.VariantFile(str(vcf_path)) as vcf:
        sample_names = list(vcf.header.samples)
        per_sample: Dict[tuple, Dict[str, int]] = {}
        for record in vcf:
            filters = filter_values(record)
            if pass_only and not ({"PASS", "MULTIALLELIC"} & filters):
                continue
            svtype = normalize_svtype(str(record.info.get("SVTYPE", "UNKNOWN")), ",".join(record.alts or ()))
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
        for label, vcf_path in ((config.vcf_a_label, config.vcf_a_path), (config.vcf_b_label, config.vcf_b_path)):
            if vcf_path is None:
                continue
            counts = collect_per_sample_counts(vcf_path, pass_only=config.pass_only)
            counts.to_csv(tables_dir / f"per_sample_counts.{label}.tsv", sep="\t", index=False)
            if counts.empty:
                continue
            site_matrix = [group["sites"].values for _, group in counts.groupby("svtype")]
            allele_matrix = [group["alleles"].values for _, group in counts.groupby("svtype")]
            labels = list(counts.groupby("svtype").groups.keys())

            fig, ax = plt.subplots(figsize=(8, 4))
            ax.boxplot(site_matrix, tick_labels=labels)
            ax.set_ylabel("Sites per sample")
            ax.set_title(f"Sites per genome by type: {label}")
            fig.savefig(output_dir / f"sites_per_genome.by_type.{label}.png", dpi=300, bbox_inches="tight")
            plt.close(fig)

            fig, ax = plt.subplots(figsize=(8, 4))
            ax.boxplot(allele_matrix, tick_labels=labels)
            ax.set_ylabel("Alleles per sample")
            ax.set_title(f"Alleles per genome by type: {label}")
            fig.savefig(output_dir / f"alleles_per_genome.by_type.{label}.png", dpi=300, bbox_inches="tight")
            plt.close(fig)
