"""Per-sample site and allele counts from the source VCFs."""

from __future__ import annotations

from pathlib import Path
from typing import Dict

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam

from ..config import AnalysisConfig
from ..dimensions import normalize_svtype
from ..vcf_format import filter_values
from .base import AnalysisModule, write_tsv_gz


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


def _plot_side_by_side_boxplots(
    counts_a: pd.DataFrame,
    counts_b: pd.DataFrame,
    value_column: str,
    output_path: Path,
    y_label: str,
    label_a: str,
    label_b: str,
) -> None:
    svtypes = sorted(set(counts_a["svtype"]).union(set(counts_b["svtype"])))
    fig, ax = plt.subplots(figsize=(max(8, 1.4 * len(svtypes)), 5))
    if not svtypes:
        ax.text(0.5, 0.5, "No sample counts", ha="center", va="center")
        ax.set_axis_off()
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close(fig)
        return
    width = 0.32
    base_positions = np.arange(len(svtypes), dtype=float)
    box_data = []
    positions = []
    colors = []
    for index, svtype in enumerate(svtypes):
        values_a = counts_a.loc[counts_a["svtype"] == svtype, value_column].to_numpy(dtype=float)
        values_b = counts_b.loc[counts_b["svtype"] == svtype, value_column].to_numpy(dtype=float)
        box_data.extend([values_a if values_a.size else np.array([np.nan]), values_b if values_b.size else np.array([np.nan])])
        positions.extend([base_positions[index] - width / 2, base_positions[index] + width / 2])
        colors.extend(["#4C72B0", "#DD8452"])
    boxplot = ax.boxplot(box_data, positions=positions, widths=width * 0.85, patch_artist=True, manage_ticks=False)
    for patch, color in zip(boxplot["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.75)
    ax.set_xticks(base_positions)
    ax.set_xticklabels(svtypes)
    ax.set_ylabel(y_label)
    ax.set_xlabel("SV type")
    ax.legend(
        handles=[
            plt.Line2D([0], [0], color="#4C72B0", lw=8, label=label_a),
            plt.Line2D([0], [0], color="#DD8452", lw=8, label=label_b),
        ],
        fontsize=8,
    )
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


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
        _plot_side_by_side_boxplots(
            counts_a,
            counts_b,
            "sites",
            output_dir / "sites_per_genome.by_type.png",
            "Sites per sample",
            config.vcf_a_label,
            config.vcf_b_label,
        )
        _plot_side_by_side_boxplots(
            counts_a,
            counts_b,
            "alleles",
            output_dir / "alleles_per_genome.by_type.png",
            "Alleles per sample",
            config.vcf_a_label,
            config.vcf_b_label,
        )
