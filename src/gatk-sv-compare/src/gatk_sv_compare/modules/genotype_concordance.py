"""SVConcordance INFO-metric summaries for shared-sample callsets."""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path
from typing import Dict

import matplotlib.pyplot as plt
import pandas as pd
import pysam

from ..aggregate import AggregatedData
from ..config import AnalysisConfig
from ..vcf_reader import _CONCORDANCE_FIELDS
from .base import AnalysisModule


def _record_identifier(record: pysam.VariantRecord) -> str:
    return record.id or f"{record.contig}:{record.pos}"


def _iter_records_for_contig(vcf: pysam.VariantFile, contig: str):
    try:
        yield from vcf.fetch(contig)
    except ValueError:
        for record in vcf:
            if record.contig == contig:
                yield record


def _safe_info_value(record: pysam.VariantRecord, key: str):
    try:
        value = record.info.get(key)
    except ValueError:
        return None
    if isinstance(value, tuple):
        return value[0] if value else None
    return value


def _extract_concordance_rows(vcf_path: Path, sites: pd.DataFrame, label: str) -> pd.DataFrame:
    target_ids_by_contig: Dict[str, set[str]] = defaultdict(set)
    for row in sites[["variant_id", "contig"]].drop_duplicates().itertuples(index=False):
        target_ids_by_contig[str(row.contig)].add(str(row.variant_id))

    site_meta = sites.set_index("variant_id")[["svtype", "size_bucket", "af_bucket", "genomic_context"]]
    rows = []
    with pysam.VariantFile(str(vcf_path)) as vcf:
        for contig, target_ids in target_ids_by_contig.items():
            for record in _iter_records_for_contig(vcf, contig):
                variant_id = _record_identifier(record)
                if variant_id not in target_ids:
                    continue
                metrics = {
                    field.lower(): _safe_info_value(record, field)
                    for field in _CONCORDANCE_FIELDS
                    if _safe_info_value(record, field) not in (None, ".")
                }
                if not metrics:
                    continue
                meta = site_meta.loc[variant_id]
                rows.append({
                    "label": label,
                    "variant_id": variant_id,
                    "svtype": meta["svtype"],
                    "size_bucket": meta["size_bucket"],
                    "af_bucket": meta["af_bucket"],
                    "genomic_context": meta["genomic_context"],
                    **{key: float(value) for key, value in metrics.items()},
                })
    return pd.DataFrame(rows)


def summarize_concordance_metrics(metrics: pd.DataFrame) -> pd.DataFrame:
    if metrics.empty:
        return pd.DataFrame(columns=["label", "metric", "svtype", "size_bucket", "af_bucket", "genomic_context", "n_sites", "mean_value", "median_value"])
    metric_columns = [column for column in metrics.columns if column in {field.lower() for field in _CONCORDANCE_FIELDS}]
    long_frame = metrics.melt(
        id_vars=["label", "variant_id", "svtype", "size_bucket", "af_bucket", "genomic_context"],
        value_vars=metric_columns,
        var_name="metric",
        value_name="value",
    ).dropna(subset=["value"])
    summary = long_frame.groupby(["label", "metric", "svtype", "size_bucket", "af_bucket", "genomic_context"], dropna=False).agg(
        n_sites=("variant_id", "count"),
        mean_value=("value", "mean"),
        median_value=("value", "median"),
    ).reset_index()
    return summary


def _plot_metric_by_type(metrics: pd.DataFrame, label: str, metric: str, output_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 4))
    label_metrics = metrics.loc[metrics["label"] == label].copy()
    if label_metrics.empty or metric not in label_metrics.columns:
        ax.text(0.5, 0.5, "No concordance metrics", ha="center", va="center")
        ax.set_axis_off()
    else:
        grouped = label_metrics.groupby("svtype", dropna=False)[metric].median().sort_index()
        ax.bar(grouped.index.astype(str), grouped.values, color="#55A868")
        ax.set_ylim(0.0, 1.0)
        ax.set_ylabel(metric.replace("_", " "))
        ax.set_xlabel("SV type")
        ax.set_title(f"{metric.replace('_', ' ')} by type: {label}")
        ax.tick_params(axis="x", rotation=30)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


class GenotypeConcordanceModule(AnalysisModule):
    @property
    def name(self) -> str:
        return "genotype_concordance"

    @property
    def requires_shared_samples(self) -> bool:
        return True

    @property
    def requires_concordance(self) -> bool:
        return True

    def run(self, data: AggregatedData, config: AnalysisConfig) -> None:
        output_dir = self.output_dir(config)
        tables_dir = output_dir / "tables"
        tables_dir.mkdir(parents=True, exist_ok=True)
        metrics_a = _extract_concordance_rows(config.vcf_a_path, data.sites_a, data.label_a)
        metrics_b = _extract_concordance_rows(config.vcf_b_path, data.sites_b, data.label_b)
        metrics = pd.concat([metrics_a, metrics_b], ignore_index=True) if not metrics_a.empty or not metrics_b.empty else pd.DataFrame()
        summary = summarize_concordance_metrics(metrics)
        summary.to_csv(tables_dir / "concordance_metrics.tsv", sep="\t", index=False)

        for label in [data.label_a, data.label_b]:
            for metric in ["genotype_concordance", "non_ref_genotype_concordance", "var_ppv", "var_sensitivity"]:
                _plot_metric_by_type(metrics, label, metric, output_dir / f"{metric}.by_type.{label}.png")
