"""SVConcordance INFO-metric summaries for shared-sample callsets."""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path
from typing import Dict

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam

from ..aggregate import AggregatedData
from ..config import AnalysisConfig
from ..vcf_reader import _CONCORDANCE_FIELDS
from .base import AnalysisModule, relabel_vcf_columns, write_tsv_gz


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


def _extract_concordance_rows(vcf_path: Path, sites: pd.DataFrame, label: str, source: str) -> pd.DataFrame:
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
                    "source": source,
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
        return pd.DataFrame(columns=["metric", "svtype", "size_bucket", "af_bucket", "genomic_context", "n_sites_a", "mean_value_a", "median_value_a", "n_sites_b", "mean_value_b", "median_value_b"])
    metric_columns = [column for column in metrics.columns if column in {field.lower() for field in _CONCORDANCE_FIELDS}]
    long_frame = metrics.melt(
        id_vars=["source", "label", "variant_id", "svtype", "size_bucket", "af_bucket", "genomic_context"],
        value_vars=metric_columns,
        var_name="metric",
        value_name="value",
    ).dropna(subset=["value"])
    summary = long_frame.groupby(["source", "metric", "svtype", "size_bucket", "af_bucket", "genomic_context"], dropna=False).agg(
        n_sites=("variant_id", "count"),
        mean_value=("value", "mean"),
        median_value=("value", "median"),
    ).reset_index()
    pivoted = summary.pivot_table(
        index=["metric", "svtype", "size_bucket", "af_bucket", "genomic_context"],
        columns="source",
        values=["n_sites", "mean_value", "median_value"],
        aggfunc="first",
    )
    if pivoted.empty:
        return pd.DataFrame(columns=["metric", "svtype", "size_bucket", "af_bucket", "genomic_context", "n_sites_a", "mean_value_a", "median_value_a", "n_sites_b", "mean_value_b", "median_value_b"])
    pivoted.columns = [f"{value}_{source}" for value, source in pivoted.columns]
    result = pivoted.reset_index()
    for column in ["n_sites_a", "mean_value_a", "median_value_a", "n_sites_b", "mean_value_b", "median_value_b"]:
        if column not in result.columns:
            result[column] = np.nan
    return result[["metric", "svtype", "size_bucket", "af_bucket", "genomic_context", "n_sites_a", "mean_value_a", "median_value_a", "n_sites_b", "mean_value_b", "median_value_b"]]


def _plot_metric_panels(metrics: pd.DataFrame, label_a: str, label_b: str, output_path: Path) -> None:
    metric_names = ["genotype_concordance", "non_ref_genotype_concordance", "var_ppv", "var_sensitivity"]
    titles = {
        "genotype_concordance": "Genotype concordance",
        "non_ref_genotype_concordance": "Non-ref genotype concordance",
        "var_ppv": "Variant PPV",
        "var_sensitivity": "Variant sensitivity",
    }
    fig, axes = plt.subplots(2, 2, figsize=(12, 10), squeeze=False, constrained_layout=True)
    flat_axes = axes.flatten()
    for ax, metric in zip(flat_axes, metric_names):
        if metrics.empty or metric not in metrics.columns:
            ax.text(0.5, 0.5, "No concordance metrics", ha="center", va="center")
            ax.set_axis_off()
            continue
        grouped_a = metrics.loc[metrics["source"] == "a"].groupby("svtype", dropna=False)[metric].median()
        grouped_b = metrics.loc[metrics["source"] == "b"].groupby("svtype", dropna=False)[metric].median()
        svtypes = sorted(set(grouped_a.index.astype(str)).union(set(grouped_b.index.astype(str))))
        if not svtypes:
            ax.text(0.5, 0.5, "No concordance metrics", ha="center", va="center")
            ax.set_axis_off()
            continue
        x_positions = np.arange(len(svtypes), dtype=float)
        values_a = [float(grouped_a.get(svtype, np.nan)) for svtype in svtypes]
        values_b = [float(grouped_b.get(svtype, np.nan)) for svtype in svtypes]
        width = 0.36
        ax.bar(x_positions - width / 2, values_a, width=width, color="#4C72B0", label=label_a, edgecolor="#E5E5E5", linewidth=0.8)
        ax.bar(x_positions + width / 2, values_b, width=width, color="#DD8452", label=label_b, edgecolor="#E5E5E5", linewidth=0.8)
        ax.set_xticks(x_positions)
        ax.set_xticklabels(svtypes, rotation=30)
        ax.set_ylim(0.0, 1.0)
        ax.set_ylabel("Median value")
        ax.set_xlabel("SV type")
        ax.set_title(titles[metric])
    handles, labels = flat_axes[0].get_legend_handles_labels()
    if handles:
        flat_axes[0].legend(handles, labels, fontsize=8)
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
        metrics_a = _extract_concordance_rows(config.vcf_a_path, data.sites_a, data.label_a, "a")
        metrics_b = _extract_concordance_rows(config.vcf_b_path, data.sites_b, data.label_b, "b")
        metrics = pd.concat([metrics_a, metrics_b], ignore_index=True) if not metrics_a.empty or not metrics_b.empty else pd.DataFrame()
        summary = relabel_vcf_columns(summarize_concordance_metrics(metrics), data.label_a, data.label_b)
        write_tsv_gz(summary, tables_dir / "concordance_metrics.tsv")
        _plot_metric_panels(metrics, data.label_a, data.label_b, output_dir / "concordance_metrics.by_type.png")
