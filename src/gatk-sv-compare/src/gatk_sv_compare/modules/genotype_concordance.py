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
from ..dimensions import af_bucket_sort_key, algorithm_sort_key, complete_genomic_context_buckets, evidence_bucket_sort_key, explode_algorithm_buckets, genomic_context_sort_key, ordered_algorithms, ordered_contexts, ordered_plot_af_buckets, ordered_plot_evidence_buckets, ordered_plot_size_buckets, ordered_svtypes, size_bucket_sort_key, svtype_sort_key
from ..plot_utils import double_column_figsize, plot_heatmap_annotated, save_figure
from ..vcf_format import safe_info_get
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
    value = safe_info_get(record, key)
    if isinstance(value, tuple):
        return value[0] if value else None
    return value


def _extract_concordance_rows(vcf_path: Path, sites: pd.DataFrame, label: str, source: str) -> pd.DataFrame:
    target_ids_by_contig: Dict[str, set[str]] = defaultdict(set)
    for row in sites[["variant_id", "contig"]].drop_duplicates().itertuples(index=False):
        target_ids_by_contig[str(row.contig)].add(str(row.variant_id))

    site_meta = sites.set_index("variant_id")[["svtype", "size_bucket", "af_bucket", "genomic_context", "algorithms", "evidence_bucket"]]
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
                    "algorithms": meta["algorithms"],
                    "evidence_bucket": meta["evidence_bucket"],
                    **{key: float(value) for key, value in metrics.items()},
                })
    return pd.DataFrame(rows)


def summarize_concordance_metrics(metrics: pd.DataFrame) -> pd.DataFrame:
    if metrics.empty:
        return pd.DataFrame(columns=["metric", "svtype", "size_bucket", "af_bucket", "genomic_context", "evidence_bucket", "algorithm", "n_sites_a", "mean_value_a", "median_value_a", "n_sites_b", "mean_value_b", "median_value_b"])
    metric_columns = [column for column in metrics.columns if column in {field.lower() for field in _CONCORDANCE_FIELDS}]
    metrics = explode_algorithm_buckets(metrics)
    cnv_mask = metrics["svtype"].astype(str) == "CNV"
    if "cnv_concordance" in metrics.columns:
        metrics.loc[cnv_mask, "genotype_concordance"] = metrics.loc[cnv_mask, "cnv_concordance"]
    long_frame = metrics.melt(
        id_vars=["source", "label", "variant_id", "svtype", "size_bucket", "af_bucket", "genomic_context", "evidence_bucket", "algorithm"],
        value_vars=metric_columns,
        var_name="metric",
        value_name="value",
    ).dropna(subset=["value"])
    cnv_metric_mask = (long_frame["svtype"].astype(str) == "CNV") & long_frame["metric"].isin(["non_ref_genotype_concordance", "var_ppv", "var_sensitivity"])
    non_cnv_metric_mask = (long_frame["svtype"].astype(str) != "CNV") & (long_frame["metric"] == "cnv_concordance")
    exclude_mask = cnv_metric_mask | non_cnv_metric_mask
    long_frame = long_frame.loc[~exclude_mask]
    summary = long_frame.groupby(["source", "metric", "svtype", "size_bucket", "af_bucket", "genomic_context", "evidence_bucket", "algorithm"], dropna=False).agg(
        n_sites=("variant_id", "count"),
        mean_value=("value", "mean"),
        median_value=("value", "median"),
    ).reset_index()
    summary = complete_genomic_context_buckets(
        summary,
        ["source", "metric", "svtype", "size_bucket", "af_bucket", "genomic_context", "evidence_bucket", "algorithm"],
        fill_values={"n_sites": 0},
    )
    summary["n_sites"] = summary["n_sites"].astype(int)
    pivoted = summary.pivot_table(
        index=["metric", "svtype", "size_bucket", "af_bucket", "genomic_context", "evidence_bucket", "algorithm"],
        columns="source",
        values=["n_sites", "mean_value", "median_value"],
        aggfunc="first",
    )
    if pivoted.empty:
        return pd.DataFrame(columns=["metric", "svtype", "size_bucket", "af_bucket", "genomic_context", "evidence_bucket", "algorithm", "n_sites_a", "mean_value_a", "median_value_a", "n_sites_b", "mean_value_b", "median_value_b"])
    pivoted.columns = [f"{value}_{source}" for value, source in pivoted.columns]
    result = pivoted.reset_index()
    for column in ["n_sites_a", "mean_value_a", "median_value_a", "n_sites_b", "mean_value_b", "median_value_b"]:
        if column not in result.columns:
            result[column] = np.nan
    result = result[["metric", "svtype", "size_bucket", "af_bucket", "genomic_context", "evidence_bucket", "algorithm", "n_sites_a", "mean_value_a", "median_value_a", "n_sites_b", "mean_value_b", "median_value_b"]]
    return result.sort_values(
        by=["metric", "svtype", "size_bucket", "af_bucket", "genomic_context", "evidence_bucket", "algorithm"],
        key=lambda series: series.map(
            lambda value: (
                svtype_sort_key(value) if series.name == "svtype" else
                size_bucket_sort_key(value) if series.name == "size_bucket" else
                af_bucket_sort_key(value) if series.name == "af_bucket" else
                genomic_context_sort_key(value) if series.name == "genomic_context" else
                evidence_bucket_sort_key(value) if series.name == "evidence_bucket" else
                algorithm_sort_key(value) if series.name == "algorithm" else
                str(value)
            )
        ),
    ).reset_index(drop=True)


def _concordance_series(metrics: pd.DataFrame, metric: str) -> pd.DataFrame:
    if metrics.empty or metric not in metrics.columns:
        return pd.DataFrame(columns=["source", "svtype", "size_bucket", "af_bucket", "genomic_context", "evidence_bucket", "algorithms", metric])
    metric_frame = metrics.copy()
    if metric == "genotype_concordance" and "cnv_concordance" in metric_frame.columns:
        cnv_mask = metric_frame["svtype"].astype(str) == "CNV"
        metric_frame.loc[cnv_mask, metric] = metric_frame.loc[cnv_mask, "cnv_concordance"]
    if metric != "genotype_concordance":
        metric_frame = metric_frame.loc[metric_frame["svtype"].astype(str) != "CNV"]
    return metric_frame[["source", "svtype", "size_bucket", "af_bucket", "genomic_context", "evidence_bucket", "algorithms", metric]].dropna(subset=[metric])


def _ordered_heatmap_rows(row_field: str, values: pd.Index) -> list[str]:
    if row_field == "af_bucket":
        return ordered_plot_af_buckets(values)
    if row_field == "algorithm":
        return ordered_algorithms(values)
    if row_field == "evidence_bucket":
        return ordered_plot_evidence_buckets(values)
    if row_field == "size_bucket":
        return ordered_plot_size_buckets(values)
    if row_field == "genomic_context":
        return ordered_contexts(values)
    return [str(value) for value in values]


def _metric_title(metric: str) -> str:
    return {
        "genotype_concordance": "Genotype concordance",
        "non_ref_genotype_concordance": "Non-ref genotype concordance",
        "var_ppv": "Variant PPV",
        "var_sensitivity": "Variant sensitivity",
    }[metric]


def _metric_colorbar_label(metric: str) -> str:
    return f"{_metric_title(metric)} (mean)"


def _metric_stem(metric: str) -> str:
    return {
        "genotype_concordance": "genotype_concordance",
        "non_ref_genotype_concordance": "non_ref_genotype_concordance",
        "var_ppv": "variant_ppv",
        "var_sensitivity": "variant_sensitivity",
    }[metric]


def _plot_metric_heatmap(
    metrics: pd.DataFrame,
    metric: str,
    row_field: str,
    output_path: Path,
    label_a: str,
    label_b: str,
) -> None:
    metric_frame = _concordance_series(metrics, metric)
    if row_field == "algorithm":
        metric_frame = explode_algorithm_buckets(metric_frame)
    fig, axes = plt.subplots(1, 2, figsize=double_column_figsize(3.2), squeeze=False)
    flat_axes = axes.flatten()
    image = None
    for ax, source, label in zip(flat_axes, ["a", "b"], [label_a, label_b]):
        source_frame = metric_frame.loc[metric_frame["source"] == source]
        if source_frame.empty:
            ax.text(0.5, 0.5, "No concordance metrics", ha="center", va="center")
            ax.set_axis_off()
            continue
        grouped = source_frame.groupby([row_field, "svtype"], dropna=False).agg(
            value=(metric, "mean"),
            n_points=(metric, "count"),
        ).reset_index()
        matrix = grouped.pivot(index=row_field, columns="svtype", values="value")
        count_matrix = grouped.pivot(index=row_field, columns="svtype", values="n_points")
        matrix = matrix.reindex(index=_ordered_heatmap_rows(row_field, matrix.index), columns=ordered_svtypes(matrix.columns))
        count_matrix = count_matrix.reindex(index=matrix.index, columns=matrix.columns)
        image = plot_heatmap_annotated(
            ax,
            matrix.values,
            list(matrix.index),
            list(matrix.columns),
            fmt="{value:.2f}",
            value_range=(0.0, 1.0),
            count_matrix=count_matrix.fillna(0.0).values,
        )
        ax.set_title(label)
        ax.set_xlabel("SV type")
        if ax is flat_axes[0]:
            ax.set_ylabel(row_field.replace("_", " "))
    if image is not None:
        colorbar = fig.colorbar(image, ax=flat_axes.tolist())
        colorbar.set_label(_metric_colorbar_label(metric))
    fig.suptitle(f"{_metric_title(metric)}: {row_field.replace('_', ' ')} × SV type", y=0.98)
    save_figure(fig, output_path)


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
        for metric in ["genotype_concordance", "non_ref_genotype_concordance", "var_ppv", "var_sensitivity"]:
            stem = _metric_stem(metric)
            _plot_metric_heatmap(metrics, metric, "af_bucket", output_dir / f"{stem}.af_bucket_x_svtype.png", data.label_a, data.label_b)
            _plot_metric_heatmap(metrics, metric, "size_bucket", output_dir / f"{stem}.size_bucket_x_svtype.png", data.label_a, data.label_b)
            _plot_metric_heatmap(metrics, metric, "genomic_context", output_dir / f"{stem}.genomic_context_x_svtype.png", data.label_a, data.label_b)
            _plot_metric_heatmap(metrics, metric, "evidence_bucket", output_dir / f"{stem}.evidence_bucket_x_svtype.png", data.label_a, data.label_b)
            _plot_metric_heatmap(metrics, metric, "algorithm", output_dir / f"{stem}.algorithm_x_svtype.png", data.label_a, data.label_b)
