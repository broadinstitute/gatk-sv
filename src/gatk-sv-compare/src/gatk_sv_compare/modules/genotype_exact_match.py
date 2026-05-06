"""Shared-sample genotype exact-match diagnostics for matched sites."""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path
from typing import Dict, Optional

import numpy as np
import pandas as pd
import pysam

from ..aggregate import AggregatedData
from ..config import AnalysisConfig
from ..dimensions import complete_genomic_context_buckets, explode_algorithm_buckets
from .base import AnalysisModule, relabel_vcf_columns, write_tsv_gz


_GT_LABELS = {0: "hom_ref", 1: "het", 2: "hom_alt"}


def _as_object_vector(values: list[object]) -> np.ndarray:
    vector = np.empty(len(values), dtype=object)
    vector[:] = values
    return vector


def _record_identifier(record: pysam.VariantRecord) -> str:
    return record.id or f"{record.contig}:{record.pos}"


def _iter_records_for_contig(vcf: pysam.VariantFile, contig: str):
    try:
        yield from vcf.fetch(contig)
    except ValueError:
        for record in vcf:
            if record.contig == contig:
                yield record


def _canonical_repcn(value: object) -> Optional[tuple[object, ...]]:
    if value in (None, "."):
        return None
    if isinstance(value, (tuple, list)):
        raw_values = list(value)
    else:
        text = str(value)
        separator = "/" if "/" in text else ","
        raw_values = text.split(separator)
    if not raw_values or any(raw_value in (None, ".", "") for raw_value in raw_values):
        return None
    try:
        return tuple(sorted(int(raw_value) for raw_value in raw_values))
    except (TypeError, ValueError):
        return tuple(sorted(str(raw_value) for raw_value in raw_values))


def _gt_code(sample: pysam.libcbcf.VariantRecordSample, svtype: str) -> Optional[object]:
    if svtype == "STR":
        return _canonical_repcn(sample.get("REPCN"))
    if svtype == "CNV":
        cn = sample.get("CN")
        ecn = sample.get("ECN")
        if cn in (None, ".") or ecn in (None, "."):
            return None
        return int(cn) - int(ecn)
    gt = sample.get("GT")
    if not gt or any(allele is None for allele in gt):
        return None
    return int(sum(1 for allele in gt if allele and allele > 0))


def _extract_genotype_vectors(
    vcf_path: Path,
    target_ids_by_contig: Dict[str, set[str]],
    sample_indices: np.ndarray,
    svtype_by_vid: Dict[str, str],
) -> Dict[str, np.ndarray]:
    vectors: Dict[str, np.ndarray] = {}
    with pysam.VariantFile(str(vcf_path)) as vcf:
        for contig, target_ids in target_ids_by_contig.items():
            for record in _iter_records_for_contig(vcf, contig):
                variant_id = _record_identifier(record)
                if variant_id not in target_ids:
                    continue
                svtype = svtype_by_vid.get(variant_id, "UNKNOWN")
                sample_values = list(record.samples.values())
                codes = [_gt_code(sample_values[int(index)], svtype) for index in sample_indices]
                vectors[variant_id] = _as_object_vector(codes)
    return vectors


def _valid_genotype_mask(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    return np.asarray([left_value is not None and right_value is not None for left_value, right_value in zip(left, right)], dtype=bool)


def _is_numeric_genotype_vector(values: np.ndarray) -> bool:
    return all(isinstance(value, (int, float, np.integer, np.floating)) for value in values)


def build_exact_match_table(data: AggregatedData, config: AnalysisConfig) -> pd.DataFrame:
    columns = [
        "variant_id_a",
        "variant_id_b",
        "svtype_a",
        "size_bucket_a",
        "af_bucket_a",
        "genomic_context_a",
        "algorithms_a",
        "evidence_bucket_a",
        "svtype_b",
        "size_bucket_b",
        "af_bucket_b",
        "genomic_context_b",
        "algorithms_b",
        "evidence_bucket_b",
        "n_compared",
        "exact_match_rate",
        "homref_to_het_rate",
        "het_to_homref_rate",
        "homref_to_homalt_rate",
        "homalt_to_homref_rate",
        "het_to_homalt_rate",
        "homalt_to_het_rate",
    ]
    if data.matched_pairs.empty or not data.shared_samples or data.sample_indices_a is None or data.sample_indices_b is None:
        return pd.DataFrame(columns=columns)

    target_ids_a: Dict[str, set[str]] = defaultdict(set)
    target_ids_b: Dict[str, set[str]] = defaultdict(set)
    svtype_by_vid: Dict[str, str] = {}
    site_meta_a = data.sites_a.set_index("variant_id")[["svtype", "size_bucket", "af_bucket", "genomic_context", "algorithms", "evidence_bucket"]]
    site_meta_b = data.sites_b.set_index("variant_id")[["svtype", "size_bucket", "af_bucket", "genomic_context", "algorithms", "evidence_bucket"]]
    for row in data.matched_pairs.itertuples(index=False):
        target_ids_a[str(row.contig_a)].add(str(row.variant_id_a))
        target_ids_b[str(row.contig_b)].add(str(row.variant_id_b))
        svtype_by_vid[str(row.variant_id_a)] = str(row.svtype_a)
        svtype_by_vid[str(row.variant_id_b)] = str(row.svtype_b)

    vectors_a = _extract_genotype_vectors(config.vcf_a_path, target_ids_a, data.sample_indices_a, svtype_by_vid)
    vectors_b = _extract_genotype_vectors(config.vcf_b_path, target_ids_b, data.sample_indices_b, svtype_by_vid)

    rows = []
    for row in data.matched_pairs.itertuples(index=False):
        vector_a = vectors_a.get(str(row.variant_id_a))
        vector_b = vectors_b.get(str(row.variant_id_b))
        if vector_a is None or vector_b is None:
            continue
        valid = _valid_genotype_mask(vector_a, vector_b)
        n_compared = int(valid.sum())
        if n_compared == 0:
            continue
        left = vector_a[valid]
        right = vector_b[valid]
        exact_match_rate = float(np.equal(left, right).mean())
        has_numeric_codes = _is_numeric_genotype_vector(left) and _is_numeric_genotype_vector(right)
        if has_numeric_codes:
            left_numeric = left.astype(int)
            right_numeric = right.astype(int)
        else:
            left_numeric = np.asarray([], dtype=int)
            right_numeric = np.asarray([], dtype=int)
        meta_a = site_meta_a.loc[str(row.variant_id_a)]
        meta_b = site_meta_b.loc[str(row.variant_id_b)]
        rows.append(
            {
                "variant_id_a": row.variant_id_a,
                "variant_id_b": row.variant_id_b,
                "svtype_a": meta_a["svtype"],
                "size_bucket_a": meta_a["size_bucket"],
                "af_bucket_a": meta_a["af_bucket"],
                "genomic_context_a": meta_a["genomic_context"],
                "algorithms_a": meta_a["algorithms"],
                "evidence_bucket_a": meta_a["evidence_bucket"],
                "svtype_b": meta_b["svtype"],
                "size_bucket_b": meta_b["size_bucket"],
                "af_bucket_b": meta_b["af_bucket"],
                "genomic_context_b": meta_b["genomic_context"],
                "algorithms_b": meta_b["algorithms"],
                "evidence_bucket_b": meta_b["evidence_bucket"],
                "n_compared": n_compared,
                "exact_match_rate": exact_match_rate,
                "homref_to_het_rate": float(((left_numeric == 0) & (right_numeric == 1)).mean()) if has_numeric_codes else 0.0,
                "het_to_homref_rate": float(((left_numeric == 1) & (right_numeric == 0)).mean()) if has_numeric_codes else 0.0,
                "homref_to_homalt_rate": float(((left_numeric == 0) & (right_numeric == 2)).mean()) if has_numeric_codes else 0.0,
                "homalt_to_homref_rate": float(((left_numeric == 2) & (right_numeric == 0)).mean()) if has_numeric_codes else 0.0,
                "het_to_homalt_rate": float(((left_numeric == 1) & (right_numeric == 2)).mean()) if has_numeric_codes else 0.0,
                "homalt_to_het_rate": float(((left_numeric == 2) & (right_numeric == 1)).mean()) if has_numeric_codes else 0.0,
            }
        )
    return pd.DataFrame(rows, columns=columns)


def summarize_exact_match(per_site: pd.DataFrame) -> pd.DataFrame:
    columns = [
        "svtype_a",
        "size_bucket_a",
        "af_bucket_a",
        "genomic_context_a",
        "evidence_bucket_a",
        "algorithm_a",
        "svtype_b",
        "size_bucket_b",
        "af_bucket_b",
        "genomic_context_b",
        "evidence_bucket_b",
        "algorithm_b",
        "n_sites",
        "mean_exact_match_rate",
        "median_exact_match_rate",
        "mean_homref_to_het_rate",
        "mean_het_to_homref_rate",
        "mean_homref_to_homalt_rate",
        "mean_homalt_to_homref_rate",
        "mean_het_to_homalt_rate",
        "mean_homalt_to_het_rate",
    ]
    if per_site.empty:
        return pd.DataFrame(columns=columns)
    summary_input = explode_algorithm_buckets(per_site, algorithms_column="algorithms_a", bucket_column="algorithm_a")
    summary_input = explode_algorithm_buckets(summary_input, algorithms_column="algorithms_b", bucket_column="algorithm_b")
    summary = summary_input.groupby(["svtype_a", "size_bucket_a", "af_bucket_a", "genomic_context_a", "evidence_bucket_a", "algorithm_a", "svtype_b", "size_bucket_b", "af_bucket_b", "genomic_context_b", "evidence_bucket_b", "algorithm_b"], dropna=False).agg(
        n_sites=("variant_id_a", "count"),
        mean_exact_match_rate=("exact_match_rate", "mean"),
        median_exact_match_rate=("exact_match_rate", "median"),
        mean_homref_to_het_rate=("homref_to_het_rate", "mean"),
        mean_het_to_homref_rate=("het_to_homref_rate", "mean"),
        mean_homref_to_homalt_rate=("homref_to_homalt_rate", "mean"),
        mean_homalt_to_homref_rate=("homalt_to_homref_rate", "mean"),
        mean_het_to_homalt_rate=("het_to_homalt_rate", "mean"),
        mean_homalt_to_het_rate=("homalt_to_het_rate", "mean"),
    ).reset_index()
    summary = complete_genomic_context_buckets(
        summary,
        ["svtype_a", "size_bucket_a", "af_bucket_a", "genomic_context_a", "evidence_bucket_a", "algorithm_a", "svtype_b", "size_bucket_b", "af_bucket_b", "genomic_context_b", "evidence_bucket_b", "algorithm_b"],
        fill_values={"n_sites": 0},
    )
    summary["n_sites"] = summary["n_sites"].astype(int)
    return summary[columns]


class GenotypeExactMatchModule(AnalysisModule):
    @property
    def name(self) -> str:
        return "genotype_exact_match"

    @property
    def requires_shared_samples(self) -> bool:
        return True

    @property
    def requires_genotype_pass(self) -> bool:
        return True

    def run(self, data: AggregatedData, config: AnalysisConfig) -> None:
        output_dir = self.output_dir(config)
        tables_dir = output_dir / "tables"
        tables_dir.mkdir(parents=True, exist_ok=True)

        per_site = build_exact_match_table(data, config)
        summary = summarize_exact_match(per_site)
        per_site_output = relabel_vcf_columns(per_site, data.label_a, data.label_b)
        summary_output = relabel_vcf_columns(summary, data.label_a, data.label_b)
        if config.enable_site_match_table:
            write_tsv_gz(per_site_output, tables_dir / "genotype_match_rates.per_site.tsv")
        write_tsv_gz(summary_output, tables_dir / "genotype_match_rates.tsv")
