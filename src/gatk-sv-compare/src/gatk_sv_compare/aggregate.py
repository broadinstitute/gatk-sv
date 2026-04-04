"""Phase 3 aggregation for compact site tables and matched pairs."""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import pysam

from .config import AnalysisConfig
from .dimensions import categorize_variant, is_filtered_pass
from .vcf_reader import SiteRecord, iter_contig

_SITE_TABLE_COLUMNS = [
    "variant_id",
    "contig",
    "pos",
    "end",
    "svtype",
    "svlen",
    "alt_allele",
    "filter",
    "af",
    "ac",
    "an",
    "n_bi_genos",
    "n_hom_ref",
    "n_het",
    "n_hom_alt",
    "n_no_call",
    "cn_nonref_freq",
    "status",
    "truth_vid",
    "genomic_context",
    "algorithms",
    "evidence_bucket",
    "size_bucket",
    "af_bucket",
    "in_filtered_pass_view",
]


@dataclass
class AggregatedData:
    """Lightweight analysis context returned by the aggregation stage."""

    sites_a: pd.DataFrame
    sites_b: pd.DataFrame
    matched_pairs: pd.DataFrame
    site_table_dir: Path
    genotype_artifacts_dir: Optional[Path]
    sample_names_a: List[str]
    sample_names_b: List[str]
    shared_samples: List[str]
    sample_indices_a: Optional[np.ndarray]
    sample_indices_b: Optional[np.ndarray]
    label_a: str
    label_b: str


def _empty_site_table() -> pd.DataFrame:
    return pd.DataFrame(columns=_SITE_TABLE_COLUMNS)


def _site_row(site: SiteRecord) -> dict:
    category = categorize_variant(site.to_dict())
    filter_values = frozenset(part for part in site.filter.split(";") if part and part != ".")
    return {
        "variant_id": site.variant_id,
        "contig": site.contig,
        "pos": site.pos,
        "end": site.end,
        "svtype": category.svtype,
        "svlen": site.svlen,
        "alt_allele": site.alt_allele,
        "filter": site.filter,
        "af": site.af,
        "ac": site.ac,
        "an": site.an,
        "n_bi_genos": site.n_bi_genos,
        "n_hom_ref": site.n_hom_ref,
        "n_het": site.n_het,
        "n_hom_alt": site.n_hom_alt,
        "n_no_call": site.n_no_call,
        "cn_nonref_freq": site.cn_nonref_freq,
        "status": site.status,
        "truth_vid": site.truth_vid,
        "genomic_context": category.genomic_context,
        "algorithms": ",".join(site.algorithms),
        "evidence_bucket": site.evidence_bucket,
        "size_bucket": category.size_bucket,
        "af_bucket": category.af_bucket,
        "in_filtered_pass_view": is_filtered_pass(filter_values),
    }


def _extract_contig(vcf_path: Path, contig: str, output_path: Path, context_overlap: float) -> Path:
    rows = [_site_row(site) for site in iter_contig(vcf_path, contig, context_overlap=context_overlap)]
    dataframe = pd.DataFrame(rows, columns=_SITE_TABLE_COLUMNS) if rows else _empty_site_table()
    dataframe.to_parquet(output_path, index=False)
    return output_path


def _sample_names(vcf_path: Path) -> List[str]:
    with pysam.VariantFile(str(vcf_path)) as vcf:
        return list(vcf.header.samples)


def _shared_samples(sample_names_a: Sequence[str], sample_names_b: Sequence[str]) -> Tuple[List[str], Optional[np.ndarray], Optional[np.ndarray]]:
    shared = [sample for sample in sample_names_a if sample in set(sample_names_b)]
    if not shared:
        return [], None, None
    index_b = {sample: idx for idx, sample in enumerate(sample_names_b)}
    indices_a = np.asarray([idx for idx, sample in enumerate(sample_names_a) if sample in set(shared)], dtype=int)
    indices_b = np.asarray([index_b[sample] for sample in shared], dtype=int)
    return shared, indices_a, indices_b


def _extract_site_tables(
    vcf_path: Path,
    contigs: Sequence[str],
    output_dir: Path,
    prefix: str,
    n_workers: int,
    context_overlap: float,
) -> List[Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    max_workers = max(1, n_workers)
    shard_paths = [output_dir / f"{prefix}.{contig}.parquet" for contig in contigs]
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        return list(
            executor.map(
                _extract_contig,
                [vcf_path] * len(contigs),
                list(contigs),
                shard_paths,
                [context_overlap] * len(contigs),
            )
        )


def _read_site_tables(paths: Sequence[Path]) -> pd.DataFrame:
    if not paths:
        return _empty_site_table()
    frames = [pd.read_parquet(path) for path in sorted(paths)]
    if not frames:
        return _empty_site_table()
    return pd.concat(frames, ignore_index=True) if frames else _empty_site_table()


def build_matched_pairs(sites_a: pd.DataFrame, sites_b: pd.DataFrame) -> pd.DataFrame:
    """Build the canonical matched-pair table from truth_vid links."""
    columns = [
        "pair_id",
        "variant_id_a",
        "variant_id_b",
        "contig_a",
        "contig_b",
        "svtype_a",
        "svtype_b",
        "af_a",
        "af_b",
        "status_a",
        "status_b",
    ]
    if sites_a.empty or sites_b.empty:
        return pd.DataFrame(columns=columns)

    matched_a = sites_a.loc[sites_a["truth_vid"].notna()].copy()
    if matched_a.empty:
        return pd.DataFrame(columns=columns)

    sites_b_lookup = sites_b[["variant_id", "contig", "svtype", "af", "status", "truth_vid"]].rename(
        columns={
            "variant_id": "variant_id_b",
            "contig": "contig_b",
            "svtype": "svtype_b",
            "af": "af_b",
            "status": "status_b",
            "truth_vid": "truth_vid_b",
        }
    )

    merged = matched_a.merge(
        sites_b_lookup,
        left_on="truth_vid",
        right_on="variant_id_b",
        how="inner",
    )
    if merged.empty:
        return pd.DataFrame(columns=columns)

    symmetric = merged.loc[(merged["truth_vid_b"].isna()) | (merged["truth_vid_b"] == merged["variant_id"])].copy()
    if symmetric.empty:
        return pd.DataFrame(columns=columns)

    symmetric["pair_id"] = symmetric["variant_id"] + "|" + symmetric["variant_id_b"]
    result = symmetric.rename(
        columns={
            "variant_id": "variant_id_a",
            "contig": "contig_a",
            "svtype": "svtype_a",
            "af": "af_a",
            "status": "status_a",
        }
    )
    return result[columns].drop_duplicates(subset=["pair_id"]).reset_index(drop=True)


def aggregate(config: AnalysisConfig) -> AggregatedData:
    """Parallel aggregation over both VCFs."""
    if config.vcf_a_path is None or config.vcf_b_path is None:
        raise ValueError("Both vcf_a_path and vcf_b_path are required for aggregation.")
    if not config.contigs:
        raise ValueError("At least one contig is required for aggregation.")

    aggregate_dir = config.output_dir / "aggregate"
    aggregate_dir.mkdir(parents=True, exist_ok=True)

    sample_names_a = _sample_names(config.vcf_a_path)
    sample_names_b = _sample_names(config.vcf_b_path)
    shared_samples, sample_indices_a, sample_indices_b = _shared_samples(sample_names_a, sample_names_b)

    site_paths_a = _extract_site_tables(
        config.vcf_a_path,
        config.contigs,
        aggregate_dir,
        "sites_a",
        config.n_workers,
        config.context_overlap,
    )
    site_paths_b = _extract_site_tables(
        config.vcf_b_path,
        config.contigs,
        aggregate_dir,
        "sites_b",
        config.n_workers,
        config.context_overlap,
    )

    sites_a = _read_site_tables(site_paths_a)
    sites_b = _read_site_tables(site_paths_b)
    if not sites_a.empty:
        sites_a.to_parquet(aggregate_dir / "sites_a.all.parquet", index=False)
    if not sites_b.empty:
        sites_b.to_parquet(aggregate_dir / "sites_b.all.parquet", index=False)

    matched_pairs = build_matched_pairs(sites_a, sites_b)
    matched_pairs.to_parquet(aggregate_dir / "matched_pairs.parquet", index=False)

    return AggregatedData(
        sites_a=sites_a,
        sites_b=sites_b,
        matched_pairs=matched_pairs,
        site_table_dir=aggregate_dir,
        genotype_artifacts_dir=None,
        sample_names_a=sample_names_a,
        sample_names_b=sample_names_b,
        shared_samples=shared_samples,
        sample_indices_a=sample_indices_a,
        sample_indices_b=sample_indices_b,
        label_a=config.vcf_a_label,
        label_b=config.vcf_b_label,
    )
