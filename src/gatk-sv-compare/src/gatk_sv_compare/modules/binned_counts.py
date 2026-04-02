"""Cross-tabulated count tables for site-level dimensions."""

from __future__ import annotations

import pandas as pd

from ..aggregate import AggregatedData
from ..config import AnalysisConfig
from ..dimensions import af_bucket_sort_key, complete_genomic_context_buckets, genomic_context_sort_key, size_bucket_sort_key, svtype_sort_key
from .base import AnalysisModule, column_safe_label, write_tsv_gz, matched_site_mask


def _filtered_sites(sites: pd.DataFrame, pass_only: bool) -> pd.DataFrame:
    if not pass_only:
        return sites.copy()
    return sites.loc[sites["in_filtered_pass_view"]].copy()


def summarize_binned_counts(sites: pd.DataFrame, pass_only: bool = False) -> pd.DataFrame:
    filtered = _filtered_sites(sites, pass_only)
    columns = ["svtype", "size_bucket", "af_bucket", "genomic_context", "n_variants", "n_matched", "n_unmatched"]
    if filtered.empty:
        return pd.DataFrame(columns=columns)
    matched = matched_site_mask(filtered)
    filtered = filtered.assign(_matched=matched.astype(int), _unmatched=(~matched).astype(int))
    grouped = filtered.groupby(["svtype", "size_bucket", "af_bucket", "genomic_context"], dropna=False)
    summary = grouped.agg(
        n_variants=("variant_id", "count"),
        n_matched=("_matched", "sum"),
        n_unmatched=("_unmatched", "sum"),
    ).reset_index()
    summary = complete_genomic_context_buckets(
        summary,
        ["svtype", "size_bucket", "af_bucket", "genomic_context"],
        fill_values={"n_variants": 0, "n_matched": 0, "n_unmatched": 0},
    )
    summary[["n_variants", "n_matched", "n_unmatched"]] = summary[["n_variants", "n_matched", "n_unmatched"]].astype(int)
    return summary[columns].sort_values(
        by=["svtype", "size_bucket", "af_bucket", "genomic_context"],
        key=lambda series: series.map(
            lambda value: (
                svtype_sort_key(value) if series.name == "svtype" else
                size_bucket_sort_key(value) if series.name == "size_bucket" else
                af_bucket_sort_key(value) if series.name == "af_bucket" else
                genomic_context_sort_key(value)
            )
        ),
    ).reset_index(drop=True)


def build_combined_binned_counts(sites_a: pd.DataFrame, sites_b: pd.DataFrame, label_a: str, label_b: str, pass_only: bool = False) -> pd.DataFrame:
    token_a = column_safe_label(label_a)
    token_b = column_safe_label(label_b)
    counts_a = summarize_binned_counts(sites_a, pass_only=pass_only).rename(
        columns={
            "n_variants": f"n_variants_{token_a}",
            "n_matched": f"n_matched_{token_a}",
            "n_unmatched": f"n_unmatched_{token_a}",
        }
    )
    counts_b = summarize_binned_counts(sites_b, pass_only=pass_only).rename(
        columns={
            "n_variants": f"n_variants_{token_b}",
            "n_matched": f"n_matched_{token_b}",
            "n_unmatched": f"n_unmatched_{token_b}",
        }
    )
    merged = counts_a.merge(counts_b, on=["svtype", "size_bucket", "af_bucket", "genomic_context"], how="outer")
    numeric_columns = [column for column in merged.columns if column.startswith("n_")]
    merged[numeric_columns] = merged[numeric_columns].fillna(0).astype(int)
    return merged.sort_values(
        by=["svtype", "size_bucket", "af_bucket", "genomic_context"],
        key=lambda series: series.map(
            lambda value: (
                svtype_sort_key(value) if series.name == "svtype" else
                size_bucket_sort_key(value) if series.name == "size_bucket" else
                af_bucket_sort_key(value) if series.name == "af_bucket" else
                genomic_context_sort_key(value)
            )
        ),
    ).reset_index(drop=True)


class BinnedCountsModule(AnalysisModule):
    @property
    def name(self) -> str:
        return "binned_counts"

    def run(self, data: AggregatedData, config: AnalysisConfig) -> None:
        output_dir = self.output_dir(config)
        counts = build_combined_binned_counts(data.sites_a, data.sites_b, data.label_a, data.label_b, pass_only=config.pass_only)
        write_tsv_gz(counts, output_dir / "counts.tsv")
        counts.to_parquet(output_dir / "counts.parquet", index=False)
