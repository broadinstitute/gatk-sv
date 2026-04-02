"""Cross-tabulated count tables for site-level dimensions."""

from __future__ import annotations

import pandas as pd

from ..aggregate import AggregatedData
from ..config import AnalysisConfig
from ..dimensions import STATUS_MATCHED, STATUS_UNMATCHED
from .base import AnalysisModule


def _filtered_sites(sites: pd.DataFrame, pass_only: bool) -> pd.DataFrame:
    if not pass_only:
        return sites.copy()
    return sites.loc[sites["in_filtered_pass_view"]].copy()


def summarize_binned_counts(sites: pd.DataFrame, pass_only: bool = False) -> pd.DataFrame:
    filtered = _filtered_sites(sites, pass_only)
    columns = ["svtype", "size_bucket", "af_bucket", "genomic_context", "n_variants", "n_matched", "n_unmatched"]
    if filtered.empty:
        return pd.DataFrame(columns=columns)
    grouped = filtered.groupby(["svtype", "size_bucket", "af_bucket", "genomic_context"], dropna=False)
    summary = grouped.agg(
        n_variants=("variant_id", "count"),
        n_matched=("status", lambda values: int((values == STATUS_MATCHED).sum())),
        n_unmatched=("status", lambda values: int((values == STATUS_UNMATCHED).sum())),
    ).reset_index()
    return summary[columns].sort_values(["svtype", "size_bucket", "af_bucket", "genomic_context"]).reset_index(drop=True)


class BinnedCountsModule(AnalysisModule):
    @property
    def name(self) -> str:
        return "binned_counts"

    def run(self, data: AggregatedData, config: AnalysisConfig) -> None:
        output_dir = self.output_dir(config)
        counts_a = summarize_binned_counts(data.sites_a, pass_only=config.pass_only)
        counts_b = summarize_binned_counts(data.sites_b, pass_only=config.pass_only)
        counts_a.to_csv(output_dir / "counts_a.tsv", sep="\t", index=False)
        counts_b.to_csv(output_dir / "counts_b.tsv", sep="\t", index=False)
        counts_a.to_parquet(output_dir / "counts_a.parquet", index=False)
        counts_b.to_parquet(output_dir / "counts_b.parquet", index=False)
