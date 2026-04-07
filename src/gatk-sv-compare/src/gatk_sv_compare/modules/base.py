"""Base analysis module interface."""

from __future__ import annotations

import re
from abc import ABC, abstractmethod
from pathlib import Path

import pandas as pd

from ..aggregate import AggregatedData
from ..config import AnalysisConfig
from ..dimensions import STATUS_MATCHED


class AnalysisModule(ABC):
    """Base class for all analysis modules."""

    @property
    @abstractmethod
    def name(self) -> str:
        """Directory name for outputs."""

    @property
    def requires_samples(self) -> bool:
        """True when the module reads sample-level FORMAT fields (GT, GQ, ECN)."""
        return False

    @property
    def requires_shared_samples(self) -> bool:
        return False

    @property
    def requires_ped_file(self) -> bool:
        return False

    @property
    def requires_gq(self) -> bool:
        return False

    @property
    def requires_concordance(self) -> bool:
        return False

    @property
    def requires_genotype_pass(self) -> bool:
        return False

    @abstractmethod
    def run(self, data: AggregatedData, config: AnalysisConfig) -> None:
        """Execute the analysis module."""

    def output_dir(self, config: AnalysisConfig) -> Path:
        directory = config.output_dir / self.name
        directory.mkdir(parents=True, exist_ok=True)
        return directory


def write_tsv_gz(dataframe: pd.DataFrame, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    resolved = path if str(path).endswith(".gz") else Path(f"{path}.gz")
    dataframe.to_csv(resolved, sep="\t", index=False, compression="gzip")
    return resolved


def matched_site_mask(sites: pd.DataFrame) -> pd.Series:
    if sites.empty:
        return pd.Series(dtype=bool)
    truth_vid_mask = sites["truth_vid"].notna() if "truth_vid" in sites.columns else pd.Series(False, index=sites.index)
    if "status" not in sites.columns:
        return truth_vid_mask.astype(bool)
    status_values = sites["status"].fillna("").astype(str).str.upper()
    return truth_vid_mask | (status_values == STATUS_MATCHED)


def column_safe_label(label: str) -> str:
    normalized = re.sub(r"\W+", "_", str(label)).strip("_")
    return normalized or "VCF"


def relabel_vcf_columns(dataframe: pd.DataFrame, label_a: str, label_b: str) -> pd.DataFrame:
    if dataframe.empty:
        return dataframe.copy()
    token_a = column_safe_label(label_a)
    token_b = column_safe_label(label_b)
    rename_map = {}
    for column in dataframe.columns:
        if column.endswith("_a"):
            rename_map[column] = f"{column[:-2]}_{token_a}"
        elif column.endswith("_b"):
            rename_map[column] = f"{column[:-2]}_{token_b}"
    return dataframe.rename(columns=rename_map)
