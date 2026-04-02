"""Base analysis module interface."""

from __future__ import annotations

from abc import ABC, abstractmethod
from pathlib import Path

from ..aggregate import AggregatedData
from ..config import AnalysisConfig


class AnalysisModule(ABC):
    """Base class for all analysis modules."""

    @property
    @abstractmethod
    def name(self) -> str:
        """Directory name for outputs."""

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
