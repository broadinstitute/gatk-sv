"""Configuration models for gatk-sv-compare."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional


@dataclass
class AnalysisConfig:
    """Top-level runtime configuration for CLI commands."""

    vcf_a_path: Optional[Path] = None
    vcf_b_path: Optional[Path] = None
    vcf_a_label: str = "VCF A"
    vcf_b_label: str = "VCF B"
    output_dir: Path = Path(".")
    reference_dict: Optional[Path] = None
    contigs: List[str] = field(default_factory=list)
    contig_lengths: Dict[str, int] = field(default_factory=dict)
    n_workers: int = 1
    modules: Optional[List[str]] = None
    pass_only: bool = False
    per_chrom: bool = False
    ped_file: Optional[Path] = None
    seg_dup_track: Optional[Path] = None
    simple_repeat_track: Optional[Path] = None
    repeatmasker_track: Optional[Path] = None
    clustering_config: Optional[Path] = None
    stratification_config: Optional[Path] = None
    gatk_path: str = "gatk"
    java_options: str = "-Xmx4g"


@dataclass
class ValidateConfig:
    """Configuration for the validate subcommand."""

    vcf_path: Path
    max_records_for_gq_scan: int = 10000
