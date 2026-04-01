"""Validate command implementation for gatk-sv-compare."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import List

import pysam

from .config import ValidateConfig
from .vcf_format import FormatIssue, PipelineStage, check_record, detect_pipeline_stage, has_precomputed_counts


@dataclass
class ValidationSummary:
    """Summary of a validation run."""

    vcf_path: Path
    stage: PipelineStage
    issues: List[FormatIssue]
    record_count: int
    precomputed_counts_available: bool

    @property
    def has_errors(self) -> bool:
        return any(issue.severity == "ERROR" for issue in self.issues)

    def severity_counts(self) -> Counter[str]:
        return Counter(issue.severity for issue in self.issues)


def validate_vcf(config: ValidateConfig) -> ValidationSummary:
    """Validate a VCF and return structured results."""
    with pysam.VariantFile(str(config.vcf_path)) as vcf:
        stage = detect_pipeline_stage(vcf)
        precomputed = has_precomputed_counts(vcf)

    issues: List[FormatIssue] = []
    record_count = 0
    with pysam.VariantFile(str(config.vcf_path)) as vcf:
        contig_lengths = {name: int(contig.length) for name, contig in vcf.header.contigs.items() if contig.length is not None}
        for record in vcf:
            record_count += 1
            issues.extend(check_record(record, contig_length=contig_lengths.get(record.contig)))
            if record_count == 1 and not precomputed:
                issues.append(
                    FormatIssue(
                        "MISSING_PRECOMPUTED_COUNTS",
                        "INFO",
                        "",
                        "Pre-computed genotype count fields are absent; validation will allow GT-based fallback",
                    )
                )
    return ValidationSummary(
        vcf_path=config.vcf_path,
        stage=stage,
        issues=issues,
        record_count=record_count,
        precomputed_counts_available=precomputed,
    )


def render_summary(summary: ValidationSummary) -> str:
    """Render a human-readable validation report."""
    counts = summary.severity_counts()
    lines = [
        f"VCF: {summary.vcf_path}",
        f"Stage: {summary.stage.value}",
        f"Records: {summary.record_count}",
        f"Issues: errors={counts.get('ERROR', 0)}, warnings={counts.get('WARN', 0)}, info={counts.get('INFO', 0)}",
    ]
    if summary.issues:
        lines.append("Details:")
        for issue in summary.issues:
            where = f" [{issue.record_id}]" if issue.record_id else ""
            lines.append(f"- {issue.severity} {issue.check_id}{where}: {issue.message}")
    else:
        lines.append("No issues found.")
    return "\n".join(lines)


def validate_and_render(vcf_path: Path) -> tuple[ValidationSummary, str]:
    """Convenience wrapper used by the CLI."""
    summary = validate_vcf(ValidateConfig(vcf_path=vcf_path))
    return summary, render_summary(summary)
