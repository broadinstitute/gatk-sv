"""Validate command implementation for gatk-sv-compare."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
import gzip
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pysam

from .config import ValidateConfig
from .vcf_format import FormatIssue, PipelineStage, check_header, check_record, detect_pipeline_stage, has_precomputed_counts

_FIXABLE_CHECK_IDS = {
    "SVLEN_SIGN",
    "SVLEN_NUMBER_DOT",
    "GQ_RANGE",
    "EMPTY_FILTER",
    "BREAKEND_NOTATION",
    "BND_END_MISMATCH",
    "INS_END_MISMATCH",
    "CPX_DDDUP_END",
    "MULTI_ALLELIC_NON_CNV",
    "CHR2_ON_NON_BND",
}


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

    def issue_type_counts(self) -> Counter[str]:
        return Counter(issue.check_id for issue in self.issues)


@dataclass
class ValidationFixResult:
    """Summary of a validate --fix run."""

    original_summary: ValidationSummary
    fixed_summary: Optional[ValidationSummary]
    out_path: Optional[Path]
    wrote_output: bool

    @property
    def has_errors(self) -> bool:
        if self.fixed_summary is not None:
            return self.fixed_summary.has_errors
        return self.original_summary.has_errors


def validate_vcf(config: ValidateConfig) -> ValidationSummary:
    """Validate a VCF and return structured results."""
    with pysam.VariantFile(str(config.vcf_path)) as vcf:
        stage = detect_pipeline_stage(vcf)
        precomputed = has_precomputed_counts(vcf)
        header_issues = check_header(vcf)

    issues: List[FormatIssue] = list(header_issues)
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


def _open_text(path: Path, mode: str):
    if path.suffix == ".gz":
        return gzip.open(path, mode + "t")
    return path.open(mode, encoding="utf-8")


def _is_bgzip_output(path: Path) -> bool:
    return path.suffix == ".gz"


def _parse_info(info_text: str) -> List[Tuple[str, Optional[str]]]:
    if info_text in {"", "."}:
        return []
    fields: List[Tuple[str, Optional[str]]] = []
    for item in info_text.split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            fields.append((key, value))
        else:
            fields.append((item, None))
    return fields


def _format_info(fields: List[Tuple[str, Optional[str]]]) -> str:
    if not fields:
        return "."
    parts = [key if value is None else f"{key}={value}" for key, value in fields]
    return ";".join(parts)


def _info_dict(fields: List[Tuple[str, Optional[str]]]) -> Dict[str, Optional[str]]:
    return {key: value for key, value in fields}


def _set_info_value(fields: List[Tuple[str, Optional[str]]], key: str, value: Optional[str]) -> List[Tuple[str, Optional[str]]]:
    updated = False
    result: List[Tuple[str, Optional[str]]] = []
    for field_key, field_value in fields:
        if field_key == key:
            if not updated:
                result.append((key, value))
                updated = True
            continue
        result.append((field_key, field_value))
    if not updated:
        result.append((key, value))
    return result


def _remove_info_key(fields: List[Tuple[str, Optional[str]]], key: str) -> List[Tuple[str, Optional[str]]]:
    return [(field_key, field_value) for field_key, field_value in fields if field_key != key]


def _rewrite_gt_for_dup_multiallelic(gt_value: str, alt_values: List[str]) -> str:
    if gt_value in {".", "./.", ".|."}:
        return gt_value
    separator = "/"
    if "|" in gt_value:
        separator = "|"
    alleles = gt_value.replace("|", "/").split("/")
    rewritten: List[str] = []
    for allele in alleles:
        if allele == ".":
            rewritten.append(".")
            continue
        allele_index = int(allele)
        if allele_index == 0:
            rewritten.append("0")
            continue
        alt_index = allele_index - 1
        cn_label = alt_values[alt_index] if 0 <= alt_index < len(alt_values) else ""
        if cn_label.startswith("<CN") and cn_label.endswith(">"):
            try:
                cn_value = int(cn_label[3:-1])
            except ValueError:
                cn_value = 2
            rewritten.append("1" if cn_value > 2 else "0")
        else:
            rewritten.append("1")
    return separator.join(rewritten)


def _fix_record_line(line: str) -> str:
    columns = line.rstrip("\n").split("\t")
    if len(columns) < 8:
        return line
    chrom, pos, record_id, ref, alt, qual, filt, info_text = columns[:8]
    format_text = columns[8] if len(columns) > 8 else None
    sample_texts = columns[9:] if len(columns) > 9 else []

    alt_values = alt.split(",") if alt not in {"", "."} else []
    info_fields = _parse_info(info_text)
    info_values = _info_dict(info_fields)
    svtype = info_values.get("SVTYPE")
    end_value = info_values.get("END")
    cpx_type = info_values.get("CPX_TYPE")

    if filt == ".":
        filt = "PASS"

    if any(bracket in alt for bracket in ("[", "]")):
        alt = "<BND>"

    if svtype == "DUP" and alt_values == ["<CN0>", "<CN1>", "<CN2>", "<CN3>"]:
        alt = "<DUP>"
        if format_text is not None:
            format_keys = format_text.split(":")
            if "GT" in format_keys:
                gt_index = format_keys.index("GT")
                rewritten_samples = []
                for sample_text in sample_texts:
                    sample_values = sample_text.split(":")
                    if gt_index < len(sample_values):
                        sample_values[gt_index] = _rewrite_gt_for_dup_multiallelic(sample_values[gt_index], alt_values)
                    rewritten_samples.append(":".join(sample_values))
                sample_texts = rewritten_samples

    svlen_text = info_values.get("SVLEN")
    if svlen_text not in (None, "."):
        try:
            svlen_value = int(svlen_text.split(",", 1)[0])
        except ValueError:
            svlen_value = None
        if svlen_value is not None and svlen_value < 0:
            info_fields = _set_info_value(info_fields, "SVLEN", str(abs(svlen_value)))
            info_values = _info_dict(info_fields)

    if svtype == "INS" and end_value not in (None, pos):
        info_fields = _set_info_value(info_fields, "END", pos)
        info_values = _info_dict(info_fields)
        end_value = pos

    if svtype in {"BND", "CTX"} and end_value not in (None, pos) and info_values.get("END2") in (None, "."):
        info_fields = _set_info_value(info_fields, "END2", str(end_value))
        info_fields = _set_info_value(info_fields, "END", pos)
        info_values = _info_dict(info_fields)

    if svtype == "CPX" and cpx_type in {"dDUP", "dDUP_iDEL"} and end_value not in (None, pos):
        info_fields = _set_info_value(info_fields, "END", pos)
        info_values = _info_dict(info_fields)

    if svtype not in {"BND", "CTX"} and "CHR2" in info_values:
        info_fields = _remove_info_key(info_fields, "CHR2")
        info_values = _info_dict(info_fields)

    if format_text is not None:
        format_keys = format_text.split(":")
        if "GQ" in format_keys:
            gq_index = format_keys.index("GQ")
            rewritten_samples = []
            for sample_text in sample_texts:
                sample_values = sample_text.split(":")
                if gq_index < len(sample_values) and sample_values[gq_index] not in {".", ""}:
                    try:
                        sample_values[gq_index] = str(int(int(sample_values[gq_index]) / 10)) if int(sample_values[gq_index]) > 99 else sample_values[gq_index]
                    except ValueError:
                        pass
                rewritten_samples.append(":".join(sample_values))
            sample_texts = rewritten_samples

    columns[4] = alt
    columns[6] = filt
    columns[7] = _format_info(info_fields)
    if format_text is not None and sample_texts:
        columns[9:] = sample_texts
    return "\t".join(columns) + "\n"


def apply_fixes(vcf_path: Path, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    temp_output_path = out_path.with_suffix("") if _is_bgzip_output(out_path) else out_path
    with _open_text(vcf_path, "r") as handle_in, temp_output_path.open("w", encoding="utf-8") as handle_out:
        for line in handle_in:
            if line.startswith("##INFO=<ID=SVLEN,") and "Number=." in line:
                handle_out.write(line.replace("Number=.", "Number=1"))
                continue
            if line.startswith("#"):
                handle_out.write(line)
                continue
            handle_out.write(_fix_record_line(line))
    if _is_bgzip_output(out_path):
        pysam.tabix_compress(str(temp_output_path), str(out_path), force=True)
        pysam.tabix_index(str(out_path), preset="vcf", force=True)
        temp_output_path.unlink()


def validate_and_fix(vcf_path: Path, out_path: Path) -> ValidationFixResult:
    original_summary = validate_vcf(ValidateConfig(vcf_path=vcf_path))
    unfixable_errors = [
        issue for issue in original_summary.issues if issue.severity == "ERROR" and issue.check_id not in _FIXABLE_CHECK_IDS
    ]
    if unfixable_errors:
        return ValidationFixResult(original_summary=original_summary, fixed_summary=None, out_path=None, wrote_output=False)
    apply_fixes(vcf_path, out_path)
    fixed_summary = validate_vcf(ValidateConfig(vcf_path=out_path))
    return ValidationFixResult(original_summary=original_summary, fixed_summary=fixed_summary, out_path=out_path, wrote_output=True)


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
        issue_type_counts = summary.issue_type_counts()
        exemplar_issues: Dict[str, FormatIssue] = {}
        for issue in summary.issues:
            exemplar_issues.setdefault(issue.check_id, issue)
        lines.append("Details:")
        for issue in exemplar_issues.values():
            where = f" [{issue.record_id}]" if issue.record_id else ""
            lines.append(f"- {issue.severity} {issue.check_id}{where}: {issue.message}")
        lines.append("Issue counts:")
        for check_id, count in sorted(issue_type_counts.items()):
            exemplar = exemplar_issues[check_id]
            lines.append(f"- {exemplar.severity} {check_id}: {count}")
    else:
        lines.append("No issues found.")
    return "\n".join(lines)


def render_fix_result(result: ValidationFixResult) -> str:
    lines = ["Original:", render_summary(result.original_summary)]
    if not result.wrote_output:
        lines.append("Fix mode aborted: unfixable validation errors remain.")
        return "\n".join(lines)
    lines.append(f"Wrote fixed VCF: {result.out_path}")
    if result.fixed_summary is not None:
        lines.append("Fixed:")
        lines.append(render_summary(result.fixed_summary))
    return "\n".join(lines)


def validate_and_render(vcf_path: Path) -> tuple[ValidationSummary, str]:
    """Convenience wrapper used by the CLI."""
    summary = validate_vcf(ValidateConfig(vcf_path=vcf_path))
    return summary, render_summary(summary)


def fix_and_render(vcf_path: Path, out_path: Path) -> tuple[ValidationFixResult, str]:
    """Convenience wrapper used by the CLI for validate --fix."""
    result = validate_and_fix(vcf_path, out_path)
    return result, render_fix_result(result)
