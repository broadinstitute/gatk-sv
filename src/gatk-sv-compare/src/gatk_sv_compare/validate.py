"""Validate command implementation for gatk-sv-compare."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
import gzip
from pathlib import Path
import re
from typing import Dict, Iterator, List, Optional, Tuple

import pysam

from .config import ValidateConfig
from .gq_utils import detect_gq_scale_factor
from .vcf_format import FormatIssue, PipelineStage, check_header, check_record, detect_pipeline_stage, has_precomputed_counts

_FIXABLE_CHECK_IDS = {
    "MISSING_SVTYPE",
    "BREAKEND_NOTATION",
    "MISSING_ECN",
    "MULTI_ALLELIC_NON_CNV",
    "MISSING_ALGORITHMS",
    "MISSING_END_HEADER",
    "MISSING_SVTYPE_HEADER",
    "MISSING_SVLEN_HEADER",
}
_BND_MATE_PATTERN = re.compile(r"[\[\]]([^:\[\]]+):(\d+)[\[\]]")
_SYMBOLIC_ALT_PATTERN = re.compile(r"^<([^>]+)>$")
_CN_ALT_PATTERN = re.compile(r"^<CN\d+>$")
_ORIGINAL_ALT_INFO_HEADER = '##INFO=<ID=ORIGINAL_ALT,Number=1,Type=String,Description="Original ALT allele before validate --fix canonicalization">'
_ECN_FORMAT_HEADER = '##FORMAT=<ID=ECN,Number=1,Type=Integer,Description="Expected copy number">'
_CHR2_INFO_HEADER = '##INFO=<ID=CHR2,Number=1,Type=String,Description="Partner chromosome">'
_END2_INFO_HEADER = '##INFO=<ID=END2,Number=1,Type=Integer,Description="Partner end">'
_ALGORITHMS_INFO_HEADER = '##INFO=<ID=ALGORITHMS,Number=.,Type=String,Description="Source algorithms">'
_END_INFO_HEADER = '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">'
_SVTYPE_INFO_HEADER = '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">'
_SVLEN_INFO_HEADER = '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">'


@dataclass
class ValidationSummary:
    """Summary of a validation run."""

    vcf_path: Path
    stage: PipelineStage
    issues: List[FormatIssue]
    record_count: int
    precomputed_counts_available: bool
    sites_only: bool = False

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
    # Errors that were truly unfixable and caused an early abort (fixed_summary=None).
    # When present, only these are reported as blocking — not all original errors.
    unfixable_errors: List[FormatIssue] = field(default_factory=list)

    @property
    def has_errors(self) -> bool:
        if self.fixed_summary is not None:
            return self.fixed_summary.has_errors
        return self.original_summary.has_errors

    def unresolved_error_counts(self) -> Counter[str]:
        if self.fixed_summary is not None:
            return Counter(issue.check_id for issue in self.fixed_summary.issues if issue.severity == "ERROR")
        if self.unfixable_errors:
            return Counter(issue.check_id for issue in self.unfixable_errors)
        return Counter(issue.check_id for issue in self.original_summary.issues if issue.severity == "ERROR")

    def unresolved_errors(self) -> List[FormatIssue]:
        if self.fixed_summary is not None:
            return [issue for issue in self.fixed_summary.issues if issue.severity == "ERROR"]
        if self.unfixable_errors:
            return list(self.unfixable_errors)
        return [issue for issue in self.original_summary.issues if issue.severity == "ERROR"]


@dataclass(frozen=True)
class MateCoordinateRepair:
    """Repair values for legacy BND/CTX mate-coordinate encoding."""

    chr2: str
    end2: int


@dataclass
class ValidationScanResult:
    """Validation results plus fix metadata collected during the same scan."""

    summary: ValidationSummary
    record_svtypes: Dict[str, Optional[str]] = field(default_factory=dict)
    mate_coordinate_repairs: Dict[str, MateCoordinateRepair] = field(default_factory=dict)


def _record_key(chrom: str, pos: str | int, record_id: str) -> str:
    if record_id and record_id != ".":
        return record_id
    return f"{chrom}:{pos}"


def _iter_record_lines(vcf_path: Path) -> Iterator[str]:
    with _open_text(vcf_path, "r") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            yield line


def _resolve_record_svtype(record: pysam.VariantRecord) -> Optional[str]:
    svtype = record.info.get("SVTYPE") if "SVTYPE" in record.info else None
    if svtype not in (None, "."):
        return str(svtype)
    alt_values = [str(alt) for alt in (record.alts or ())]
    return _infer_svtype(record.contig, alt_values)


def _extract_mate_coordinate_repair(record: pysam.VariantRecord, raw_line: str, svtype: Optional[str]) -> Optional[MateCoordinateRepair]:
    if svtype not in {"BND", "CTX"}:
        return None
    columns = raw_line.rstrip("\n").split("\t")
    if len(columns) < 8:
        return None
    raw_pos = columns[1]
    raw_info_values = _info_dict(_parse_info(columns[7]))

    current_chr2 = record.info.get("CHR2") if "CHR2" in record.info else None
    current_end2 = record.info.get("END2") if "END2" in record.info else None
    missing_chr2 = current_chr2 in (None, ".")
    missing_end2 = current_end2 in (None, ".")
    if not missing_chr2 and not missing_end2:
        return None

    resolved_chr2 = str(current_chr2) if not missing_chr2 else raw_info_values.get("CHR2")
    if resolved_chr2 in (None, "."):
        return None

    if not missing_end2:
        resolved_end2 = int(current_end2)
    else:
        raw_end = raw_info_values.get("END")
        if raw_end in (None, "."):
            resolved_end2 = int(raw_pos)
        else:
            try:
                resolved_end2 = int(raw_end)
            except ValueError:
                return None
    return MateCoordinateRepair(chr2=str(resolved_chr2), end2=resolved_end2)


def _scan_vcf(config: ValidateConfig, *, collect_fix_metadata: bool = False) -> ValidationScanResult:
    """Validate a VCF and optionally collect repair metadata during the same record scan."""
    with pysam.VariantFile(str(config.vcf_path)) as vcf:
        stage = detect_pipeline_stage(vcf)
        precomputed = has_precomputed_counts(vcf)
        header_issues = check_header(vcf)
        sites_only = len(vcf.header.samples) == 0

    issues: List[FormatIssue] = list(header_issues)
    if sites_only:
        issues.append(
            FormatIssue(
                "SITES_ONLY",
                "INFO",
                "",
                "VCF is sites-only (no sample columns); sample-level checks are skipped",
            )
        )
    record_count = 0
    record_svtypes: Dict[str, Optional[str]] = {}
    mate_coordinate_repairs: Dict[str, MateCoordinateRepair] = {}
    record_lines = _iter_record_lines(config.vcf_path) if collect_fix_metadata else None
    with pysam.VariantFile(str(config.vcf_path)) as vcf:
        contig_lengths = {name: int(contig.length) for name, contig in vcf.header.contigs.items() if contig.length is not None}
        for record in vcf:
            raw_line = next(record_lines) if record_lines is not None else None
            record_count += 1
            issues.extend(check_record(record, contig_length=contig_lengths.get(record.contig), sites_only=sites_only))
            if collect_fix_metadata and raw_line is not None:
                record_key = _record_key(record.contig, record.pos, record.id)
                svtype = _resolve_record_svtype(record)
                record_svtypes[record_key] = svtype
                mate_coordinate_repair = _extract_mate_coordinate_repair(record, raw_line, svtype)
                if mate_coordinate_repair is not None:
                    mate_coordinate_repairs[record_key] = mate_coordinate_repair
            if record_count == 1 and not precomputed:
                issues.append(
                    FormatIssue(
                        "MISSING_PRECOMPUTED_COUNTS",
                        "INFO",
                        "",
                        "Pre-computed genotype count fields are absent; validation will allow GT-based fallback",
                    )
                )
    return ValidationScanResult(
        summary=ValidationSummary(
            vcf_path=config.vcf_path,
            stage=stage,
            issues=issues,
            record_count=record_count,
            precomputed_counts_available=precomputed,
            sites_only=sites_only,
        ),
        record_svtypes=record_svtypes,
        mate_coordinate_repairs=mate_coordinate_repairs,
    )


def validate_vcf(config: ValidateConfig) -> ValidationSummary:
    """Validate a VCF and return structured results."""
    return _scan_vcf(config).summary


def _open_text(path: Path, mode: str):
    if path.suffix == ".gz":
        return gzip.open(path, mode + "t")
    return path.open(mode, encoding="utf-8")


def _is_bgzip_output(path: Path) -> bool:
    return path.suffix == ".gz"


def _normalize_bgzip_output_path(path: Path) -> Path:
    if path.name.endswith(".vcf.gz"):
        return path
    if path.suffix == ".vcf":
        return path.with_name(path.name + ".gz")
    return path.with_name(path.name + ".vcf.gz")


def _parse_ploidy_table(path: Path) -> Dict[str, Dict[str, int]]:
    ploidy_dict: Dict[str, Dict[str, int]] = {}
    with path.open("r", encoding="utf-8") as handle:
        header = handle.readline().strip().split("\t")
        for line in handle:
            tokens = line.strip().split("\t")
            if not tokens or len(tokens) < 2:
                continue
            sample = tokens[0]
            ploidy_dict[sample] = {header[index]: int(tokens[index]) for index in range(1, len(header))}
    return ploidy_dict


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


def _parse_breakend_alt(alt: str) -> Optional[Tuple[str, str]]:
    match = _BND_MATE_PATTERN.search(alt)
    if match is None:
        return None
    mate_chrom, mate_pos = match.groups()
    return mate_chrom, mate_pos


def _infer_svtype(chrom: str, alt_values: List[str]) -> Optional[str]:
    if not alt_values:
        return None
    if len(alt_values) > 1 and all(_CN_ALT_PATTERN.fullmatch(alt_value) for alt_value in alt_values):
        return "CNV"
    if len(alt_values) != 1:
        return None
    alt_value = alt_values[0]
    if "[" in alt_value or "]" in alt_value:
        mate = _parse_breakend_alt(alt_value)
        if mate is None:
            return None
        mate_chrom, _ = mate
        return "CTX" if mate_chrom != chrom else "BND"
    symbolic_match = _SYMBOLIC_ALT_PATTERN.fullmatch(alt_value)
    if symbolic_match is None:
        return None
    symbolic_type = symbolic_match.group(1)
    if symbolic_type in {"DEL", "DUP", "INS", "INV", "BND", "CTX", "CPX", "CNV"}:
        return symbolic_type
    if symbolic_type.startswith("INS:"):
        return "INS"
    if symbolic_type.startswith("CN") and symbolic_type[2:].isdigit():
        return "CNV"
    return None


def _is_cn_alt_list(alt_values: List[str]) -> bool:
    return len(alt_values) > 1 and all(_CN_ALT_PATTERN.fullmatch(alt_value) for alt_value in alt_values)


def _clear_sample_gt(sample_text: str, format_keys: List[str]) -> str:
    if not format_keys:
        return sample_text
    sample_values = sample_text.split(":") if sample_text else []
    while len(sample_values) < len(format_keys):
        sample_values.append(".")
    gt_index = format_keys.index("GT") if "GT" in format_keys else None
    if gt_index is not None:
        sample_values[gt_index] = "./."
    return ":".join(sample_values)


def _normalize_sample_gq(sample_text: str, format_keys: List[str], gq_scale_factor: float) -> str:
    if gq_scale_factor == 1.0 or not format_keys or "GQ" not in format_keys:
        return sample_text
    sample_values = sample_text.split(":") if sample_text else []
    while len(sample_values) < len(format_keys):
        sample_values.append(".")
    gq_index = format_keys.index("GQ")
    gq_value = sample_values[gq_index]
    if gq_value not in {"", "."}:
        scaled = min(99, int(float(gq_value) / gq_scale_factor))
        sample_values[gq_index] = str(scaled)
    return ":".join(sample_values)


def _has_bracket_alt(alt_values: List[str]) -> bool:
    return len(alt_values) == 1 and any(bracket in alt_values[0] for bracket in ("[", "]"))


def _missing_mate_coords(info_values: Dict[str, Optional[str]]) -> bool:
    return info_values.get("CHR2") in (None, ".") or info_values.get("END2") in (None, ".")


def _fix_record_line(
    line: str,
    ploidy_dict: Optional[Dict[str, Dict[str, int]]] = None,
    sample_names: Optional[List[str]] = None,
    drop_bnd: bool = False,
    drop_ctx: bool = False,
    mate_coordinate_repairs: Optional[Dict[str, MateCoordinateRepair]] = None,
    gq_scale_factor: float = 1.0,
) -> Optional[str]:
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
    record_key = _record_key(chrom, pos, record_id)

    if svtype in (None, "."):
        inferred_svtype = _infer_svtype(chrom, alt_values)
        if inferred_svtype is not None:
            info_fields = _set_info_value(info_fields, "SVTYPE", inferred_svtype)
            info_values = _info_dict(info_fields)
            svtype = inferred_svtype

    if info_values.get("ALGORITHMS") in (None, "."):
        info_fields = _set_info_value(info_fields, "ALGORITHMS", "UNKNOWN")
        info_values = _info_dict(info_fields)

    applied_legacy_mate_repair = False
    mate_coordinate_repair = mate_coordinate_repairs.get(record_key) if mate_coordinate_repairs is not None else None
    if svtype in {"BND", "CTX"} and mate_coordinate_repair is not None:
        if info_values.get("CHR2") in (None, "."):
            info_fields = _set_info_value(info_fields, "CHR2", mate_coordinate_repair.chr2)
            applied_legacy_mate_repair = True
        if info_values.get("END2") in (None, "."):
            info_fields = _set_info_value(info_fields, "END2", str(mate_coordinate_repair.end2))
            applied_legacy_mate_repair = True
        if applied_legacy_mate_repair:
            info_fields = _set_info_value(info_fields, "END", pos)
            info_values = _info_dict(info_fields)

    if svtype == "BND" and drop_bnd:
        return None

    if svtype == "CTX" and drop_ctx:
        return None

    if _is_cn_alt_list(alt_values):
        alt = "<CNV>"
        info_fields = _set_info_value(info_fields, "SVTYPE", "CNV")
        info_values = _info_dict(info_fields)
        svtype = "CNV"
        if format_text is not None and sample_texts:
            format_keys = format_text.split(":")
            sample_texts = [_clear_sample_gt(sample_text, format_keys) for sample_text in sample_texts]

    if svtype == "CTX" and _missing_mate_coords(info_values):
        if drop_ctx:
            return None

    if _has_bracket_alt(alt_values):
        breakend_alt = alt_values[0]
        mate = _parse_breakend_alt(breakend_alt)
        if mate is not None:
            mate_chrom, mate_pos = mate
            breakend_svtype = svtype if svtype in {"BND", "CTX"} else ("CTX" if mate_chrom != chrom else "BND")
            info_fields = _set_info_value(info_fields, "SVTYPE", breakend_svtype)
            info_fields = _set_info_value(info_fields, "END", pos)
            info_fields = _set_info_value(info_fields, "CHR2", mate_chrom)
            info_fields = _set_info_value(info_fields, "END2", mate_pos)
            info_fields = _set_info_value(info_fields, "ORIGINAL_ALT", breakend_alt)
            info_values = _info_dict(info_fields)
            alt = f"<{breakend_svtype}>"

    if ploidy_dict is not None and format_text is not None and "ECN" not in format_text.split(":"):
        if sample_names is None or len(sample_names) != len(sample_texts):
            raise ValueError("Cannot repair ECN without aligned sample names")
        format_keys = format_text.split(":") if format_text else []
        format_keys.append("ECN")
        rewritten_samples: List[str] = []
        for sample_name, sample_text in zip(sample_names, sample_texts):
            contig_ploidies = ploidy_dict.get(sample_name)
            if contig_ploidies is None or chrom not in contig_ploidies:
                raise ValueError(f"Missing ploidy entry for sample {sample_name} on contig {chrom}")
            sample_values = sample_text.split(":") if sample_text else []
            sample_values.append(str(contig_ploidies[chrom]))
            rewritten_samples.append(":".join(sample_values))
        format_text = ":".join(format_keys)
        sample_texts = rewritten_samples

    if format_text is not None and sample_texts:
        format_keys = format_text.split(":")
        sample_texts = [_normalize_sample_gq(sample_text, format_keys, gq_scale_factor) for sample_text in sample_texts]

    columns[4] = alt
    columns[7] = _format_info(info_fields)
    if format_text is not None:
        columns[8] = format_text
    if sample_texts:
        columns[9:] = sample_texts
    return "\t".join(columns) + "\n"


def _write_augmented_header(header_lines: List[str], handle_out) -> None:
    has_end = any(line.startswith("##INFO=<ID=END,") for line in header_lines)
    has_svtype = any(line.startswith("##INFO=<ID=SVTYPE,") for line in header_lines)
    has_svlen = any(line.startswith("##INFO=<ID=SVLEN,") for line in header_lines)
    has_chr2 = any(line.startswith("##INFO=<ID=CHR2,") for line in header_lines)
    has_end2 = any(line.startswith("##INFO=<ID=END2,") for line in header_lines)
    has_original_alt = any(line.startswith("##INFO=<ID=ORIGINAL_ALT,") for line in header_lines)
    has_ecn = any(line.startswith("##FORMAT=<ID=ECN,") for line in header_lines)
    has_algorithms = any(line.startswith("##INFO=<ID=ALGORITHMS,") for line in header_lines)

    for line in header_lines[:-1]:
        if line.startswith("##INFO=<ID=SVLEN,") and "Number=." in line:
            handle_out.write(line.replace("Number=.", "Number=1"))
        else:
            handle_out.write(line)
    if not has_end:
        handle_out.write(_END_INFO_HEADER + "\n")
    if not has_svtype:
        handle_out.write(_SVTYPE_INFO_HEADER + "\n")
    if not has_svlen:
        handle_out.write(_SVLEN_INFO_HEADER + "\n")
    if not has_chr2:
        handle_out.write(_CHR2_INFO_HEADER + "\n")
    if not has_end2:
        handle_out.write(_END2_INFO_HEADER + "\n")
    if not has_original_alt:
        handle_out.write(_ORIGINAL_ALT_INFO_HEADER + "\n")
    if not has_ecn:
        handle_out.write(_ECN_FORMAT_HEADER + "\n")
    if not has_algorithms:
        handle_out.write(_ALGORITHMS_INFO_HEADER + "\n")
    handle_out.write(header_lines[-1])


def apply_fixes(
    vcf_path: Path,
    out_path: Path,
    ploidy_table_path: Optional[Path] = None,
    drop_bnd: bool = False,
    drop_ctx: bool = False,
    mate_coordinate_repairs: Optional[Dict[str, MateCoordinateRepair]] = None,
) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    temp_output_path = out_path.with_suffix("") if _is_bgzip_output(out_path) else out_path
    ploidy_dict = _parse_ploidy_table(ploidy_table_path) if ploidy_table_path is not None else None
    with pysam.VariantFile(str(vcf_path)) as _vcf:
        _sites_only = len(_vcf.header.samples) == 0
    gq_scale_factor = 1.0 if _sites_only else detect_gq_scale_factor(vcf_path)
    sample_names: List[str] | None = None
    with _open_text(vcf_path, "r") as handle_in, temp_output_path.open("w", encoding="utf-8") as handle_out:
        header_lines: List[str] = []
        for line in handle_in:
            if line.startswith("##"):
                header_lines.append(line)
                continue
            if line.startswith("#"):
                header_lines.append(line)
                sample_names = line.rstrip("\n").split("\t")[9:]
                _write_augmented_header(header_lines, handle_out)
                continue
            fixed_line = _fix_record_line(
                line,
                ploidy_dict=ploidy_dict,
                sample_names=sample_names,
                drop_bnd=drop_bnd,
                drop_ctx=drop_ctx,
                mate_coordinate_repairs=mate_coordinate_repairs,
                gq_scale_factor=gq_scale_factor,
            )
            if fixed_line is None:
                continue
            handle_out.write(fixed_line)
    if _is_bgzip_output(out_path):
        pysam.tabix_compress(str(temp_output_path), str(out_path), force=True)
        pysam.tabix_index(str(out_path), preset="vcf", force=True)
        temp_output_path.unlink()


def _is_fixable_error(
    issue: FormatIssue,
    ploidy_table_path: Optional[Path],
    drop_bnd: bool,
    drop_ctx: bool,
    mate_coordinate_repairs: Optional[Dict[str, MateCoordinateRepair]] = None,
) -> bool:
    if issue.check_id == "MISSING_ECN":
        return ploidy_table_path is not None
    if issue.check_id == "MISSING_BND_COORDS":
        return drop_bnd or (issue.record_id in (mate_coordinate_repairs or {}))
    if issue.check_id == "MISSING_CTX_COORDS":
        return drop_ctx or (issue.record_id in (mate_coordinate_repairs or {}))
    return issue.check_id in _FIXABLE_CHECK_IDS


def validate_and_fix(
    vcf_path: Path,
    out_path: Path,
    ploidy_table_path: Optional[Path] = None,
    drop_bnd: bool = False,
    drop_ctx: bool = False,
) -> ValidationFixResult:
    resolved_out_path = _normalize_bgzip_output_path(out_path)
    scan_result = _scan_vcf(ValidateConfig(vcf_path=vcf_path), collect_fix_metadata=True)
    original_summary = scan_result.summary
    dropped_record_ids = {
        record_id
        for record_id, svtype in scan_result.record_svtypes.items()
        if (drop_bnd and svtype == "BND") or (drop_ctx and svtype == "CTX")
    }
    unfixable_errors = []
    for issue in original_summary.issues:
        if issue.severity != "ERROR":
            continue
        if issue.record_id and issue.record_id in dropped_record_ids:
            continue
        if _is_fixable_error(
            issue,
            ploidy_table_path,
            drop_bnd,
            drop_ctx,
            mate_coordinate_repairs=scan_result.mate_coordinate_repairs,
        ):
            continue
        unfixable_errors.append(issue)
    if unfixable_errors:
        return ValidationFixResult(original_summary=original_summary, fixed_summary=None, out_path=None, wrote_output=False, unfixable_errors=unfixable_errors)
    apply_fixes(
        vcf_path,
        resolved_out_path,
        ploidy_table_path=ploidy_table_path,
        drop_bnd=drop_bnd,
        drop_ctx=drop_ctx,
        mate_coordinate_repairs=scan_result.mate_coordinate_repairs,
    )
    fixed_summary = validate_vcf(ValidateConfig(vcf_path=resolved_out_path))
    if fixed_summary.has_errors:
        resolved_out_path.unlink(missing_ok=True)
        Path(str(resolved_out_path) + ".tbi").unlink(missing_ok=True)
        return ValidationFixResult(original_summary=original_summary, fixed_summary=fixed_summary, out_path=None, wrote_output=False)
    return ValidationFixResult(original_summary=original_summary, fixed_summary=fixed_summary, out_path=resolved_out_path, wrote_output=True)


def render_summary(summary: ValidationSummary) -> str:
    """Render a human-readable validation report."""
    counts = summary.severity_counts()
    lines = [
        f"VCF: {summary.vcf_path}",
        f"Stage: {summary.stage.value}",
        f"Records: {summary.record_count}",
        f"Sites-only: {summary.sites_only}",
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
        unresolved_errors = result.unresolved_error_counts()
        if unresolved_errors:
            offending_checks = ", ".join(
                f"{check_id} ({count})" if count > 1 else check_id
                for check_id, count in sorted(unresolved_errors.items())
            )
            lines.append(f"Fix mode aborted: unresolved critical checks remain: {offending_checks}")
            if "MISSING_ECN" in unresolved_errors:
                lines.append("Hint: rerun with --ploidy-table to repair missing ECN fields.")
            if "MISSING_BND_COORDS" in unresolved_errors:
                lines.append("Hint: rerun with --drop-bnd to omit all BND records from the comparison.")
            if "MISSING_CTX_COORDS" in unresolved_errors:
                lines.append("Hint: rerun with --drop-ctx to omit all CTX records from the comparison.")
            exemplar_errors: Dict[str, FormatIssue] = {}
            for issue in result.unresolved_errors():
                exemplar_errors.setdefault(issue.check_id, issue)
            lines.append("Examples:")
            for check_id, issue in sorted(exemplar_errors.items()):
                where = f" [{issue.record_id}]" if issue.record_id else ""
                lines.append(f"- {check_id}{where}: {issue.message}")
        else:
            lines.append("Fix mode aborted: unresolved critical checks remain.")
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


def fix_and_render(
    vcf_path: Path,
    out_path: Path,
    ploidy_table_path: Optional[Path] = None,
    drop_bnd: bool = False,
    drop_ctx: bool = False,
) -> tuple[ValidationFixResult, str]:
    """Convenience wrapper used by the CLI for validate --fix."""
    result = validate_and_fix(
        vcf_path,
        out_path,
        ploidy_table_path=ploidy_table_path,
        drop_bnd=drop_bnd,
        drop_ctx=drop_ctx,
    )
    return result, render_fix_result(result)
