"""VCF format inspection helpers for gatk-sv-compare."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum, unique
from typing import Iterable, List, Optional, Set, TypeVar, overload

import pysam

from .dimensions import SVTYPES_CORE

_ALLOWED_SVTYPES = set(SVTYPES_CORE)
_DefaultT = TypeVar("_DefaultT")


@dataclass(frozen=True)
class FormatIssue:
    """A validation issue emitted while inspecting a VCF."""

    check_id: str
    severity: str
    record_id: str
    message: str


@unique
class PipelineStage(Enum):
    FINAL_ANNOTATED = "final_annotated"
    CLEANED = "cleaned"
    REGENOTYPED = "regenotyped"
    GENOTYPED = "genotyped"
    CLUSTER_BATCH = "cluster_batch"
    UNKNOWN = "unknown"


def _record_identifier(record: pysam.VariantRecord) -> str:
    return record.id or f"{record.contig}:{record.pos}"


def filter_values(record: pysam.VariantRecord) -> Set[str]:
    """Return the FILTER values for a record as a normalized set."""
    values = set(record.filter.keys())
    if not values and record.filter is not None:
        raw = str(record.filter)
        if raw and raw != ".":
            values.add(raw)
    return values


def _iter_alt_alleles(record: pysam.VariantRecord) -> Iterable[str]:
    return tuple(record.alts or ())


@overload
def safe_info_get(record: pysam.VariantRecord, key: str) -> object | None:
    ...


@overload
def safe_info_get(record: pysam.VariantRecord, key: str, default: _DefaultT) -> object | _DefaultT:
    ...


def safe_info_get(record: pysam.VariantRecord, key: str, default: _DefaultT | None = None) -> object | _DefaultT | None:
    """Safely retrieve an INFO value even when the field is absent from the header."""
    if key not in record.header.info:
        return default
    try:
        return record.info.get(key, default)
    except ValueError:
        return default


def is_mei(alt_allele: str) -> bool:
    return "<INS:ME:" in alt_allele or alt_allele == "<INS:ME>"


def is_cnv(record: pysam.VariantRecord) -> bool:
    return str(safe_info_get(record, "SVTYPE", "")) == "CNV"


def has_precomputed_counts(vcf: pysam.VariantFile) -> bool:
    info = vcf.header.info
    return all(field in info for field in ("N_BI_GENOS", "N_HOMREF", "N_HET", "N_HOMALT"))


def detect_pipeline_stage(vcf: pysam.VariantFile) -> PipelineStage:
    """Infer the pipeline stage from header content."""
    info = vcf.header.info
    formats = vcf.header.formats

    has_ecn = "ECN" in formats
    has_members = "MEMBERS" in info
    has_vargq = "varGQ" in info
    has_precomputed = has_precomputed_counts(vcf)
    has_ogq = "OGQ" in formats
    has_sl = "SL" in formats
    has_gnomad = any(str(key).startswith("gnomad_") for key in info.keys())
    has_gq = "GQ" in formats

    if has_precomputed and has_ogq and has_sl and has_gnomad:
        return PipelineStage.FINAL_ANNOTATED
    if has_ecn and not has_precomputed and not has_members:
        return PipelineStage.CLEANED
    if has_vargq and not has_ecn:
        for record in vcf:
            alts = set(_iter_alt_alleles(record))
            if {"<CN0>", "<CN1>", "<CN2>", "<CN3>"}.issubset(alts):
                return PipelineStage.REGENOTYPED
            break
    if has_members and has_vargq and not has_ecn:
        return PipelineStage.GENOTYPED
    if has_members and not has_gq:
        return PipelineStage.CLUSTER_BATCH
    return PipelineStage.UNKNOWN


def check_header(vcf: pysam.VariantFile) -> List[FormatIssue]:
    """Inspect VCF header-level format issues."""
    issues: list[FormatIssue] = []
    if "SVLEN" in vcf.header.info and vcf.header.info["SVLEN"].number != 1:
        issues.append(
            FormatIssue(
                "SVLEN_NUMBER_DOT",
                "WARN",
                "",
                "SVLEN header Number is not 1",
            )
        )
    return issues


def check_record(record: pysam.VariantRecord, contig_length: Optional[int] = None) -> List[FormatIssue]:
    """Inspect a single record and return any format issues."""
    issues: list[FormatIssue] = []
    record_id = _record_identifier(record)
    svtype = safe_info_get(record, "SVTYPE")
    alts = tuple(_iter_alt_alleles(record))
    filters = filter_values(record)

    if len(alts) > 1 and svtype != "CNV":
        issues.append(FormatIssue("MULTI_ALLELIC_NON_CNV", "ERROR", record_id, "Non-CNV record has more than one ALT allele"))

    if svtype is None:
        issues.append(FormatIssue("MISSING_SVTYPE", "ERROR", record_id, "SVTYPE INFO field is missing"))
        return issues

    svtype_text = str(svtype)
    if svtype_text not in _ALLOWED_SVTYPES:
        issues.append(FormatIssue("UNKNOWN_SVTYPE", "WARN", record_id, f"Unexpected SVTYPE: {svtype_text}"))

    has_bracket_alt = any("]" in alt or "[" in alt for alt in alts)
    if has_bracket_alt:
        issues.append(FormatIssue("BREAKEND_NOTATION", "ERROR", record_id, "Breakend notation ALT allele is not supported"))

    has_chr2 = "CHR2" in record.info
    has_end2 = "END2" in record.info
    if svtype_text == "BND" and not has_bracket_alt and (not has_chr2 or not has_end2):
        issues.append(
            FormatIssue(
                "MISSING_BND_COORDS",
                "ERROR",
                record_id,
                "BND record requires CHR2 and END2 unless ALT uses bracket notation",
            )
        )
    if svtype_text == "CTX" and (not has_chr2 or not has_end2):
        issues.append(
            FormatIssue(
                "MISSING_CTX_COORDS",
                "ERROR",
                record_id,
                "CTX record requires CHR2 and END2",
            )
        )

    svlen = safe_info_get(record, "SVLEN")
    if svtype_text in {"DEL", "DUP", "INS", "INV", "CPX"} and svlen in (None, "."):
        issues.append(FormatIssue("MISSING_SVLEN", "WARN", record_id, f"SVLEN missing for {svtype_text}"))
    if svlen not in (None, "."):
        if isinstance(svlen, tuple):
            svlen_value = int(svlen[0])
        else:
            svlen_value = int(svlen)
        if svlen_value < 0:
            issues.append(FormatIssue("SVLEN_SIGN", "WARN", record_id, "SVLEN is negative"))

    end = record.stop
    if svtype_text == "INS" and end != record.pos:
        issues.append(FormatIssue("INS_END_MISMATCH", "WARN", record_id, "INS END does not equal POS"))
    if svtype_text in {"BND", "CTX"} and end != record.pos and "END2" not in record.info:
        issues.append(FormatIssue("BND_END_MISMATCH", "WARN", record_id, "BND/CTX END does not equal POS and END2 is missing"))
    cpx_type = str(safe_info_get(record, "CPX_TYPE", ""))
    if svtype_text == "CPX" and cpx_type in {"dDUP", "dDUP_iDEL"} and end != record.pos:
        issues.append(FormatIssue("CPX_DDDUP_END", "WARN", record_id, "Dispersed duplication CPX END does not equal POS"))

    if "GT" not in record.format:
        issues.append(FormatIssue("MISSING_GT", "ERROR", record_id, "GT FORMAT field is missing"))
    if "ECN" not in record.format:
        issues.append(FormatIssue("MISSING_ECN", "ERROR", record_id, "ECN FORMAT field is missing"))
    if "GQ" in record.format:
        for sample in record.samples.values():
            gq = sample.get("GQ")
            if gq not in (None, ".") and int(gq) > 99:
                issues.append(FormatIssue("GQ_RANGE", "WARN", record_id, "GQ exceeds expected 0-99 range"))
                break
    if not filters:
        issues.append(FormatIssue("EMPTY_FILTER", "WARN", record_id, "FILTER is empty"))
    if "CHR2" in record.info and svtype_text not in {"BND", "CTX"}:
        issues.append(FormatIssue("CHR2_ON_NON_BND", "INFO", record_id, "CHR2 present on non-BND/CTX record"))
    if "MEMBERS" in record.info:
        issues.append(FormatIssue("MEMBERS_PRESENT", "INFO", record_id, "MEMBERS INFO field is present"))
    if "varGQ" in record.info:
        issues.append(FormatIssue("VARGQ_PRESENT", "INFO", record_id, "varGQ INFO field is present"))
    if "MULTIALLELIC" in record.info:
        issues.append(FormatIssue("MULTIALLELIC_INFO_FLAG", "INFO", record_id, "MULTIALLELIC is encoded as INFO rather than FILTER"))
    return issues
