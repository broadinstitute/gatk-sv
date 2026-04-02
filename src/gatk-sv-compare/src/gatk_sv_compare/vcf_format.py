"""VCF format inspection helpers for gatk-sv-compare."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum, unique
from typing import Iterable, List, Optional, Set

import pysam

from .dimensions import SVTYPES_CORE

_ALLOWED_SVTYPES = set(SVTYPES_CORE)


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


def is_mei(alt_allele: str) -> bool:
    return "<INS:ME:" in alt_allele or alt_allele == "<INS:ME>"


def is_cnv(record: pysam.VariantRecord) -> bool:
    return str(record.info.get("SVTYPE", "")) == "CNV"


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


def check_record(record: pysam.VariantRecord, contig_length: Optional[int] = None) -> List[FormatIssue]:
    """Inspect a single record and return any format issues."""
    issues: list[FormatIssue] = []
    record_id = _record_identifier(record)
    svtype = record.info.get("SVTYPE")
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

    if any("]" in alt or "[" in alt for alt in alts):
        issues.append(FormatIssue("BREAKEND_NOTATION", "ERROR", record_id, "Breakend notation ALT allele is not supported"))

    svlen = record.info.get("SVLEN")
    if svtype_text in {"DEL", "DUP", "INS", "INV", "CPX"} and svlen in (None, "."):
        issues.append(FormatIssue("MISSING_SVLEN", "WARN", record_id, f"SVLEN missing for {svtype_text}"))
    if svlen not in (None, "."):
        if isinstance(svlen, tuple):
            svlen_value = int(svlen[0])
        else:
            svlen_value = int(svlen)
        if svlen_value < 0:
            issues.append(FormatIssue("SVLEN_SIGN", "WARN", record_id, "SVLEN is negative"))
        if contig_length is not None and abs(svlen_value) > contig_length * 0.20:
            issues.append(FormatIssue("IMPLAUSIBLE_SVLEN", "WARN", record_id, "SVLEN exceeds 20% of chromosome length and may be artifactual"))

    end = record.stop
    if svtype_text == "INS" and end != record.pos:
        issues.append(FormatIssue("INS_END_MISMATCH", "WARN", record_id, "INS END does not equal POS"))
    if svtype_text in {"BND", "CTX"} and end != record.pos and "END2" not in record.info:
        issues.append(FormatIssue("BND_END_MISMATCH", "WARN", record_id, "BND/CTX END does not equal POS and END2 is missing"))
    cpx_type = str(record.info.get("CPX_TYPE", ""))
    if svtype_text == "CPX" and cpx_type in {"dDUP", "dDUP_iDEL"} and end != record.pos:
        issues.append(FormatIssue("CPX_DDDUP_END", "WARN", record_id, "Dispersed duplication CPX END does not equal POS"))

    if "GT" not in record.format:
        issues.append(FormatIssue("MISSING_GT", "ERROR", record_id, "GT FORMAT field is missing"))
    if "ECN" not in record.format:
        issues.append(FormatIssue("MISSING_ECN", "WARN", record_id, "ECN FORMAT field is missing"))
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
