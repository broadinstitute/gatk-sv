"""
VCF utility helpers shared by multiple subcommands.

These replace the scattered helpers that lived in ``svtk.utils``,
``overlap_breakpoint_filter.py``, ``postCPX_cleanup.py``, etc.
"""

from __future__ import annotations

import logging
from typing import List, Set, Tuple

import pysam

from gatk_sv_cpx.constants import (
    ALGORITHMS,
    CHR2,
    CPX_INTERVALS,
    END2,
    EVIDENCE,
    SOURCE,
    STRANDS,
    SVLEN,
    SVTYPE,
    VCF_HEADER_LINES,
)

log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# VCF header helpers
# ---------------------------------------------------------------------------

def ensure_header_lines(header: pysam.VariantHeader) -> None:
    """Add any missing INFO / FILTER / ALT lines that the pipeline expects."""
    for key, line in VCF_HEADER_LINES.items():
        try:
            header.add_line(line)
        except (ValueError, OSError):
            pass  # already present


# ---------------------------------------------------------------------------
# Basic record accessors
# ---------------------------------------------------------------------------

def get_svtype(record: pysam.VariantRecord) -> str:
    return record.info[SVTYPE]


def get_svlen(record: pysam.VariantRecord) -> int:
    return record.info.get(SVLEN, 0)


def get_strands(record: pysam.VariantRecord) -> str:
    return record.info.get(STRANDS, "")


def get_evidence(record: pysam.VariantRecord) -> Set[str]:
    ev = record.info.get(EVIDENCE, ())
    return set(ev) if ev else set()


def get_algorithms(record: pysam.VariantRecord) -> Set[str]:
    alg = record.info.get(ALGORITHMS, ())
    return set(alg) if alg else set()


def get_called_samples(record: pysam.VariantRecord) -> List[str]:
    """Return sample IDs with non-ref genotype (any allele > 0)."""
    called = []
    for sample in record.samples:
        gt = record.samples[sample]["GT"]
        if gt is not None and any(a is not None and a > 0 for a in gt):
            called.append(sample)
    return called


def get_end(record: pysam.VariantRecord) -> int:
    """Canonical end position, respecting END2 for BND."""
    if get_svtype(record) == "BND":
        return record.info.get(END2, record.stop)
    return record.stop


# ---------------------------------------------------------------------------
# Breakpoint keys (used by overlap filter)
# ---------------------------------------------------------------------------

def breakend_key(chrom: str, pos: int, strand: str) -> str:
    return f"{chrom}_{pos}_{strand}"


def record_breakend_keys(record: pysam.VariantRecord) -> Tuple[str, str]:
    """Return (bnd1_key, bnd2_key) for a variant record."""
    strands = get_strands(record)
    end = get_end(record)
    chr2 = record.info.get(CHR2, record.chrom)
    bnd1 = breakend_key(record.chrom, record.pos, strands[0] if strands else "+")
    bnd2 = breakend_key(chr2, end, strands[1] if len(strands) > 1 else "+")
    return bnd1, bnd2


# ---------------------------------------------------------------------------
# CPX interval parsing
# ---------------------------------------------------------------------------

def parse_cpx_interval(interval_str: str) -> Tuple[str, str, int, int]:
    """
    Parse a CPX_INTERVALS token like ``DUP_chrY:3125606-3125667``.

    Returns
    -------
    (cnv_type, chrom, start, end)
    """
    cnv_type, region = interval_str.split("_", 1)
    chrom, coords = region.split(":")
    start_s, end_s = coords.split("-")
    return cnv_type, chrom, int(start_s), int(end_s)


def parse_source(source_str: str) -> Tuple[str, str, int, int]:
    """
    Parse a SOURCE field like ``DUP_chrY:3125606-3125667``.

    Returns same tuple as :func:`parse_cpx_interval`.
    """
    return parse_cpx_interval(source_str)


def extract_cnv_segments(
    record: pysam.VariantRecord,
    min_size: int = 0,
) -> List[Tuple[str, str, int, int, str]]:
    """
    Extract DEL/DUP CNV segments from CPX_INTERVALS and SOURCE.

    Returns list of (cnv_type, chrom, start, end, variant_id).
    """
    segments = []
    vid = record.id

    cpx_intervals = record.info.get(CPX_INTERVALS, ())
    if cpx_intervals:
        for tok in cpx_intervals:
            if "DEL_" in tok or "DUP_" in tok:
                cnv_type, chrom, start, end = parse_cpx_interval(tok)
                if end - start >= min_size:
                    segments.append((cnv_type, chrom, start, end, vid))

    source = record.info.get(SOURCE, None)
    if source and ("DEL_" in source or "DUP_" in source):
        cnv_type, chrom, start, end = parse_cpx_interval(source)
        if end - start >= min_size:
            segments.append((cnv_type, chrom, start, end, vid))

    return segments


# ---------------------------------------------------------------------------
# Genomic overlap helpers
# ---------------------------------------------------------------------------

def get_overlap(
    chrom1: str, start1: int, end1: int,
    chrom2: str, start2: int, end2: int,
) -> int:
    if chrom1 != chrom2:
        return 0
    if not (start1 <= end2 and start2 <= end1):
        return 0
    return min(end1, end2) - max(start1, start2)


def has_reciprocal_overlap(
    chrom1: str, start1: int, end1: int,
    chrom2: str, start2: int, end2: int,
    min_ro: float,
) -> bool:
    overlap = get_overlap(chrom1, start1, end1, chrom2, start2, end2)
    size_1 = end1 - start1
    size_2 = end2 - start2
    return overlap / max(size_1, size_2, 1) >= min_ro


# ---------------------------------------------------------------------------
# Sex-chromosome ploidy helpers
# ---------------------------------------------------------------------------

def load_ped_sex(
    ped_path: str,
) -> Tuple[Set[str], Set[str]]:
    """
    Read a PED file and return (male_samples, female_samples).

    PED sex column: 1 = male, 2 = female.
    """
    males: Set[str] = set()
    females: Set[str] = set()
    with open(ped_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            tokens = line.strip().split("\t")
            if len(tokens) < 5:
                continue
            sample = tokens[1]
            sex = tokens[4]
            if sex == "1":
                males.add(sample)
            elif sex == "2":
                females.add(sample)
    return males, females
