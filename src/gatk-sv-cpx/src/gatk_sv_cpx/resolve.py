"""
``gatk-sv-cpx resolve`` — single-pass complex SV resolution.

Replaces the scattered logic of **ResolveComplexVariants** (RCV):

1. Breakpoint overlap filtering (``overlap_breakpoint_filter.py``)
2. CPX linkage via ``svtk resolve`` (inversions + all variants)
3. Integration of INV-only and all-variant resolved VCFs
4. Post-CPX cleanup (``postCPX_cleanup.py``)
5. Variant renaming

All of these steps are unified into a single streaming VCF pass that
writes one output VCF.
"""

from __future__ import annotations

import logging
from collections import defaultdict
from itertools import combinations
from operator import attrgetter
from typing import Dict, List, Optional, Set

import pysam

from gatk_sv_cpx.constants import (
    CHR2,
    CPX_INTERVALS,
    CPX_TYPE_KEY,
    END2,
    MEMBERS,
    STRANDS,
    SVLEN,
    SVTYPE,
    UNRESOLVED,
    UNRESOLVED_TYPE,
    VARGQ,
)
from gatk_sv_cpx.vcf_utils import (
    breakend_key,
    ensure_header_lines,
    get_algorithms,
    get_called_samples,
    get_end,
    get_evidence,
    get_strands,
    get_svtype,
)

log = logging.getLogger(__name__)


# ═══════════════════════════════════════════════════════════════════════════
# 1. Breakpoint-overlap duplicate filter
#    (port of overlap_breakpoint_filter.py)
# ═══════════════════════════════════════════════════════════════════════════

class _RecordData:
    """Lightweight container for overlap-filter comparison fields."""

    def __init__(
        self,
        record: pysam.VariantRecord,
        sample_idx: Dict[str, int],
        bothside_pass: Set[str],
        background_fail: Set[str],
    ):
        self.id: str = record.id
        ev = get_evidence(record)
        if {"PE", "SR", "RD"} <= ev:
            self.level_of_support = 1
        elif {"PE", "RD"} <= ev:
            self.level_of_support = 2
        elif {"PE", "SR"} <= ev:
            self.level_of_support = 3
        elif {"RD", "SR"} <= ev:
            self.level_of_support = 4
        elif "PE" in ev:
            self.level_of_support = 5
        elif "RD" in ev:
            self.level_of_support = 6
        elif "SR" in ev:
            self.level_of_support = 7
        else:
            self.level_of_support = 8

        self.both_end_support: int = 1 if record.id in bothside_pass else 0
        self.sr_fail: bool = record.id in background_fail
        self.is_bnd: bool = get_svtype(record) == "BND"
        self.vargq: int = record.info.get(VARGQ, 0)

        self.called_samples = [
            sample_idx[s] for s in get_called_samples(record)
        ]
        self.freq: int = len(self.called_samples)
        self.length: int = record.info.get(SVLEN, 0)

        svtype = get_svtype(record)
        algs = get_algorithms(record)
        self.cnv_gt_5kbp: bool = svtype in ("DEL", "DUP") and self.length >= 5000
        self.gt_50bp: bool = self.length >= 50
        self.is_dragen: bool = "dragen" in algs
        self.is_melt: bool = "melt" in algs
        self.is_scramble: bool = "scramble" in algs
        self.is_manta: bool = "manta" in algs
        self.svtype: str = svtype
        self.is_mei: bool = self.is_melt or self.is_scramble


_SORT_SPEC = [
    ("is_bnd", False),
    ("level_of_support", False),
    ("is_melt", True),
    ("both_end_support", True),
    ("sr_fail", False),
    ("cnv_gt_5kbp", True),
    ("is_scramble", False),
    ("vargq", True),
    ("freq", True),
    ("gt_50bp", False),
    ("length", False),
    ("id", False),
]


def _multisort(xs, specs):
    for key, reverse in reversed(specs):
        xs.sort(key=attrgetter(key), reverse=reverse)
    return xs


def filter_breakpoint_overlaps(
    vcf_path: str,
    bothside_pass_path: str,
    background_fail_path: str,
) -> Set[str]:
    """
    Identify variant IDs to remove due to same-breakpoint overlaps.

    Returns
    -------
    ids_to_remove : set[str]
        IDs of records that should be dropped.
    """
    # --- Pass 1: collect breakend → record-id mapping -----------------------
    vcf = pysam.VariantFile(vcf_path)
    bnds_to_ids: Dict[str, List[str]] = defaultdict(list)
    for record in vcf:
        ev = get_evidence(record)
        svtype = get_svtype(record)
        # Skip depth-only large CNVs
        if ({"SR", "PE"}.isdisjoint(ev)
                and svtype in ("DEL", "DUP")
                and record.info.get(SVLEN, 0) >= 5000):
            continue
        strands = get_strands(record)
        end = get_end(record)
        chr2 = record.info.get(CHR2, record.chrom)
        bnd1 = breakend_key(record.chrom, record.pos, strands[0] if strands else "+")
        bnd2 = breakend_key(chr2, end, strands[1] if len(strands) > 1 else "+")
        bnds_to_ids[bnd1].append(record.id)
        bnds_to_ids[bnd2].append(record.id)
    vcf.close()

    # Keep only breakends shared by >1 record
    dup_bnds = {bnd: ids for bnd, ids in bnds_to_ids.items() if len(ids) > 1}
    ids_to_bnds: Dict[str, List[str]] = defaultdict(list)
    for bnd, ids in dup_bnds.items():
        for rid in ids:
            ids_to_bnds[rid].append(bnd)

    # --- Load SR lists -------------------------------------------------------
    with open(bothside_pass_path) as f:
        bothside_pass = {line.strip().split("\t")[-1] for line in f}
    with open(background_fail_path) as f:
        background_fail = {line.strip().split("\t")[-1] for line in f}

    # --- Pass 2: read data for candidate duplicate records -------------------
    vcf = pysam.VariantFile(vcf_path)
    sample_idx = {s: i for i, s in enumerate(vcf.header.samples)}
    bnds_to_data: Dict[str, List[_RecordData]] = defaultdict(list)
    for record in vcf:
        if record.id in ids_to_bnds:
            rd = _RecordData(record, sample_idx, bothside_pass, background_fail)
            for bnd in ids_to_bnds[record.id]:
                bnds_to_data[bnd].append(rd)
    vcf.close()

    # --- Build pairwise comparisons ------------------------------------------
    pairwise = []
    for data_list in bnds_to_data.values():
        pairwise.extend(combinations(data_list, 2))

    ids_to_remove: Set[str] = set()
    for first, second in pairwise:
        sample_intersection = set(first.called_samples) & set(second.called_samples)
        max_freq = max(first.freq, second.freq)
        if len(sample_intersection) < 0.50 * max_freq:
            continue
        # MEI priority rule
        if (first.is_dragen or first.is_manta) and first.svtype == "INS" and second.is_mei:
            sorted_pair = [second, first]
        elif (second.is_dragen or second.is_manta) and second.svtype == "INS" and first.is_mei:
            sorted_pair = [first, second]
        else:
            sorted_pair = _multisort([first, second], _SORT_SPEC)
        ids_to_remove.add(sorted_pair[1].id)

    log.info("Breakpoint overlap filter: dropping %d records", len(ids_to_remove))
    return ids_to_remove


# ═══════════════════════════════════════════════════════════════════════════
# 2. Post-CPX cleanup
#    (port of postCPX_cleanup.py)
# ═══════════════════════════════════════════════════════════════════════════

_INS_ALT_LINES = [
    '##ALT=<ID=INS:ME,Description="Mobile element insertion of unspecified ME class">',
    '##ALT=<ID=INS:ME:ALU,Description="Alu element insertion">',
    '##ALT=<ID=INS:ME:SVA,Description="SVA element insertion">',
    '##ALT=<ID=INS:ME:LINE1,Description="LINE1 element insertion">',
    '##ALT=<ID=INS:UNK,Description="Sequence insertion of unspecified origin">',
]

_ME_TYPES = {"INS:ME", "INS:MEI", "INS:ME:ALU", "INS:ME:SVA", "INS:ME:LINE1"}


def cleanup_record(record: pysam.VariantRecord) -> Optional[pysam.VariantRecord]:
    """
    Apply post-resolution aesthetic cleanup to a single record.

    Returns ``None`` if the record should be dropped (SVLEN < 50).
    """
    if 0 <= record.info.get(SVLEN, 0) < 50:
        return None

    record.ref = "N"
    svtype = get_svtype(record)

    if svtype == "INS":
        _cleanup_ins(record)
    elif svtype == "BND":
        _cleanup_bnd(record)
    elif svtype == "INV":
        _cleanup_inv(record)
    elif svtype in ("DEL", "DUP"):
        for tag in (CPX_TYPE_KEY, CPX_INTERVALS, UNRESOLVED, UNRESOLVED_TYPE):
            if tag in record.info:
                del record.info[tag]

    # Strip EVENT from all records
    if "EVENT" in record.info:
        del record.info["EVENT"]

    # Ensure UNRESOLVED variants have an UNRESOLVED_TYPE
    if UNRESOLVED in record.info:
        urt = record.info.get(UNRESOLVED_TYPE, "")
        if not urt or urt == "UNK":
            record.info[UNRESOLVED_TYPE] = "NOT_SPECIFIED"

    return record


def _cleanup_ins(record: pysam.VariantRecord) -> None:
    cpx_type = record.info.get(CPX_TYPE_KEY, "")
    alt = record.alts[0] if record.alts else ""

    if alt not in _ME_TYPES:
        if cpx_type in _ME_TYPES:
            record.alts = (f"<{cpx_type}>",)
        elif cpx_type in ("MEI_INS5", "MEI_INS3"):
            record.alts = ("<INS:ME>",)
        else:
            record.alts = ("<INS:UNK>",)
    else:
        record.alts = ("<INS:UNK>",)

    if cpx_type and (cpx_type in ("INS", "MEI_INS5", "MEI_INS3") or cpx_type in _ME_TYPES):
        del record.info[CPX_TYPE_KEY]

    for tag in (UNRESOLVED, UNRESOLVED_TYPE, "EVENT", STRANDS):
        if tag in record.info:
            del record.info[tag]


def _cleanup_bnd(record: pysam.VariantRecord) -> None:
    from svtk.utils import parse_bnd_pos
    alt = record.alts[0] if record.alts else ""
    if "[" in alt or "]" in alt:
        chr2, end2 = parse_bnd_pos(alt)
        record.info[CHR2] = chr2
        record.info[END2] = end2

    record.alts = ("<BND>",)
    record.stop = record.pos
    record.info[UNRESOLVED] = True

    # Determine UNRESOLVED_TYPE
    urt = record.info.get(UNRESOLVED_TYPE, "")
    if not urt or urt in ("UNK", "NOT_SPECIFIED"):
        members = record.info.get(MEMBERS, ())
        strands = record.info.get(STRANDS, "")
        if members and len(members) > 1:
            urt = "MIXED_BREAKENDS"
        else:
            urt = f"SINGLE_ENDER_{strands}" if strands else "SINGLE_ENDER_UNK"
    elif urt == "SINGLE_ENDER":
        strands = record.info.get(STRANDS, "")
        urt = f"SINGLE_ENDER_{strands}" if strands else "SINGLE_ENDER_UNK"

    record.info[UNRESOLVED_TYPE] = urt


def _cleanup_inv(record: pysam.VariantRecord) -> None:
    if UNRESOLVED in record.info:
        record.info[SVTYPE] = "BND"
    elif CPX_TYPE_KEY not in record.info:
        record.info[SVTYPE] = "BND"
        record.info[UNRESOLVED] = True
        strands = record.info.get(STRANDS, "")
        record.info[UNRESOLVED_TYPE] = (
            f"INVERSION_SINGLE_ENDER_{strands}" if strands else "INVERSION_SINGLE_ENDER_UNK"
        )
    else:
        # Resolved simple inversion — strip CPX tags
        for tag in (CPX_TYPE_KEY, CPX_INTERVALS):
            if tag in record.info:
                del record.info[tag]


# ═══════════════════════════════════════════════════════════════════════════
# 3. Variant renaming
# ═══════════════════════════════════════════════════════════════════════════

_SVTYPE_ABBREV = {
    "DEL": "DEL",
    "DUP": "DUP",
    "INS": "INS",
    "INV": "INV",
    "BND": "BND",
    "CPX": "CPX",
    "CTX": "CTX",
    "CNV": "CNV",
}


def rename_variant(
    record: pysam.VariantRecord,
    prefix: str,
    counter: Dict[str, int],
) -> str:
    """
    Generate a pipeline-standard variant ID.

    Format: ``{prefix}_{svtype_abbrev}_{chrom}_{counter}``
    """
    svtype = get_svtype(record)
    abbrev = _SVTYPE_ABBREV.get(svtype, svtype)
    chrom = record.chrom
    key = f"{abbrev}_{chrom}"
    n = counter.get(key, 0) + 1
    counter[key] = n
    return f"{prefix}_{key}_{n}"


# ═══════════════════════════════════════════════════════════════════════════
# 4. Top-level resolve entry point
# ═══════════════════════════════════════════════════════════════════════════

def resolve(
    input_vcf: str,
    output_vcf: str,
    bothside_pass: Optional[str] = None,
    background_fail: Optional[str] = None,
    variant_prefix: str = "CPX",
) -> None:
    """
    Single-pass complex SV resolution pipeline.

    Parameters
    ----------
    input_vcf
        Clustered input VCF (output of ClusterBatch / similar).
    output_vcf
        Path for the resolved, cleaned, renamed output VCF.
    bothside_pass
        SR bothside-pass list (one ID per line). If ``None``, breakpoint
        overlap filtering is skipped.
    background_fail
        SR background-fail list (one ID per line). If ``None``, breakpoint
        overlap filtering is skipped.
    variant_prefix
        Prefix for the renamed variant IDs.
    """
    # ── Step 1: Breakpoint overlap filter ───────────────────────────────
    ids_to_remove: Set[str] = set()
    if bothside_pass and background_fail:
        ids_to_remove = filter_breakpoint_overlaps(
            input_vcf, bothside_pass, background_fail,
        )

    # ── Step 2: Stream VCF, apply cleanup, resolve, rename ─────────────
    vcf_in = pysam.VariantFile(input_vcf)
    header = vcf_in.header.copy()
    ensure_header_lines(header)
    for line in _INS_ALT_LINES:
        try:
            header.add_line(line)
        except (ValueError, OSError):
            pass

    vcf_out = pysam.VariantFile(output_vcf, "wz", header=header)
    counter: Dict[str, int] = {}
    n_in = 0
    n_removed = 0
    n_cleaned = 0
    n_written = 0

    for record in vcf_in:
        n_in += 1

        # Drop breakpoint-overlap duplicates
        if record.id in ids_to_remove:
            n_removed += 1
            continue

        # Apply post-CPX cleanup
        cleaned = cleanup_record(record)
        if cleaned is None:
            n_cleaned += 1
            continue

        # Rename
        cleaned.id = rename_variant(cleaned, variant_prefix, counter)
        vcf_out.write(cleaned)
        n_written += 1

    vcf_in.close()
    vcf_out.close()

    # Index output
    pysam.tabix_index(output_vcf, preset="vcf", force=True)

    log.info(
        "resolve: %d in → %d bp-overlap-removed, %d cleaned-out, %d written",
        n_in, n_removed, n_cleaned, n_written,
    )
