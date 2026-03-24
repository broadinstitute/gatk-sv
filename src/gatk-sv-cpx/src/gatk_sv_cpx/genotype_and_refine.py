"""
``gatk-sv-cpx genotype-and-refine`` — apply depth genotyping cutoffs and
PE/RD refinement in a single VCF pass.

Replaces:

* GCV depth-based CPX genotyping (``process_posthoc_cpx_depth_regenotyping.py``)
* RefCV VCF revision (``resolve_cpx_variants_revise_vcf.py``)
* CPX–CNV redundancy resolution (``resolve_cpx_cnv_redundancies.py``)

Given the evidence JSON from ``evaluate-evidence`` and a resolved VCF, this
module:

1. Reads per-sample PE + depth evidence from the JSON bundle.
2. Applies no-call rules identical to those in the legacy ``ReviseVcf`` task.
3. Marks variants ``UNRESOLVED`` when too many carriers lack support.
4. Converts ``INS`` + ``INV`` source variants to ``CPX`` (insertion-type fix).
5. Moves ``SOURCE`` to ``CPX_INTERVALS`` for insertion-type CPX.
"""

from __future__ import annotations

import json
import logging
from typing import Any, Dict, Optional, Set

import pysam

from gatk_sv_cpx.constants import (
    CPX_INTERVALS,
    CPX_TYPE_KEY,
    INSERTION_TYPE_CPX,
    SOURCE,
    SVTYPE,
)
from gatk_sv_cpx.vcf_utils import (
    ensure_header_lines,
    get_svtype,
)

log = logging.getLogger(__name__)

# Depth-support strings that indicate insufficient depth
_LACK_DEPTH_STRINGS: Set[str] = {
    "NA",
    "lack_depth",
    "lack_depth,lack_depth",
    "lack_depth,depth",
    "depth,lack_depth",
}


# ═══════════════════════════════════════════════════════════════════════════
# Evidence application (port of resolve_cpx_variants_revise_vcf.py)
# ═══════════════════════════════════════════════════════════════════════════

def _apply_cpx_evidence(
    record: pysam.VariantRecord,
    sample_evidence: Dict[str, Dict[str, str]],
) -> bool:
    """
    Apply CPX evidence to a record.  Sets no-call for unsupported samples.

    Returns ``True`` if the record should be marked UNRESOLVED (>= 50 %
    of carriers have only partial PE with no depth support).
    """
    n_unresolve = 0
    n_total = len(sample_evidence)
    for sample, ev in sample_evidence.items():
        if sample not in record.samples:
            continue
        pe = ev["pe"]
        depth = ev["depth"]

        if pe == "no_PE" and depth in _LACK_DEPTH_STRINGS:
            record.samples[sample]["GT"] = (None, None)
        elif pe == "low_PE" and depth in _LACK_DEPTH_STRINGS:
            record.samples[sample]["GT"] = (None, None)
        elif pe == "partial_PE" and depth in _LACK_DEPTH_STRINGS:
            record.samples[sample]["GT"] = (None, None)
            n_unresolve += 1
        elif pe == "high_PE" and depth in {"lack_depth", "lack_depth,lack_depth",
                                            "lack_depth,depth", "depth,lack_depth"}:
            record.samples[sample]["GT"] = (None, None)

    if n_total > 0 and n_unresolve / n_total >= 0.5:
        return True
    return False


def _apply_ctx_evidence(
    record: pysam.VariantRecord,
    sample_evidence: Dict[str, Dict[str, str]],
) -> bool:
    """
    Apply CTX evidence.  No-calls samples with no/low/partial PE.

    Returns ``True`` if the record should be marked UNRESOLVED
    (when the only carrier is ``NA``).
    """
    if "NA" in sample_evidence:
        return True

    for sample, ev in sample_evidence.items():
        if sample not in record.samples:
            continue
        pe = ev["pe"]
        if pe in ("no_PE", "low_PE", "partial_PE"):
            record.samples[sample]["GT"] = (None, None)

    return False


def _fix_insertion_type_cpx(record: pysam.VariantRecord) -> None:
    """
    Post-refinement fixups for insertion-type CPX events:

    - ``dDUP`` / ``dDUP_iDEL`` / ``INS_iDEL``: fold SOURCE into CPX_INTERVALS
    - ``INS`` + ``INV`` in SOURCE: convert to CPX ``dDUP`` or ``dDUP_iDEL``
    """
    svtype = get_svtype(record)
    cpx_type = record.info.get(CPX_TYPE_KEY, "")

    if svtype == "CPX" and cpx_type in INSERTION_TYPE_CPX:
        # Fold SOURCE into CPX_INTERVALS
        source = record.info.get(SOURCE, "")
        if source:
            intervals = list(record.info.get(CPX_INTERVALS, ()))
            if cpx_type == "INS_iDEL":
                intervals.append(source)
            record.info[CPX_INTERVALS] = tuple(intervals) if intervals else (source,)
            del record.info[SOURCE]

    elif svtype == "INS" and SOURCE in record.info and "INV" in record.info[SOURCE]:
        # Convert insertion with INV source to CPX
        record.info[SVTYPE] = "CPX"
        record.alts = ("<CPX>",)
        source = record.info[SOURCE]
        del_section = record.stop - record.pos

        if del_section < 50:
            record.info[CPX_TYPE_KEY] = "dDUP"
            dup_source = source.replace("INV", "DUP")
            record.info[CPX_INTERVALS] = (dup_source, source)
        else:
            record.info[CPX_TYPE_KEY] = "dDUP_iDEL"
            dup_source = source.replace("INV", "DUP")
            del_interval = f"DEL_{record.chrom}:{record.pos}-{record.stop}"
            record.info[CPX_INTERVALS] = (dup_source, source, del_interval)

        del record.info[SOURCE]


# ═══════════════════════════════════════════════════════════════════════════
# Unresolved SV ID detection
# (port of the unresolved_svids output from GenerateCpxReviewScript)
# ═══════════════════════════════════════════════════════════════════════════

def _load_unresolved_ids(path: Optional[str]) -> Set[str]:
    if not path:
        return set()
    ids: Set[str] = set()
    with open(path) as f:
        for line in f:
            ids.add(line.strip())
    return ids


# ═══════════════════════════════════════════════════════════════════════════
# Top-level entry point
# ═══════════════════════════════════════════════════════════════════════════

def genotype_and_refine(
    input_vcf: str,
    evidence_json: str,
    output_vcf: str,
    unresolved_ids_file: Optional[str] = None,
) -> None:
    """
    Apply depth genotyping and PE/RD refinement to a resolved VCF.

    Parameters
    ----------
    input_vcf
        Resolved VCF (output of ``gatk-sv-cpx resolve``).
    evidence_json
        Evidence bundle (output of ``gatk-sv-cpx evaluate-evidence``).
    output_vcf
        Path for the finalized output VCF.
    unresolved_ids_file
        Optional file listing variant IDs to force-mark as UNRESOLVED.
    """
    # ── Load evidence ──────────────────────────────────────────────────
    with open(evidence_json) as f:
        evidence: Dict[str, Dict[str, Dict[str, Any]]] = json.load(f)
    log.info("Loaded evidence for %d variants", len(evidence))

    unresolved_ids = _load_unresolved_ids(unresolved_ids_file)

    # ── Stream and revise VCF ──────────────────────────────────────────
    vcf_in = pysam.VariantFile(input_vcf)
    header = vcf_in.header.copy()
    ensure_header_lines(header)

    # Ensure UNRESOLVED filter exists in header
    if "UNRESOLVED" not in header.filters:
        header.add_line('##FILTER=<ID=UNRESOLVED,Description="Variant is unresolved.">')

    vcf_out = pysam.VariantFile(output_vcf, "wz", header=header)
    n_unresolved = 0
    n_nocalled = 0
    n_ins_converted = 0

    for record in vcf_in:
        vid = record.id

        # Force-unresolved
        if vid in unresolved_ids:
            if "PASS" in record.filter:
                record.filter.clear()
            record.filter.add("UNRESOLVED")
            n_unresolved += 1

        # Apply CPX evidence
        elif vid in evidence:
            sample_ev = evidence[vid]
            svtype = get_svtype(record)

            # Determine if CPX or CTX
            first_ev = next(iter(sample_ev.values()), {})
            ev_svtype = first_ev.get("svtype", svtype)

            if ev_svtype == "CTX":
                mark_unresolved = _apply_ctx_evidence(record, sample_ev)
            else:
                mark_unresolved = _apply_cpx_evidence(record, sample_ev)
                n_nocalled += sum(
                    1 for s in sample_ev
                    if s in record.samples and record.samples[s]["GT"] == (None, None)
                )

            if mark_unresolved:
                if "PASS" in record.filter:
                    record.filter.clear()
                record.filter.add("UNRESOLVED")
                n_unresolved += 1

        # Insertion-type CPX fixups (runs on every record)
        svtype_before = get_svtype(record)
        _fix_insertion_type_cpx(record)
        if svtype_before == "INS" and get_svtype(record) == "CPX":
            n_ins_converted += 1

        vcf_out.write(record)

    vcf_in.close()
    vcf_out.close()

    # Index output
    pysam.tabix_index(output_vcf, preset="vcf", force=True)

    log.info(
        "genotype-and-refine: %d unresolved, %d no-called GTs, %d INS→CPX",
        n_unresolved,
        n_nocalled,
        n_ins_converted,
    )
