"""
``gatk-sv-cpx evaluate-evidence`` — gather PE and RD evidence for CPX/CTX
variants in a single pass.

Replaces the fragmented evidence gathering spread across
**GenotypeComplexVariants** (RD) and **RefineComplexVariants** (PE + depth).

The module reads a resolved VCF and the relevant evidence files (bincov,
PE metrics, depth DEL/DUP beds) and emits a JSON evidence bundle consumed
by ``genotype-and-refine``.
"""

from __future__ import annotations

import gzip
import json
import logging
import os
from typing import Any, Dict, List, Optional, Tuple

import pysam

from gatk_sv_cpx.constants import (
    CPX_TYPE_KEY,
    DEFAULT_DEPTH_FRACTION_THRESHOLD,
    DEFAULT_MIN_PE_CPX,
    DEFAULT_MIN_PE_CTX,
    SOURCE,
)
from gatk_sv_cpx.vcf_utils import (
    get_called_samples,
    get_svtype,
)

log = logging.getLogger(__name__)


# ═══════════════════════════════════════════════════════════════════════════
# PE evidence collection
# (port of calculate_cpx_evidences.py, calculate_ctx_evidences.py, and
#  reformat_CPX_bed_and_generate_script.py)
# ═══════════════════════════════════════════════════════════════════════════

def _read_pe_metrics(pe_path: str) -> Dict[str, Dict[str, List[List]]]:
    """
    Read tabix-indexed PE metrics file.

    Returns nested dict: ``{variant_id: {sample: [[fwd_count, strand, pos], [rev_count, strand, pos]]}}``
    """
    out: Dict[str, Dict[str, List[List]]] = {}
    opener = gzip.open if pe_path.endswith(".gz") else open
    with opener(pe_path, "rt") as f:
        for line in f:
            tokens = line.strip().split()
            if len(tokens) < 5:
                continue
            svid = tokens[4]
            sample = tokens[3]
            if svid not in out:
                out[svid] = {}
            if sample not in out[svid]:
                out[svid][sample] = [[], []]
            count = int(tokens[0])
            strand = tokens[1]
            if strand == "+":
                out[svid][sample][0] = [count, strand, tokens[2]]
            elif strand == "-":
                out[svid][sample][1] = [count, strand, tokens[2]]
    return out


def classify_pe_support(
    fwd_rev: List[List],
    min_pe: int,
) -> str:
    """
    Classify PE support for a single sample.

    Categories: ``no_PE``, ``partial_PE``, ``low_PE``, ``high_PE``.
    """
    fwd = fwd_rev[0]
    rev = fwd_rev[1]
    fwd_count = fwd[0] if fwd else 0
    rev_count = rev[0] if rev else 0

    if fwd_count <= 0 and rev_count <= 0:
        return "no_PE"
    if fwd_count > 0 and rev_count > 0:
        if fwd_count >= min_pe and rev_count >= min_pe:
            return "high_PE"
        return "low_PE"
    return "partial_PE"


# ═══════════════════════════════════════════════════════════════════════════
# Depth evidence collection
# (port of generate_cnv_segment_from_cpx.py + depth-support logic)
# ═══════════════════════════════════════════════════════════════════════════

def _read_depth_bed(bed_path: str) -> Dict[str, Dict[str, List[float]]]:
    """
    Read a depth support BED (output of SeekDepthSupportForCPX or similar).

    Returns ``{variant_id: {sample: [frac1, frac2, ...]}}``.
    """
    out: Dict[str, Dict[str, List[float]]] = {}
    opener = gzip.open if bed_path.endswith(".gz") else open
    with opener(bed_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            tokens = line.strip().split()
            if len(tokens) < 8:
                continue
            svid = tokens[4]
            sample = tokens[5]
            frac = float(tokens[7])
            if svid not in out:
                out[svid] = {}
            if sample not in out[svid]:
                out[svid][sample] = []
            out[svid][sample].append(frac)
    return out


def classify_depth_support(
    fractions: List[float],
    threshold: float = DEFAULT_DEPTH_FRACTION_THRESHOLD,
) -> str:
    """
    Classify depth support for a single sample.

    Returns a comma-separated string of ``depth`` / ``lack_depth`` per
    segment, or ``NA`` if no data.
    """
    if not fractions:
        return "NA"
    labels = ["depth" if f > threshold else "lack_depth" for f in fractions]
    return ",".join(labels)


# ═══════════════════════════════════════════════════════════════════════════
# Combined evidence evaluation
# ═══════════════════════════════════════════════════════════════════════════

def _collect_variant_carriers(
    vcf_path: str,
) -> Tuple[
    Dict[str, List[str]],  # svid → [samples]
    Dict[str, str],        # svid → svtype
    Dict[str, str],        # svid → cpx_type
]:
    """Scan VCF to collect carrier samples and type info for CPX/CTX."""
    svid_samples: Dict[str, List[str]] = {}
    svid_svtype: Dict[str, str] = {}
    svid_cpx_type: Dict[str, str] = {}

    vcf = pysam.VariantFile(vcf_path)
    for record in vcf:
        svtype = get_svtype(record)
        if svtype not in ("CPX", "CTX", "INS"):
            continue
        # INS with INV in SOURCE get treated as CPX during refinement
        if svtype == "INS":
            source = record.info.get(SOURCE, "")
            if "INV" not in source:
                continue
        vid = record.id
        carriers = get_called_samples(record)
        if not carriers:
            continue
        svid_samples[vid] = carriers
        svid_svtype[vid] = svtype
        svid_cpx_type[vid] = record.info.get(CPX_TYPE_KEY, "")
    vcf.close()
    return svid_samples, svid_svtype, svid_cpx_type


def evaluate_evidence(
    input_vcf: str,
    output_json: str,
    pe_metrics: Optional[str] = None,
    depth_support: Optional[str] = None,
    min_pe_cpx: int = DEFAULT_MIN_PE_CPX,
    min_pe_ctx: int = DEFAULT_MIN_PE_CTX,
    depth_threshold: float = DEFAULT_DEPTH_FRACTION_THRESHOLD,
) -> None:
    """
    Gather PE and RD evidence for CPX/CTX variants and write a JSON
    evidence bundle.

    Parameters
    ----------
    input_vcf
        Resolved VCF (output of ``gatk-sv-cpx resolve``).
    output_json
        Output evidence file (JSON).
    pe_metrics
        Concatenated/merged PE metrics file (optional).
    depth_support
        Depth support BED file (optional).
    min_pe_cpx
        Minimum PE read count for CPX ``high_PE`` classification.
    min_pe_ctx
        Minimum PE read count for CTX ``high_PE`` classification.
    depth_threshold
        Fraction threshold for depth support.
    """
    # ── Collect carriers from VCF ──────────────────────────────────────
    svid_samples, svid_svtype, svid_cpx_type = _collect_variant_carriers(input_vcf)
    log.info("Collected %d CPX/CTX/INS variants with carriers", len(svid_samples))

    # ── Read PE evidence ───────────────────────────────────────────────
    pe_data: Dict[str, Dict[str, List[List]]] = {}
    if pe_metrics and os.path.isfile(pe_metrics):
        pe_data = _read_pe_metrics(pe_metrics)
        log.info("Read PE metrics for %d variants", len(pe_data))

    # ── Read depth evidence ────────────────────────────────────────────
    depth_data: Dict[str, Dict[str, List[float]]] = {}
    if depth_support and os.path.isfile(depth_support):
        depth_data = _read_depth_bed(depth_support)
        log.info("Read depth support for %d variants", len(depth_data))

    # ── Build per-sample evidence records ──────────────────────────────
    evidence: Dict[str, Dict[str, Dict[str, Any]]] = {}  # svid → sample → {pe, depth}

    for svid, samples in svid_samples.items():
        svtype = svid_svtype[svid]
        min_pe = min_pe_ctx if svtype == "CTX" else min_pe_cpx
        evidence[svid] = {}

        for sample in samples:
            # PE classification
            if svid in pe_data and sample in pe_data[svid]:
                pe_class = classify_pe_support(pe_data[svid][sample], min_pe)
            else:
                pe_class = "no_PE"

            # Depth classification
            if svid in depth_data and sample in depth_data[svid]:
                depth_class = classify_depth_support(
                    depth_data[svid][sample], depth_threshold,
                )
            else:
                depth_class = "NA"

            evidence[svid][sample] = {
                "pe": pe_class,
                "depth": depth_class,
                "svtype": svtype,
                "cpx_type": svid_cpx_type.get(svid, ""),
            }

    # ── Write output ───────────────────────────────────────────────────
    with open(output_json, "w") as f:
        json.dump(evidence, f, indent=1)

    log.info(
        "evaluate-evidence: wrote evidence for %d variants to %s",
        len(evidence),
        output_json,
    )
