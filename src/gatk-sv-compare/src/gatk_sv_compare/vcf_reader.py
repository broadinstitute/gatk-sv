"""Streaming VCF extraction utilities for Phase 3 aggregation."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Sequence, Tuple

import numpy as np
import pysam

from .dimensions import GENOMIC_CONTEXTS, normalize_algorithms, normalize_evidence_bucket, normalize_svtype
from .vcf_format import filter_values, has_precomputed_counts, safe_info_get

_GENOMIC_CONTEXT_TOKENS = {
    "segdup": ("segdup", "segmentaldup", "segmentalduplication", "segduplication"),
    "simple_repeat": ("simplerepeat", "simple_repeat", "tandemrepeat", "trf"),
    "repeatmasker": ("repeatmasker", "repeat_masker", "rmsk"),
}
# TODO: some CPX subtypes are not span-based
_SPAN_BASED_CONTEXT_SVTYPES = frozenset({"DEL", "DUP", "CNV", "CPX", "INV"})

_CONCORDANCE_FIELDS = (
    "VAR_SENSITIVITY",
    "VAR_PPV",
    "VAR_SPECIFICITY",
    "HET_SENSITIVITY",
    "HET_PPV",
    "HOMVAR_SENSITIVITY",
    "HOMVAR_PPV",
    "GENOTYPE_CONCORDANCE",
    "NON_REF_GENOTYPE_CONCORDANCE",
    "CNV_CONCORDANCE",
)


@dataclass
class SiteRecord:
    """Flat representation of one VCF record."""

    variant_id: str
    contig: str
    pos: int
    end: int
    svtype: str
    svlen: Optional[int]
    alt_allele: str
    filter: str
    af: Optional[float]
    ac: Optional[int]
    an: Optional[int]
    n_bi_genos: int
    n_hom_ref: int
    n_het: int
    n_hom_alt: int
    n_no_call: int
    cn_nonref_freq: Optional[float]
    status: Optional[str]
    truth_vid: Optional[str]
    genomic_context: str
    algorithms: tuple[str, ...]
    evidence_bucket: str
    gq_array: Optional[np.ndarray]
    concordance_metrics: Optional[Dict[str, float]]

    def to_dict(self) -> Dict[str, object]:
        """Convert to a dict suitable for DataFrame construction."""
        return asdict(self)


def _record_identifier(record: pysam.VariantRecord) -> str:
    return record.id or f"{record.contig}:{record.pos}"


def _info_value(record: pysam.VariantRecord, key: str) -> object:
    value = safe_info_get(record, key)
    if isinstance(value, tuple):
        return value[0] if value else None
    return value


def _algorithms(record: pysam.VariantRecord) -> tuple[str, ...]:
    return normalize_algorithms(safe_info_get(record, "ALGORITHMS"))


def _evidence_bucket(record: pysam.VariantRecord) -> str:
    return normalize_evidence_bucket(safe_info_get(record, "EVIDENCE"))


def _int_or_none(value: object) -> Optional[int]:
    if value in (None, "."):
        return None
    return int(value)


def _float_or_none(value: object) -> Optional[float]:
    if value in (None, "."):
        return None
    return float(value)


def _normalize_info_key(key: str) -> str:
    return key.lower().replace("_", "").replace("-", "")


def _context_metrics(
    record: pysam.VariantRecord,
    normalized_keys: Dict[str, str],
    context: str,
) -> Tuple[Optional[float], Optional[int], bool]:
    tokens = _GENOMIC_CONTEXT_TOKENS.get(context, (_normalize_info_key(context),))
    overlap_fraction: Optional[float] = None
    end_overlap_count: Optional[int] = None
    legacy_match = False

    for normalized_key, raw_key in normalized_keys.items():
        if not any(token in normalized_key for token in tokens):
            continue
        if "overlapfrac" in normalized_key:
            value = _float_or_none(_info_value(record, raw_key))
            if value is not None:
                overlap_fraction = value if overlap_fraction is None else max(overlap_fraction, value)
            continue
        if "numendoverlaps" in normalized_key:
            value = _int_or_none(_info_value(record, raw_key))
            if value is not None:
                end_overlap_count = value if end_overlap_count is None else max(end_overlap_count, value)
            continue

    if overlap_fraction is None and end_overlap_count is None:
        for token in tokens:
            raw_key = normalized_keys.get(token)
            if raw_key is None:
                continue
            value = _info_value(record, raw_key)
            if value not in (None, False, 0, "."):
                legacy_match = True
                break

    return overlap_fraction, end_overlap_count, legacy_match


def _context_passes(
    *,
    svtype: str,
    overlap_fraction: Optional[float],
    end_overlap_count: Optional[int],
    legacy_match: bool,
    context_overlap: float,
) -> bool:
    if svtype in _SPAN_BASED_CONTEXT_SVTYPES:
        return overlap_fraction is not None and overlap_fraction >= context_overlap
    if overlap_fraction is not None and overlap_fraction > 0:
        return True
    if end_overlap_count is not None and end_overlap_count > 0:
        return True
    return legacy_match


def _alt_allele(record: pysam.VariantRecord) -> str:
    alts = tuple(record.alts or ())
    return ",".join(alts) if alts else ""


def _filter_text(record: pysam.VariantRecord) -> str:
    values = filter_values(record)
    if not values:
        return "."
    return ";".join(sorted(values))


def _compute_gt_counts(record: pysam.VariantRecord) -> Dict[str, int]:
    hom_ref = 0
    het = 0
    hom_alt = 0
    no_call = 0
    for sample in record.samples.values():
        gt = sample.get("GT")
        if not gt or any(allele is None for allele in gt):
            no_call += 1
            continue
        alt_allele_count = sum(1 for allele in gt if allele and allele > 0)
        if alt_allele_count == 0:
            hom_ref += 1
        elif alt_allele_count == 1:
            het += 1
        else:
            hom_alt += 1
    return {
        "n_bi_genos": hom_ref + het + hom_alt,
        "n_hom_ref": hom_ref,
        "n_het": het,
        "n_hom_alt": hom_alt,
        "n_no_call": no_call,
    }


def _extract_gq_array(record: pysam.VariantRecord, sample_indices: Optional[np.ndarray]) -> Optional[np.ndarray]:
    if "GQ" not in record.format:
        return None
    sample_values = list(record.samples.values())
    indices: Sequence[int]
    if sample_indices is None:
        indices = range(len(sample_values))
    else:
        indices = [int(index) for index in sample_indices]
    values: List[float] = []
    for index in indices:
        gq = sample_values[index].get("GQ")
        values.append(np.nan if gq in (None, ".") else float(gq))
    return np.asarray(values, dtype=float)


def _extract_concordance_metrics(record: pysam.VariantRecord) -> Optional[Dict[str, float]]:
    metrics: Dict[str, float] = {}
    for field in _CONCORDANCE_FIELDS:
        value = _float_or_none(_info_value(record, field))
        if value is not None:
            metrics[field] = value
    return metrics or None


def _genomic_context(record: pysam.VariantRecord, svtype: str, context_overlap: float = 0.5) -> str:
    normalized_keys = {_normalize_info_key(str(key)): str(key) for key in record.info.keys()}
    candidates: List[Tuple[float, int, int, str]] = []
    for context_index, context in enumerate(GENOMIC_CONTEXTS):
        if context == "none":
            continue
        overlap_fraction, end_overlap_count, legacy_match = _context_metrics(record, normalized_keys, context)
        if not _context_passes(
            svtype=svtype,
            overlap_fraction=overlap_fraction,
            end_overlap_count=end_overlap_count,
            legacy_match=legacy_match,
            context_overlap=context_overlap,
        ):
            continue
        candidates.append(
            (
                overlap_fraction if overlap_fraction is not None else -1.0,
                end_overlap_count if end_overlap_count is not None else -1,
                -context_index,
                context,
            )
        )
    if not candidates:
        return "none"
    candidates.sort(reverse=True)
    return candidates[0][3]


def _iter_records_for_contig(vcf: pysam.VariantFile, contig: str) -> Iterator[pysam.VariantRecord]:
    try:
        yield from vcf.fetch(contig)
    except ValueError:
        for record in vcf:
            if record.contig == contig:
                yield record


def iter_contig(
    vcf_path: Path,
    contig: str,
    extract_gq: bool = False,
    extract_concordance: bool = False,
    sample_indices: Optional[np.ndarray] = None,
    context_overlap: float = 0.5,
) -> Iterator[SiteRecord]:
    """Yield SiteRecord values for one contig."""
    with pysam.VariantFile(str(vcf_path)) as vcf:
        use_precomputed_counts = has_precomputed_counts(vcf)
        for record in _iter_records_for_contig(vcf, contig):
            svtype_raw = str(_info_value(record, "SVTYPE") or "UNKNOWN")
            alt_allele = _alt_allele(record)
            svtype = normalize_svtype(svtype_raw, alt_allele)
            svlen = _int_or_none(_info_value(record, "SVLEN"))
            cn_nonref_freq = _float_or_none(_info_value(record, "CN_NONREF_FREQ"))
            af = cn_nonref_freq if svtype == "CNV" else _float_or_none(_info_value(record, "AF"))
            ac = _int_or_none(_info_value(record, "AC"))
            an = _int_or_none(_info_value(record, "AN"))

            if use_precomputed_counts:
                n_bi_genos = int(_info_value(record, "N_BI_GENOS") or 0)
                n_hom_ref = int(_info_value(record, "N_HOMREF") or 0)
                n_het = int(_info_value(record, "N_HET") or 0)
                n_hom_alt = int(_info_value(record, "N_HOMALT") or 0)
                n_no_call = int(_info_value(record, "NCN") or max(len(record.samples) - n_bi_genos, 0))
            elif svtype == "CNV":
                n_bi_genos = int(_info_value(record, "CN_NUMBER") or len(record.samples))
                n_hom_ref = 0
                n_het = 0
                n_hom_alt = 0
                n_no_call = 0
            else:
                gt_counts = _compute_gt_counts(record)
                n_bi_genos = gt_counts["n_bi_genos"]
                n_hom_ref = gt_counts["n_hom_ref"]
                n_het = gt_counts["n_het"]
                n_hom_alt = gt_counts["n_hom_alt"]
                n_no_call = gt_counts["n_no_call"]

            yield SiteRecord(
                variant_id=_record_identifier(record),
                contig=record.contig,
                pos=int(record.pos),
                end=int(record.stop),
                svtype=svtype,
                svlen=svlen,
                alt_allele=alt_allele,
                filter=_filter_text(record),
                af=af,
                ac=ac,
                an=an,
                n_bi_genos=n_bi_genos,
                n_hom_ref=n_hom_ref,
                n_het=n_het,
                n_hom_alt=n_hom_alt,
                n_no_call=n_no_call,
                cn_nonref_freq=cn_nonref_freq,
                status=str(_info_value(record, "STATUS")) if _info_value(record, "STATUS") not in (None, ".") else None,
                truth_vid=str(_info_value(record, "TRUTH_VID")) if _info_value(record, "TRUTH_VID") not in (None, ".") else None,
                genomic_context=_genomic_context(record, svtype, context_overlap=context_overlap),
                algorithms=_algorithms(record),
                evidence_bucket=_evidence_bucket(record),
                gq_array=_extract_gq_array(record, sample_indices) if extract_gq else None,
                concordance_metrics=_extract_concordance_metrics(record) if extract_concordance else None,
            )
