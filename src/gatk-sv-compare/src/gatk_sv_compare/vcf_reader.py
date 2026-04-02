"""Streaming VCF extraction utilities for Phase 3 aggregation."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Sequence

import numpy as np
import pysam

from .dimensions import GENOMIC_CONTEXTS, normalize_svtype
from .vcf_format import filter_values, has_precomputed_counts, safe_info_get

_GENOMIC_CONTEXT_ALIASES = {
    "segdup": ("segdup", "SEG_DUP", "SEGDUP", "segmental_duplication", "segmentalduplication"),
    "simple_repeat": ("simple_repeat", "SIMPLE_REPEAT", "simplerepeat", "tandem_repeat", "trf"),
    "repeatmasker": ("repeatmasker", "REPEATMASKER", "repeat_masker", "rmsk"),
}
_GENOMIC_CONTEXT_SUBSTRINGS = {
    "segdup": ("segdup", "segdup", "segmentaldup", "seg_dup"),
    "simple_repeat": ("simplerepeat", "simple_repeat", "tandemrepeat", "trf"),
    "repeatmasker": ("repeatmasker", "repeat_masker", "rmsk"),
}

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


def _int_or_none(value: object) -> Optional[int]:
    if value in (None, "."):
        return None
    return int(value)


def _float_or_none(value: object) -> Optional[float]:
    if value in (None, "."):
        return None
    return float(value)


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


def _genomic_context(record: pysam.VariantRecord) -> str:
    available_keys = {str(key).lower().replace("_", "").replace("-", ""): str(key) for key in record.info.keys()}
    for context in GENOMIC_CONTEXTS:
        if context == "none":
            continue
        for alias in _GENOMIC_CONTEXT_ALIASES.get(context, (context,)):
            key = available_keys.get(alias.lower().replace("_", "").replace("-", ""), alias)
            value = _info_value(record, key)
            if value in (None, False, 0, "."):
                continue
            return context
        for normalized_key, raw_key in available_keys.items():
            if any(fragment in normalized_key for fragment in _GENOMIC_CONTEXT_SUBSTRINGS.get(context, ())):
                value = _info_value(record, raw_key)
                if value not in (None, False, 0, "."):
                    return context
    return "none"


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
                genomic_context=_genomic_context(record),
                gq_array=_extract_gq_array(record, sample_indices) if extract_gq else None,
                concordance_metrics=_extract_concordance_metrics(record) if extract_concordance else None,
            )
