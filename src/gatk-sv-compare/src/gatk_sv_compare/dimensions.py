"""Shared bucketing logic and constants for gatk-sv-compare."""

from __future__ import annotations

from dataclasses import dataclass
from typing import FrozenSet, Iterable, Mapping, Optional, Set, Union

SVTYPE_ORDER = ["DEL", "DUP", "CNV", "INS", "INS:MEI", "INV", "CPX", "CTX", "BND"]
SVTYPES_CORE = ["DEL", "DUP", "INS", "INV", "BND", "CPX", "CTX", "CNV"]
SVTYPE_INS_MEI = "INS:MEI"
MEI_ALT_PATTERNS = frozenset({"<INS:ME:", "<INS:ME>"})
INS_NON_MEI_ALTS = frozenset({"<INS>", "<INS:UNK>"})

AF_BUCKETS = [
    ("AC=1", None),
    ("AC>1,AF<1%", (0.0, 0.01)),
    ("1-10%", (0.01, 0.10)),
    ("10-50%", (0.10, 0.50)),
    (">50%", (0.50, 1.00)),
]

SIZE_BUCKETS = [
    ("<100bp", (0, 100)),
    ("100-500bp", (100, 500)),
    ("500bp-2.5kb", (500, 2500)),
    ("2.5-10kb", (2500, 10000)),
    ("10-50kb", (10000, 50000)),
    (">50kb", (50000, float("inf"))),
]

GENOMIC_CONTEXTS = ["simple_repeat", "segdup", "repeatmasker", "none"]
_AF_BUCKET_ORDER = [label for label, _ in AF_BUCKETS] + ["unknown"]
_SIZE_BUCKET_ORDER = [label for label, _ in SIZE_BUCKETS] + ["unknown", "N/A"]

FILTER_PASS = "PASS"
FILTER_MULTIALLELIC = "MULTIALLELIC"

STATUS_MATCHED = "MATCHED"
STATUS_UNMATCHED = "UNMATCHED"


@dataclass(frozen=True)
class VariantCategory:
    """Categorization key used by downstream analysis modules."""

    svtype: str
    size_bucket: str
    af_bucket: str
    genomic_context: str


def is_mei_alt(alt_allele: Optional[str]) -> bool:
    """Return True when the ALT allele corresponds to a mobile element insertion."""
    if not alt_allele:
        return False
    return any(pattern in alt_allele for pattern in MEI_ALT_PATTERNS)


def normalize_svtype(svtype: Optional[str], alt_allele: Optional[str] = None) -> str:
    """Collapse raw ALT/SVTYPE combinations into a stable analysis label."""
    if svtype == "INS" and is_mei_alt(alt_allele):
        return SVTYPE_INS_MEI
    if svtype:
        return svtype
    return "UNKNOWN"


def is_filtered_pass(filter_values: Union[Set[str], FrozenSet[str]]) -> bool:
    """Return True when a record belongs to the filtered-pass analysis view."""
    return FILTER_PASS in filter_values or FILTER_MULTIALLELIC in filter_values


def bucket_size(svtype: str, svlen: Optional[int]) -> str:
    """Assign an SV to a size bucket."""
    if svtype in {"BND", "CTX"}:
        return "N/A"
    if svlen is None or svlen <= 0:
        return "unknown"
    for label, (lower, upper) in SIZE_BUCKETS:
        if lower < svlen <= upper:
            return label
    return "unknown"


def bucket_af(ac: Optional[int], af: Optional[float]) -> str:
    """Assign an SV to an AF bucket."""
    if ac == 1:
        return "AC=1"
    if af is None:
        return "unknown"
    for label, bounds in AF_BUCKETS:
        if bounds is None:
            continue
        lower, upper = bounds
        if lower < af <= upper:
            return label
    return "unknown"


def normalize_context(genomic_context: Optional[str]) -> str:
    """Normalize genomic context values to the supported set."""
    if genomic_context in GENOMIC_CONTEXTS:
        return genomic_context
    return "none"


def svtype_sort_key(value: object) -> tuple[int, str]:
    normalized = str(value)
    try:
        return (SVTYPE_ORDER.index(normalized), normalized)
    except ValueError:
        return (len(SVTYPE_ORDER), normalized)


def size_bucket_sort_key(value: object) -> tuple[int, str]:
    normalized = str(value)
    try:
        return (_SIZE_BUCKET_ORDER.index(normalized), normalized)
    except ValueError:
        return (len(_SIZE_BUCKET_ORDER), normalized)


def af_bucket_sort_key(value: object) -> tuple[int, str]:
    normalized = str(value)
    try:
        return (_AF_BUCKET_ORDER.index(normalized), normalized)
    except ValueError:
        return (len(_AF_BUCKET_ORDER), normalized)


def genomic_context_sort_key(value: object) -> tuple[int, str]:
    normalized = str(value)
    try:
        return (GENOMIC_CONTEXTS.index(normalized), normalized)
    except ValueError:
        return (len(GENOMIC_CONTEXTS), normalized)


def ordered_svtypes(values: Iterable[object]) -> list[str]:
    return [str(value) for value in sorted({str(value) for value in values if value is not None}, key=svtype_sort_key)]


def ordered_size_buckets(values: Iterable[object]) -> list[str]:
    return [str(value) for value in sorted({str(value) for value in values if value is not None}, key=size_bucket_sort_key)]


def ordered_af_buckets(values: Iterable[object]) -> list[str]:
    return [str(value) for value in sorted({str(value) for value in values if value is not None}, key=af_bucket_sort_key)]


def ordered_contexts(values: Iterable[object]) -> list[str]:
    return [str(value) for value in sorted({str(value) for value in values if value is not None}, key=genomic_context_sort_key)]


def categorize_variant(record_info: Mapping[str, object]) -> VariantCategory:
    """Return the analysis category for an extracted record."""
    raw_svtype = str(record_info.get("svtype") or "UNKNOWN")
    alt_allele = record_info.get("alt_allele")
    svtype = normalize_svtype(raw_svtype, str(alt_allele) if alt_allele is not None else None)
    is_cnv = svtype == "CNV"
    af_value = record_info.get("cn_nonref_freq") if is_cnv else record_info.get("af")
    af = float(af_value) if af_value not in (None, ".") else None
    ac_value = record_info.get("ac")
    ac = int(ac_value) if ac_value not in (None, ".") else None
    svlen_value = record_info.get("svlen")
    svlen = int(svlen_value) if svlen_value not in (None, ".") else None
    context = normalize_context(str(record_info.get("genomic_context") or "none"))
    return VariantCategory(
        svtype=svtype,
        size_bucket=bucket_size(svtype, svlen),
        af_bucket=bucket_af(ac=ac, af=af),
        genomic_context=context,
    )
