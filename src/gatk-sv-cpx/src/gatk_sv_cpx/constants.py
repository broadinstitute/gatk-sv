"""
Shared constants, CPX type definitions, and INFO-field keys used throughout
the ``gatk-sv-cpx`` package.
"""

from enum import Enum
from typing import Dict, FrozenSet

# ---------------------------------------------------------------------------
# CPX type taxonomy
# ---------------------------------------------------------------------------

class CpxType(str, Enum):
    """Complex SV subtypes recognised by the GATK-SV pipeline."""

    # Inversion-flanked
    delINV = "delINV"
    INVdel = "INVdel"
    dupINV = "dupINV"
    INVdup = "INVdup"
    delINVdel = "delINVdel"
    delINVdup = "delINVdup"
    dupINVdel = "dupINVdel"
    dupINVdup = "dupINVdup"

    # Dispersed duplication
    dDUP = "dDUP"
    dDUP_iDEL = "dDUP_iDEL"

    # Insertion with internal deletion
    INS_iDEL = "INS_iDEL"

    # Palindromic duplications
    piDUP_RF = "piDUP_RF"
    piDUP_FR = "piDUP_FR"

    # Translocations
    CTX_PQ = "CTX_PQ/QP"
    CTX_PP = "CTX_PP/QQ"

    # Simple resolved inversion
    INV = "INV"


# CPX types that contain DEL/DUP CNV sub-intervals
CPX_TYPES_WITH_CNV: FrozenSet[str] = frozenset({
    "delINV", "INVdel", "dupINV", "INVdup",
    "delINVdel", "delINVdup", "dupINVdel", "dupINVdup",
    "dDUP", "dDUP_iDEL", "INS_iDEL",
    "piDUP_RF", "piDUP_FR",
})

# CPX types for which SOURCE field is relevant
CPX_TYPES_WITH_SOURCE: FrozenSet[str] = frozenset({
    "dDUP", "dDUP_iDEL", "INS_iDEL",
})

# Insertion-type CPX that get special handling during refinement
INSERTION_TYPE_CPX: FrozenSet[str] = frozenset({
    "dDUP", "dDUP_iDEL", "INS_iDEL",
})

# ---------------------------------------------------------------------------
# VCF INFO field keys
# ---------------------------------------------------------------------------

SVTYPE = "SVTYPE"
CPX_TYPE_KEY = "CPX_TYPE"
CPX_INTERVALS = "CPX_INTERVALS"
SOURCE = "SOURCE"
MEMBERS = "MEMBERS"
EVIDENCE = "EVIDENCE"
ALGORITHMS = "ALGORITHMS"
STRANDS = "STRANDS"
SVLEN = "SVLEN"
CHR2 = "CHR2"
END2 = "END2"
VARGQ = "varGQ"

# Filter tags
UNRESOLVED = "UNRESOLVED"
UNRESOLVED_TYPE = "UNRESOLVED_TYPE"

# ---------------------------------------------------------------------------
# PE support classification levels
# ---------------------------------------------------------------------------

class PeSupport(str, Enum):
    """Per-sample PE support categories."""
    NO_PE = "no_PE"
    PARTIAL_PE = "partial_PE"
    LOW_PE = "low_PE"
    HIGH_PE = "high_PE"


class DepthSupport(str, Enum):
    """Per-sample depth support categories."""
    DEPTH = "depth"
    LACK_DEPTH = "lack_depth"
    NA = "NA"


# ---------------------------------------------------------------------------
# Header lines to inject
# ---------------------------------------------------------------------------

VCF_HEADER_LINES: Dict[str, str] = {
    "UNRESOLVED_FILTER":
        '##FILTER=<ID=UNRESOLVED,Description="Variant is unresolved.">',
    "CPX_TYPE":
        '##INFO=<ID=CPX_TYPE,Number=1,Type=String,Description="Class of complex variant.">',
    "CPX_INTERVALS":
        '##INFO=<ID=CPX_INTERVALS,Number=.,Type=String,Description="Intervals constituting complex variant.">',
    "SOURCE":
        '##INFO=<ID=SOURCE,Number=1,Type=String,Description="Source interval for dispersed dup / insertion.">',
    "UNRESOLVED_FLAG":
        '##INFO=<ID=UNRESOLVED,Number=0,Type=Flag,Description="Variant is unresolved.">',
    "UNRESOLVED_TYPE_INFO":
        '##INFO=<ID=UNRESOLVED_TYPE,Number=1,Type=String,Description="Class of unresolved variant.">',
    "END2":
        '##INFO=<ID=END2,Number=1,Type=Integer,Description="Position of breakpoint on CHR2.">',
    "CHR2_INFO":
        '##INFO=<ID=CHR2,Number=1,Type=String,Description="Chromosome for END2 coordinate.">',
    "BPID":
        '##INFO=<ID=BPID,Number=.,Type=String,'
        'Description="ID of retained variant from breakpoint overlap filtering">',
}

# ---------------------------------------------------------------------------
# Default thresholds
# ---------------------------------------------------------------------------

DEFAULT_MIN_PE_CPX = 3
DEFAULT_MIN_PE_CTX = 3
DEFAULT_MIN_DEPTH_SIZE = 5000
DEFAULT_DEPTH_FRACTION_THRESHOLD = 0.5
DEFAULT_MIN_SIZE = 1000
DEFAULT_MIN_DIFF = 0.4
DEFAULT_MIN_SIZE_IDEL = 150
DEFAULT_MIN_SVLEN = 50
DEFAULT_MIN_RECIPROCAL_OVERLAP = 0.5
