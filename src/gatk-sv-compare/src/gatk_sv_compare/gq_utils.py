"""Helpers for detecting and normalizing genotype-quality scales."""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path

import numpy as np
import pysam

_DEFAULT_GQ_SCAN_RECORD_LIMIT = 10_000


@lru_cache(maxsize=None)
def detect_gq_scale_factor(vcf_path: Path, scan_record_limit: int = _DEFAULT_GQ_SCAN_RECORD_LIMIT) -> float:
    """Return the divisor needed to normalize GQ values to the 0-99 range.

    Supported inputs are standard 0-99 GQs and expanded 0-999 GQs. Detection
    scans up to ``scan_record_limit`` records and switches to a 10x divisor as
    soon as any finite GQ exceeds 99.
    """
    with pysam.VariantFile(str(vcf_path)) as vcf:
        for record_index, record in enumerate(vcf, start=1):
            if "GQ" not in record.format:
                continue
            for sample in record.samples.values():
                gq = sample.get("GQ")
                if gq in (None, "."):
                    continue
                gq_value = float(gq)
                if np.isfinite(gq_value) and gq_value > 99.0:
                    return 10.0
            if record_index >= scan_record_limit:
                break
    return 1.0


def normalize_gq_value(gq: object, scale_factor: float) -> float:
    """Normalize one GQ value to the canonical 0-99 scale."""
    return float(gq) / scale_factor
