#!/usr/bin/env python3
"""Report BED variant pairs with matching SVTYPE, reciprocal interval overlap, and sample overlap.

Expected input columns (at minimum):
1. chrom
2. start
3. end
4. name
5. svtype
6. samples (comma-separated)
"""

import argparse
import gzip
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Set, TextIO, Tuple


@dataclass
class Variant:
    chrom: str
    start: int
    end: int
    name: str
    svtype: str
    samples: Set[str]


def open_text(path: str, mode: str) -> TextIO:
    if path.endswith(".gz"):
        return gzip.open(path, mode)  # type: ignore[return-value]
    return open(path, mode)


def parse_line(line: str) -> Optional[Variant]:
    s = line.strip()
    if not s or s.startswith("#"):
        return None

    parts = s.split("\t") if "\t" in s else s.split()
    if len(parts) < 6:
        return None

    chrom = parts[0]
    try:
        start = int(parts[1])
        end = int(parts[2])
    except ValueError:
        return None

    if end < start:
        start, end = end, start

    name = parts[3]
    svtype = parts[4]
    sample_field = parts[5]
    samples = {x.strip() for x in sample_field.split(",") if x.strip() and x.strip() != "."}

    return Variant(chrom=chrom, start=start, end=end, name=name, svtype=svtype, samples=samples)


def interval_length(start: int, end: int) -> int:
    return max(0, end - start)


def reciprocal_overlap(v1: Variant, v2: Variant) -> Tuple[float, float, int]:
    overlap = max(0, min(v1.end, v2.end) - max(v1.start, v2.start))
    len1 = interval_length(v1.start, v1.end)
    len2 = interval_length(v2.start, v2.end)

    ro1 = overlap / len1 if len1 > 0 else 0.0
    ro2 = overlap / len2 if len2 > 0 else 0.0
    return ro1, ro2, overlap


def reciprocal_sample_overlap(v1: Variant, v2: Variant) -> Tuple[float, float, int]:
    if not v1.samples or not v2.samples:
        return 0.0, 0.0, 0

    inter = v1.samples & v2.samples
    count = len(inter)
    s1 = count / len(v1.samples) if v1.samples else 0.0
    s2 = count / len(v2.samples) if v2.samples else 0.0
    return s1, s2, count


def load_variants(path: str) -> List[Variant]:
    variants: List[Variant] = []
    with open_text(path, "rt") as fin:
        for raw in fin:
            v = parse_line(raw)
            if v is not None:
                variants.append(v)
    return variants


def group_by_chrom_svtype(variants: Sequence[Variant]) -> Dict[Tuple[str, str], List[Variant]]:
    grouped: Dict[Tuple[str, str], List[Variant]] = {}
    for v in variants:
        key = (v.chrom, v.svtype)
        grouped.setdefault(key, []).append(v)

    for key in grouped:
        grouped[key].sort(key=lambda x: (x.start, x.end, x.name))

    return grouped


def find_pairs(
    grouped: Dict[Tuple[str, str], List[Variant]],
    min_genomic_recip: float,
    min_sample_recip: float,
) -> Iterable[str]:
    header = [
        "chrom",
        "svtype",
        "name1",
        "start1",
        "end1",
        "name2",
        "start2",
        "end2",
        "genomic_overlap_bp",
        "genomic_recip_1",
        "genomic_recip_2",
        "shared_samples_count",
        "samples_count_1",
        "samples_count_2",
        "sample_recip_1",
        "sample_recip_2",
    ]
    yield "\t".join(header)

    for (chrom, svtype), items in grouped.items():
        active: List[Variant] = []

        for v in items:
            active = [a for a in active if a.end > v.start]

            for a in active:
                gro1, gro2, gbp = reciprocal_overlap(a, v)
                if gro1 < min_genomic_recip or gro2 < min_genomic_recip:
                    continue

                sro1, sro2, shared = reciprocal_sample_overlap(a, v)
                if sro1 < min_sample_recip or sro2 < min_sample_recip:
                    continue

                row = [
                    chrom,
                    svtype,
                    a.name,
                    str(a.start),
                    str(a.end),
                    v.name,
                    str(v.start),
                    str(v.end),
                    str(gbp),
                    f"{gro1:.4f}",
                    f"{gro2:.4f}",
                    str(shared),
                    str(len(a.samples)),
                    str(len(v.samples)),
                    f"{sro1:.4f}",
                    f"{sro2:.4f}",
                ]
                yield "\t".join(row)

            active.append(v)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Report variant pairs in a BED-like file that share SVTYPE, have at least "
            "a minimum reciprocal genomic overlap, and a minimum reciprocal sample overlap."
        )
    )
    parser.add_argument("--input", required=True, help="Input BED-like file path")
    parser.add_argument("--output", required=True, help="Output TSV path")
    parser.add_argument(
        "--min-genomic-recip",
        type=float,
        default=0.5,
        help="Minimum reciprocal genomic overlap for both variants (default: 0.5)",
    )
    parser.add_argument(
        "--min-sample-recip",
        type=float,
        default=0.5,
        help="Minimum reciprocal sample overlap for both variants (default: 0.5)",
    )
    args = parser.parse_args()

    variants = load_variants(args.input)
    grouped = group_by_chrom_svtype(variants)

    n_written = 0
    with open_text(args.output, "wt") as fout:
        for line in find_pairs(grouped, args.min_genomic_recip, args.min_sample_recip):
            fout.write(line + "\n")
            n_written += 1

    print(
        {
            "input_variants": len(variants),
            "groups": len(grouped),
            "output": args.output,
            "output_lines_including_header": n_written,
            "pair_count": max(0, n_written - 1),
        }
    )


if __name__ == "__main__":
    main()
