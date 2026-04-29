#!/usr/bin/env python3
"""Extract columns 1-6 and remap comma-separated sample names in column 6.

Mapping file format:
- Column 1: replacement value
- Column 2: sample name key to match against column 6 entries
"""

import argparse
import gzip
from typing import Dict, TextIO


def open_text(path: str, mode: str) -> TextIO:
    if path.endswith(".gz"):
        return gzip.open(path, mode)  # type: ignore[return-value]
    return open(path, mode)


def parse_first_six(line: str):
    stripped = line.rstrip("\n")
    if "\t" in stripped:
        parts = stripped.split("\t")
    else:
        parts = stripped.split()
    if len(parts) < 6:
        return None
    return parts[:6]


def load_mapping(mapping_path: str) -> Dict[str, str]:
    mapping: Dict[str, str] = {}
    with open_text(mapping_path, "rt") as fin:
        for raw in fin:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue

            if "\t" in line:
                cols = line.split("\t")
            else:
                cols = line.split()

            if len(cols) < 2:
                continue

            replacement = cols[0].strip()
            sample_name = cols[1].strip()
            if sample_name:
                mapping[sample_name] = replacement

    return mapping


def remap_sample_list(sample_list: str, mapping: Dict[str, str]) -> str:
    names = [s.strip() for s in sample_list.split(",")]
    remapped = [mapping.get(name, name) for name in names if name]
    return ",".join(remapped)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Extract columns 1-6 from input and remap comma-separated sample names "
            "in column 6 using mapping file col2->col1."
        )
    )
    parser.add_argument("--input", required=True, help="Input file (e.g. test.variants.chr8_15_x)")
    parser.add_argument(
        "--mapping",
        required=True,
        help="Mapping file (e.g. gnomAD_SV_v3.final_retained.releasable.sample_batch)",
    )
    parser.add_argument("--output", required=True, help="Output file path")
    args = parser.parse_args()

    mapping = load_mapping(args.mapping)

    processed = 0
    skipped = 0

    with open_text(args.input, "rt") as fin, open_text(args.output, "wt") as fout:
        for raw in fin:
            if not raw.strip():
                continue

            first6 = parse_first_six(raw)
            if first6 is None:
                skipped += 1
                continue

            first6[5] = remap_sample_list(first6[5], mapping)
            fout.write("\t".join(first6) + "\n")
            processed += 1

    print(
        {
            "processed_rows": processed,
            "skipped_rows_lt_6_cols": skipped,
            "mapped_keys": len(mapping),
            "output": args.output,
        }
    )


if __name__ == "__main__":
    main()
