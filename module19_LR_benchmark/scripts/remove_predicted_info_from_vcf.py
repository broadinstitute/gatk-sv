#!/usr/bin/env python3
"""Remove INFO fields starting with PREDICTED_ from a VCF.

This script updates both:
1) Header INFO declarations (##INFO=<ID=PREDICTED_...>)
2) Per-record INFO entries where key starts with PREDICTED_

Supports plain text and .gz input/output based on file extension.
"""

import argparse
import gzip
import re
from typing import List, TextIO


def open_text(path: str, mode: str) -> TextIO:
    if path.endswith(".gz"):
        return gzip.open(path, mode)  # type: ignore[return-value]
    return open(path, mode)


REMOVE_PREFIXES = ("PREDICTED_",)
REMOVE_SUBSTRINGS = ("_AC", "_AF", "_AN", "_N_", "_FREQ_")


def should_remove_info_key(key: str) -> bool:
    if any(key.startswith(p) for p in REMOVE_PREFIXES):
        return True
    if any(sub in key for sub in REMOVE_SUBSTRINGS):
        return True
    return False


def is_predicted_info_header(line: str) -> bool:
    if not line.startswith("##INFO=<ID="):
        return False
    m = re.match(r"##INFO=<ID=([^,>]+)", line)
    if not m:
        return False
    return should_remove_info_key(m.group(1))


def filter_info_field(info: str) -> str:
    if info == "." or info == "":
        return "."

    kept: List[str] = []
    for entry in info.split(";"):
        if not entry:
            continue
        key = entry.split("=", 1)[0]
        if should_remove_info_key(key):
            continue
        kept.append(entry)

    return ";".join(kept) if kept else "."


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Remove INFO fields with names starting with PREDICTED_ from a VCF"
    )
    parser.add_argument("--input", required=True, help="Input VCF or VCF.GZ")
    parser.add_argument("--output", required=True, help="Output VCF or VCF.GZ")
    args = parser.parse_args()

    total_records = 0
    updated_records = 0
    removed_header_info = 0

    with open_text(args.input, "rt") as fin, open_text(args.output, "wt") as fout:
        for line in fin:
            if line.startswith("##"):
                if is_predicted_info_header(line):
                    removed_header_info += 1
                    continue
                fout.write(line)
                continue

            if line.startswith("#"):
                fout.write(line)
                continue

            total_records += 1
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                fout.write(line)
                continue

            old_info = fields[7]
            new_info = filter_info_field(old_info)
            if new_info != old_info:
                updated_records += 1
            fields[7] = new_info

            fout.write("\t".join(fields) + "\n")

    print(
        {
            "input": args.input,
            "output": args.output,
            "removed_header_info_definitions": removed_header_info,
            "total_records": total_records,
            "updated_records": updated_records,
        }
    )


if __name__ == "__main__":
    main()
