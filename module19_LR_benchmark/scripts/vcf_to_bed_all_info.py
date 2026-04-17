#!/usr/bin/env python3

"""Convert VCF/VCF.GZ to BED with all INFO fields as columns.

Output columns start with:
    chrom, start, end, id, SVLEN, SVTYPE, AC, AN, SVANN

Then all remaining INFO IDs from the VCF header are appended, preserving
header order. REF/ALT sequences are intentionally not included.

Notes:
- BED start is 0-based.
- BED end is half-open; uses INFO/END when present, otherwise POS + len(REF) - 1.
"""

import argparse
import gzip
import sys
from typing import Dict, Iterable, List, TextIO


PRIORITY_INFO_KEYS = ["SVLEN", "SVTYPE", "AC", "AN", "SVANN"]


def open_text(path: str, mode: str) -> TextIO:
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def parse_info_field(info_str: str) -> Dict[str, str]:
    info: Dict[str, str] = {}
    if info_str == "." or not info_str:
        return info

    for item in info_str.split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            info[key] = value
        else:
            # INFO flags are represented as present/absent in VCF.
            info[item] = "1"
    return info


def _first_value(raw: str) -> str:
    if "," in raw:
        return raw.split(",", 1)[0]
    return raw


def _abs_int_str(raw: str) -> str:
    first = _first_value(raw).strip()
    if first == "":
        return "."
    try:
        return str(abs(int(float(first))))
    except ValueError:
        return "."


def get_priority_value(info_map: Dict[str, str], key: str, missing_value: str) -> str:
    if key == "SVLEN":
        if "SVLEN" in info_map:
            v = _abs_int_str(info_map["SVLEN"])
            return v if v != "." else missing_value
        if "allele_length" in info_map:
            v = _abs_int_str(info_map["allele_length"])
            return v if v != "." else missing_value
        return missing_value

    if key == "SVTYPE":
        if "SVTYPE" in info_map and info_map["SVTYPE"]:
            return _first_value(info_map["SVTYPE"])
        if "allele_type" in info_map and info_map["allele_type"]:
            return _first_value(info_map["allele_type"]).upper()
        return missing_value

    if key in info_map and info_map[key] != "":
        return _first_value(info_map[key])

    return missing_value


def extract_info_ids_from_header(vcf_path: str) -> List[str]:
    info_ids: List[str] = []
    seen = set()

    with open_text(vcf_path, "rt") as fh:
        for line in fh:
            if line.startswith("##INFO=<ID="):
                payload = line.strip()[len("##INFO=<") :]
                payload = payload[:-1] if payload.endswith(">") else payload
                for entry in payload.split(","):
                    if entry.startswith("ID="):
                        info_id = entry[3:]
                        if info_id and info_id not in seen:
                            seen.add(info_id)
                            info_ids.append(info_id)
                        break
            elif line.startswith("#CHROM"):
                break

    return info_ids


def to_bed_end(pos_1based: int, ref: str, info_map: Dict[str, str]) -> int:
    end_raw = info_map.get("END")
    if end_raw is not None:
        try:
            end_val = int(end_raw)
            if end_val >= pos_1based:
                return end_val
        except ValueError:
            pass

    ref_len = max(1, len(ref))
    # For [POS, POS+ref_len-1] in 1-based closed coordinates,
    # BED half-open end is POS+ref_len-1.
    return pos_1based + ref_len - 1


def convert(vcf_path: str, out_bed: str, missing_value: str) -> None:
    info_ids = extract_info_ids_from_header(vcf_path)
    remaining_info_ids = [k for k in info_ids if k not in PRIORITY_INFO_KEYS]

    with open_text(vcf_path, "rt") as vcf, open_text(out_bed, "wt") as out:
        header = ["chrom", "start", "end", "id"] + PRIORITY_INFO_KEYS + remaining_info_ids
        out.write("\t".join(header) + "\n")

        for line in vcf:
            if not line or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue

            chrom = fields[0]
            pos_1based = int(fields[1])
            var_id = fields[2]
            ref = fields[3]
            info_map = parse_info_field(fields[7])

            start_0based = pos_1based - 1
            end_1based = to_bed_end(pos_1based, ref, info_map)

            row = [chrom, str(start_0based), str(end_1based), var_id]
            row.extend(get_priority_value(info_map, k, missing_value) for k in PRIORITY_INFO_KEYS)
            row.extend(info_map.get(k, missing_value) for k in remaining_info_ids)
            out.write("\t".join(row) + "\n")


def parse_args(argv: Iterable[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert VCF to BED and include all INFO fields as columns."
    )
    parser.add_argument("--input", required=True, help="Input VCF or VCF.GZ")
    parser.add_argument("--output", required=True, help="Output BED (or BED.GZ)")
    parser.add_argument(
        "--missing",
        default=".",
        help="Placeholder for missing INFO values (default: .)",
    )
    return parser.parse_args(list(argv))


def main(argv: Iterable[str]) -> int:
    args = parse_args(argv)
    convert(args.input, args.output, args.missing)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
