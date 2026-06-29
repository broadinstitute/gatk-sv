#!/usr/bin/env python3
"""Compare two VCF/BCF files with relaxed formatting rules.

Allowed differences:
- Header lines may differ.
- INFO field order may differ.
- FORMAT field order may differ.

Disallowed differences:
- Record order/content differences.
- INFO key/value differences.
- Genotype/value differences per sample.

For each mismatch, print the variant ID. For genotype mismatches, also print
how many samples mismatch and one example sample ID.
"""

from __future__ import annotations

import argparse
from typing import Any, Iterable, Tuple

import pysam

CANONICAL_FMT_ORDER = (
    "GT",
    "EV",
    "GQ",
    "PE_GQ",
    "PE_GT",
    "RD_CN",
    "RD_GQ",
    "SR_GQ",
    "SR_GT",
    "CN",
    "CNQ",
    "ECN",
    "DP",
    "AD",
    "PL",
    "FT",
)
FMT_RANK = {key: index for index, key in enumerate(CANONICAL_FMT_ORDER)}
CANONICAL_FMT_ORDER_BYTES = tuple(key.encode("ascii") for key in CANONICAL_FMT_ORDER)
FMT_RANK_BYTES = {key: index for index, key in enumerate(CANONICAL_FMT_ORDER_BYTES)}


def open_variant_lines(path: str):
    if path.endswith(".gz") or path.endswith(".bgz"):
        return pysam.BGZFile(path, "r")
    return open(path, "rb")


def next_record_line(handle) -> bytes | None:
    for line in handle:
        if line[:1] != b"#":
            return line.rstrip(b"\n")
    return None


def read_sample_names(path: str) -> Tuple[str, ...] | None:
    with open_variant_lines(path) as handle:
        for raw_line in handle:
            if raw_line.startswith(b"#CHROM"):
                header_fields = raw_line.rstrip(b"\n").split(b"\t")
                return tuple(field.decode("utf-8") for field in header_fields[9:])
    return None


def ordered_format_keys_bytes(keys: Iterable[bytes]) -> Tuple[bytes, ...]:
    return tuple(sorted(keys, key=lambda key: (FMT_RANK_BYTES.get(key, len(FMT_RANK_BYTES)), key)))


def normalize_info_bytes(info: bytes) -> Tuple[bytes, ...]:
    if info in (b"", b"."):
        return ()
    return tuple(sorted(info.split(b";")))


def normalize_filter_bytes(value: bytes) -> Tuple[bytes, ...]:
    if value in (b"", b".", b"PASS"):
        return ()
    return tuple(sorted(value.split(b";")))


def bytes_to_text(value: bytes) -> str:
    return value.decode("utf-8", errors="replace")


def describe_site_differences(left_fields: list[bytes], right_fields: list[bytes]) -> list[str]:
    diffs: list[str] = []

    field_names = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL"]
    for index, name in enumerate(field_names):
        left_value = left_fields[index] if index < len(left_fields) else b""
        right_value = right_fields[index] if index < len(right_fields) else b""
        if left_value != right_value:
            diffs.append(
                f"{name}(left={bytes_to_text(left_value)},right={bytes_to_text(right_value)})"
            )

    left_filter = left_fields[6] if len(left_fields) > 6 else b""
    right_filter = right_fields[6] if len(right_fields) > 6 else b""
    left_filter_norm = normalize_filter_bytes(left_filter)
    right_filter_norm = normalize_filter_bytes(right_filter)
    if left_filter_norm != right_filter_norm:
        diffs.append(
            "FILTER("
            f"left={bytes_to_text(left_filter)},"
            f"right={bytes_to_text(right_filter)},"
            f"left_norm={','.join(bytes_to_text(x) for x in left_filter_norm)},"
            f"right_norm={','.join(bytes_to_text(x) for x in right_filter_norm)}"
            ")"
        )

    left_info = left_fields[7] if len(left_fields) > 7 else b""
    right_info = right_fields[7] if len(right_fields) > 7 else b""
    left_info_norm = normalize_info_bytes(left_info)
    right_info_norm = normalize_info_bytes(right_info)
    if left_info_norm != right_info_norm:
        diffs.append(
            "INFO("
            f"left={bytes_to_text(left_info)},"
            f"right={bytes_to_text(right_info)},"
            f"left_norm={';'.join(bytes_to_text(x) for x in left_info_norm)},"
            f"right_norm={';'.join(bytes_to_text(x) for x in right_info_norm)}"
            ")"
        )

    return diffs


def sample_value_map(fmt_keys: Tuple[bytes, ...], sample_value: bytes) -> dict[bytes, bytes]:
    values = sample_value.split(b":") if sample_value else []
    return {
        key: values[index] if index < len(values) and values[index] else b"."
        for index, key in enumerate(fmt_keys)
    }


def normalized_sample_bytes(fmt_keys: Tuple[bytes, ...], sample_value: bytes) -> bytes:
    value_map = sample_value_map(fmt_keys, sample_value)
    ordered_keys = ordered_format_keys_bytes(value_map)
    return b":".join(value_map[key] for key in ordered_keys)


def normalized_format_signature(fmt_keys: Tuple[bytes, ...]) -> Tuple[bytes, ...]:
    return ordered_format_keys_bytes(fmt_keys)


def ordered_format_keys(keys: Iterable[str]) -> Tuple[str, ...]:
    return tuple(sorted(keys, key=lambda key: (FMT_RANK.get(key, len(FMT_RANK)), key)))


def normalize_atom(value: Any) -> str:
    if value is None:
        return "."
    if isinstance(value, bytes):
        return value.decode("utf-8")
    if isinstance(value, float):
        if value != value:
            return "."
        return format(value, ".12g")
    return str(value)


def normalize_info_value(value: Any) -> Tuple[str, ...] | str:
    if isinstance(value, (list, tuple)):
        if not value:
            return "."
        normalized = tuple(normalize_atom(item) for item in value)
        if all(item == "." for item in normalized):
            return "."
        return normalized
    return normalize_atom(value)


def normalize_format_value(key: str, value: Any) -> str:
    if value is None:
        return "."

    if key == "GT":
        if isinstance(value, (list, tuple)):
            if not value:
                return "."
            alleles = tuple("." if item is None else str(item) for item in value)
            if all(item == "." for item in alleles):
                return "."
            return "/".join(alleles)
        return normalize_atom(value)

    if isinstance(value, (list, tuple)):
        if not value:
            return "."
        normalized = tuple(normalize_atom(item) for item in value)
        if all(item == "." for item in normalized):
            return "."
        return ",".join(normalized)

    return normalize_atom(value)


def variant_label(record: pysam.libcbcf.VariantRecord | None) -> str:
    if record is None:
        return "<none>"
    if record.id and record.id != ".":
        return record.id
    alts = ",".join(record.alts or ())
    return f"{record.chrom}:{record.pos}:{record.ref}:{alts}"


def site_equal(left: pysam.libcbcf.VariantRecord, right: pysam.libcbcf.VariantRecord) -> bool:
    if left.chrom != right.chrom or left.pos != right.pos:
        return False
    if (left.id or ".") != (right.id or "."):
        return False
    if left.ref != right.ref:
        return False
    if (left.alts or ()) != (right.alts or ()):
        return False
    if left.qual != right.qual:
        return False
    if set(left.filter.keys()) != set(right.filter.keys()):
        return False

    if len(left.info) != len(right.info):
        return False

    for key, left_value in left.info.items():
        if key not in right.info:
            return False
        right_value = right.info[key]
        if left_value == right_value:
            continue
        if normalize_info_value(left_value) != normalize_info_value(right_value):
            return False

    return True


def genotype_mismatch_count(
    left: pysam.libcbcf.VariantRecord,
    right: pysam.libcbcf.VariantRecord,
    sample_names: Iterable[str],
) -> Tuple[int, str | None]:
    left_keys = tuple(left.format.keys())
    right_keys = tuple(right.format.keys())

    if not left_keys and not right_keys:
        return 0, None
    if set(left_keys) != set(right_keys):
        sample_list = tuple(sample_names)
        return len(sample_list), sample_list[0] if sample_list else None

    keys = left_keys if left_keys == right_keys else ordered_format_keys(left_keys)

    mismatch_count = 0
    first_sample = None

    for sample_name in sample_names:
        left_sample = left.samples[sample_name]
        right_sample = right.samples[sample_name]

        sample_equal = True
        for key in keys:
            left_raw = left_sample.get(key)
            right_raw = right_sample.get(key)

            if left_raw == right_raw:
                continue

            if normalize_format_value(key, left_raw) != normalize_format_value(key, right_raw):
                sample_equal = False
                break

        if not sample_equal:
            mismatch_count += 1
            if first_sample is None:
                first_sample = sample_name

    return mismatch_count, first_sample


def compare_vcfs(left_path: str, right_path: str) -> int:
    mismatch_count = 0
    compared_records = 0

    left_samples = read_sample_names(left_path)
    right_samples = read_sample_names(right_path)

    if left_samples is None or right_samples is None:
        print("missing_header")
        return 2

    if left_samples != right_samples:
        print("sample_mismatch")
        print(
            f"left_sample_count={len(left_samples)} "
            f"right_sample_count={len(right_samples)}"
        )
        return 2

    with open_variant_lines(left_path) as left_handle, open_variant_lines(right_path) as right_handle:
        while True:
            left_line = next_record_line(left_handle)
            right_line = next_record_line(right_handle)

            if left_line is None and right_line is None:
                break

            compared_records += 1

            if left_line is None or right_line is None:
                mismatch_count += 1
                mismatched_fields = (left_line or right_line or b"").split(b"\t")
                label = "<none>"
                if mismatched_fields and mismatched_fields[0]:
                    label = (
                        mismatched_fields[2].decode("utf-8")
                        if len(mismatched_fields) > 2 and mismatched_fields[2] != b"."
                        else f"{mismatched_fields[0].decode('utf-8')}:{mismatched_fields[1].decode('utf-8')}"
                    )
                print(f"variant={label}\ttype=record_count_mismatch")
                continue

            if left_line == right_line:
                continue

            left_fields = left_line.split(b"\t")
            right_fields = right_line.split(b"\t")

            record_site_equal = (
                left_fields[0] == right_fields[0]
                and left_fields[1] == right_fields[1]
                and left_fields[2] == right_fields[2]
                and left_fields[3] == right_fields[3]
                and left_fields[4] == right_fields[4]
                and left_fields[5] == right_fields[5]
                and normalize_filter_bytes(left_fields[6]) == normalize_filter_bytes(right_fields[6])
                and normalize_info_bytes(left_fields[7]) == normalize_info_bytes(right_fields[7])
            )

            gt_mismatch_count = 0
            gt_example_sample = None
            if len(left_fields) > 8 or len(right_fields) > 8:
                left_fmt = tuple(left_fields[8].split(b":")) if len(left_fields) > 8 and left_fields[8] else ()
                right_fmt = tuple(right_fields[8].split(b":")) if len(right_fields) > 8 and right_fields[8] else ()

                if set(left_fmt) != set(right_fmt):
                    gt_mismatch_count = len(left_samples)
                    gt_example_sample = left_samples[0] if left_samples else None
                else:
                    left_signature = normalized_format_signature(left_fmt)
                    right_signature = normalized_format_signature(right_fmt)
                    sample_count = len(left_samples)
                    for sample_index in range(sample_count):
                        left_sample = left_fields[9 + sample_index] if 9 + sample_index < len(left_fields) else b""
                        right_sample = right_fields[9 + sample_index] if 9 + sample_index < len(right_fields) else b""

                        if left_fmt == right_fmt and left_sample == right_sample:
                            continue

                        if left_signature == right_signature:
                            if normalized_sample_bytes(left_fmt, left_sample) == normalized_sample_bytes(right_fmt, right_sample):
                                continue

                        gt_mismatch_count += 1
                        if gt_example_sample is None:
                            gt_example_sample = left_samples[sample_index]

            if not record_site_equal or gt_mismatch_count > 0:
                mismatch_count += 1
                if left_fields[2] != b".":
                    label = left_fields[2].decode("utf-8")
                else:
                    label = (
                        f"{left_fields[0].decode('utf-8')}:"
                        f"{left_fields[1].decode('utf-8')}:"
                        f"{left_fields[3].decode('utf-8')}:"
                        f"{left_fields[4].decode('utf-8')}"
                    )

                site_differences = describe_site_differences(left_fields, right_fields)

                if gt_mismatch_count > 0:
                    print(
                        f"variant={label}\t"
                        f"site_equal={record_site_equal}\t"
                        f"site_diff_count={len(site_differences)}\t"
                        f"gt_mismatch_count={gt_mismatch_count}\t"
                        f"example_sample={gt_example_sample}"
                    )
                else:
                    print(
                        f"variant={label}\t"
                        "site_equal=false\t"
                        f"site_diff_count={len(site_differences)}\t"
                        "gt_mismatch_count=0"
                    )

                for diff in site_differences:
                    print(f"variant={label}\tsite_diff={diff}")

    print(f"records_compared {compared_records}")
    print(f"mismatch_count {mismatch_count}")
    return 0 if mismatch_count == 0 else 1


def main() -> int:
    parser = argparse.ArgumentParser(description="Compare two VCFs with relaxed formatting rules.")
    parser.add_argument("left_vcf", help="Path to first VCF/BCF")
    parser.add_argument("right_vcf", help="Path to second VCF/BCF")
    args = parser.parse_args()

    return compare_vcfs(args.left_vcf, args.right_vcf)


if __name__ == "__main__":
    main()
