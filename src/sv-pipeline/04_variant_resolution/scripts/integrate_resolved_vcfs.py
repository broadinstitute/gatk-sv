#!/bin/env python

import argparse
import sys
from typing import Optional, List, Text, Dict, Set, Callable

import pysam

MEMBERS_KEY = "MEMBERS"
UNRESOLVED_KEY = "UNRESOLVED"


def update_header(header: pysam.VariantHeader) -> None:
    header.add_line('##INFO=<ID=UNRESOLVED,Number=0,Type=Flag,Description="Variant is unresolved.">')
    header.add_line('##INFO=<ID=UNRESOLVED_TYPE,Number=1,Type=String,Description=\"Class of unresolved variant.\">')


def is_unresolved(record: pysam.VariantRecord):
    return record.info.get(UNRESOLVED_KEY, None)


def is_resolved(record: pysam.VariantRecord):
    return not is_unresolved(record)


def get_members(record: pysam.VariantRecord):
    return list(record.info[MEMBERS_KEY]) if isinstance(record.info[MEMBERS_KEY], tuple) \
        else [record.info[MEMBERS_KEY]] if record.info[MEMBERS_KEY] is not None \
        else list()


def get_vids_and_members_sets(vcf: pysam.VariantFile,
                              predicate: Callable) -> Dict:
    unresolved_vids_set = set()
    unresolved_members_set = set()
    for r in vcf:
        if predicate(r):
            unresolved_vids_set.add(r.id)
            unresolved_members_set.update(get_members(r))
    vcf.reset()
    return unresolved_vids_set, unresolved_members_set


def write_vcf(header: pysam.VariantHeader,
              all_vcf: pysam.VariantFile,
              inv_vcf: pysam.VariantFile,
              inv_resolved_vids_set: Set,
              inv_resolved_members_set: Set,
              all_unresolved_vids_set: Set,
              all_unresolved_members_set: Set) -> None:
    sys.stdout.write(str(header))
    for r in all_vcf:
        if r.id not in all_unresolved_vids_set or r.id not in inv_resolved_members_set:
            # Resolved in ALL vcf, or unresolved in both VCFs
            sys.stdout.write(str(r))
    for r in inv_vcf:
        if r.id in inv_resolved_vids_set:
            # Resolved variant in the INV vcf
            members = get_members(r)
            if all((m in all_unresolved_members_set) for m in members):
                # Resolved in the INV vcf and every member unresolved in the ALL vcf
                sys.stdout.write(str(r))


def __parse_arguments(argv: List[Text]) -> argparse.Namespace:
    # noinspection PyTypeChecker
    parser = argparse.ArgumentParser(
        description="Integrates inversion-only and all-SV VCFs from the complex resolve module. "
                    "Unsorted output is written to stdout.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--all-vcf", type=str, required=True,
                        help="Complex-resolved VCF containing all SVs")
    parser.add_argument("--inv-only-vcf", type=str, required=True,
                        help="Complex-resolved VCF containing only inversions")
    if len(argv) <= 1:
        parser.parse_args(["--help"])
        sys.exit(0)
    parsed_arguments = parser.parse_args(argv[1:])
    return parsed_arguments


def main(argv: Optional[List[Text]] = None):
    if argv is None:
        argv = sys.argv
    arguments = __parse_arguments(argv)
    with pysam.VariantFile(arguments.all_vcf) as all_vcf, \
            pysam.VariantFile(arguments.inv_only_vcf) as inv_vcf:
        header = all_vcf.header
        update_header(header)
        inv_resolved_vids_set, inv_resolved_members_set = get_vids_and_members_sets(inv_vcf, is_resolved)
        all_unresolved_vids_set, all_unresolved_members_set = get_vids_and_members_sets(all_vcf, is_unresolved)
        write_vcf(header=header,
                  all_vcf=all_vcf,
                  inv_vcf=inv_vcf,
                  inv_resolved_vids_set=inv_resolved_vids_set,
                  inv_resolved_members_set=inv_resolved_members_set,
                  all_unresolved_vids_set=all_unresolved_vids_set,
                  all_unresolved_members_set=all_unresolved_members_set)


if __name__ == "__main__":
    main()
