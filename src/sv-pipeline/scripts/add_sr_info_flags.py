#!/bin/env python

import argparse
import pysam
import sys
from typing import Optional, List, Text, Set


def read_variant_ids(path: Text) -> Set[Text]:
    """
    Reads a list of variant IDs from a file (one per line).

    Parameters
    ----------
    path: Text
        path to file containing variant IDs

    Returns
    -------
    ids: Set[Text]
        set of variant IDs
    """
    with open(path, 'r') as f:
        return set(line.strip() for line in f if line.strip())


def process_record(record: pysam.VariantRecord,
                   vcf_out: pysam.VariantFile,
                   bothsides_support_ids: Set[Text],
                   high_sr_background_ids: Set[Text]) -> pysam.VariantRecord:
    """
    Processes a record and sets INFO flags based on variant ID membership in the provided sets.

    Parameters
    ----------
    record: pysam.VariantRecord
        input record
    vcf_out: pysam.VariantFile
        output vcf file
    bothsides_support_ids: Set[Text]
        set of variant IDs that should have BOTHSIDES_SUPPORT flag set
    high_sr_background_ids: Set[Text]
        set of variant IDs that should have HIGH_SR_BACKGROUND flag set

    Returns
    -------
    record: pysam.VariantRecord
        record with flags set appropriately
    """

    # Set BOTHSIDES_SUPPORT flag if variant ID is in the set
    if record.id in bothsides_support_ids:
        record.info['BOTHSIDES_SUPPORT'] = True
    else:
        record.info['BOTHSIDES_SUPPORT'] = False

    # Set HIGH_SR_BACKGROUND flag if variant ID is in the set
    if record.id in high_sr_background_ids:
        record.info['HIGH_SR_BACKGROUND'] = True
    else:
        record.info['HIGH_SR_BACKGROUND'] = False

    return record


def __parse_arguments(argv: List[Text]) -> argparse.Namespace:
    # noinspection PyTypeChecker
    parser = argparse.ArgumentParser(
        description="Add BOTHSIDES_SUPPORT and HIGH_SR_BACKGROUND INFO flags to variants based on ID lists.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--vcf", type=str, required=True,
                        help="Input VCF file")
    parser.add_argument("--out", type=str, required=True,
                        help="Output VCF file")
    parser.add_argument("--bothsides-support-list", type=str, required=False,
                        help="File containing variant IDs that should have BOTHSIDES_SUPPORT flag set (one per line)")
    parser.add_argument("--high-sr-background-list", type=str, required=False,
                        help="File containing variant IDs that should have HIGH_SR_BACKGROUND flag set (one per line)")
    if len(argv) <= 1:
        parser.parse_args(["--help"])
        sys.exit(0)
    parsed_arguments = parser.parse_args(argv[1:])
    return parsed_arguments


def main(argv: Optional[List[Text]] = None):
    if argv is None:
        argv = sys.argv
    arguments = __parse_arguments(argv)

    # Read variant ID lists
    bothsides_support_ids = set()
    if arguments.bothsides_support_list:
        bothsides_support_ids = read_variant_ids(arguments.bothsides_support_list)

    high_sr_background_ids = set()
    if arguments.high_sr_background_list:
        high_sr_background_ids = read_variant_ids(arguments.high_sr_background_list)

    # Process VCF
    with pysam.VariantFile(arguments.vcf) as vcf_in:
        with pysam.VariantFile(arguments.out, mode='w', header=vcf_in.header) as vcf_out:
            for record in vcf_in:
                new_record = process_record(
                    record=record,
                    vcf_out=vcf_out,
                    bothsides_support_ids=bothsides_support_ids,
                    high_sr_background_ids=high_sr_background_ids
                )
                vcf_out.write(new_record)


if __name__ == "__main__":
    main()
