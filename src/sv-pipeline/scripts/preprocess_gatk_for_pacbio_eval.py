#!/bin/env python

import argparse
import pysam
import sys
from typing import List, Text, Optional


def _parse_arguments(argv: List[Text]) -> argparse.Namespace:
    # noinspection PyTypeChecker
    parser = argparse.ArgumentParser(
        description="Prepare call set for evaluation against PacBio variants, including subset to DEL/DUP/INS, "
                    "converting DUP to INS, and filtering variants over 5kbp",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("vcf", type=str,
                        help="Input VCF. Usually this will be the cleaned vcf.")
    parser.add_argument("--max-variant-size", type=int, default=5000,
                        help="Maximum size of variants to include in the output VCF.")

    if len(argv) <= 1:
        parser.parse_args(["--help"])
        sys.exit(0)
    parsed_arguments = parser.parse_args(argv[1:])
    return parsed_arguments


def add_header_lines(header: pysam.VariantHeader) -> pysam.VariantHeader:
    """
    Ingests the given header, removes specified fields, and adds necessary fields.

    Parameters
    ----------
    header: pysam.VariantHeader
        input header
    """
    # new fields
    header.add_line('##INFO=<ID=ORIGINAL_SVTYPE,Number=1,Type=String,Description="Original type of structural variant">')
    header.add_line('##INFO=<ID=ORIGINAL_ALT,Number=1,Type=String,Description="Original ALT allele">')


def process(vcf: pysam.VariantFile) -> None:
    add_header_lines(vcf.header)
    sys.stdout.write(str(vcf.header))
    allowed_svtypes = set(['DEL', 'DUP', 'INS', 'INV'])
    for record in vcf:
        svtype = record.info['SVTYPE']
        if svtype not in allowed_svtypes:
            continue
        svlen = record.info.get('SVLEN', record.stop - record.pos)
        if svlen > max_variant_size:
            continue
        record.info['ORIGINAL_SVTYPE'] = svtype
        record.info['ORIGINAL_ALT'] = record.alts[0]
        if svtype != 'DUP':
            sys.stdout.write(str(record))
        else:
            record.info['SVLEN'] = record.stop - record.pos
            record.info['SVTYPE'] = 'INS'
            record.stop = record.pos + 1
            sys.stdout.write(str(record).replace("\t<DUP>", "\t<INS>"))


def main(argv: Optional[List[Text]] = None):
    if argv is None:
        argv = sys.argv
    arguments = _parse_arguments(argv)

    # convert vcf header and records
    with pysam.VariantFile(arguments.vcf) as vcf:
        process(vcf, arguments.max_variant_size)


if __name__ == "__main__":
    main()
