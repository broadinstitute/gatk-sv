#!/bin/env python

import argparse
import pysam
import sys
from typing import Any, List, Text, Set, Dict, Optional


def _parse_arguments(argv: List[Text]) -> argparse.Namespace:
    # noinspection PyTypeChecker
    parser = argparse.ArgumentParser(
        description="Prepare call set for evaluation against PacBio variants, including subset to DEL/DUP/INS, "
                    "converting DUP to INS, and filtering variants over 5kbp",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("vcf", type=str,
                        help="Input VCF. Usually this will be the cleaned vcf.")

    if len(argv) <= 1:
        parser.parse_args(["--help"])
        sys.exit(0)
    parsed_arguments = parser.parse_args(argv[1:])
    return parsed_arguments


def process(vcf: pysam.VariantFile) -> None:
    sys.stdout.write(str(vcf.header))
    allowed_svtypes = set(['DEL', 'DUP', 'INS', 'INV'])
    for record in vcf:
        svtype = record.info['SVTYPE']
        if svtype not in allowed_svtypes:
            continue
        svlen = record.info.get('SVLEN', record.stop - record.pos)
        if svlen > 5000:
            continue
        if svtype != 'DUP':
            sys.stdout.write(str(record))
        else:
            record.info['SVLEN'] = record.stop - record.pos
            record.info['SVTYPE'] = 'INS'
            record.stop = record.pos + 1
            sys.stdout.write(str(record).replace("<DUP>", "<INS>"))


def main(argv: Optional[List[Text]] = None):
    if argv is None:
        argv = sys.argv
    arguments = _parse_arguments(argv)

    # convert vcf header and records
    with pysam.VariantFile(arguments.vcf) as vcf:
        process(vcf)


if __name__ == "__main__":
    main()
