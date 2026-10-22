#!/bin/env python

import argparse
import pysam
import sys
from typing import Optional, List, Text


def __parse_arguments(argv: List[Text]) -> argparse.Namespace:
    # noinspection PyTypeChecker
    parser = argparse.ArgumentParser(
        description="Extract FORMAT field values from VCF to TSV table.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--vcf", type=str, required=True,
                        help="Input VCF file")
    parser.add_argument("--out", type=str, required=True,
                        help="Output TSV file")
    parser.add_argument("--format-field", type=str, required=True,
                        help="FORMAT field ID to extract (e.g., GT)")
    parser.add_argument("--id-column", type=str, default="ID",
                        help="Name of the variant ID column")
    if len(argv) <= 1:
        parser.parse_args(["--help"])
        sys.exit(0)
    parsed_arguments = parser.parse_args(argv[1:])
    return parsed_arguments


def main(argv: Optional[List[Text]] = None):
    if argv is None:
        argv = sys.argv
    arguments = __parse_arguments(argv)

    with pysam.VariantFile(arguments.vcf) as vcf_in:
        samples = list(vcf_in.header.samples)

        with open(arguments.out, 'w') as out:
            # Write header
            header_fields = ["chr", "start", "end", arguments.id_column] + samples
            out.write("\t".join(header_fields) + "\n")

            # Process records
            for record in vcf_in:
                row = [
                    record.contig,
                    str(record.start),  # pysam already uses 0-based start
                    str(record.stop),
                    record.id if record.id else "."
                ]

                # Extract FORMAT field for each sample
                for sample in samples:
                    value = record.samples[sample].get(arguments.format_field, None)
                    if value is None:
                        row.append(".")
                    else:
                        row.append(str(value))

                out.write("\t".join(row) + "\n")


if __name__ == "__main__":
    main()
