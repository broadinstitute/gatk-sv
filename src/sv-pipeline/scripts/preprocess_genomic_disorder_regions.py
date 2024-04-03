#!/bin/python

import argparse
import sys
import math
from typing import List, Text, Optional

import pysam


# Delimiter suffix appended to the end of interval IDs before the index, e.g. "intervalA__0", "intervalA__1", ...
INDEX_DELIMITER = "__"
LEFT_INDEX_PREFIX = "L"
MIDDLE_INDEX_PREFIX = "M"
RIGHT_INDEX_PREFIX = "R"


def write_subdivisions(fout, ref, n_subdivisions, subdivision_size, start, name, chrom, extra_fields, index_prefix):
    current_pos = start
    for i in range(n_subdivisions):
        current_end = current_pos + subdivision_size
        current_name = f"{name}{INDEX_DELIMITER}{index_prefix}{i}"
        fields = [str(x) for x in [chrom, current_pos, current_end, current_name] + extra_fields]
        if 0 <= current_pos < ref.get_reference_length(chrom) - 1:
            fout.write("\t".join(fields) + "\n")
        current_pos = current_end


def subdivide_intervals(input_path, output_path, reference_path, default_subdivisions, min_size, padding_subdivisions):
    with open(input_path) as fin, open(output_path, "w") as fout, pysam.FastaFile(reference_path) as ref:
        for line in fin:
            record = line.strip().split('\t')
            if record[0].startswith("#"):
                continue
            chrom = record[0]
            start = int(record[1])
            stop = int(record[2])
            name = record[3]
            if INDEX_DELIMITER in name:
                raise ValueError(f"Interval ID cannot contain \"{INDEX_DELIMITER}\": {name}")
            length = stop - start
            subdivision_size = math.ceil(length / default_subdivisions)
            if subdivision_size < min_size:
                n_subdivisions = max(1, length // min_size)
                subdivision_size = round(length / n_subdivisions)
            else:
                n_subdivisions = default_subdivisions
            extra_fields = record[4:] if len(record) > 4 else list()
            # Left padding
            left_start_pos = start - n_subdivisions * subdivision_size
            write_subdivisions(
                fout=fout, ref=ref, n_subdivisions=n_subdivisions, subdivision_size=subdivision_size,
                start=left_start_pos, name=name, chrom=chrom, extra_fields=extra_fields,
                index_prefix=LEFT_INDEX_PREFIX
            )
            # Middle, over the actual region
            middle_start_pos = start
            write_subdivisions(
                fout=fout, ref=ref, n_subdivisions=n_subdivisions, subdivision_size=subdivision_size,
                start=middle_start_pos, name=name, chrom=chrom, extra_fields=extra_fields,
                index_prefix=MIDDLE_INDEX_PREFIX
            )
            # Right padding
            right_start_pos = stop
            write_subdivisions(
                fout=fout, ref=ref, n_subdivisions=n_subdivisions, subdivision_size=subdivision_size,
                start=right_start_pos, name=name, chrom=chrom, extra_fields=extra_fields,
                index_prefix=RIGHT_INDEX_PREFIX
            )


def _parse_arguments(argv: List[Text]) -> argparse.Namespace:
    # noinspection PyTypeChecker
    parser = argparse.ArgumentParser(
        description="Subdivides genomic disorder region intervals",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--input', required=True, type=str, help='Input path with columns CHROM,POS,END,ID (.bed)')
    parser.add_argument('--out', required=True, type=str, help='Output path (.bed)')
    parser.add_argument('--reference', required=True, type=str, help='Indexed reference path (.fasta)')
    parser.add_argument('--subdivisions', type=int, default=10, help='Subdivisions per interval')
    parser.add_argument('--min-size', type=int, default=10000, help='Minimum size for resulting intervals')
    parser.add_argument('--padding-subdivisions', type=float, default=1.0,
                        help='Number of subdivisions to pad on each side')
    if len(argv) <= 1:
        parser.parse_args(["--help"])
        sys.exit(0)
    parsed_arguments = parser.parse_args(argv[1:])
    return parsed_arguments


def main(argv: Optional[List[Text]] = None):
    if argv is None:
        argv = sys.argv
    args = _parse_arguments(argv)
    subdivide_intervals(input_path=args.input, output_path=args.out,
                        reference_path=args.reference, default_subdivisions=args.subdivisions,
                        min_size=args.min_size, padding_subdivisions=args.padding_subdivisions)


if __name__ == "__main__":
    main()
