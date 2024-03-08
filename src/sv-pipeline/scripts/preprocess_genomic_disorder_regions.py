#!/bin/python

import argparse
import sys
import math
from typing import List, Text, Optional


# Delimiter suffix appended to the end of interval IDs before the index, e.g. "intervalA__0", "intervalA__1", ...
INDEX_DELIMITER = "__"


def subdivide_intervals(input_path, output_path, default_subdivisions, min_size):
    with open(input_path) as fin, open(output_path, "w") as fout:
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
                subdivisions = math.ceil(length / min_size)
                subdivision_size = round(length / subdivisions)
            else:
                subdivisions = default_subdivisions
            current_pos = start
            extra_fields = record[4:] if len(record) > 4 else list()
            for i in range(subdivisions):
                current_end = current_pos + subdivision_size
                current_name = f"{name}{INDEX_DELIMITER}{i}"
                fields = [str(x) for x in [chrom, current_pos, current_end, current_name] + extra_fields]
                fout.write("\t".join(fields) + "\n")
                current_pos = current_end


def _parse_arguments(argv: List[Text]) -> argparse.Namespace:
    # noinspection PyTypeChecker
    parser = argparse.ArgumentParser(
        description="Subdivides genomic disorder region intervals",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--input', type=str, help='Input path with columns CHROM,POS,END,ID (.bed)')
    parser.add_argument('--out', type=str, help='Output path (.bed)')
    parser.add_argument('--subdivisions', type=int, default=10, help='Subdivisions per interval')
    parser.add_argument('--min-size', type=int, default=10000, help='Minimum size for resulting intervals')
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
                        default_subdivisions=args.subdivisions, min_size=args.min_size)


if __name__ == "__main__":
    main()
