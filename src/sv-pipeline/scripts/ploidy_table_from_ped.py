#!/bin/env python

import argparse
import sys
from collections import defaultdict
from typing import Optional, List, Text


def create_header(contigs: List[Text]) -> Text:
    """
    Creates header for the table

    Parameters
    ----------
    contigs: List[Text]
        ordered list of contigs

    Returns
    -------
    Text
        header line
    """
    return '\t'.join(['SAMPLE'] + contigs)


def convert_ped_record(ped_record: Text,
                       contigs: List[Text],
                       chr_x: Text = 'chrX',
                       chr_y: Text = 'chrY') -> Text:
    """
    Converts a ped file record to a table record.

    Parameters
    ----------
    ped_record: Text
        ped file record
    contigs: List[Text]
        ordered list of contigs
    chr_x: Text = 'chrX'
        chromosome X name
    chr_y: Text = 'chrY'
        chromosome Y name

    Returns
    -------
    Text
        ploidy table record
    """
    tokens = ped_record.strip().split('\t')
    sample = tokens[1]
    ploidy = defaultdict(lambda: 2)
    if tokens[4] == "1":
        ploidy[chr_x] = 1
        ploidy[chr_y] = 1
    elif tokens[4] == "2":
        ploidy[chr_x] = 2
        ploidy[chr_y] = 0
    else:
        ploidy[chr_x] = 0
        ploidy[chr_y] = 0
    return "\t".join([sample] + [str(ploidy[c]) for c in contigs])


def __read_contigs(path: Text) -> List[Text]:
    with open(path, 'r') as f:
        return [line.strip().split('\t')[0] for line in f]


def __parse_arguments(argv: List[Text]) -> argparse.Namespace:
    # noinspection PyTypeChecker
    parser = argparse.ArgumentParser(
        description="Create a ploidy table from a PED file",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--ped", type=str, required=True,
                        help="PED file")
    parser.add_argument("--contigs", type=str, required=True,
                        help="Ordered list of contigs")
    parser.add_argument("--out", type=str, required=True,
                        help="Output VCF")
    parser.add_argument("--chr-x", type=str, default="chrX",
                        help="Chromosome X name")
    parser.add_argument("--chr-y", type=str, default="chrY",
                        help="Chromosome Y name")
    if len(argv) <= 1:
        parser.parse_args(["--help"])
        sys.exit(0)
    parsed_arguments = parser.parse_args(argv[1:])
    return parsed_arguments


def main(argv: Optional[List[Text]] = None):
    if argv is None:
        argv = sys.argv
    arguments = __parse_arguments(argv)
    contigs = __read_contigs(arguments.contigs)
    with open(arguments.ped, 'r') as ped, open(arguments.out, 'w') as out:
        out.write(create_header(contigs=contigs) + "\n")
        for line in ped:
            if line.startswith('#'):
                # skip comments / headers
                continue
            out.write(convert_ped_record(
                ped_record=line,
                contigs=contigs,
                chr_x=arguments.chr_x,
                chr_y=arguments.chr_y
            ) + "\n")


if __name__ == "__main__":
    main()
