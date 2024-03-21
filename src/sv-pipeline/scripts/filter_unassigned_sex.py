#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

import argparse
import sys
import pysam


def _parse_arguments(argv):
    # noinspection PyTypeChecker
    parser = argparse.ArgumentParser(
        description="Sets depth-only sex chromosome genotypes to no-call for samples with unassigned sex",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--vcf', type=str, required=True, help='Input vcf')
    parser.add_argument('--ped-file', type=str, required=True, help='Ped file')
    parser.add_argument('--out', type=str, required=True, help='Output vcf')
    parser.add_argument('--chr-x', type=str, default="chrX", help='Chromosome X name')
    parser.add_argument('--chr-y', type=str, default="chrY", help='Chromosome Y name')
    if len(argv) <= 1:
        parser.parse_args(["--help"])
        sys.exit(0)
    parsed_arguments = parser.parse_args(argv[1:])
    return parsed_arguments


def main(argv=None):
    if argv is None:
        argv = sys.argv
    args = _parse_arguments(argv)

    # Parse ped file
    with open(args.ped_file) as f:
        unassigned_sex_samples = set()
        for line in f:
            tokens = line.strip().split("\t")
            if tokens[4] != "1" and tokens[4] != "2":
                unassigned_sex_samples.add(tokens[1])

    # Filter genotypes
    with pysam.VariantFile(args.vcf) as fin, pysam.VariantFile(args.out, mode="w", header=fin.header) as fout:
        samples_to_reset = list(unassigned_sex_samples.intersection(set(fin.header.samples)))
        for record in fin:
            if record.chrom == args.chr_x or record.chrom == args.chr_y:
                for s in samples_to_reset:
                    gt = record.samples[s]
                    gt["GT"] = tuple(None for _ in gt["GT"])
                    gt["GQ"] = None
                    gt["EV"] = None
                    gt["RD_CN"] = None
                    gt["RD_GQ"] = None
            fout.write(record)


if __name__ == '__main__':
    main()
