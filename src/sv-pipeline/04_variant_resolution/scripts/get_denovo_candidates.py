#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""

"""


import argparse
import sys
import pysam
import svtk.utils as svu
from svtk.famfile import parse_famfile


def get_denovo_candidates(record, fam, max_parents=10):
    """
    Obtain list of samples which are putatively called de novo
    """
    called = svu.get_called_samples(record)
    parents = [s for s in called if fam.samples[s].is_parent]

    if len(parents) > max_parents:
        return []

    denovo = []
    for ID in called:
        sample = fam.samples[ID]
        if sample.has_parents:
            if sample.mother not in called and sample.father not in called:
                denovo.append(sample.ID)

    return denovo


def filter_denovo_records(vcf, fam, max_parents=10):
    for record in vcf:
        # Skip records without any Mendelian violations
        candidates = get_denovo_candidates(record, fam, max_parents)
        if len(candidates) == 0:
            continue

        # Restrict to rare (parental VF<0.1%) variants
        called = svu.get_called_samples(record)
        parents = [s for s in called if fam.samples[s].is_parent]
        if len(parents) > max_parents:
            continue

        # Skip non-stranded (wham)
        if 'STRANDS' in record.info.keys():
            if record.info['STRANDS'] not in '+- -+ ++ --'.split():
                continue

        yield record


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf')
    parser.add_argument('famfile', type=argparse.FileType('r'))
    parser.add_argument('fout')
    parser.add_argument('--max-parents', type=int, default=10)
    args = parser.parse_args()

    if args.vcf in '- stdin'.split():
        vcf = pysam.VariantFile(sys.stdin)
    else:
        vcf = pysam.VariantFile(args.vcf)

    if args.fout in '- stdout'.split():
        fout = pysam.VariantFile(sys.stdout, 'w', header=vcf.header)
    else:
        fout = pysam.VariantFile(args.fout, 'w', header=vcf.header)

    fam = parse_famfile(args.famfile)

    for record in filter_denovo_records(vcf, fam, args.max_parents):
        fout.write(record)


if __name__ == '__main__':
    main()
