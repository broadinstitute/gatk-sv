#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""
Clean records in VCF based on de novo filter results
"""

import argparse
import sys
from collections import defaultdict
import pysam
import svtk.utils as svu
from svtk.famfile import parse_famfile


def set_null(record, sample):
    dat = record.samples[sample].items()

    # Set genotype to no-call
    record.samples[sample]['GT'] = (0, 0)

    for fmt, value in dat:
        if fmt == 'GT':
            continue

        # Get type and count of FORMAT
        n = record.format[fmt].number
        dtype = record.format[fmt].type

        # Set null value based on FORMAT type
        if dtype in 'Integer Float'.split():
            null_val = 0
        elif dtype in 'String Character'.split():
            null_val = ''
        else:
            raise ValueError('Invalid VCF FORMAT type: {0}'.format(dtype))

        # Set the appropriate count of values
        if n == 1:
            record.samples[sample][fmt] = null_val
        elif n == '.':
            record.samples[sample][fmt] = (null_val, )
        else:
            record.samples[sample][fmt] = tuple(null_val for i in range(n))


def parse_filtered(filterfile, fam):
    """Get dictionaries of added/removed samples from de novo filter"""

    add = defaultdict(list)
    remove = defaultdict(list)

    for line in filterfile:
        name, family, called, support = line.strip().split()
        called = called.split(',')

        samples = [s.ID for s in fam.families[family]]

        for sample in samples:
            if sample in called:
                add[name].append(sample)
            else:
                remove[name].append(sample)

    return add, remove


def filter_dn_variants(vcf, filterfile, fam, fout):
    """
    Add parent false negatives and remove child false positives from dn filter
    Arguments
    ---------
    vcf : pysam.VariantFile
    filterfile : file
    fout : pysam.VariantFile
    """

    # Get dictionaries of samples to add and remove from each variant
    add, remove = parse_filtered(filterfile, fam)

    for record in vcf:
        # Write records unaltered if they weren't included in the de novo check
        if record.id not in add.keys() and record.id not in remove.keys():
            fout.write(record)

        # Otherwise set samples appropriately
        else:
            for sample in add.get(record.id, []):
                record.samples[sample]['GT'] = (0, 1)

            for sample in remove.get(record.id, []):
                set_null(record, sample)

            # Only report record if any samples made it through de novo check
            if len(svu.get_called_samples(record)) > 0:
                fout.write(record)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf')
    parser.add_argument('filtered', type=argparse.FileType('r'),
                        help='De novo filter results')
    parser.add_argument('famfile', type=argparse.FileType('r'))
    parser.add_argument('fout', help='Filtered VCF')
    args = parser.parse_args()

    vcf = pysam.VariantFile(args.vcf)
    fam = parse_famfile(args.famfile)

    fout = sys.stdout if args.fout in 'stdout -'.split() else args.fout
    fout = pysam.VariantFile(fout, 'w', header=vcf.header)

    filter_dn_variants(vcf, args.filtered, fam, fout)


if __name__ == '__main__':
    main()
