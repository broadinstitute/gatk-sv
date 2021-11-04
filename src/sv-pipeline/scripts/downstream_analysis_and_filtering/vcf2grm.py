#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""
Creates a genetic relatedness matrix (GRM) from biallelic sites in an input VCF
"""

import argparse
import sys
import pysam


def get_allele_dosage(record):
    ads = []
    for s in record.samples:
        GT = list(record.samples[s]['GT'])
        if len([a for a in GT if a is None]) > 0:
            ad = 'NA'
        else:
            ad = str(sum([a for a in GT]))
        ads.append(ad)

    grmline = '\t'.join(map(str, ads))

    return grmline


def _is_multiallelic(record):
    status = False

    if 'MULTIALLELIC' in record.filter \
            or 'MULTIALLELIC' in record.info.keys() \
            or len(record.alleles) - 1 > 1:
        status = True

    return status


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf', help='Input vcf (supports "stdin").')
    parser.add_argument('fout', help='Output file (supports "stdout").')

    args = parser.parse_args()

    if args.vcf in '- stdin'.split():
        vcf = pysam.VariantFile(sys.stdin)
    else:
        vcf = pysam.VariantFile(args.vcf)

    if args.fout in '- stdout'.split():
        fout = sys.stdout
    else:
        fout = open(args.fout, 'w')

    # Write header to outfile
    samples = [s for s in vcf.header.samples]
    outheader = '{0}\t{1}\n'.format('#VID', '\t'.join(samples))
    fout.write(outheader)

    # Iterate over records in vcf and count allele dosage per sample
    for record in vcf.fetch():
        # Do not process multiallelic variants, unless optioned
        if not _is_multiallelic(record):
            ad = get_allele_dosage(record)
            newline = '{0}\t{1}\n'.format(record.id, ad)
            fout.write(newline)

    fout.close()


if __name__ == '__main__':
    main()
