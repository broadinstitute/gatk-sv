#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Annotate sample-level Boost scores directly into an input VCF
"""

import csv
import numpy as np
import argparse
import pysam
import sys


def load_boost_scores(inpath):
    """
    Loads Boost scores per variant for a single sample from .tsv
    """

    scores = {}

    with open(inpath) as scores_tsv:
        for vid, score in csv.reader(scores_tsv, delimiter='\t'):
            scores[vid] = np.round(float(score), 3)

    return scores


def main():
    """
    Command-line main block
    """

    # Parse command line arguments and options
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--vcf', required=True, help='VCF to be annotated. Required.')
    parser.add_argument('--boost-tsv', required=True, help='two-column .tsv of ' +
                                                           'sample ID and path to Boost scores for that sample. ' +
                                                           'Required.')
    parser.add_argument('-o', '--outfile', default='stdout', help='Path to output ' +
                                                                  'VCF. Default: stdout.')
    args = parser.parse_args()

    # Open connection to input VCF
    if args.vcf in '- stdin'.split():
        vcf = pysam.VariantFile(sys.stdin)
    else:
        vcf = pysam.VariantFile(args.vcf)

    # Add new FORMAT line to VCF header
    vcf.header.add_meta('FORMAT',
                        items=[('ID', "BS"), ('Number', "1"), ('Type', "Float"),
                               ('Description', "lgBoost score")])

    # Load dict of each sample's boost scores
    boost_lookup = {}
    with open(args.boost_tsv) as boost_tsv:
        for sample, scores_path in csv.reader(boost_tsv, delimiter='\t'):
            boost_lookup[sample] = load_boost_scores(scores_path)

    # Open connection to output VCF
    if args.outfile in '- stdout'.split():
        outvcf = pysam.VariantFile(sys.stdout, 'w', header=vcf.header)
    else:
        outvcf = pysam.VariantFile(args.outfile, 'w', header=vcf.header)

    # Annotate each record
    for record in vcf:
        for sample in record.samples:
            bs = boost_lookup.get(sample, {}).get(record.id, None)
            record.samples[sample]['BS'] = bs
        outvcf.write(record)
    outvcf.close()


if __name__ == '__main__':
    main()
