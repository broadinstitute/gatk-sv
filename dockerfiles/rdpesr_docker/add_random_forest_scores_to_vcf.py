#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Contact: Xuefang Zhao <xzhao@broadinstitute.org> 

"""
Annotate sample-level Boost scores directly into an input VCF
"""


import csv
import numpy as np
import argparse
import pysam
import sys


def load_RF_scores(inpath):
    """
    Loads Boost scores per variant for a single sample from .tsv
    """

    scores_pbsv = {}
    scores_vapor = {}
    with open(inpath) as scores_tsv:
        for SVID, PBSV_neg, PBSV_pos, VaPoR_neg, VaPoR_pos in csv.reader(scores_tsv, delimiter='\t'):
            if not SVID=="SVID":
                scores_pbsv[SVID] = np.round(float(PBSV_pos), 3)
                scores_vapor[SVID] = np.round(float(VaPoR_pos), 3)

    return [scores_pbsv, scores_vapor]


def main():
    """
    Command-line main block
    """

    # Parse command line arguments and options
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--vcf', required=True, help='VCF to be annotated. Required.')
    parser.add_argument('--RF_tsv', required=True, help='two-column .tsv of ' +
                        'sample ID and path to RF scores for that sample. ' +
                        'Required.')
    parser.add_argument('-o', '--outfile', default='stdout', help='Path to output ' +
                        'VCF. Default: stdout.')
    args = parser.parse_args()

    #Open connection to input VCF
    if args.vcf in '- stdin'.split():
        vcf = pysam.VariantFile(sys.stdin) 
    else:
        vcf = pysam.VariantFile(args.vcf)

    # Add new FORMAT line to VCF header
    vcf.header.add_meta('FORMAT', 
                        items=[('ID', "PBRF"), ('Number', "1"), ('Type', "Float"), 
                               ('Description', "RF score trained from Pacbio SV callset (Ebert et al)")])
    vcf.header.add_meta('FORMAT', 
                        items=[('ID', "VaRF"), ('Number', "1"), ('Type', "Float"), 
                               ('Description', "RF score trained from VaPoR evaluation")])


    # Load dict of each sample's boost scores
    RF_pbsv_lookup = {}
    RF_vapor_lookup = {}
    with open(args.RF_tsv) as RF_tsv:
        for sample, scores_path in csv.reader(RF_tsv, delimiter='\t'):
            score_lookups = load_RF_scores(scores_path)
            RF_pbsv_lookup[sample] = score_lookups[0]
            RF_vapor_lookup[sample] = score_lookups[1]
    RF_samples = set(RF_pbsv_lookup.keys())

    #Open connection to output VCF
    if args.outfile in '- stdout'.split():
        outvcf = pysam.VariantFile(sys.stdout, 'w', header=vcf.header)
    else:
        outvcf = pysam.VariantFile(args.outfile, 'w', header=vcf.header)

    # Annotate each record
    for record in vcf:
        for sample in record.samples:
            bs_1 = RF_pbsv_lookup.get(sample, {}).get(record.id, None)
            record.samples[sample]['PBRF'] = bs_1
            bs_2 = RF_vapor_lookup.get(sample, {}).get(record.id, None)
            record.samples[sample]['VaRF'] = bs_2

        outvcf.write(record)
    outvcf.close()


if __name__ == '__main__':
    main()