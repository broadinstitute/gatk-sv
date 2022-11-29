#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2022 The Broad Institute of M.I.T. and Harvard
# Distributed under terms of the MIT license.
# Contact: Ryan Collins <rlcollins@g.harvard.edu>

"""
Annotate no-call GT rates for all records in a VCF
"""


import numpy as np
import csv
import argparse
import pysam
import sys


def is_mcnv(record):
    """
    Checks if a record is multiallelic
    """

    if 'MULTIALLELIC' in record.filter.keys() \
    or '<CNV>' in record.alleles \
    or len(record.alleles) > 2 \
    or record.info['SVTYPE'] == 'CNV':
        return True
    else:
        return False


def add_ncrs(record, subsets):
    """
    Add no-call GT rates to record INFO field
    """

    nocalls = 0
    subset_nocalls = {prefix : 0 for prefix in subsets.keys()}
    
    # Process each sample in serial
    for sample in record.samples.keys():
        if record.samples[sample]['GT'] == (None, None):
            nocalls += 1
            for prefix, subdat in subsets.items():
                if sample in subdat['samples']:
                    subset_nocalls[prefix] += 1

    # Annotate no-call rate (only for biallelic variants)
    if not is_mcnv(record):
        record.info['NCR'] = nocalls / len(record.samples.keys())
        for prefix, sub_nc in subset_nocalls.items():
            sub_ns = subsets[prefix]['N']
            record.info['{}_NCR'.format(prefix)] = sub_nc / sub_ns

    return record


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf_in', help='VCF to annotate.')
    parser.add_argument('vcf_out', help='Path to output VCF. Accepts "-" ' +
                        'and "stdout". Default: stdout', default='stdout')
    parser.add_argument('--sample-subsets', help='Two-column .tsv with sample ' +
                        'subsets for which NCR should be annotated separately. ' +
                        'First column: subset prefix. Second column: path to .txt ' +
                        'file with all sample IDs in subset.')
    parser.add_argument('--stats-tsv', help='Optional path to output .tsv file ' +
                        'with NCR stats for all records in vcf_in.')
    args = parser.parse_args()

    # Open connection to input VCF
    vcf = pysam.VariantFile(args.vcf_in)
    allsamps = set(vcf.header.samples)

    # Read sample subsets, if any
    subsets = {}
    if args.sample_subsets is not None:
        with open(args.sample_subsets) as ssin:
            for prefix, slpath in csv.reader(ssin, delimiter='\t'):
                samps = set([s.rstrip() for s in open(slpath).readlines()])
                samps = samps.intersection(allsamps)
                subsets[prefix] = {'samples' : samps, 'N' : len(samps)}

    # Update VCF header as needed
    vcf.header.add_meta('INFO',
                        items=[('ID', "NCR"), ('Number', "1"), ('Type', "Float"),
                               ('Description', "Proportion of no-call GTs")])
    for prefix in subsets.keys():
        descrip = "Proportion of no-call GTs in {} samples".format(prefix)
        vcf.header.add_meta('INFO',
                            items=[('ID', "{}_NCR".format(prefix)), 
                                   ('Number', "1"), ('Type', "Float"),
                                   ('Description', descrip)])

    # Open connection to output VCF
    if args.vcf_out in '- stdout'.split():
        outvcf = pysam.VariantFile(sys.stdout, 'w', header=vcf.header)
    else:
        outvcf = pysam.VariantFile(args.vcf_out, 'w', header=vcf.header)

    # Open connection to output .tsv, if optioned
    if args.stats_tsv is not None:
        stats_out = open(args.stats_tsv, 'w')

    # Process records
    for record in vcf:
        record = add_ncrs(record, subsets)
        outvcf.write(record)
        if args.stats_tsv is not None:
            stats = [record.id]
            for key in ['NCR'] + [x + '_NCR' for x in subsets.keys()]:
                if key in record.info.keys():
                    stats.append(str(np.round(record.info[key], 8)))
                else:
                    stats.append('.')
            stats_out.write('\t'.join(stats) + '\n')

    # Close connection to output VCF and (optionally) stats .tsv
    outvcf.close()
    if args.stats_tsv is not None:
        stats_out.close()


if __name__ == '__main__':
    main()

