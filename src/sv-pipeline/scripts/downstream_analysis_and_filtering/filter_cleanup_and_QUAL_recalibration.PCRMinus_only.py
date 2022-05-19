#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""
Apply final FILTER cleanup and QUAL score recalibration
"""


import argparse
import sys
import pysam
import csv
from numpy import median
from svtk.utils import is_biallelic


# Define global variables
filts_for_info = 'PESR_GT_OVERDISPERSION HIGH_SR_BACKGROUND BOTHSIDES_SUPPORT VARIABLE_ACROSS_BATCHES'.split(
    ' ')
filts_to_remove = 'HIGH_PCRPLUS_NOCALL_RATE HIGH_PCRMINUS_NOCALL_RATE'.split(
    ' ')
filts_to_remove = filts_to_remove + filts_for_info
NULL_GTs = [(None, None), (None, )]
REF_GTs = [(0, 0), (0, ), (None, 2)]
NULL_and_REF_GTs = NULL_GTs + REF_GTs
HET_GTs = [(0, 1), (None, 1), (None, 3)]


def import_callrates(table_in):
    """
    Import table of variant callrates
    """

    callrates = {}

    with open(table_in) as tsvfile:
        reader = csv.reader(tsvfile, delimiter='\t')
        for vid, callrate in reader:
            if vid not in callrates.keys():
                callrates[vid] = float(callrate)

    return callrates


# def get_call_rate(record, samples):
#     """
#     Get fraction of samples with non-null genotypes
#     """
#     total_s = [s for s in record.samples if s in samples]
#     total = len(total_s)
#     nocall_s = [s for s in total_s if record.samples[s]['GT'] in NULL_GTs]
#     nocall = len(nocall_s)
#     callrate = 1 - ( nocall / total )
#     return callrate


def recal_qual_score(record):
    """
    Recalibrate quality score for a single variant
    """
    quals = []
    for s in [s for s in record.samples]:
        GT = record.samples[s]['GT']
        if GT in NULL_and_REF_GTs:
            continue
        elif GT in HET_GTs:
            quals.append(record.samples[s]['GQ'])
        else:
            quals.append(99)

    if len(quals) > 0:
        return int(median(quals))


def cleanup_vcf(vcf, fout, callrates, min_callrate_global=0.85,
                min_callrate_smallDels=0.95):

    # minus_samples = [s for s in vcf.header.samples if s not in plus_samples]
    # male_minus_samples = [s for s in minus_samples if s not in male_samples]

    for record in vcf:
        # Move several filters from FILTER to INFO
        for filt in filts_for_info:
            if filt in record.filter:
                record.info[filt] = True

        # Move HIGH_SR_BACKGROUND

        # Remove all HIGH_NOCALL_RATE and HIGH_SR_BACKGROUND tags from FILTER column
        newfilts = [
            filt for filt in record.filter if filt not in filts_to_remove]
        record.filter.clear()
        for filt in newfilts:
            record.filter.add(filt)
        if len(record.filter) == 0:
            record.filter.add('PASS')

        # #Mark sites with low PCR+ call rate
        # plus_callrate = get_call_rate(record, plus_samples)
        # if plus_callrate < min_callrate:
        #     if 'LOW_PCRPLUS_CALL_RATE' not in record.info.keys():
        #         record.info.keys().append('LOW_PCRPLUS_CALL_RATE')
        #     record.info['LOW_PCRPLUS_CALL_RATE'] = True

        # Mark sites with low PCR- call rate
        if record.id in callrates.keys():
            callrate = callrates[record.id]
            # Mark small (300bp-1kb) deletions with stricter 5% null gt rate,
            # and mark all other variants at specified null gt rate
            if record.info['SVTYPE'] == 'DEL' \
                    and record.info['SVLEN'] < 1000 \
                    and record.info['SVLEN'] > 300:
                if callrate < min_callrate_smallDels:
                    record.filter.add('LOW_CALL_RATE')
            else:
                if callrate < min_callrate_global:
                    record.filter.add('LOW_CALL_RATE')

        # Recalibrate QUAL score for biallelic variants
        if is_biallelic(record):
            newQUAL = recal_qual_score(record)
            if newQUAL is not None:
                record.qual = newQUAL

        # Only check for non-empty GTs for biallelic variants
        if is_biallelic(record):
            for s in record.samples:
                if record.samples[s]['GT'] not in NULL_and_REF_GTs:
                    fout.write(record)
                    break
        else:
            fout.write(record)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf', help='Input vcf (supports "stdin").')
    # parser.add_argument('PCRPLUS_samples', help='List of PCRPLUS sample IDs.')
    # parser.add_argument('male_samples', help='List of male sample IDs.')
    parser.add_argument('fout', help='Output file (supports "stdout").')
    parser.add_argument('--callrate-table', help='TSV of variant IDs and ' +
                        'their corresponding callrates.', required=True)
    parser.add_argument('--min-callrate-global', type=float, help='Minimum fraction ' +
                        'of samples required to have non-missing genotypes for ' +
                        'all variants.', default=0.85)
    parser.add_argument('--min-callrate-smallDels', type=float, help='Minimum fraction ' +
                        'of samples required to have non-missing genotypes for ' +
                        'DEL variants between 300bp-1kb.', default=0.95)

    args = parser.parse_args()

    # Open connection to input VCF
    if args.vcf in '- stdin'.split():
        vcf = pysam.VariantFile(sys.stdin)
    else:
        vcf = pysam.VariantFile(args.vcf)

    # Add new FILTER lines to VCF header
    NEW_FILTERS = ['##FILTER=<ID=LOW_CALL_RATE,Description="Site does not meet ' +
                   'minimum requirements for fraction of PCR- samples with non-null ' +
                   'genotypes. Flags sites more prone to false discoveries.">']
    header = vcf.header
    for filt in NEW_FILTERS:
        header.add_line(filt)

    # Remove unused FILTER lines from VCF header
    for filt in filts_to_remove:
        if filt in header.filters:
            header.filters.remove_header(filt)

    # Add new INFO lines to VCF header
    NEW_INFOS = ['##INFO=<ID=PESR_GT_OVERDISPERSION,Number=0,Type=Flag,Description=' +
                 '"PESR genotyping data is overdispersed. Flags sites where genotypes' +
                 ' are likely noisier.">',
                 '##INFO=<ID=HIGH_SR_BACKGROUND,Number=0,Type=Flag,Description=' +
                 '"Suspicious accumulation of split reads in predicted non-carrier ' +
                 'samples. Flags sites more prone to false discoveries and where ' +
                 'breakpoint precision is reduced.">',
                 '##INFO=<ID=BOTHSIDES_SUPPORT,Number=0,Type=Flag,Description=' +
                 '"Variant has read-level support for both sides of breakpoint. ' +
                 'Indicates higher-confidence variants.">',
                 '##INFO=<ID=VARIABLE_ACROSS_BATCHES,Number=0,Type=Flag,Description=' +
                 '"Site appears at variable frequencies across batches. Accuracy ' +
                 'of allele frequency estimates for these sites may be reduced.">']
    for info in NEW_INFOS:
        header.add_line(info)

    # #Read list of PCR+ samples
    # f_plus_samples = open(args.PCRPLUS_samples, 'r')
    # plus_samples = f_plus_samples.read().splitlines()
    # f_plus_samples.close()

    # #Read list of male samples
    # f_male_samples = open(args.male_samples, 'r')
    # male_samples = f_male_samples.read().splitlines()
    # f_male_samples.close()

    # Read callrates
    callrates = import_callrates(args.callrate_table)

    # Open connection to output VCF
    if args.fout in '- stdout'.split():
        fout = pysam.VariantFile(sys.stdout, 'w', header=vcf.header)
    else:
        fout = pysam.VariantFile(args.fout, 'w', header=vcf.header)

    # Cleanup VCF
    cleanup_vcf(vcf, fout, callrates,
                min_callrate_global=args.min_callrate_global,
                min_callrate_smallDels=args.min_callrate_smallDels,)

    fout.close()


if __name__ == '__main__':
    main()
