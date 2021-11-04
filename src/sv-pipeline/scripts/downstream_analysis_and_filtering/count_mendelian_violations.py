#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""
Tabulate Mendelian violation rate (MVR) per SV & append other info
"""


import argparse
import sys
import csv
import pysam
import numpy as np


# Convert famfile into list of dicts for complete trios present in vcf
def parse_famfile(famfile, vcf):
    # Get list of samples present in VCF
    samples_in_vcf = [s for s in vcf.header.samples]

    # Iterate over famfile and check for presence of all three members
    # If all three members exist in VCF, add them as a dict to the trios list
    trios = []
    with open(famfile) as f:
        reader = csv.reader(f, delimiter='\t')

        for famID, proband, father, mother, sex, pheno in reader:
            # Skip header line, if any
            if "#" in famID:
                continue

            # Check for presence of all three family members in VCF
            # If all are present, add to list of trios
            famIDs = [proband, father, mother]
            matches = len([s for s in famIDs if s in samples_in_vcf])
            if matches == 3:
                newfam = {"proband": proband,
                          "father": father,
                          "mother": mother}
                trios.append(newfam)

    return trios


# Build lookup dictionary of samples by PCR+ and PCR- status
def get_pcrdict(vcf, PCRPLUS_samples_file):
    # Get list of all samples in VCF
    samples = [s for s in vcf.header.samples]

    # Read list of all PCRPLUS samples
    PCRPLUS_samples = [line.rstrip('\n') for line in PCRPLUS_samples_file]

    # Make dictionary of samples with PCR statuses
    pcrdict = {}
    for s in samples:
        if s not in pcrdict.keys():
            if s in PCRPLUS_samples:
                pcrdict[s] = 'PCRPLUS'
            else:
                pcrdict[s] = 'PCRMINUS'

    return pcrdict


# Get allele count for a given sample
def get_AC(record, ID):
    GTs_to_skip = './. None/None 0/None None/0'.split()
    GT = '/'.join([str(i) for i in record.samples[ID]['GT']])
    if GT not in GTs_to_skip:
        AC = sum([int(a) for a in GT.split('/')])
        return AC


# Classify single trio allele count
def classify_trio_AC(ACs):
    n_missing = len([c for c in ACs if c is None])
    n_homref = len([c for c in ACs if c == 0])
    pro_homref = ACs[0] == 0
    n_homref_parents = len([c for c in ACs[1:3] if c == 0])
    # n_het = len([c for c in ACs if c == 1])
    # pro_het = ACs[0] == 1
    # n_het_parents = len([c for c in ACs[1:3] if c == 1])
    # n_homalt = len([c for c in ACs if c == 2])
    pro_homalt = ACs[0] == 2
    n_homalt_parents = len([c for c in ACs[1:3] if c == 2])

    if n_missing > 0:
        return 'INCOMPLETE'
    elif n_homref == 3:
        return 'NO_VARIANT'
    elif pro_homalt and n_homref_parents == 2:
        return 'APPARENT_DE_NOVO'
    elif pro_homref and n_homalt_parents > 0:
        return 'UNTRANSMITTED_HOMOZYGOTE'
    elif pro_homalt and n_homref_parents > 0:
        return 'SPONTANEOUS_HOMOZYGOTE'
    else:
        return 'MENDELIAN'


# Count Mendelian violation categories for all trios
def count_MVs(record, trios, pcrdict):
    mendel = {}
    cats = 'INCOMPLETE NO_VARIANT APPARENT_DE_NOVO UNTRANSMITTED_HOMOZYGOTE ' + \
           'SPONTANEOUS_HOMOZYGOTE MENDELIAN TRIOS_CONSIDERED'
    for PCR in 'PCRPLUS PCRMINUS'.split(' '):
        for k in cats.split(' '):
            newkey = '{0}_{1}'.format(PCR, k)
            if newkey not in mendel.keys():
                mendel[newkey] = 0

    for fam in trios:
        propcr = pcrdict.get(fam['proband'], 'PCRMINUS')
        ACs = [get_AC(record, s) for s in fam.values()]
        classification = classify_trio_AC(ACs)
        mendel['{0}_{1}'.format(propcr, classification)] += 1
        skips = '{0}_INCOMPLETE {0}_NO_VARIANT {1}_INCOMPLETE {1}_NO_VARIANT'.format(
            'PCRPLUS', 'PCRMINUS')
        if classification not in skips.split(' '):
            mendel['{0}_TRIOS_CONSIDERED'.format(propcr)] += 1

    return mendel


# Iterate over a vcf and gather MVR data & other info
def gather_info(vcf, trios, pcrdict, fout, no_header=False):
    sex_chroms = 'X Y chrX chrY'.split()

    # Write header to output file
    if not no_header:
        header = '#VID\tSVLEN\tSVTYPE\tFILTER\tQUAL\t' + \
                 '{0}_NULL_GTs\t{0}_REF_GTs\t{0}_NONREF_GTs\t' + \
                 '{1}_NULL_GTs\t{1}_REF_GTs\t{1}_NONREF_GTs\t' + \
                 '{0}_REF_MIN_GQ\t{0}_REF_Q1_GQ\t{0}_REF_MEDIAN_GQ\t' + \
                 '{0}_REF_Q3_GQ\t{0}_REF_MAX_GQ\t' + \
                 '{0}_NONREF_MIN_GQ\t{0}_NONREF_Q1_GQ\t{0}_NONREF_MEDIAN_GQ\t' + \
                 '{0}_NONREF_Q3_GQ\t{0}_NONREF_MAX_GQ\t' + \
                 '{1}_REF_MIN_GQ\t{1}_REF_Q1_GQ\t{1}_REF_MEDIAN_GQ\t' + \
                 '{1}_REF_Q3_GQ\t{1}_REF_MAX_GQ\t' + \
                 '{1}_NONREF_MIN_GQ\t{1}_NONREF_Q1_GQ\t{1}_NONREF_MEDIAN_GQ\t' + \
                 '{1}_NONREF_Q3_GQ\t{1}_NONREF_MAX_GQ\t' + \
                 '{0}_NONREF_RD\t{0}_NONREF_PE\t{0}_NONREF_SR\t' + \
                 '{0}_NONREF_RDPE\t{0}_NONREF_RDSR\t{0}_NONREF_PESR\t' + \
                 '{0}_NONREF_RDPESR\t' + \
                 '{1}_NONREF_RD\t{1}_NONREF_PE\t{1}_NONREF_SR\t' + \
                 '{1}_NONREF_RDPE\t{1}_NONREF_RDSR\t{1}_NONREF_PESR\t' + \
                 '{1}_NONREF_RDPESR\t' + \
                 '{0}_TRIOS_CONSIDERED\t{0}_MENDELIAN\t{0}_APPARENT_DE_NOVO\t' + \
                 '{0}_UNTRANSMITTED_PARENT_HOMOZYGOTE\t{0}_SPONTANEOUS_CHILD_HOMOZYGOTE\t' + \
                 '{0}_NO_VARIANT_TRIOS\t{0}_INCOMPLETE_GENOTYPE_TRIOS\t' + \
                 '{1}_TRIOS_CONSIDERED\t{1}_MENDELIAN\t{1}_APPARENT_DE_NOVO\t' + \
                 '{1}_UNTRANSMITTED_PARENT_HOMOZYGOTE\t{1}_SPONTANEOUS_CHILD_HOMOZYGOTE\t' + \
                 '{1}_NO_VARIANT_TRIOS\t{1}_INCOMPLETE_GENOTYPE_TRIOS\n'
        header = header.format('PCRPLUS', 'PCRMINUS')
        fout.write(header)

    # Iterate over VCF
    for record in vcf:
        # Do not include variants from sex chromosomes
        if record.chrom in sex_chroms:
            continue

        # Do not include multiallelic variants
        if 'MULTIALLELIC' in record.info.keys() \
                or 'MULTIALLELIC' in record.filter \
                or len(record.alts) > 1:
            continue

        # Iterate over samples and get count of genotypes by PCR status,
        # lists of GQs based on genotype, and counts of non-ref evidence
        GT_counts = {'PCRPLUS_NULL': 0,
                     'PCRPLUS_REF': 0,
                     'PCRPLUS_NONREF': 0,
                     'PCRMINUS_NULL': 0,
                     'PCRMINUS_REF': 0,
                     'PCRMINUS_NONREF': 0}
        GQ_dict = {'PCRPLUS_REF_GQs': [],
                   'PCRPLUS_NONREF_GQs': [],
                   'PCRMINUS_REF_GQs': [],
                   'PCRMINUS_NONREF_GQs': []}
        EV_counts = {'PCRPLUS_RD': 0, 'PCRPLUS_PE': 0, 'PCRPLUS_SR': 0,
                     'PCRPLUS_RDPE': 0, 'PCRPLUS_RDSR': 0, 'PCRPLUS_PESR': 0,
                     'PCRPLUS_RDPESR': 0,
                     'PCRMINUS_RD': 0, 'PCRMINUS_PE': 0, 'PCRMINUS_SR': 0,
                     'PCRMINUS_RDPE': 0, 'PCRMINUS_RDSR': 0, 'PCRMINUS_PESR': 0,
                     'PCRMINUS_RDPESR': 0}
        for s in vcf.header.samples:
            sgt = get_AC(record, s)
            spcr = pcrdict.get(s, 'PCRMINUS')
            if sgt is None:
                GT_counts['{0}_NULL'.format(spcr)] += 1
            elif sgt == 0:
                GT_counts['{0}_REF'.format(spcr)] += 1
                GQ_dict['{0}_REF_GQs'.format(spcr)].append(
                    record.samples[s]['GQ'])
            else:
                GT_counts['{0}_NONREF'.format(spcr)] += 1
                GQ_dict['{0}_NONREF_GQs'.format(spcr)].append(
                    record.samples[s]['GQ'])
                sev = record.samples[s]['EV']
                if isinstance(sev, tuple):
                    sev = ''.join(list(sev))
                EV_counts['{0}_{1}'.format(spcr, sev)] += 1

        # Summarize GQ stats
        GQ_stats = {}
        for PCR in 'PCRPLUS PCRMINUS'.split(' '):
            for GT in 'REF NONREF'.split(' '):
                GQs = GQ_dict['{0}_{1}_GQs'.format(PCR, GT)]
                if len(GQs) > 0:
                    gmin = min(GQs)
                    g1q = int(np.percentile(GQs, 0.25))
                    gmed = int(np.percentile(GQs, 0.50))
                    g3q = int(np.percentile(GQs, 0.75))
                    gmax = max(GQs)
                    gstats = '\t'.join([str(i)
                                        for i in [gmin, g1q, gmed, g3q, gmax]])
                else:
                    gstats = '\t'.join(['NA'] * 5)
                GQ_stats['{0}_{1}'.format(PCR, GT)] = gstats

        # Compute MVR stats across all trios & format into string for output
        mendel = count_MVs(record, trios, pcrdict)
        mendel_info = ''
        info_order = 'TRIOS_CONSIDERED MENDELIAN APPARENT_DE_NOVO ' + \
                     'UNTRANSMITTED_HOMOZYGOTE SPONTANEOUS_HOMOZYGOTE ' + \
                     'NO_VARIANT INCOMPLETE'
        for PCR in 'PCRPLUS PCRMINUS'.split(' '):
            for info in info_order.split(' '):
                ikey = '{0}_{1}'.format(PCR, info)
                mendel_info = mendel_info + str(mendel.get(ikey, 'NA')) + '\t'
        mendel_info = mendel_info.rstrip('\t')

        # Get minimal variant info
        vid = record.id
        size = str(record.info['SVLEN'])
        svtype = record.info['SVTYPE']
        filt = ','.join([f for f in record.filter])
        VQ = record.qual
        GT_count_string = '\t'.join(str(i) for i in GT_counts.values())
        EV_count_string = '\t'.join(str(i) for i in EV_counts.values())
        vinfo = '\t'.join([str(i) for i in [vid, size, svtype, filt, VQ,
                                            GT_count_string,
                                            GQ_stats['PCRPLUS_REF'],
                                            GQ_stats['PCRPLUS_NONREF'],
                                            GQ_stats['PCRMINUS_REF'],
                                            GQ_stats['PCRMINUS_NONREF'],
                                            EV_count_string]])

        # Write record to file
        newline = '{0}\t{1}'.format(vinfo, mendel_info)
        fout.write(newline + '\n')


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf', help='Input vcf (supports "stdin").')
    parser.add_argument('famfile', help='Input famfile.')
    parser.add_argument('PCRPLUS_samples', help='List of PCRPLUS sample IDs.')
    parser.add_argument('fout', help='Output file (supports "stdout").')
    parser.add_argument('--no-header', help='Do not write header line.',
                        action='store_true', default=False)

    args = parser.parse_args()

    # Open input files
    if args.vcf in '- stdin'.split():
        vcf = pysam.VariantFile(sys.stdin)
    else:
        vcf = pysam.VariantFile(args.vcf)

    trios = parse_famfile(args.famfile, vcf)

    PCRPLUS_samples_file = open(args.PCRPLUS_samples, 'r')
    pcrdict = get_pcrdict(vcf, PCRPLUS_samples_file)

    # Open connection to output vcf
    if args.fout in '- stdout'.split():
        fout = sys.stdout
    else:
        fout = open(args.fout, 'w')

    # Gather info on all records in VCF
    gather_info(vcf, trios, pcrdict, fout, args.no_header)

    fout.close()


if __name__ == '__main__':
    main()
