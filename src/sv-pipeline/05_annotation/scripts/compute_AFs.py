#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2018 Ryan Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

"""
Compute allele frequencies for sex & population combinations given an SV VCF
"""


import sys
import argparse
import pysam
import numpy as np
from svtk import utils as svu
from collections import Counter


def create_pop_dict(popfile):
    """
    Makes dictionary of sample-population pairs
    """

    pop_dict = {}

    for sample in popfile:
        pop_dict[sample.split('\t')[0]] = sample.split('\t')[1]

    return pop_dict


def _is_biallelic(record):
    """
    Check if record is biallelic
    """
    if 'MULTIALLELIC' not in record.filter \
    and len(record.alleles) <= 2 \
    and record.info['SVTYPE'] not in 'CNV MCNV'.split():
        return True
    else:
        return False


def gather_allele_freqs(record, all_samples, males, females, pop_dict, pops, no_combos = False):
    """
    Wrapper to compute allele frequencies for all sex & population pairings
    """

    #Get allele frequencies
    calc_allele_freq(record, all_samples)
    if len(males) > 0:
        calc_allele_freq(record, males, prefix = 'MALE')
    if len(females) > 0:
        calc_allele_freq(record, females, prefix = 'FEMALE')
    if len(pops) > 0:
        for pop in pops:
            pop_samps = [s for s in all_samples if pop_dict.get(s, None) == pop]
            calc_allele_freq(record, pop_samps, prefix = pop)
            if len(males) > 0 and not no_combos:
                calc_allele_freq(record, [s for s in pop_samps if s in males], 
                                 prefix = pop + '_MALE')
            if len(females) > 0 and not no_combos:
                calc_allele_freq(record, [s for s in pop_samps if s in females], 
                                 prefix = pop + '_FEMALE')

        #Get POPMAX AF biallelic sites only
        if _is_biallelic(record):
            AFs = [record.info['{0}_AF'.format(pop)][0] for pop in pops]
            popmax = max(AFs)
            record.info['POPMAX_AF'] = popmax

    return record


def calc_allele_freq(record, samples, prefix = None):
    """
    Computes allele frequencies for a single record based on a list of samples
    """

    #Treat biallelic and multiallelic records differently
    #For biallelic sites, count number of non-ref, non-no-call GTs
    if _is_biallelic(record):

        # Get all sample GTs
        GTs = [s['GT'] for s in record.samples.values() if s.name in samples]
        
        #Count alleles & genotypes
        AC = 0
        AN = 0
        nhomref = 0
        nhet = 0
        nhomalt = 0
        for GT in GTs:
            AN += len([allele for allele in GT if allele != '.' and allele is not None])
            AC += len([allele for allele in GT if allele != '.' and allele != 0 and allele is not None])
            if GT == (0, 0):
                nhomref += 1
            if len([allele for allele in GT if allele == 0 and allele != '.' and allele is not None]) == 1 \
            and len([allele for allele in GT if allele != 0 and allele != '.' and allele is not None]) == 1:
                nhet += 1
            if len([allele for allele in GT if allele != 0 and allele != '.' and allele is not None]) == 2:
                nhomalt += 1

        #Calculate allele frequency
        if AN > 0:
            AF = AC / AN
            AF = round(AF, 6)
        else:
            AF = 0

        #Add AN, AC, and AF to INFO field
        record.info[(prefix + '_' if prefix else '') + 'AN'] = AN
        record.info[(prefix + '_' if prefix else '') + 'AC'] = AC
        record.info[(prefix + '_' if prefix else '') + 'AF'] = AF

        #Calculate genotype frequencies
        n_bi_genos = nhomref + nhet + nhomalt
        if n_bi_genos > 0:
            freq_homref = nhomref / n_bi_genos
            freq_het = nhet / n_bi_genos
            freq_homalt = nhomalt / n_bi_genos
        else:
            freq_homref = 0
            freq_het = 0
            freq_homalt = 0

        #Add N_BI_GENOS, N_HOMREF, N_HET, N_HOMALT, FREQ_HOMREF, FREQ_HET, and FREQ_HOMALT to INFO field
        record.info[(prefix + '_' if prefix else '') + 'N_BI_GENOS'] = n_bi_genos
        record.info[(prefix + '_' if prefix else '') + 'N_HOMREF'] = nhomref
        record.info[(prefix + '_' if prefix else '') + 'N_HET'] = nhet
        record.info[(prefix + '_' if prefix else '') + 'N_HOMALT'] = nhomalt
        record.info[(prefix + '_' if prefix else '') + 'FREQ_HOMREF'] = freq_homref
        record.info[(prefix + '_' if prefix else '') + 'FREQ_HET'] = freq_het
        record.info[(prefix + '_' if prefix else '') + 'FREQ_HOMALT'] = freq_homalt

    #Multiallelic sites should reference FORMAT:CN rather than GT
    #Compute CN_NUMBER, CN_NONREF_COUNT, CN_NONREF_FREQ, and CN_COUNT/CN_FREQ for each copy state
    else:

        # Get all sample CNs
        CNs = [s['CN'] for s in record.samples.values() if s.name in samples]

        # Count number of samples per CN and total CNs observed
        CN_counts = dict(Counter(CNs))
        nonnull_CNs = sum([int(v) for k, v in CN_counts.items() if k not in '. NA'.split()])

        # Get max observed CN and enumerate counts/frequencies per CN as list starting from CN=0
        max_CN = max([int(k) for k, v in CN_counts.items()])
        CN_dist = [int(CN_counts.get(k, 0)) for k in range(max_CN+1)]
        CN_freqs = [round(v / nonnull_CNs, 6) for v in CN_dist]

        # Get total non-reference CN counts and freq
        # Note: assumes reference is diploid (this will not be true for sex chromosomes,
        # but can be clarified when passing --FAMFILE and annotating per sex, 
        # and is a more internally consistent solution than passing a BED file with 
        # PAR regions etc. specifying assumed ploidy genome-wide)
        nonref_CN_count = sum([int(CN_counts.get(k, 0)) for k in range(max_CN+1) if k != 2])
        nonref_CN_freq = round(nonref_CN_count / nonnull_CNs, 6)

        #Add values to INFO field
        record.info[(prefix + '_' if prefix else '') + 'CN_NUMBER'] = nonnull_CNs
        record.info[(prefix + '_' if prefix else '') + 'CN_COUNT'] = tuple(CN_dist)
        record.info[(prefix + '_' if prefix else '') + 'CN_FREQ'] = tuple(CN_freqs)
        record.info[(prefix + '_' if prefix else '') + 'CN_NONREF_COUNT'] = nonref_CN_count
        record.info[(prefix + '_' if prefix else '') + 'CN_NONREF_FREQ'] = nonref_CN_freq

    return record


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf', help='Input vcf. Also accepts "stdin" and "-".')
    parser.add_argument('-p','--popfile', help='Two-column file of samples & ' +
                        'their population assignments. A "." denotes no assignment.',
                        default = None)
    parser.add_argument('-f','--famfile', help='Input .fam file (used for sex-specific AFs).',
                        default = None)
    parser.add_argument('--no-combos', help='Do not compute combinations of populations ' +
                        'and sexes.', action = 'store_true', default = False)
    parser.add_argument('fout', help='Output vcf. Also accepts "stdout" and "-".')
    args = parser.parse_args()

    #Open connections to input VCF
    if args.vcf in '- stdin'.split():
        vcf = pysam.VariantFile(sys.stdin) 
    else:
        vcf = pysam.VariantFile(args.vcf)

    #Get list of all samples in vcf
    all_samples = list(vcf.header.samples)

    #Get lists of males and femailes
    if args.famfile is not None:
        famfile = [line.rstrip('\n') for line in open(args.famfile)]
        males = [line.split('\t')[1] for line in famfile if line.split('\t')[4] == '1']
        females = [line.split('\t')[1] for line in famfile if line.split('\t')[4] == '2']
        sexes = 'MALE FEMALE'.split()
    else:
        males = []
        females = []
        sexes = []


    #Get dictionary of populations
    if args.popfile is not None:
        popfile = [line.rstrip('\n') for line in open(args.popfile)]
        pop_dict = create_pop_dict(popfile)
        pops = list(set(pop_dict.values()))
        pops = sorted([p for p in pops if p != "."])
    else:
        pop_dict = {}
        pops = []
    
    #Add relevant fields to header
    INFO_ADD = [
    '##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles genotyped (biallelic sites only).">',
    '##INFO=<ID=AC,Number=A,Type=Integer,Description="Number of non-reference alleles observed (biallelic sites only).">',
    '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency (biallelic sites only).">',
    '##INFO=<ID=N_BI_GENOS,Number=1,Type=Integer,Description="Total number of individuals with complete genotypes (biallelic sites only).">',
    '##INFO=<ID=N_HOMREF,Number=1,Type=Integer,Description="Number of individuals with homozygous reference genotypes (biallelic sites only).">',
    '##INFO=<ID=N_HET,Number=1,Type=Integer,Description="Number of individuals with heterozygous genotypes (biallelic sites only).">',
    '##INFO=<ID=N_HOMALT,Number=1,Type=Integer,Description="Number of individuals with homozygous alternate genotypes (biallelic sites only).">',
    '##INFO=<ID=FREQ_HOMREF,Number=1,Type=Float,Description="Homozygous reference genotype frequency (biallelic sites only).">',
    '##INFO=<ID=FREQ_HET,Number=1,Type=Float,Description="Heterozygous genotype frequency (biallelic sites only).">',
    '##INFO=<ID=FREQ_HOMALT,Number=1,Type=Float,Description="Homozygous alternate genotype frequency (biallelic sites only).">',
    '##INFO=<ID=CN_NUMBER,Number=1,Type=Integer,Description="Total number of individuals with estimated copy numbers (multiallelic CNVs only).">',
    '##INFO=<ID=CN_COUNT,Number=.,Type=Integer,Description="Number of individuals observed at each copy state, starting from CN=0 (multiallelic CNVs only).">',
    '##INFO=<ID=CN_FREQ,Number=.,Type=Float,Description="Frequency of individuals observed at each copy state, starting from CN=0 (multiallelic CNVs only).">',
    '##INFO=<ID=CN_NONREF_COUNT,Number=1,Type=Integer,Description="Sum of all individuals with non-reference copy states (multiallelic CNVs only).">',
    '##INFO=<ID=CN_NONREF_FREQ,Number=1,Type=Float,Description="Total frequency of all individuals across all non-reference copy states (multiallelic CNVs only).">'
    ]
    if len(sexes) > 0:
        for sex in sexes:
            INFO_ADD.append('##INFO=<ID=%s_AN,Number=1,Type=Integer,Description="Total number of %s alleles genotyped (biallelic sites only).">' % (sex, sex))
            INFO_ADD.append('##INFO=<ID=%s_AC,Number=A,Type=Integer,Description="Number of non-reference %s alleles observed (biallelic sites only).">' % (sex, sex))
            INFO_ADD.append('##INFO=<ID=%s_AF,Number=A,Type=Float,Description="%s allele frequency (biallelic sites only).">' % (sex, sex))
            INFO_ADD.append('##INFO=<ID=%s_N_BI_GENOS,Number=1,Type=Integer,Description="Total number of %s individuals with complete genotypes (biallelic sites only).">' % (sex, sex))
            INFO_ADD.append('##INFO=<ID=%s_N_HOMREF,Number=1,Type=Integer,Description="Number of %s individuals with homozygous reference genotypes (biallelic sites only).">' % (sex, sex))
            INFO_ADD.append('##INFO=<ID=%s_N_HET,Number=1,Type=Integer,Description="Number of %s individuals with heterozygous genotypes (biallelic sites only).">' % (sex, sex))
            INFO_ADD.append('##INFO=<ID=%s_N_HOMALT,Number=1,Type=Integer,Description="Number of %s individuals with homozygous alternate genotypes (biallelic sites only).">' % (sex, sex))
            INFO_ADD.append('##INFO=<ID=%s_FREQ_HOMREF,Number=1,Type=Float,Description="%s homozygous reference genotype frequency (biallelic sites only).">' % (sex, sex))
            INFO_ADD.append('##INFO=<ID=%s_FREQ_HET,Number=1,Type=Float,Description="%s heterozygous genotype frequency (biallelic sites only).">' % (sex, sex))
            INFO_ADD.append('##INFO=<ID=%s_FREQ_HOMALT,Number=1,Type=Float,Description="%s homozygous alternate genotype frequency (biallelic sites only).">' % (sex, sex))
            INFO_ADD.append('##INFO=<ID=%s_CN_NUMBER,Number=1,Type=Integer,Description="Total number of %s individuals with estimated copy numbers (multiallelic CNVs only).">' % (sex, sex))
            INFO_ADD.append('##INFO=<ID=%s_CN_COUNT,Number=.,Type=Integer,Description="Number of %s individuals observed at each copy state, starting from CN=0 (multiallelic CNVs only).">' % (sex, sex))
            INFO_ADD.append('##INFO=<ID=%s_CN_FREQ,Number=.,Type=Float,Description="Frequency of %s individuals observed at each copy state, starting from CN=0 (multiallelic CNVs only).">' % (sex, sex))
            INFO_ADD.append('##INFO=<ID=%s_CN_NONREF_COUNT,Number=1,Type=Integer,Description="Sum of all %s individuals with non-reference copy states (multiallelic CNVs only).">' % (sex, sex))
            INFO_ADD.append('##INFO=<ID=%s_CN_NONREF_FREQ,Number=1,Type=Float,Description="Total frequency of all %s individuals across all non-reference copy states (multiallelic CNVs only).">' % (sex, sex))
    if len(pops) > 0:
        INFO_ADD.append('##INFO=<ID=POPMAX_AF,Number=1,Type=Float,Description="Maximum allele frequency across any population (biallelic sites only).">')
        for pop in pops:
            INFO_ADD.append('##INFO=<ID=%s_AN,Number=1,Type=Integer,Description="Total number of %s alleles genotyped (biallelic sites only).">' % (pop, pop))
            INFO_ADD.append('##INFO=<ID=%s_AC,Number=A,Type=Integer,Description="Number of non-reference %s alleles observed (biallelic sites only).">' % (pop, pop))
            INFO_ADD.append('##INFO=<ID=%s_AF,Number=A,Type=Float,Description="%s allele frequency (biallelic sites only).">' % (pop, pop))
            INFO_ADD.append('##INFO=<ID=%s_N_BI_GENOS,Number=1,Type=Integer,Description="Total number of %s individuals with complete genotypes (biallelic sites only).">' % (pop, pop))
            INFO_ADD.append('##INFO=<ID=%s_N_HOMREF,Number=1,Type=Integer,Description="Number of %s individuals with homozygous reference genotypes (biallelic sites only).">' % (pop, pop))
            INFO_ADD.append('##INFO=<ID=%s_N_HET,Number=1,Type=Integer,Description="Number of %s individuals with heterozygous genotypes (biallelic sites only).">' % (pop, pop))
            INFO_ADD.append('##INFO=<ID=%s_N_HOMALT,Number=1,Type=Integer,Description="Number of %s individuals with homozygous alternate genotypes (biallelic sites only).">' % (pop, pop))
            INFO_ADD.append('##INFO=<ID=%s_FREQ_HOMREF,Number=1,Type=Float,Description="%s homozygous reference genotype frequency (biallelic sites only).">' % (pop, pop))
            INFO_ADD.append('##INFO=<ID=%s_FREQ_HET,Number=1,Type=Float,Description="%s heterozygous genotype frequency (biallelic sites only).">' % (pop, pop))
            INFO_ADD.append('##INFO=<ID=%s_FREQ_HOMALT,Number=1,Type=Float,Description="%s homozygous alternate genotype frequency (biallelic sites only).">' % (pop, pop))
            INFO_ADD.append('##INFO=<ID=%s_CN_NUMBER,Number=1,Type=Integer,Description="Total number of %s individuals with estimated copy numbers (multiallelic CNVs only).">' % (pop, pop))
            INFO_ADD.append('##INFO=<ID=%s_CN_COUNT,Number=.,Type=Integer,Description="Number of %s individuals observed at each copy state, starting from CN=0 (multiallelic CNVs only).">' % (pop, pop))
            INFO_ADD.append('##INFO=<ID=%s_CN_FREQ,Number=.,Type=Float,Description="Frequency of %s individuals observed at each copy state, starting from CN=0 (multiallelic CNVs only).">' % (pop, pop))
            INFO_ADD.append('##INFO=<ID=%s_CN_NONREF_COUNT,Number=1,Type=Integer,Description="Sum of all %s individuals with non-reference copy states (multiallelic CNVs only).">' % (pop, pop))
            INFO_ADD.append('##INFO=<ID=%s_CN_NONREF_FREQ,Number=1,Type=Float,Description="Total frequency of all %s individuals across all non-reference copy states (multiallelic CNVs only).">' % (pop, pop))
            if len(sexes) > 0 and not args.no_combos:
                for sex in sexes:
                    INFO_ADD.append('##INFO=<ID=%s_AN,Number=1,Type=Integer,Description="Total number of %s alleles genotyped (biallelic sites only).">' % ('_'.join((pop, sex)), ' '.join((pop, sex))))
                    INFO_ADD.append('##INFO=<ID=%s_AC,Number=A,Type=Integer,Description="Number of non-reference %s alleles observed (biallelic sites only).">' % ('_'.join((pop, sex)), ' '.join((pop, sex))))
                    INFO_ADD.append('##INFO=<ID=%s_AF,Number=A,Type=Float,Description="%s allele frequency (biallelic sites only).">' % ('_'.join((pop, sex)), ' '.join((pop, sex))))
                    INFO_ADD.append('##INFO=<ID=%s_N_BI_GENOS,Number=1,Type=Integer,Description="Total number of %s individuals with complete genotypes (biallelic sites only).">' % ('_'.join((pop, sex)), ' '.join((pop, sex))))
                    INFO_ADD.append('##INFO=<ID=%s_N_HOMREF,Number=1,Type=Integer,Description="Number of %s individuals with homozygous reference genotypes (biallelic sites only).">' % ('_'.join((pop, sex)), ' '.join((pop, sex))))
                    INFO_ADD.append('##INFO=<ID=%s_N_HET,Number=1,Type=Integer,Description="Number of %s individuals with heterozygous genotypes (biallelic sites only).">' % ('_'.join((pop, sex)), ' '.join((pop, sex))))
                    INFO_ADD.append('##INFO=<ID=%s_N_HOMALT,Number=1,Type=Integer,Description="Number of %s individuals with homozygous alternate genotypes (biallelic sites only).">' % ('_'.join((pop, sex)), ' '.join((pop, sex))))
                    INFO_ADD.append('##INFO=<ID=%s_FREQ_HOMREF,Number=1,Type=Float,Description="%s homozygous reference genotype frequency (biallelic sites only).">' % ('_'.join((pop, sex)), ' '.join((pop, sex))))
                    INFO_ADD.append('##INFO=<ID=%s_FREQ_HET,Number=1,Type=Float,Description="%s heterozygous genotype frequency (biallelic sites only).">' % ('_'.join((pop, sex)), ' '.join((pop, sex))))
                    INFO_ADD.append('##INFO=<ID=%s_FREQ_HOMALT,Number=1,Type=Float,Description="%s homozygous alternate genotype frequency (biallelic sites only).">' % ('_'.join((pop, sex)), ' '.join((pop, sex))))
                    INFO_ADD.append('##INFO=<ID=%s_CN_NUMBER,Number=1,Type=Integer,Description="Total number of %s individuals with estimated copy numbers (multiallelic CNVs only).">' % ('_'.join((pop, sex)), ' '.join((pop, sex))))
                    INFO_ADD.append('##INFO=<ID=%s_CN_COUNT,Number=.,Type=Integer,Description="Number of %s individuals observed at each copy state, starting from CN=0 (multiallelic CNVs only).">' % ('_'.join((pop, sex)), ' '.join((pop, sex))))
                    INFO_ADD.append('##INFO=<ID=%s_CN_FREQ,Number=.,Type=Float,Description="Frequency of %s individuals observed at each copy state, starting from CN=0 (multiallelic CNVs only).">' % ('_'.join((pop, sex)), ' '.join((pop, sex))))
                    INFO_ADD.append('##INFO=<ID=%s_CN_NONREF_COUNT,Number=1,Type=Integer,Description="Sum of all %s individuals with non-reference copy states (multiallelic CNVs only).">' % ('_'.join((pop, sex)), ' '.join((pop, sex))))
                    INFO_ADD.append('##INFO=<ID=%s_CN_NONREF_FREQ,Number=1,Type=Float,Description="Total frequency of all %s individuals across all non-reference copy states (multiallelic CNVs only).">' % ('_'.join((pop, sex)), ' '.join((pop, sex))))
    for line in INFO_ADD:
        vcf.header.add_line(line)
    
    #Prep output VCF
    if args.fout in '- stdout'.split():
        fout = pysam.VariantFile(sys.stdout, 'w', header=vcf.header)
    else:
        fout = pysam.VariantFile(args.fout, 'w', header=vcf.header)

    #Get allele frequencies for each record & write to new VCF
    for r in vcf.fetch():
        newrec = gather_allele_freqs(r, all_samples, males, females, pop_dict, pops, args.no_combos)
        fout.write(newrec)

    fout.close()

if __name__ == '__main__':
    main()
