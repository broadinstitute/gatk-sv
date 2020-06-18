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


def create_pop_dict(popfile):
    """
    Makes dictionary of sample-population pairs
    """

    pop_dict = {}

    for sample in popfile:
        pop_dict[sample.split('\t')[0]] = sample.split('\t')[1]

    return pop_dict


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

        #Get POPMAX AF for biallelic sites
        if 'MULTIALLELIC' not in record.filter and len(record.alleles) <= 2:
            AFs = [record.info['{0}_AF'.format(pop)][0] for pop in pops]
            popmax = max(AFs)
            record.info['POPMAX_AF'] = popmax

    return record


def calc_allele_freq(record, samples, prefix = None):
    """
    Computes allele frequencies for a single record based on a list of samples
    """

    GTs = [s['GT'] for s in record.samples.values() if s.name in samples]

    #Treat biallelic and multiallelic records differently
    #For biallelic sites, count number of non-ref, non-no-call GTs
    if 'MULTIALLELIC' not in record.filter and len(record.alleles) <= 2:
        
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

    #For multiallelic sites, count AN, AC, and AF for each non-ref copy state
    else:

        #Get list of all allele indexes to evaluate
        alleles_to_eval = list(np.arange(0, len(record.alleles)-1))

        #Count number of observed alleles per copy state
        AC_list = []
        for cn in alleles_to_eval:
            
            AC = 0
            for GT in GTs:
                AC += len([allele for allele in GT if allele == cn])

            AC_list.append(AC)

        #Calculate AN & AF
        AN = sum(AC_list)
        if AN > 0:
            AF_list = [round(c / AN, 6) for c in AC_list]
        else:
            AF_list = [0] * len(AC_list)

        #Add AN, AC, and AF to INFO field
        record.info[(prefix + '_' if prefix else '') + 'AN'] = AN
        record.info[(prefix + '_' if prefix else '') + 'AC'] = tuple(AC_list)
        record.info[(prefix + '_' if prefix else '') + 'AF'] = tuple(AF_list)        

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
    '##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles genotyped (for biallelic sites) or individuals with copy-state estimates (for multiallelic sites).">',
    '##INFO=<ID=AC,Number=A,Type=Integer,Description="Number of non-reference alleles observed (for biallelic sites) or individuals at each copy state (for multiallelic sites).">',
    '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency (for biallelic sites) or copy-state frequency (for multiallelic sites).">',
    '##INFO=<ID=N_BI_GENOS,Number=1,Type=Integer,Description="Total number of individuals with complete genotypes (biallelic sites only).">',
    '##INFO=<ID=N_HOMREF,Number=1,Type=Integer,Description="Number of individuals with homozygous reference genotypes (biallelic sites only).">',
    '##INFO=<ID=N_HET,Number=1,Type=Integer,Description="Number of individuals with heterozygous genotypes (biallelic sites only).">',
    '##INFO=<ID=N_HOMALT,Number=1,Type=Integer,Description="Number of individuals with homozygous alternate genotypes (biallelic sites only).">',
    '##INFO=<ID=FREQ_HOMREF,Number=1,Type=Float,Description="Homozygous reference genotype frequency (biallelic sites only).">',
    '##INFO=<ID=FREQ_HET,Number=1,Type=Float,Description="Heterozygous genotype frequency (biallelic sites only).">',
    '##INFO=<ID=FREQ_HOMALT,Number=1,Type=Float,Description="Homozygous alternate genotype frequency (biallelic sites only).">'
    ]
    if len(sexes) > 0:
        for sex in sexes:
            INFO_ADD.append('##INFO=<ID=%s_AN,Number=1,Type=Integer,Description="Total number of %s alleles genotyped (for biallelic sites) or %s individuals with copy-state estimates (for multiallelic sites).">' % (sex, sex, sex))
            INFO_ADD.append('##INFO=<ID=%s_AC,Number=A,Type=Integer,Description="Number of non-reference %s alleles observed (for biallelic sites) or %s individuals at each copy state (for multiallelic sites).">' % (sex, sex, sex))
            INFO_ADD.append('##INFO=<ID=%s_AF,Number=A,Type=Float,Description="%s allele frequency (for biallelic sites) or %s copy-state frequency (for multiallelic sites).">' % (sex, sex, sex))
            INFO_ADD.append('##INFO=<ID=%s_N_BI_GENOS,Number=1,Type=Integer,Description="Total number of %s individuals with complete genotypes (biallelic sites only).">' % (sex, sex))
            INFO_ADD.append('##INFO=<ID=%s_N_HOMREF,Number=1,Type=Integer,Description="Number of %s individuals with homozygous reference genotypes (biallelic sites only).">' % (sex, sex))
            INFO_ADD.append('##INFO=<ID=%s_N_HET,Number=1,Type=Integer,Description="Number of %s individuals with heterozygous genotypes (biallelic sites only).">' % (sex, sex))
            INFO_ADD.append('##INFO=<ID=%s_N_HOMALT,Number=1,Type=Integer,Description="Number of %s individuals with homozygous alternate genotypes (biallelic sites only).">' % (sex, sex))
            INFO_ADD.append('##INFO=<ID=%s_FREQ_HOMREF,Number=1,Type=Float,Description="%s homozygous reference genotype frequency (biallelic sites only).">' % (sex, sex))
            INFO_ADD.append('##INFO=<ID=%s_FREQ_HET,Number=1,Type=Float,Description="%s heterozygous genotype frequency (biallelic sites only).">' % (sex, sex))
            INFO_ADD.append('##INFO=<ID=%s_FREQ_HOMALT,Number=1,Type=Float,Description="%s homozygous alternate genotype frequency (biallelic sites only).">' % (sex, sex))
    if len(pops) > 0:
        INFO_ADD.append('##INFO=<ID=POPMAX_AF,Number=1,Type=Float,Description="Maximum allele frequency across any population (biallelic sites only).">')
        for pop in pops:
            INFO_ADD.append('##INFO=<ID=%s_AN,Number=1,Type=Integer,Description="Total number of %s alleles genotyped (for biallelic sites) or %s individuals with copy-state estimates (for multiallelic sites).">' % (pop, pop, pop))
            INFO_ADD.append('##INFO=<ID=%s_AC,Number=A,Type=Integer,Description="Number of non-reference %s alleles observed (for biallelic sites) or %s individuals at each copy state (for multiallelic sites).">' % (pop, pop, pop))
            INFO_ADD.append('##INFO=<ID=%s_AF,Number=A,Type=Float,Description="%s allele frequency (for biallelic sites) or %s copy-state frequency (for multiallelic sites).">' % (pop, pop, pop))
            INFO_ADD.append('##INFO=<ID=%s_N_BI_GENOS,Number=1,Type=Integer,Description="Total number of %s individuals with complete genotypes (biallelic sites only).">' % (pop, pop))
            INFO_ADD.append('##INFO=<ID=%s_N_HOMREF,Number=1,Type=Integer,Description="Number of %s individuals with homozygous reference genotypes (biallelic sites only).">' % (pop, pop))
            INFO_ADD.append('##INFO=<ID=%s_N_HET,Number=1,Type=Integer,Description="Number of %s individuals with heterozygous genotypes (biallelic sites only).">' % (pop, pop))
            INFO_ADD.append('##INFO=<ID=%s_N_HOMALT,Number=1,Type=Integer,Description="Number of %s individuals with homozygous alternate genotypes (biallelic sites only).">' % (pop, pop))
            INFO_ADD.append('##INFO=<ID=%s_FREQ_HOMREF,Number=1,Type=Float,Description="%s homozygous reference genotype frequency (biallelic sites only).">' % (pop, pop))
            INFO_ADD.append('##INFO=<ID=%s_FREQ_HET,Number=1,Type=Float,Description="%s heterozygous genotype frequency (biallelic sites only).">' % (pop, pop))
            INFO_ADD.append('##INFO=<ID=%s_FREQ_HOMALT,Number=1,Type=Float,Description="%s homozygous alternate genotype frequency (biallelic sites only).">' % (pop, pop))
            if len(sexes) > 0 and not args.no_combos:
                for sex in sexes:
                    INFO_ADD.append('##INFO=<ID=%s_AN,Number=1,Type=Integer,Description="Total number of %s alleles genotyped (for biallelic sites) or %s individuals with copy-state estimates (for multiallelic sites).">' % ('_'.join((pop, sex)), ' '.join((pop, sex)), ' '.join((pop, sex))))
                    INFO_ADD.append('##INFO=<ID=%s_AC,Number=A,Type=Integer,Description="Number of non-reference %s alleles observed (for biallelic sites) or %s individuals at each copy state (for multiallelic sites).">' % ('_'.join((pop, sex)), ' '.join((pop, sex)), ' '.join((pop, sex))))
                    INFO_ADD.append('##INFO=<ID=%s_AF,Number=A,Type=Float,Description="%s allele frequency (for biallelic sites) or %s copy-state frequency (for multiallelic sites).">' % ('_'.join((pop, sex)), ' '.join((pop, sex)), ' '.join((pop, sex))))
                    INFO_ADD.append('##INFO=<ID=%s_N_BI_GENOS,Number=1,Type=Integer,Description="Total number of %s individuals with complete genotypes (biallelic sites only).">' % ('_'.join((pop, sex)), ' '.join((pop, sex))))
                    INFO_ADD.append('##INFO=<ID=%s_N_HOMREF,Number=1,Type=Integer,Description="Number of %s individuals with homozygous reference genotypes (biallelic sites only).">' % ('_'.join((pop, sex)), ' '.join((pop, sex))))
                    INFO_ADD.append('##INFO=<ID=%s_N_HET,Number=1,Type=Integer,Description="Number of %s individuals with heterozygous genotypes (biallelic sites only).">' % ('_'.join((pop, sex)), ' '.join((pop, sex))))
                    INFO_ADD.append('##INFO=<ID=%s_N_HOMALT,Number=1,Type=Integer,Description="Number of %s individuals with homozygous alternate genotypes (biallelic sites only).">' % ('_'.join((pop, sex)), ' '.join((pop, sex))))
                    INFO_ADD.append('##INFO=<ID=%s_FREQ_HOMREF,Number=1,Type=Float,Description="%s homozygous reference genotype frequency (biallelic sites only).">' % ('_'.join((pop, sex)), ' '.join((pop, sex))))
                    INFO_ADD.append('##INFO=<ID=%s_FREQ_HET,Number=1,Type=Float,Description="%s heterozygous genotype frequency (biallelic sites only).">' % ('_'.join((pop, sex)), ' '.join((pop, sex))))
                    INFO_ADD.append('##INFO=<ID=%s_FREQ_HOMALT,Number=1,Type=Float,Description="%s homozygous alternate genotype frequency (biallelic sites only).">' % ('_'.join((pop, sex)), ' '.join((pop, sex))))
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
