#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""
Compute allele frequencies for sex & population combinations given an SV VCF
"""


import sys
import argparse
import pysam
from svtk import utils as svu
from collections import Counter
import pybedtools as pbt


ALLOWED_PARTS = set([
    'XX','XY','afr','ami',
    'amr','asj','eas','fin',
    'mid','nfe','remaining','sas'
])


def create_pop_dict(popfile):
    """
    Makes dictionary of sample-population pairs
    """

    pop_dict = {}

    for sample in popfile:
        pop_dict[sample.split('\t')[0]] = sample.split('\t')[1]

    return pop_dict


def in_par(record, parbt):
    """
    Check if variant overlaps pseudoautosomal region
    """

    # Sort start & end to handle edge cases and ensure end > start for pbt functionality
    sstart, send = [str(x) for x in sorted([record.start, record.stop])]
    svbt_str = '\t'.join([record.chrom, sstart, send]) + '\n'
    svbt = pbt.BedTool(svbt_str, from_string=True)
    if len(svbt.intersect(parbt)) > 0:
        return True
    else:
        return False


def update_sex_freqs(record, pop=None):
    """
    Recompute allele frequencies for variants on sex chromosomes outside of PARs
    """

    if pop is not None:
        m_prefix = '_'.join([pop, 'XY'])
        f_prefix = '_'.join([pop, 'XX'])
    else:
        m_prefix = 'XY'
        f_prefix = 'XX'

    m_an = record.info.get('AN_' + m_prefix, 0)
    m_ac = sum(record.info.get('AC_' + m_prefix, 0))
    # m_af = sum(record.info.get('AF_' + m_prefix, 0))

    f_an = record.info.get('AN_' + f_prefix , 0)
    f_ac = sum(record.info.get('AC_' + f_prefix , 0))
    # f_af = sum(record.info.get('AF_' + f_prefix , 0))

    adj_an = m_an + f_an
    adj_ac = m_ac + f_ac
    if adj_an > 0:
        adj_af = adj_ac / adj_an
    else:
        adj_af = 0

    if pop is None:
        record.info['AN'] = adj_an
        record.info['AC'] = (adj_ac, )
        record.info['AF'] = (adj_af, )
    else:
        record.info['AN_' + pop ] = adj_an
        record.info['AC_' + pop ] = (adj_ac, )
        record.info['AF_' + pop ] = (adj_af, )

    return record


def gather_allele_freqs(record, samples, males_set, females_set, parbt, pop_dict, pops,
                        sex_chroms, no_combos=False):
    """
    Wrapper to compute allele frequencies for all sex & population pairings
    """

    # Add PAR annotation to record (if optioned)
    if record.chrom in sex_chroms and len(parbt) > 0:
        if in_par(record, parbt):
            rec_in_par = True
            record.info['par'] = True
        else:
            rec_in_par = False
    else:
        rec_in_par = False

    # Get allele frequencies for all populations
    calc_allele_freq(record, samples)
    if len(males_set) > 0:
        if record.chrom in sex_chroms and not rec_in_par:
            calc_allele_freq(record, males_set, prefix='XY', hemi=True)
        else:
            calc_allele_freq(record, males_set, prefix='XY')
    if len(females_set) > 0:
        calc_allele_freq(record, females_set, prefix='XX')

    # Adjust global allele frequencies on sex chromosomes, if famfile provided
    if record.chrom in sex_chroms and not rec_in_par \
            and svu.is_biallelic(record) and len(males_set) + len(females_set) > 0:
        update_sex_freqs(record)
    
    # Get allele frequencies per population
    if len(pops) > 0:
        for pop in pops:
            pop_samps = [
                s for s in samples if pop_dict.get(s, None) == pop]
            calc_allele_freq(record, pop_samps, prefix=pop)
            if len(males_set) > 0 and not no_combos:
                if record.chrom in sex_chroms and not rec_in_par:
                    calc_allele_freq(record, list([s for s in pop_samps if s in males_set]),
                                     prefix=pop + '_XY', hemi=True)
                else:
                    calc_allele_freq(record, list([s for s in pop_samps if s in males_set]),
                                     prefix=pop + '_XY')
            if len(females_set) > 0 and not no_combos:
                calc_allele_freq(record, list([s for s in pop_samps if s in females_set]),
                                 prefix=pop + '_XX')

            # Adjust per-pop allele frequencies on sex chromosomes, if famfile provided
            if record.chrom in sex_chroms and not rec_in_par \
                    and svu.is_biallelic(record) and len(males_set) + len(females_set) > 0:
                update_sex_freqs(record, pop=pop)

        # Get grpmax AF biallelic sites only
        if svu.is_biallelic(record):
            AFs = [record.info['AF_{0}'.format(pop)][0] for pop in pops]
            grpmax = max(AFs)
            record.info['AF_grpmax'] = grpmax

            grp_label = pops[AFs.index(grpmax)]
            record.info['grpmax'] = grp_label

            record.info['AC_grpmax'] = int(record.info['AC_{0}'.format(grp_label)][0])
            record.info['AN_grpmax'] = int(record.info['AN_{0}'.format(grp_label)])
            record.info['nhomalt_grpmax'] = int(record.info['nhomalt_{0}'.format(grp_label)])

        return record


def calc_allele_freq(record, samples, prefix=None, hemi=False):
    """
    Computes allele frequencies for a single record based on a list of samples
    """

    # Treat biallelic and multiallelic records differently
    # For biallelic sites, count number of non-ref, non-no-call GTs
    if svu.is_biallelic(record):

        # Get all sample GTs
        GTs = [record.samples[s]['GT'] for s in samples]

        # Count alleles & genotypes
        AC = 0
        AN = 0
        n_alt_count_0 = 0
        n_alt_count_1 = 0
        n_alt_count_2 = 0
        n_gts_with_gt_0_alts = 0  # Used specifically for hemizygous sites
        for GT in GTs:
            AN += len([allele for allele in GT if allele !=
                       '.' and allele is not None])
            AC += len([allele for allele in GT if allele !=
                       '.' and allele != 0 and allele is not None])
            if GT == (0, 0):
                n_alt_count_0 += 1
            if len([allele for allele in GT if allele == 0 and allele != '.' and allele is not None]) == 1 \
                    and len([allele for allele in GT if allele != 0 and allele != '.' and allele is not None]) == 1:
                n_alt_count_1 += 1
                n_gts_with_gt_0_alts += 1
            if len([allele for allele in GT if allele != 0 and allele != '.' and allele is not None]) == 2:
                n_alt_count_2 += 1
                n_gts_with_gt_0_alts += 1

        # Adjust hemizygous allele number and allele count, if optioned
        if hemi:
            AN = AN / 2
            # For hemizygous sites, AC must be the sum of all non-reference *genotypes*, not alleles
            AC = n_gts_with_gt_0_alts

        # Calculate allele frequency
        if AN > 0:
            AF = AC / AN
            AF = round(AF, 6)
        else:
            AF = 0

        # Add AN, AC, and AF to INFO field
        record.info[('AN' if prefix is None else 'AN_' + prefix)] = AN
        record.info[('AC' if prefix is None else 'AC_' + prefix)] = AC
        record.info[('AF' if prefix is None else 'AF_' + prefix)] = AF

        # Calculate genotype frequencies
        n_bi_genos = n_alt_count_0 + n_alt_count_1 + n_alt_count_2
        if n_bi_genos > 0:
            freq_homref = n_alt_count_0 / n_bi_genos
            freq_het = n_alt_count_1 / n_bi_genos
            freq_homalt = n_alt_count_2 / n_bi_genos
        else:
            freq_homref = 0
            freq_het = 0
            freq_homalt = 0
        if hemi:
            freq_hemialt = freq_het + freq_homalt

        # Add n_bi_genos, n_homref, nhet, nhomalt, freq_homref, freq_het, and freq_homalt to INFO field
        record.info['n_bi_genos' + ('_' + prefix if prefix else '')] = n_bi_genos
        if hemi:
            record.info['nhemiref' + ('_' + prefix if prefix else '')] = n_alt_count_0
            record.info['nhemialt' + ('_' + prefix if prefix else '')] = n_gts_with_gt_0_alts
            record.info['freq_hemiref' + ('_' + prefix if prefix else '')] = freq_homref
            record.info['freq_hemialt' + ('_' + prefix if prefix else '')] = freq_hemialt
        record.info['nhomref'    + ('_' + prefix if prefix else '')] = n_alt_count_0
        record.info['nhet'       + ('_' + prefix if prefix else '')] = n_alt_count_1
        record.info['nhomalt'    + ('_' + prefix if prefix else '')] = n_alt_count_2
        record.info['freq_homref' + ('_' + prefix if prefix else '')] = freq_homref
        record.info['freq_het'    + ('_' + prefix if prefix else '')] = freq_het
        record.info['freq_homalt' + ('_' + prefix if prefix else '')] = freq_homalt

    # Multiallelic sites should reference FORMAT:CN rather than GT
    # Compute CN_NUMBER, CN_NONREF_COUNT, CN_NONREF_FREQ, and CN_COUNT/CN_FREQ for each copy state
    else:

        # Get all sample CNs and remove Nones
        CNs_wNones = [record.samples[s]['CN'] for s in samples]
        CNs = [c for c in CNs_wNones if c is not None and c not in '. NA'.split()]

        if len(CNs) == 0:
            nonnull_CNs, nonref_CN_count, nonref_CN_freq = [0] * 3
            CN_dist = (0, )
            CN_freqs = (0, )
            CN_status = (0, )
        else:
            # Count number of samples per CN and total CNs observed
            CN_counts = dict(Counter(CNs))
            nonnull_CNs = len(CNs)

            # Get max observed CN and enumerate counts/frequencies per CN as list starting from CN=0
            max_CN = max([int(k) for k, v in CN_counts.items()])
            CN_dist = [int(CN_counts.get(k, 0)) for k in range(max_CN + 1)]
            CN_freqs = [round(v / nonnull_CNs, 6) for v in CN_dist]
            CN_status = [s for s in range(max_CN + 1)]

            # Get total non-reference CN counts and freq
            if hemi:
                ref_CN = 1
            else:
                ref_CN = 2
            nonref_CN_count = sum([int(CN_counts.get(k, 0)) for k in range(max_CN + 1) if k != ref_CN])
            nonref_CN_freq = round(nonref_CN_count / nonnull_CNs, 6)

        # Add values to INFO field
        record.info['CN_number' + ('_' + prefix if prefix else '')] = nonnull_CNs
        record.info['CN_count' + ('_' + prefix if prefix else '')] = tuple(CN_dist)
        record.info['CN_freq' + ('_' + prefix if prefix else '')] = tuple(CN_freqs)
        record.info['CN_status' + ('_' + prefix if prefix else '')] = tuple(CN_status)
        record.info['CN_nnonref' + ('_' + prefix if prefix else '')] = nonref_CN_count
        record.info['cn_freq_nonref' + ('_' + prefix if prefix else '')] = nonref_CN_freq

    return record


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf', help='Input vcf. Also accepts "stdin" and "-".')
    parser.add_argument('-p', '--popfile', help='Two-column file of samples & ' +
                        'their population assignments. A "." denotes no assignment.',
                        default=None)
    parser.add_argument('-f', '--famfile', help='Input .fam file (used for sex-specific AFs).',
                        default=None)
    parser.add_argument('--no-combos', help='Do not compute combinations of populations ' +
                        'and sexes.', action='store_true', default=False)
    parser.add_argument('--allosomes-list', help='TSV of sex chromosomes (used for ' +
                        'sex-specific AFs).', default=None)
    parser.add_argument('--par', help='BED file of pseudoautosomal regions (used ' +
                        'for sex-specific AFs).', default=None)
    parser.add_argument(
        'fout', help='Output vcf. Also accepts "stdout" and "-".')
    args = parser.parse_args()

    # Open connections to input VCF
    if args.vcf in '- stdin'.split():
        vcf = pysam.VariantFile(sys.stdin)
    else:
        vcf = pysam.VariantFile(args.vcf)

    # Get list of all samples in vcf
    samples_list = list(vcf.header.samples)

    # Get lists of males and females
    parbt = pbt.BedTool('', from_string=True)
    if args.famfile is not None:
        famfile = [line.rstrip('\n') for line in open(args.famfile)]
        males_set = set([line.split('\t')[1]
                 for line in famfile if line.split('\t')[4] == '1'])
        males_set = set(s for s in samples_list if s in males_set)
        females_set = set([line.split('\t')[1]
                   for line in famfile if line.split('\t')[4] == '2'])
        females_set = set(s for s in samples_list if s in females_set)
        sexes = 'XY XX'.split()
        if args.par is not None:
            parbt = pbt.BedTool(args.par)

    else:
        males_set = set()
        females_set = set()
        sexes = list()

    # Get dictionary of populations
    if args.popfile is not None:
        popfile = [line.rstrip('\n') for line in open(args.popfile)]
        pop_dict = create_pop_dict(popfile)
        pops = list(set(pop_dict.values()))
        pops = sorted([p for p in pops if p != "."])
        for label in pops:
            if label not in ALLOWED_PARTS:
                raise ValueError(f"Invalid label: '{label}'.")
    else:
        pop_dict = {}
        pops = []


    # Get list of sex chromosomes, if optioned
    if args.allosomes_list is not None:
        sex_chroms = [
            l.split('\t')[0]
            for l in open(args.allosomes_list).readlines()
        ]
    else:
        sex_chroms = 'X Y chrX chrY'.split()

    # Add relevant fields to header
    INFO_ADD = [
        '##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles genotyped (biallelic sites only).">',
        '##INFO=<ID=AC,Number=A,Type=Integer,Description="Number of non-reference alleles observed (biallelic sites only).">',
        '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency (biallelic sites only).">',
        '##INFO=<ID=n_bi_genos,Number=1,Type=Integer,Description="Total number of samples with complete genotypes (biallelic sites only).">',
        '##INFO=<ID=nhomref,Number=1,Type=Integer,Description="Number of samples with homozygous reference genotypes (biallelic sites only).">',
        '##INFO=<ID=nhet,Number=1,Type=Integer,Description="Number of samples with heterozygous genotypes (biallelic sites only).">',
        '##INFO=<ID=nhomalt,Number=1,Type=Integer,Description="Number of samples with homozygous alternate genotypes (biallelic sites only).">',
        '##INFO=<ID=freq_homref,Number=1,Type=Float,Description="Homozygous reference genotype frequency (biallelic sites only).">',
        '##INFO=<ID=freq_het,Number=1,Type=Float,Description="Heterozygous genotype frequency (biallelic sites only).">',
        '##INFO=<ID=freq_homalt,Number=1,Type=Float,Description="Homozygous alternate genotype frequency (biallelic sites only).">',
        '##INFO=<ID=CN_number,Number=1,Type=Integer,Description="Total number of samples with estimated copy numbers (multiallelic CNVs only).">',
        '##INFO=<ID=CN_count,Number=.,Type=Integer,Description="Number of samples observed at each copy state, starting from CN=0 (multiallelic CNVs only).">',
        '##INFO=<ID=CN_status,Number=.,Type=Integer,Description="Copy states corresponding to CN_count, CN_freq: 0,1,...,maximum observed copy state (multiallelic CNVs only).">',
        '##INFO=<ID=CN_freq,Number=.,Type=Float,Description="Frequency of samples observed at each copy state, starting from CN=0 (multiallelic CNVs only).">',
        '##INFO=<ID=CN_nnonref,Number=1,Type=Integer,Description="Number of samples with non-reference copy states (multiallelic CNVs only).">',
        '##INFO=<ID=CN_freq_nonref,Number=1,Type=Float,Description="Frequency of samples with non-reference copy states (multiallelic CNVs only).">'
    ]
    if len(sexes) > 0:
        for sex in sexes:
            INFO_ADD.append('##INFO=<ID=AN_%s,Number=1,Type=Integer,Description="Total number of %s alleles genotyped (biallelic sites only).">' % (sex, sex))
            INFO_ADD.append('##INFO=<ID=AC_%s,Number=A,Type=Integer,Description="Number of non-reference %s alleles observed (biallelic sites only).">' % (sex, sex))
            INFO_ADD.append('##INFO=<ID=AF_%s,Number=A,Type=Float,Description="%s allele frequency (biallelic sites only).">' % (sex, sex))
            INFO_ADD.append('##INFO=<ID=n_bi_genos_%s,Number=1,Type=Integer,Description="Total number of %s samples with complete genotypes (biallelic sites only).">' % (sex, sex))
            INFO_ADD.append('##INFO=<ID=nhomref_%s,Number=1,Type=Integer,Description="Number of %s samples with homozygous reference genotypes (biallelic sites only).">' % (sex, sex))
            INFO_ADD.append('##INFO=<ID=nhet_%s,Number=1,Type=Integer,Description="Number of %s samples with heterozygous genotypes (biallelic sites only).">' % (sex, sex))
            INFO_ADD.append('##INFO=<ID=nhomalt_%s,Number=1,Type=Integer,Description="Number of %s samples with homozygous alternate genotypes (biallelic sites only).">' % (sex, sex))
            INFO_ADD.append('##INFO=<ID=freq_homref_%s,Number=1,Type=Float,Description="%s homozygous reference genotype frequency (biallelic sites only).">' % (sex, sex))
            INFO_ADD.append('##INFO=<ID=freq_het_%s,Number=1,Type=Float,Description="%s heterozygous genotype frequency (biallelic sites only).">' % (sex, sex))
            INFO_ADD.append('##INFO=<ID=freq_homalt_%s,Number=1,Type=Float,Description="%s homozygous alternate genotype frequency (biallelic sites only).">' % (sex, sex))
            INFO_ADD.append('##INFO=<ID=CN_number_%s,Number=1,Type=Integer,Description="Total number of %s samples with estimated copy numbers (multiallelic CNVs only).">' % (sex, sex))
            INFO_ADD.append('##INFO=<ID=CN_count_%s,Number=.,Type=Integer,Description="Number of %s samples observed at each copy state, starting from CN=0 (multiallelic CNVs only).">' % (sex, sex))
            INFO_ADD.append('##INFO=<ID=CN_status_%s,Number=.,Type=Integer,Description="Copy states corresponding to CN_COUNT_%s, CN_FREQ_%s: 0,1,...,maximum observed copy state (multiallelic CNVs only).">' % (sex, sex, sex))
            INFO_ADD.append('##INFO=<ID=CN_freq_%s,Number=.,Type=Float,Description="Frequency of %s samples observed at each copy state, starting from CN=0 (multiallelic CNVs only).">' % (sex, sex))
            INFO_ADD.append('##INFO=<ID=CN_nnonref_%s,Number=1,Type=Integer,Description="Number of %s samples with non-reference copy states (multiallelic CNVs only).">' % (sex, sex))
            INFO_ADD.append('##INFO=<ID=CN_freq_nonref_%s,Number=1,Type=Float,Description="Frequency of %s samples with non-reference copy states (multiallelic CNVs only).">' % (sex, sex))
            if sex == 'XY':
                INFO_ADD.append('##INFO=<ID=nhemiref_%s,Number=1,Type=Integer,Description="Number of %s samples with hemizygous reference genotypes (biallelic sites only).">' % (sex, sex))
                INFO_ADD.append('##INFO=<ID=nhemialt_%s,Number=1,Type=Integer,Description="Number of %s samples with hemizygous alternate genotypes (biallelic sites only).">' % (sex, sex))
                INFO_ADD.append('##INFO=<ID=freq_hemiref_%s,Number=1,Type=Float,Description="%s hemizygous reference genotype frequency (biallelic sites only).">' % (sex, sex))
                INFO_ADD.append('##INFO=<ID=freq_hemialt_%s,Number=1,Type=Float,Description="%s hemizygous alternate genotype frequency (biallelic sites only).">' % (sex, sex))
                if len(parbt) > 0:
                    INFO_ADD.append('##INFO=<ID=par,Number=0,Type=Flag,Description="Variant overlaps pseudoautosomal region.">')
    if len(pops) > 0:
        INFO_ADD.append('##INFO=<ID=AF_grpmax,Number=1,Type=Float,Description="Maximum allele frequency across any population (biallelic sites only).">')
        INFO_ADD.append('##INFO=<ID=grpmax,Number=1,Type=String,Description="Population label with maximum allele frequency (biallelic sites only).">')
        INFO_ADD.append('##INFO=<ID=AC_grpmax,Number=1,Type=Integer,Description="Allele count for the grpmax population (biallelic sites only).">')
        INFO_ADD.append('##INFO=<ID=AN_grpmax,Number=1,Type=Integer,Description="Allele number for the grpmax population (biallelic sites only).">')
        INFO_ADD.append('##INFO=<ID=nhomalt_grpmax,Number=1,Type=Integer,Description="Number of homozygous-alternate genotypes in the grpmax population (biallelic sites only).">')
        for pop in pops:
            INFO_ADD.append('##INFO=<ID=AN_%s,Number=1,Type=Integer,Description="Total number of %s alleles genotyped (biallelic sites only).">' % (pop, pop))
            INFO_ADD.append('##INFO=<ID=AC_%s,Number=A,Type=Integer,Description="Number of non-reference %s alleles observed (biallelic sites only).">' % (pop, pop))
            INFO_ADD.append('##INFO=<ID=AF_%s,Number=A,Type=Float,Description="%s allele frequency (biallelic sites only).">' % (pop, pop))
            INFO_ADD.append('##INFO=<ID=n_bi_genos_%s,Number=1,Type=Integer,Description="Total number of %s samples with complete genotypes (biallelic sites only).">' % (pop, pop))
            INFO_ADD.append('##INFO=<ID=nhomref_%s,Number=1,Type=Integer,Description="Number of %s samples with homozygous reference genotypes (biallelic sites only).">' % (pop, pop))
            INFO_ADD.append('##INFO=<ID=nhet_%s,Number=1,Type=Integer,Description="Number of %s samples with heterozygous genotypes (biallelic sites only).">' % (pop, pop))
            INFO_ADD.append('##INFO=<ID=nhomalt_%s,Number=1,Type=Integer,Description="Number of %s samples with homozygous alternate genotypes (biallelic sites only).">' % (pop, pop))
            INFO_ADD.append('##INFO=<ID=freq_homref_%s,Number=1,Type=Float,Description="%s homozygous reference genotype frequency (biallelic sites only).">' % (pop, pop))
            INFO_ADD.append('##INFO=<ID=freq_het_%s,Number=1,Type=Float,Description="%s heterozygous genotype frequency (biallelic sites only).">' % (pop, pop))
            INFO_ADD.append('##INFO=<ID=freq_homalt_%s,Number=1,Type=Float,Description="%s homozygous alternate genotype frequency (biallelic sites only).">' % (pop, pop))
            INFO_ADD.append('##INFO=<ID=CN_number_%s,Number=1,Type=Integer,Description="Total number of %s samples with estimated copy numbers (multiallelic CNVs only).">' % (pop, pop))
            INFO_ADD.append('##INFO=<ID=CN_count_%s,Number=.,Type=Integer,Description="Number of %s samples observed at each copy state, starting from CN=0 (multiallelic CNVs only).">' % (pop, pop))
            INFO_ADD.append('##INFO=<ID=CN_status_%s,Number=.,Type=Integer,Description="Copy states corresponding to CN_COUNT_%s, CN_FREQ_%s: 0,1,...,maximum observed copy state (multiallelic CNVs only).">' % (pop, pop, pop))
            INFO_ADD.append('##INFO=<ID=CN_freq_%s,Number=.,Type=Float,Description="Frequency of %s samples observed at each copy state, starting from CN=0 (multiallelic CNVs only).">' % (pop, pop))
            INFO_ADD.append('##INFO=<ID=CN_nnonref_%s,Number=1,Type=Integer,Description="Number of %s samples with non-reference copy states (multiallelic CNVs only).">' % (pop, pop))
            INFO_ADD.append('##INFO=<ID=CN_freq_nonref_%s,Number=1,Type=Float,Description="Frequency of %s samples with non-reference copy states (multiallelic CNVs only).">' % (pop, pop))
            if len(sexes) > 0 and not args.no_combos:
                for sex in sexes:
                    INFO_ADD.append('##INFO=<ID=AN_%s,Number=1,Type=Integer,Description="Total number of %s alleles genotyped (biallelic sites only).">' % ('_'.join((pop, sex)), ' '.join((pop, sex))))
                    INFO_ADD.append('##INFO=<ID=AC_%s,Number=A,Type=Integer,Description="Number of non-reference %s alleles observed (biallelic sites only).">' % ('_'.join((pop, sex)), ' '.join((pop, sex))))
                    INFO_ADD.append('##INFO=<ID=AF_%s,Number=A,Type=Float,Description="%s allele frequency (biallelic sites only).">' % ('_'.join((pop, sex)), ' '.join((pop, sex))))
                    INFO_ADD.append('##INFO=<ID=n_bi_genos_%s,Number=1,Type=Integer,Description="Total number of %s samples with complete genotypes (biallelic sites only).">' % ('_'.join((pop, sex)), ' '.join((pop, sex))))
                    INFO_ADD.append('##INFO=<ID=nhomref_%s,Number=1,Type=Integer,Description="Number of %s samples with homozygous reference genotypes (biallelic sites only).">' % ('_'.join((pop, sex)), ' '.join((pop, sex))))
                    INFO_ADD.append('##INFO=<ID=nhet_%s,Number=1,Type=Integer,Description="Number of %s samples with heterozygous genotypes (biallelic sites only).">' % ('_'.join((pop, sex)), ' '.join((pop, sex))))
                    INFO_ADD.append('##INFO=<ID=nhomalt_%s,Number=1,Type=Integer,Description="Number of %s samples with homozygous alternate genotypes (biallelic sites only).">' % ('_'.join((pop, sex)), ' '.join((pop, sex))))
                    INFO_ADD.append('##INFO=<ID=freq_homref_%s,Number=1,Type=Float,Description="%s homozygous reference genotype frequency (biallelic sites only).">' % ('_'.join((pop, sex)), ' '.join((pop, sex))))
                    INFO_ADD.append('##INFO=<ID=freq_het_%s,Number=1,Type=Float,Description="%s heterozygous genotype frequency (biallelic sites only).">' % ('_'.join((pop, sex)), ' '.join((pop, sex))))
                    INFO_ADD.append('##INFO=<ID=freq_homalt_%s,Number=1,Type=Float,Description="%s homozygous alternate genotype frequency (biallelic sites only).">' % ('_'.join((pop, sex)), ' '.join((pop, sex))))
                    INFO_ADD.append('##INFO=<ID=CN_number_%s,Number=1,Type=Integer,Description="Total number of %s samples with estimated copy numbers (multiallelic CNVs only).">' % ('_'.join((pop, sex)), ' '.join((pop, sex))))
                    INFO_ADD.append('##INFO=<ID=CN_count_%s,Number=.,Type=Integer,Description="Number of %s samples observed at each copy state, starting from CN=0 (multiallelic CNVs only).">' % ('_'.join((pop, sex)), ' '.join((pop, sex))))
                    INFO_ADD.append('##INFO=<ID=CN_status_%s,Number=.,Type=Integer,Description="Copy states corresponding to CN_COUNT_%s, CN_FREQ_%s: 0,1,...,maximum observed copy state (multiallelic CNVs only).">' % ('_'.join((pop, sex)), '_'.join((pop, sex)), '_'.join((pop, sex))))
                    INFO_ADD.append('##INFO=<ID=CN_freq_%s,Number=.,Type=Float,Description="Frequency of %s samples observed at each copy state, starting from CN=0 (multiallelic CNVs only).">' % ('_'.join((pop, sex)), ' '.join((pop, sex))))
                    INFO_ADD.append('##INFO=<ID=CN_nnonref_%s,Number=1,Type=Integer,Description="Number of %s samples with non-reference copy states (multiallelic CNVs only).">' % ('_'.join((pop, sex)), ' '.join((pop, sex))))
                    INFO_ADD.append('##INFO=<ID=CN_freq_nonref_%s,Number=1,Type=Float,Description="Frequency of %s samples with non-reference copy states (multiallelic CNVs only).">' % ('_'.join((pop, sex)), ' '.join((pop, sex))))
                    if sex == 'XY':
                        INFO_ADD.append('##INFO=<ID=nhemiref_%s,Number=1,Type=Integer,Description="Number of %s samples with hemizygous reference genotypes (biallelic sites only).">' % ('_'.join((pop, sex)), ' '.join((pop, sex))))
                        INFO_ADD.append('##INFO=<ID=nhemialt_%s,Number=1,Type=Integer,Description="Number of %s samples with hemizygous alternate genotypes (biallelic sites only).">' % ('_'.join((pop, sex)), ' '.join((pop, sex))))
                        INFO_ADD.append('##INFO=<ID=freq_hemiref_%s,Number=1,Type=Float,Description="%s hemizygous reference genotype frequency (biallelic sites only).">' % ('_'.join((pop, sex)), ' '.join((pop, sex))))
                        INFO_ADD.append('##INFO=<ID=freq_hemialt_%s,Number=1,Type=Float,Description="%s hemizygous alternate genotype frequency (biallelic sites only).">' % ('_'.join((pop, sex)), ' '.join((pop, sex))))

    for line in INFO_ADD:
        vcf.header.add_line(line)

    # Prep output VCF
    if args.fout in '- stdout'.split():
        fout = pysam.VariantFile(sys.stdout, 'w', header=vcf.header)
    else:
        fout = pysam.VariantFile(args.fout, 'w', header=vcf.header)

    # Get allele frequencies for each record & write to new VCF
    for r in vcf.fetch():
        newrec = gather_allele_freqs(r, samples_list, males_set, females_set, parbt, pop_dict,
                                     pops, sex_chroms, args.no_combos)
        fout.write(newrec)

    fout.close()


if __name__ == '__main__':
    main()
