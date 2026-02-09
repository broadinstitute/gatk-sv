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


ALLOWED_POPS = set(['afr', 'ami', 'amr', 'asj', 'eas', 'fin', 'mid', 'nfe', 'sas', 'oth'])


def check_motifs_canonically_unique(motifs):
    """
    Check if motifs are canonically unique.
    Two motifs are NOT canonically unique if they have the same character counts
    and one can be created from the other by rotation.
    """

    for i in range(len(motifs)):
        for j in range(i + 1, len(motifs)):
            motif1 = motifs[i]
            motif2 = motifs[j]
            if len(motif1) != len(motif2):
                continue
            if Counter(motif1) != Counter(motif2):
                continue
            doubled = motif1 + motif1
            if motif2 in doubled:
                return False
    return True


def create_pop_dict(popfile):
    """
    Makes dictionary of sample-population pairs
    """

    pop_dict = {}
    for sample in popfile:
        pop_dict[sample.split('\t')[0]] = sample.split('\t')[1]
    return pop_dict


def load_lps_dict(lpsfile):
    """
    Load LPS TSV file into a nested dictionary structure:
    {variant_id: {sample_name: [val1, val2, ...]}}
    """
    lps_dict = {}
    lines = [line.rstrip('\n').rstrip('\r') for line in open(lpsfile) if line.strip()]
    
    if len(lines) == 0:
        return lps_dict
    
    header = lines[-1].split('\t')
    sample_names = header[2:]
    
    for line in lines[:-1]:
        fields = line.split('\t')
        if len(fields) < 2:
            continue
        
        variant_id = fields[0]
        sample_values = fields[2:]
        
        lps_dict[variant_id] = {}
        for i, sample_name in enumerate(sample_names):
            if i < len(sample_values):
                val_str = sample_values[i].strip()
                if val_str == '.' or val_str == '':
                    lps_dict[variant_id][sample_name] = None
                else:
                    lps_dict[variant_id][sample_name] = [int(v.strip()) for v in val_str.split(',')]
            else:
                lps_dict[variant_id][sample_name] = None
    
    return lps_dict


def in_par(record, parbt):
    """
    Check if variant overlaps pseudoautosomal region
    """

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
    m_ac = sum(record.info.get('AC_' + m_prefix, (0, )))
    # m_af = sum(record.info.get('AF_' + m_prefix, 0))

    f_an = record.info.get('AN_' + f_prefix , 0)
    f_ac = sum(record.info.get('AC_' + f_prefix , (0, )))
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
                        sex_chroms, lps_dict, no_combos=False):
    """
    Wrapper to compute allele frequencies for all sex & population pairings
    """

    # Add PAR annotation to record if applicable
    if record.chrom in sex_chroms and len(parbt) > 0:
        if in_par(record, parbt):
            rec_in_par = True
            record.info['par'] = True
        else:
            rec_in_par = False
    else:
        rec_in_par = False

    # Get allele frequencies for all populations
    calc_allele_freq(record, samples, lps_dict=lps_dict)
    if len(males_set) > 0:
        if record.chrom in sex_chroms and not rec_in_par:
            calc_allele_freq(record, males_set, prefix='XY', hemi=True, lps_dict=lps_dict)
        else:
            calc_allele_freq(record, males_set, prefix='XY', lps_dict=lps_dict)
    if len(females_set) > 0:
        calc_allele_freq(record, females_set, prefix='XX', lps_dict=lps_dict)

    # Adjust global allele frequencies on sex chromosomes, if famfile provided
    if record.chrom in sex_chroms and not rec_in_par \
            and svu.is_biallelic(record) and len(males_set) + len(females_set) > 0:
        update_sex_freqs(record)
    
    # Get allele frequencies per population
    if len(pops) > 0:
        for pop in pops:
            pop_samps = [
                s for s in samples if pop_dict.get(s, None) == pop]
            calc_allele_freq(record, pop_samps, prefix=pop, lps_dict=lps_dict)
            if len(males_set) > 0 and not no_combos:
                if record.chrom in sex_chroms and not rec_in_par:
                    calc_allele_freq(record, list([s for s in pop_samps if s in males_set]),
                                     prefix=pop + '_XY', hemi=True, lps_dict=lps_dict)
                else:
                    calc_allele_freq(record, list([s for s in pop_samps if s in males_set]),
                                     prefix=pop + '_XY', lps_dict=lps_dict)
            if len(females_set) > 0 and not no_combos:
                calc_allele_freq(record, list([s for s in pop_samps if s in females_set]),
                                 prefix=pop + '_XX', lps_dict=lps_dict)

            # Adjust per-pop allele frequencies on sex chromosomes, if famfile provided
            if record.chrom in sex_chroms and not rec_in_par \
                    and svu.is_biallelic(record) and len(males_set) + len(females_set) > 0:
                update_sex_freqs(record, pop=pop)

        # Get grpmax AF biallelic sites only
        if svu.is_biallelic(record) and (record.alts and len(record.alts) == 1):
            AFs = [record.info['AF_{0}'.format(pop)][0] for pop in pops]
            grpmax = max(AFs)
            record.info['AF_grpmax'] = grpmax

            grp_label = pops[AFs.index(grpmax)]
            record.info['grpmax'] = grp_label

            record.info['AC_grpmax'] = int(record.info['AC_{0}'.format(grp_label)][0])
            record.info['AN_grpmax'] = int(record.info['AN_{0}'.format(grp_label)])
            record.info['nhomalt_grpmax'] = int(record.info['nhomalt_{0}'.format(grp_label)])

        return record


def calc_allele_freq(record, samples, prefix=None, hemi=False, lps_dict=None):
    """
    Computes allele frequencies for a single record based on a list of samples
    """

    # Treat monoallelic, biallelic and multiallelic records differently
    #  - For monoallelic sites, do nothing
    #  - For biallelic sites, count number of non-ref, non-no-call GTs
    #  - For multiallelic sites, compute allele counts & frequencies for all alleles
    if not record.alts or len(record.alts) == 0:
        return record

    elif svu.is_biallelic(record):
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
        record.info[('AC' if prefix is None else 'AC_' + prefix)] = (AC, )
        record.info[('AF' if prefix is None else 'AF_' + prefix)] = (AF, )

        # Add NCR to INFO field
        n_samples = len(samples)
        if n_samples > 0:
            if hemi:
                ncr = 1 - (AN / n_samples)
            else:
                ncr = 1 - (AN / (2 * n_samples))
            ncr = round(ncr, 6)
        else:
            ncr = 1.0
        record.info[('NCR' if prefix is None else 'NCR_' + prefix)] = ncr

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

    else:
        # Get all sample GTs and filter out missing genotypes
        GTs_raw = [(s, record.samples[s]['GT']) for s in samples]
        valid_GTs = [(s, gt) for s, gt in GTs_raw if gt is not None and any(a is not None for a in gt)]

        # Flatten all alleles into a single list of indices
        all_alleles = [int(a) for s, gt in valid_GTs for a in gt if a is not None and a != '.']
        allele_counts = Counter(all_alleles)
        
        # Calculate AN, AC and AF for each allele
        AN = len(all_alleles)
        num_alleles = len(record.alts) + 1
        allele_indices = list(range(num_alleles))
        ac_allele = [int(allele_counts.get(k, 0)) for k in allele_indices]
        af_allele = [round(v / AN, 6) if AN > 0 else 0 for v in ac_allele]

        # Add AN, AC, and AF to INFO field
        s_suffix = ('_' + prefix if prefix else '')
        record.info[('AN' if prefix is None else 'AN_' + prefix)] = AN
        record.info[('AC' if prefix is None else 'AC_' + prefix)] = tuple(ac_allele[1:])
        record.info[('AF' if prefix is None else 'AF_' + prefix)] = tuple(af_allele[1:])

        # Add NCR to INFO field
        n_samples = len(samples)
        if n_samples > 0:
            ncr = 1 - (AN / (2 * n_samples))
            ncr = round(ncr, 6)
        else:
            ncr = 1.0
        record.info[('NCR' if prefix is None else 'NCR_' + prefix)] = ncr
        
        # Compute AP_allele, MC_allele and LPS_allele only for all samples
        if prefix is None:
            # Compute AP_allele
            ap_allele = []
            for allele_idx in allele_indices:
                ap_val = None
                for sample_name, gt in valid_GTs:
                    if allele_idx in gt:
                        gt_list = list(gt)
                        if allele_idx in gt_list:
                            pos_in_gt = gt_list.index(allele_idx)
                            ap_field = record.samples[sample_name].get('AP')
                            if ap_field is not None:
                                ap_values = ap_field.split(',') if isinstance(ap_field, str) else ap_field
                                if pos_in_gt < len(ap_values) and ap_values[pos_in_gt] is not None:
                                    ap_val = float(ap_values[pos_in_gt])
                            break
                ap_allele.append(ap_val if ap_val is not None else 0)
            record.info['AP_allele'] = tuple(ap_allele)
            
            # Compute MC_allele and LPS_allele for single motif sites
            motifs_field = record.info.get('MOTIFS')
            motifs = list(motifs_field) if isinstance(motifs_field, tuple) or isinstance(motifs_field, list) else motifs_field.split(',')
            if len(motifs) == 1 and check_motifs_canonically_unique(motifs):
                mc_allele = []
                for allele_idx in allele_indices:
                    mc_val = None
                    for sample_name, gt in valid_GTs:
                        if allele_idx in gt:
                            gt_list = list(gt)
                            if allele_idx in gt_list:
                                pos_in_gt = gt_list.index(allele_idx)
                                mc_field = record.samples[sample_name].get('MC')
                                if mc_field is not None:
                                    mc_values = list(mc_field) if isinstance(mc_field, tuple) or isinstance(mc_field, list) else mc_field.split(',')
                                    if pos_in_gt < len(mc_values) and mc_values[pos_in_gt] is not None:
                                        mc_val = int(mc_values[pos_in_gt])
                                break
                    mc_allele.append(mc_val if mc_val is not None else 0)
                record.info['MC_allele'] = tuple(mc_allele)

                if lps_dict is not None and record.id in lps_dict:
                    lps_allele = []
                    for allele_idx in allele_indices:
                        lps_val = None
                        for sample_name, gt in valid_GTs:
                            if allele_idx in gt:
                                gt_list = list(gt)
                                if allele_idx in gt_list:
                                    pos_in_gt = gt_list.index(allele_idx)
                                    lps_field = lps_dict[record.id].get(sample_name)
                                    if lps_field is not None and pos_in_gt < len(lps_field):
                                        lps_val = lps_field[pos_in_gt]
                                    break
                        lps_allele.append(lps_val if lps_val is not None else 0)
                    record.info['LPS_allele'] = tuple(lps_allele)

    return record


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('vcf', help='Input vcf. Also accepts "stdin" and "-".')
    parser.add_argument('fout', help='Output vcf. Also accepts "stdout" and "-".')
    parser.add_argument('-p', '--popfile', help='Two-column file of samples & population assignments', default=None)
    parser.add_argument('-f', '--famfile', help='Input .fam file (used for sex-specific AFs).', default=None)
    parser.add_argument('-l', '--lpsfile', help='TSV file of LPS values per sample (for multiallelic sites).', default=None)
    parser.add_argument('--no-combos', help='Do not compute combinations of populations and sexes.', action='store_true', default=False)
    parser.add_argument('--allosomes-list', help='TSV of sex chromosomes (used for sex-specific AFs).', default=None)
    parser.add_argument('--par', help='BED file of pseudoautosomal regions (used for sex-specific AFs).', default=None)
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
            if label not in ALLOWED_POPS:
                raise ValueError(f"Invalid label: '{label}'.")
    else:
        pop_dict = {}
        pops = []

    # Load LPS dictionary if provided
    if args.lpsfile is not None:
        lps_dict = load_lps_dict(args.lpsfile)
    else:
        lps_dict = None


    # Get list of sex chromosomes, if optioned
    if args.allosomes_list is not None:
        sex_chroms = [
            l.split('\t')[0]
            for l in open(args.allosomes_list).readlines()
        ]
    else:
        sex_chroms = 'X Y chrX chrY'.split()

    # Add base fields to header
    INFO_ADD = [
        '##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles genotyped.">',
        '##INFO=<ID=AC,Number=A,Type=Integer,Description="Number of alleles observed.">',
        '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency.">',
        '##INFO=<ID=NCR,Number=1,Type=Float,Description="No-call rate.">',
        '##INFO=<ID=n_bi_genos,Number=1,Type=Integer,Description="Total number of samples with complete genotypes (biallelic sites only).">',
        '##INFO=<ID=nhomref,Number=1,Type=Integer,Description="Number of samples with homozygous reference genotypes (biallelic sites only).">',
        '##INFO=<ID=nhet,Number=1,Type=Integer,Description="Number of samples with heterozygous genotypes (biallelic sites only).">',
        '##INFO=<ID=nhomalt,Number=1,Type=Integer,Description="Number of samples with homozygous alternate genotypes (biallelic sites only).">',
        '##INFO=<ID=freq_homref,Number=1,Type=Float,Description="Homozygous reference genotype frequency (biallelic sites only).">',
        '##INFO=<ID=freq_het,Number=1,Type=Float,Description="Heterozygous genotype frequency (biallelic sites only).">',
        '##INFO=<ID=freq_homalt,Number=1,Type=Float,Description="Homozygous alternate genotype frequency (biallelic sites only).">',
        '##INFO=<ID=AP_allele,Number=.,Type=Float,Description="Allele purity for each allele index (multiallelic sites only).">',
        '##INFO=<ID=MC_allele,Number=.,Type=Integer,Description="Motif count for each allele index (multiallelic sites only).">',
        '##INFO=<ID=LPS_allele,Number=.,Type=Integer,Description="Longest polymer sequence for each allele index (multiallelic sites only).">'
    ]

    # Add sex fields to header
    if len(sexes) > 0:
        for sex in sexes:
            INFO_ADD.append('##INFO=<ID=AN_%s,Number=1,Type=Integer,Description="Total number of %s alleles genotyped.">' % (sex, sex))
            INFO_ADD.append('##INFO=<ID=AC_%s,Number=A,Type=Integer,Description="Number of %s alleles observed.">' % (sex, sex))
            INFO_ADD.append('##INFO=<ID=AF_%s,Number=A,Type=Float,Description="%s allele frequency.">' % (sex, sex))
            INFO_ADD.append('##INFO=<ID=NCR_%s,Number=1,Type=Float,Description="%s no-call rate.">' % (sex, sex))
            INFO_ADD.append('##INFO=<ID=n_bi_genos_%s,Number=1,Type=Integer,Description="Total number of %s samples with complete genotypes (biallelic sites only).">' % (sex, sex))
            INFO_ADD.append('##INFO=<ID=nhomref_%s,Number=1,Type=Integer,Description="Number of %s samples with homozygous reference genotypes (biallelic sites only).">' % (sex, sex))
            INFO_ADD.append('##INFO=<ID=nhet_%s,Number=1,Type=Integer,Description="Number of %s samples with heterozygous genotypes (biallelic sites only).">' % (sex, sex))
            INFO_ADD.append('##INFO=<ID=nhomalt_%s,Number=1,Type=Integer,Description="Number of %s samples with homozygous alternate genotypes (biallelic sites only).">' % (sex, sex))
            INFO_ADD.append('##INFO=<ID=freq_homref_%s,Number=1,Type=Float,Description="%s homozygous reference genotype frequency (biallelic sites only).">' % (sex, sex))
            INFO_ADD.append('##INFO=<ID=freq_het_%s,Number=1,Type=Float,Description="%s heterozygous genotype frequency (biallelic sites only).">' % (sex, sex))
            INFO_ADD.append('##INFO=<ID=freq_homalt_%s,Number=1,Type=Float,Description="%s homozygous alternate genotype frequency (biallelic sites only).">' % (sex, sex))
            if sex == 'XY':
                INFO_ADD.append('##INFO=<ID=nhemiref_%s,Number=1,Type=Integer,Description="Number of %s samples with hemizygous reference genotypes (biallelic sites only).">' % (sex, sex))
                INFO_ADD.append('##INFO=<ID=nhemialt_%s,Number=1,Type=Integer,Description="Number of %s samples with hemizygous alternate genotypes (biallelic sites only).">' % (sex, sex))
                INFO_ADD.append('##INFO=<ID=freq_hemiref_%s,Number=1,Type=Float,Description="%s hemizygous reference genotype frequency (biallelic sites only).">' % (sex, sex))
                INFO_ADD.append('##INFO=<ID=freq_hemialt_%s,Number=1,Type=Float,Description="%s hemizygous alternate genotype frequency (biallelic sites only).">' % (sex, sex))
                if len(parbt) > 0:
                    INFO_ADD.append('##INFO=<ID=par,Number=0,Type=Flag,Description="Variant overlaps pseudoautosomal region.">')
    
    # Add pop fields to header
    if len(pops) > 0:
        INFO_ADD.append('##INFO=<ID=AN_grpmax,Number=1,Type=Integer,Description="Allele number for the grpmax population (biallelic sites only.">')
        INFO_ADD.append('##INFO=<ID=AC_grpmax,Number=1,Type=Integer,Description="Allele count for the grpmax population (biallelic sites only.">')
        INFO_ADD.append('##INFO=<ID=AF_grpmax,Number=1,Type=Float,Description="Maximum allele frequency across any population (biallelic sites only.">')
        INFO_ADD.append('##INFO=<ID=grpmax,Number=1,Type=String,Description="Population label with maximum allele frequency (biallelic sites only).">')
        INFO_ADD.append('##INFO=<ID=nhomalt_grpmax,Number=1,Type=Integer,Description="Number of homozygous-alternate genotypes in the grpmax population (biallelic sites only).">')
        for pop in pops:
            INFO_ADD.append('##INFO=<ID=AN_%s,Number=1,Type=Integer,Description="Total number of %s alleles genotyped.">' % (pop, pop))
            INFO_ADD.append('##INFO=<ID=AC_%s,Number=A,Type=Integer,Description="Number of non-reference %s alleles observed.">' % (pop, pop))
            INFO_ADD.append('##INFO=<ID=AF_%s,Number=A,Type=Float,Description="%s allele frequency.">' % (pop, pop))
            INFO_ADD.append('##INFO=<ID=NCR_%s,Number=1,Type=Float,Description="%s no-call rate.">' % (pop, pop))
            INFO_ADD.append('##INFO=<ID=n_bi_genos_%s,Number=1,Type=Integer,Description="Total number of %s samples with complete genotypes (biallelic sites only).">' % (pop, pop))
            INFO_ADD.append('##INFO=<ID=nhomref_%s,Number=1,Type=Integer,Description="Number of %s samples with homozygous reference genotypes (biallelic sites only).">' % (pop, pop))
            INFO_ADD.append('##INFO=<ID=nhet_%s,Number=1,Type=Integer,Description="Number of %s samples with heterozygous genotypes (biallelic sites only).">' % (pop, pop))
            INFO_ADD.append('##INFO=<ID=nhomalt_%s,Number=1,Type=Integer,Description="Number of %s samples with homozygous alternate genotypes (biallelic sites only).">' % (pop, pop))
            INFO_ADD.append('##INFO=<ID=freq_homref_%s,Number=1,Type=Float,Description="%s homozygous reference genotype frequency (biallelic sites only).">' % (pop, pop))
            INFO_ADD.append('##INFO=<ID=freq_het_%s,Number=1,Type=Float,Description="%s heterozygous genotype frequency (biallelic sites only).">' % (pop, pop))
            INFO_ADD.append('##INFO=<ID=freq_homalt_%s,Number=1,Type=Float,Description="%s homozygous alternate genotype frequency (biallelic sites only).">' % (pop, pop))
            if len(sexes) > 0 and not args.no_combos:
                for sex in sexes:
                    INFO_ADD.append('##INFO=<ID=AN_%s,Number=1,Type=Integer,Description="Total number of %s alleles genotyped.">' % ('_'.join((pop, sex)), ' '.join((pop, sex))))
                    INFO_ADD.append('##INFO=<ID=AC_%s,Number=A,Type=Integer,Description="Number of non-reference %s alleles observed.">' % ('_'.join((pop, sex)), ' '.join((pop, sex))))
                    INFO_ADD.append('##INFO=<ID=AF_%s,Number=A,Type=Float,Description="%s allele frequency.">' % ('_'.join((pop, sex)), ' '.join((pop, sex))))
                    INFO_ADD.append('##INFO=<ID=NCR_%s,Number=1,Type=Float,Description="%s no-call rate.">' % ('_'.join((pop, sex)), ' '.join((pop, sex))))
                    INFO_ADD.append('##INFO=<ID=n_bi_genos_%s,Number=1,Type=Integer,Description="Total number of %s samples with complete genotypes (biallelic sites only).">' % ('_'.join((pop, sex)), ' '.join((pop, sex))))
                    INFO_ADD.append('##INFO=<ID=nhomref_%s,Number=1,Type=Integer,Description="Number of %s samples with homozygous reference genotypes (biallelic sites only).">' % ('_'.join((pop, sex)), ' '.join((pop, sex))))
                    INFO_ADD.append('##INFO=<ID=nhet_%s,Number=1,Type=Integer,Description="Number of %s samples with heterozygous genotypes (biallelic sites only).">' % ('_'.join((pop, sex)), ' '.join((pop, sex))))
                    INFO_ADD.append('##INFO=<ID=nhomalt_%s,Number=1,Type=Integer,Description="Number of %s samples with homozygous alternate genotypes (biallelic sites only).">' % ('_'.join((pop, sex)), ' '.join((pop, sex))))
                    INFO_ADD.append('##INFO=<ID=freq_homref_%s,Number=1,Type=Float,Description="%s homozygous reference genotype frequency (biallelic sites only).">' % ('_'.join((pop, sex)), ' '.join((pop, sex))))
                    INFO_ADD.append('##INFO=<ID=freq_het_%s,Number=1,Type=Float,Description="%s heterozygous genotype frequency (biallelic sites only).">' % ('_'.join((pop, sex)), ' '.join((pop, sex))))
                    INFO_ADD.append('##INFO=<ID=freq_homalt_%s,Number=1,Type=Float,Description="%s homozygous alternate genotype frequency (biallelic sites only).">' % ('_'.join((pop, sex)), ' '.join((pop, sex))))
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
                                     pops, sex_chroms, lps_dict, args.no_combos)
        fout.write(newrec)

    fout.close()


if __name__ == '__main__':
    main()
