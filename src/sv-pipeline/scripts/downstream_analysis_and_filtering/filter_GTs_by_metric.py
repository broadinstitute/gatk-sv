#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2022 Ryan Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

"""
Apply per-sample filters to one or more metrics in FORMAT field in an input VCF
"""


import argparse
import sys
import pandas as pd
from collections import defaultdict
import csv
import numpy as np
import pysam
import svtk.utils as svu


def load_conditions(minMetricTable):
    """
    Loads conditions as a pandas DataFrame and a dict of FORMAT : value pairs keyed by condition ID
    """

    # Load table of conditions
    conditions_table = pd.read_csv(minMetricTable, sep='\t')
    conditions_table.rename(columns={'#condition' : 'condition'}, inplace=True)

    # Load dictionary of FORMAT : value pairs
    conditions_dict = {row.condition : {row.metric : row.min_metric} \
                       for row in conditions_table.itertuples()}

    # Reformat table of conditions for quick lookups
    conditions_table.set_index(keys='condition', drop=True, inplace=True)
    set_cols = 'includeSVTYPE excludeSVTYPE includeFILTER excludeFILTER includeEV excludeEV'.split()
    for col in set_cols:
        conditions_table.loc[:, col] = \
            conditions_table.loc[:, col].apply(lambda x: set(x.split(',')))
    conditions_table.drop(columns='metric min_metric source'.split(), inplace=True)

    return {'table' : conditions_table, 'dict' : conditions_dict}


def get_variant_conditions(record, conditions):
    """
    Find conditions that apply to a single VCF record
    """

    # Get variant info
    SVTYPE = set([record.info.get('SVTYPE')])
    SVLEN = int(record.info.get('SVLEN', 0))
    AF = float(record.info.get('AF', [None])[0])
    FILTER = set(record.filter.keys())
    n_FILTER = len(FILTER)
    
    # Find conditions that apply to the record
    SVTYPE_hits = np.logical_and(conditions['table'].includeSVTYPE.apply(lambda x: len(x.intersection(SVTYPE)) > 0), 
                                 conditions['table'].excludeSVTYPE.apply(lambda x: len(x.intersection(SVTYPE)) == 0 or x == {'None'}))
    SVLEN_hits = np.logical_and(SVLEN >= conditions['table'].minSVLEN, SVLEN <= conditions['table'].maxSVLEN)
    AF_hits = np.logical_and(AF >= conditions['table'].minAF, AF <= conditions['table'].maxAF)
    FILTER_hits = np.logical_and(conditions['table'].includeFILTER.apply(lambda x: len(x.intersection(FILTER)) == n_FILTER or x == {'None'}),
                                 conditions['table'].excludeFILTER.apply(lambda x: len(x.intersection(FILTER)) == 0 or x == {'None'}))
    all_hits = (SVTYPE_hits & SVLEN_hits & AF_hits & FILTER_hits)
    cond_ids = set(list(conditions['table'].index[all_hits]))
    
    return cond_ids


def get_sample_conditions(record, sample, conditions):
    """
    Find conditions that apply to a single sample's GT for a given VCF record
    """

    # Get variant info
    EV = set(record.samples[sample].get('EV', []))
    n_EV = len(EV)
    
    # Find conditions that apply to the record
    EV_hits = np.logical_and(conditions['table'].includeEV.apply(lambda x: len(x.intersection(EV)) == n_EV or x == {'None'}),
                             conditions['table'].excludeEV.apply(lambda x: len(x.intersection(EV)) == 0 or x == {'None'}))
    all_hits = EV_hits
    cond_ids = set(list(conditions['table'].index[all_hits]))
    
    return cond_ids


def apply_conditions(record, conditions, filter_homref, filter_homalt, 
                     fail_missing_scores, require_all_criteria, annotate_ncr, 
                     max_ncr, ncr_prefix="COHORT"):
    """
    Apply GT filtering conditions to a single record
    """

    # Get conditions that apply to all samples for this variant
    variant_conditions = get_variant_conditions(record, conditions)

    n_samples = len(record.samples)

    bl = 0
    if len(variant_conditions) > 0:
        for s in record.samples:
            # Handle certain genotypes differently based on inputs
            samp_gt = record.samples[s]['GT']

            # Don't need to process GTs already set to ./., but should be counted towards NCR
            if samp_gt == (None, None):
                bl += 1
                continue

            # Only filter homref GTs if optioned
            if samp_gt == (0, 0) and not filter_homref:
                continue

            # Only filter homalt GTs if optioned
            if samp_gt == (1, 1) and not filter_homalt:
                continue

            # Get conditions to apply to this particular sample
            sample_conditions = get_sample_conditions(record, s, conditions)
            relevant_conditions = \
                [conditions['dict'][cid] for cid in \
                 variant_conditions.intersection(sample_conditions)]


            # Check metric cutoff(s) versus each condition
            samp_pass = []
            for cond in relevant_conditions:
                metric, cutoff = list(cond.items())[0]
                samp_val = record.samples[s].get(metric)
                if samp_val is None:
                    if fail_missing_scores:
                        samp_pass.append(False)
                    else:
                        samp_pass.append(True)
                elif samp_val < cutoff:
                    samp_pass.append(False)
                else:
                    samp_pass.append(True)

            # Rewrite sample GT based on filter pass/fail
            if require_all_criteria:
                gt_pass = all(samp_pass)
            else:
                if len(samp_pass) > 0:
                    gt_pass = any(samp_pass)
                else:
                    gt_pass = True
            if not gt_pass:
                record.samples[s]['GT'] = (None, None)
                bl += 1

    if annotate_ncr:
        if n_samples > 0:
            ncr = bl / n_samples
            record.info['{}_NCR'.format(ncr_prefix)] = ncr
            if ncr > max_ncr:
                record.filter.add('HIGH_{0}_NOCALL_RATE'.format(ncr_prefix))


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
    parser.add_argument('minMetricTable', help='Tab-delimited filtering lookup ' + 
                        'table in the same format as generated by ' +
                        'create_minGQ_lookup_table.R.')
    parser.add_argument('fout', help='Output file (supports "stdout").')
    parser.add_argument('--multiallelics', default=False, action='store_true',
                        help='Also apply filtering to multiallelic sites ' + 
                        '(default: do not filter multiallelics).')
    parser.add_argument('--filter-homref', default=False, action='store_true',
                        help='Apply filters to homozygous reference GTs [default: ' +
                        'do not filter any homref GTs]')
    parser.add_argument('--filter-homalt', default=False, action='store_true',
                        help='Apply filters to homozygous alternate GTs [default: ' +
                        'do not filter any homalt GTs]')
    parser.add_argument('--fail-missing-scores', default=False, action='store_true',
                        help='Treat GTs with missing scores as failures [default: ' +
                        'do not filter GTs with missing scores]')
    parser.add_argument('--dropEmpties', default=False, action='store_true',
                        help='After GT reassignments, drop any SV with no remaining ' + 
                        ' non-ref samples (default: keep all SV).')
    parser.add_argument('--require-all-criteria', default=False, action='store_true',
                        help='Require that each GT must pass all criteria When ' +
                        'multiple criteria are specified in minMetricTable for the ' +
                        'same variant. [default: require only one passing metric]')
    parser.add_argument('--simplify-INS-SVTYPEs', default=False, action='store_true',
                        help='Resets the SVTYPE of all INS variants, including MEIs, ' + 
                        'to be SVTYPE=INS (default: keep original SVTYPEs).')
    parser.add_argument('--annotate-ncr', default=False, action='store_true',
                        help='Annotate output VCF with no-call rates.')
    parser.add_argument('--max-ncr', help='Max no-call rate among all ' + 
                        'samples before adding a flag to the record\'s FILTER field' + 
                        ' (default: 0.05)', type=float, default=0.05)
    parser.add_argument('--cleanAFinfo', help='Remove all AF-related terms from ' + 
                        ' the INFO field and VCF header (default: keep all terms).', 
                        default=False, action='store_true')
    parser.add_argument('--prefix', help='Cohort label to append to NCR FILTER.', 
                        default='COHORT', dest='prefix')

    args = parser.parse_args()

    if args.vcf in '- stdin'.split():
        vcf = pysam.VariantFile(sys.stdin) 
    else:
        vcf = pysam.VariantFile(args.vcf)

    #Add HIGH_NOCALL_RATE filter to vcf header
    if args.annotate_ncr:
        if '{}_NCR'.format(args.prefix) not in vcf.header.info.keys():
            NEW_INFO = '##INFO=<ID={0}_NCR,Number=1,Type=Float,Description="' + \
                       'Fraction of {0} sample GTs that are no-calls.">'
            vcf.header.add_line(NEW_INFO.format(args.prefix))
        NEW_FILTER = '##FILTER=<ID=HIGH_{0}_NOCALL_RATE,Description="More than '.format(args.prefix) + \
                     '{:.2%}'.format(args.max_ncr) + ' of {0} sample GTs were '.format(args.prefix) + \
                     'masked as no-call GTs. Indicates a possibly noisy locus ' + \
                     'in {0} samples.>'.format(args.prefix)
        vcf.header.add_line(NEW_FILTER)

    # Check to ensure either AF or AC & AN are provided in input VCF header
    if 'AF' in vcf.header.info.keys():
        pass
    elif ('AC' in vcf.header.info.keys() and 'AN' in vcf.header.info.keys()):
        vcf.header.add_line('##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">')
    else:
        exit('Input VCF must have either AF or AC & AN defined.')

    if args.fout in '- stdout'.split():
        fout = pysam.VariantFile(sys.stdout, 'w', header=vcf.header)
    else:
        fout = pysam.VariantFile(args.fout, 'w', header=vcf.header)

    # Load filter table as pd.DataFrame
    conditions = load_conditions(args.minMetricTable)

    #Iterate over records in vcf and apply filter
    for record in vcf.fetch():

        #Do not process multiallelic variants, unless optioned
        if args.multiallelics or \
        (not args.multiallelics and 
         not _is_multiallelic(record)):

            # Infer record's AF if AC & AN (but not AF) are provided
            if 'AF' in record.info.keys():
                pass
            else:
                AF = tuple([ac / record.info['AN'] for ac in record.info['AC']])
                record.info['AF'] = AF
            apply_conditions(record, conditions, args.filter_homref, 
                             args.filter_homalt, args.fail_missing_scores, 
                             args.require_all_criteria, args.annotate_ncr, 
                             args.max_ncr, ncr_prefix=args.prefix)

        if args.cleanAFinfo:
            # Clean biallelic AF annotation
            for key in 'AN AC AF N_BI_GENOS N_HOMREF N_HET N_HOMALT FREQ_HOMREF FREQ_HET FREQ_HOMALT'.split(' '):
                if key in record.info.keys():
                    record.info.pop(key)
            # Clean CN frequency annotation
            for key in 'CN_NUMBER CN_COUNT CN_FREQ CN_NONREF_COUNT CN_NONREF_FREQ'.split():
                if key in record.info.keys():
                    record.info.pop(key)

        # Standardize SVTYPE for all INS variants, if optioned
        if any([keyword in record.info['SVTYPE'].split(':') for keyword in 'INS MEI'.split()]):
            if args.simplify_INS_SVTYPEs:
                record.info['SVTYPE'] = 'INS'

        if args.dropEmpties:
            samps = svu.get_called_samples(record, include_null=False)
            if len(samps) > 0:
                fout.write(record)
        else: 
            fout.write(record)

    fout.close()


if __name__ == '__main__':
    main()
