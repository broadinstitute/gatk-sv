#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2022 The Broad Institute of M.I.T. and Harvard
# Distributed under terms of the MIT license.
# Contact: Ryan Collins <rlcollins@g.harvard.edu>

"""
Integrate VCFs filtered by various different models
"""


import numpy as np
import csv
from itertools import combinations
from collections import defaultdict
from copy import deepcopy
import argparse
import pysam
import sys
from time import time
from datetime import timedelta


def load_filtered_vcfs(info_tsv, homref_exclude=[]):
    """
    Load connections to filtered VCFs
    """

    filtered_vcfs = {}

    with open(info_tsv) as fin:
        for prefix, vcf in csv.reader(fin, delimiter='\t'):
            prefix = prefix.lower()
            filtered_vcfs[prefix] = pysam.VariantFile(vcf)

    filter_models = list(filtered_vcfs.keys())
    homref_filter_models = [m for m in filter_models if m not in homref_exclude]

    return filtered_vcfs, filter_models, homref_filter_models


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


def tokenize_numeric(x):
    """
    Translates log10-scaled numeric values (SVLEN, AF, etc.) into string keys for dicts
    """

    return str(int(np.floor(x)))


def tokenize_EV(EV):
    """
    Translates EV values into a sorted token for keying into dicts
    """

    if isinstance(EV, str):
        token = ','.join(sorted(EV.split(',')))
    elif isinstance(EV, tuple):
        token = ','.join(sorted([str(x) for x in EV]))

    return token


def calc_af(record):
    """
    Computes a simple AF for a single record
    """

    an = 0
    ac = 0

    for sample, sgt in record.samples.items():
        alleles = sgt['GT']
        an += len([a for a in alleles if a is not None])
        ac += len([a for a in alleles if a is not None and str(a) != "0"])

    if an > 0:
        af = ac / an
    else:
        af = np.NaN

    return af


def load_rules(rules_tsv, filter_models, strict=False):
    """
    Loads rules .tsv as a nested series of defaultdicts, keyed on:
    1. SVTYPE
    2. SVLEN, keyed by tokenize_numeric(SVLEN)
    3. AF before any filtering, keyed by tokenize_numeric(AF)
    4. EV, keyed by tokenize_EV(EV)
    """

    # Set baseline: keep any model if no information is provided unless --strict option specified
    if strict:
        default_pass = [set(filter_models)]
    else:
        default_pass = [set(x) for x in combinations(filter_models, 1)] + \
                       [set(x) for x in combinations(filter_models, 2)] + \
                       [set(x) for x in combinations(filter_models, 3)]
    rules = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: default_pass))))

    if rules_tsv is not None:
        with open(rules_tsv) as tsvin:
            reader = csv.reader(tsvin, delimiter='\t')
            for svtype, min_log10_SVLEN, min_log10_AF, EV, combos_str in reader:
                if svtype.startswith('#'):
                    continue
                combos = [set(cstr.split(';')) for cstr in combos_str.split('|')]
                if svtype not in rules.keys():
                    rules[svtype] = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: default_pass)))
                svlen = tokenize_numeric(float(min_log10_SVLEN))
                if svlen not in rules[svtype].keys():
                    rules[svtype][svlen] = defaultdict(lambda: defaultdict(lambda: default_pass))
                svaf = tokenize_numeric(float(min_log10_AF))
                if svaf not in rules[svtype][svlen].keys():
                    rules[svtype][svlen][svaf] = defaultdict(lambda: default_pass)
                ev = tokenize_EV(EV)
                rules[svtype][svlen][svaf].update({ev : combos})

    return rules


def modify_homref_rules(rules, homref_filter_models=[], strict=False):
    """
    Modifies a rules dictionary (see load_rules()) for integrating homozygous ref GTs
    Don't consider models specified in homref_filter_models for homref GTs (e.g., because Boost only scores non-ref GTs)
    """

    homref_rules = deepcopy(rules)
    
    # Update all defaults
    if strict:
        default_pass = [set(homref_filter_models)]
    else:
        default_pass = [set(x) for x in combinations(homref_filter_models, 1)] + \
                       [set(x) for x in combinations(homref_filter_models, 2)]
    homref_rules.default_factory = lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: default_pass)))

    # Update all rules
    for svtype in rules.keys():
        homref_rules[svtype].default_factory = lambda: defaultdict(lambda: defaultdict(lambda: default_pass))
        for svlen in rules[svtype].keys():
            homref_rules[svtype][svlen].default_factory = lambda: defaultdict(lambda: default_pass)
            for svaf in rules[svtype][svlen].keys():
                homref_rules[svtype][svlen][svaf].default_factory = lambda: default_pass
                for svev in rules[svtype][svlen][svaf].keys():
                    old = rules[svtype][svlen][svaf][svev]
                    new = [x - set(homref_filter_models) for x in old]
                    homref_rules[svtype][svlen][svaf][svev] = \
                        [x for x in new if len(x) > 0]

    return homref_rules


def unify_records(record, record_matches, rules, homref_rules, 
                  homref_filter_models=[], bs_key=None, gqr_key=None):
    """
    Unify information across all versions of the same variant record
    """

    # Get descriptive information about record
    svtype = record.info['SVTYPE']
    svsize = tokenize_numeric(np.log10(np.nanmax([1, record.info.get('SVLEN', 1)])))
    multiallelic = is_mcnv(record)
    # Compute simple AF on the fly if AF is not present in record INFO
    if not multiallelic:
        if 'AF' in record.info.keys():
            af = record.info['AF']
        else:
            af = calc_af(record)
        if af > 0:
            svaf = tokenize_numeric(np.log10(af))
        else:
            svaf = 'AF0'
    nocalls = 0
    nonref = 0

    # Process each sample in serial
    for sample in record.samples.keys():

        # Add BS from boost
        if bs_key is not None:
            if record_matches[bs_key] is not None:
                BS = record_matches[bs_key].samples[sample].get('BS', None)
            else:
                BS = None
            record.samples[sample]['BS'] = BS

        # Add GQR and SL from gqrecal, if provided
        if gqr_key is not None:
            if record_matches[gqr_key] is not None:
                GQR = record_matches[gqr_key].samples[sample].get('GQ', None)
                record.samples[sample]['GQR'] = GQR
                SL = record_matches[gqr_key].samples[sample].get('SL', None)
            else:
                SL = None
            record.samples[sample]['SL'] = SL

        # Do not modify GT for multiallelic variants
        if multiallelic:
            continue

        # Get list of filtering methods where this genotype passes
        pass_filts = set()
        fail_filts = set()
        for prefix, matching_record in record_matches.items():
            if matching_record is not None:
                GT = matching_record.samples[sample]['GT']

                # Don't consider certain models for homref GTs (e.g., Boost only scores non-ref GTs)
                if GT == (0, 0) and prefix not in homref_filter_models:
                    continue

                # If any GT other than ./. is observed, consider this GT passing
                if GT == (None, None):
                    fail_filts.add(prefix)
                else:
                    pass_filts.add(prefix)

        # Get list of filter combinations eligible for this GT
        EV = tokenize_EV(record.samples[sample]['EV'])
        if record.samples[sample]['GT'] == (0, 0):
            relevant_rules = homref_rules
        else:
            relevant_rules = rules
        elig_combos = relevant_rules[svtype][svsize][svaf][EV]

        # Check if GT passes any combination of eligible filters
        passing = any([len(combo.intersection(pass_filts)) == len(combo) for combo in elig_combos])

        # If not passing, set GT to ./.
        if not passing:
            record.samples[sample]['GT'] = (None, None)
            nocalls += 1
        else:
            if record.samples[sample]['GT'] != (None, None) \
            and record.samples[sample]['GT'] != (0, 0):
                nonref += 1

        # Update sample FT field based on failed filters
        if len(fail_filts) > 0:
            record.samples[sample]['FT'] = tuple(fail_filts)
        else:
            record.samples[sample]['FT'] = tuple(['PASS'])

    # Annotate no-call rate (only for biallelic variants)
    if not multiallelic:
        record.info['NCR'] = nocalls / len(record.samples.keys())
    else:
        nonref = len(record.samples.keys())

    return record, nonref


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--unfiltered-vcf', help='Original vcf prior to filtering',
                        required=True)
    parser.add_argument('--filtered-vcfs', help='Two-column .tsv listing all ' +
                        'filtered VCF prefixes and paths', required=True)
    parser.add_argument('--rules', help='.tsv of rules for integrating genotypes. ' +
                        'If provided, must have the following five columns: SVTYPE, ' +
                        'min_log10_SVLEN, min_log10_AF, EV, passing_filter_combos. ' +
                        'If omitted, will keep all informative GTs over no-calls.')
    parser.add_argument('--strict', default=False, action='store_true',
                        help='Require intersection of all three filters for any SVs ' +
                        'not otherwise specified in --rules. [default: keep any ' +
                        'passing GT]')
    parser.add_argument('--homref-exclude', action='append', help='Specify filtered ' +
                        'VCF prefix(es) to exclude when filtering homozygous ' +
                        'reference GTs.')
    parser.add_argument('--keep-empty-records', default=False, action='store_true',
                        help='Retain records where 100\% of non-reference GTs are ' +
                        'removed during filtering. Default: drop records with AC=0.')
    parser.add_argument('-o','--vcf-out', help='Path to output VCF. Accepts "-" ' +
                        'and "stdout". Default: stdout', default='stdout')
    parser.add_argument('-v', '--verbose', default=False, action='store_true',
                        help='Print logging messages.')
    args = parser.parse_args()

    # Open connections to all input VCFs
    raw_vcf = pysam.VariantFile(args.unfiltered_vcf)
    filtered_vcfs, filter_models, homref_filter_models = \
        load_filtered_vcfs(args.filtered_vcfs, args.homref_exclude)

    # Update VCF header as needed
    bs_prefix_hits = [p for p in filtered_vcfs.keys() if 'boost' in p]
    if len(bs_prefix_hits) > 0:
        bs_prefix = bs_prefix_hits[0]
        raw_vcf.header.add_meta('FORMAT', 
                                items=[('ID', "BS"), ('Number', "1"), ('Type', "Float"), 
                                       ('Description', "lgBoost score")])
    else:
        bs_prefix = None
    gqr_prefix_hits = [p for p in filtered_vcfs.keys() if 'gqrecal' in p]
    if len(gqr_prefix_hits) > 0:
        gqr_prefix = gqr_prefix_hits[0]
        raw_vcf.header.add_meta('FORMAT', 
                                items=[('ID', "GQR"), ('Number', "1"), ('Type', "Integer"), 
                                       ('Description', "Recalibrated genotype " + \
                                        "quality from GATK GQRecalibrator")])
        raw_vcf.header.add_meta('FORMAT', 
                                items=[('ID', "SL"), ('Number', "1"), ('Type', "Integer"), 
                                       ('Description', "250 times the logits that " + \
                                                       "the genotype is correct")])
    else:
        gqr_prefix = None
    raw_vcf.header.add_meta('INFO',
                            items=[('ID', "NCR"), ('Number', "1"), ('Type', "Float"),
                                   ('Description', "Proportion of no-call GTs")])
    raw_vcf.header.add_meta('FORMAT',
                            items=[('ID', "FT"), ('Number', "."), ('Type', "String"),
                                   ('Description', "GT filters")])

    # Load integration rules
    rules = load_rules(args.rules, filter_models, args.strict)
    homref_rules = modify_homref_rules(rules, homref_filter_models)

    # Open connection to output VCF
    if args.vcf_out in '- stdout'.split():
        outvcf = pysam.VariantFile(sys.stdout, 'w', header=raw_vcf.header)
    else:
        outvcf = pysam.VariantFile(args.vcf_out, 'w', header=raw_vcf.header)

    # Load first records from all filtered VCFs into memory
    next_records = {k : v.__next__() for k, v in filtered_vcfs.items()}

    # Iterate over records
    k = 0
    k_out = 0
    seen_vids = set()
    if args.verbose:
        start_time = time()
    for record in raw_vcf:
        vid = record.id

        # Check to ensure no duplicate variant IDs
        if vid in seen_vids:
            raise Exception('Encountered record ID in unfiltered VCF that ' + \
                            'was already processed: ' + vid + '\n')

        # Check for matching records
        record_matches = {}
        for prefix, vcf in filtered_vcfs.items():
            record_matches[prefix] = None
            if next_records[prefix] is not None:
                
                # Check to ensure no duplicate variant IDs
                if next_records[prefix].id in seen_vids:
                    raise Exception('Encountered record ID in ' + prefix + \
                                    ' VCF that was already processed: ' + \
                                    next_records[prefix].id + '\n')

                if next_records[prefix].id == vid:
                    record_matches[prefix] = next_records[prefix].copy()
                    try:
                        next_records[prefix] = filtered_vcfs[prefix].__next__()
                    except:
                        next_records[prefix] = None

        # Unify information from all versions of record
        record, n_nonref = \
            unify_records(record, record_matches, rules, homref_rules, 
                          homref_filter_models, bs_prefix, gqr_prefix)

        # Update list of seen VIDs
        seen_vids.add(vid)

        # Write updated record to output VCF if at least one non-reference GT was
        # retained (or unless overridden by --keep-empty-records)
        if n_nonref > 0 or args.keep_empty_records:
            outvcf.write(record)
            k_out += 1

        # Write logging message if --verbose is specified
        k += 1
        if k % 10 == 0 and args.verbose:
            msg = 'Processed {:,} records ({:,} written to VCF). Elapsed time: {} ({} seconds per record).'
            elapsed = time() - start_time
            print(msg.format(k, k_out, timedelta(seconds=np.round(elapsed)),
                             np.round(elapsed / k, 3)))

    # Close connection to output VCF
    outvcf.close()


if __name__ == '__main__':
    main()

