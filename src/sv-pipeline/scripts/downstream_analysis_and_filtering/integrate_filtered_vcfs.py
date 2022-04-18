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


filter_models = ['mingq', 'boost', 'gqrecalibrator']
filter_models_noBoost = [x for x in filter_models if x != 'boost']


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


def tokenize_numeric(svlen):
    """
    Translates log10-scaled numeric values (SVLEN, AF, etc.) into string keys for dicts
    """

    return str(int(np.floor(svlen)))


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


def load_rules(rules_tsv, strict=False):
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


def modify_homref_rules(rules, strict=False):
    """
    Modifies a rules dictionary (see load_rules()) for integrating homozygous ref GTs
    Don't consider Boost for homref GTs because Boost only scores non-ref GTs
    """

    homref_rules = deepcopy(rules)
    
    # Update all defaults
    if strict:
        default_pass = [set(filter_models_noBoost)]
    else:
        default_pass = [set(x) for x in combinations(filter_models_noBoost, 1)] + \
                       [set(x) for x in combinations(filter_models_noBoost, 2)]
    homref_rules.default_factory = lambda: default_pass

    # Update all rules
    for svtype in rules.keys():
        homref_rules[svtype].default_factory = lambda: default_pass
        for svlen in rules[svtype].keys():
            homref_rules[svtype][svlen].default_factory = lambda: default_pass
            for svaf in rules[svtype][svlen].keys():
                homref_rules[svtype][svlen][svaf].default_factory = lambda: default_pass
                for svev in rules[svtype][svlen][svaf].keys():
                    old = rules[svtype][svlen][svaf][svev]
                    new = [x - set(['boost']) for x in old]
                    homref_rules[svtype][svlen][svaf][svev] = \
                        [x for x in new if len(x) > 0]

    return homref_rules


def unify_records(record, mingq_r, boost_r, gqrecal_r, rules, homref_rules):
    """
    Unify information across all versions of the same variant record
    """

    # Get descriptive information about record
    svtype = record.info['SVTYPE']
    svsize = tokenize_numeric(np.log10(record.info.get('SVLEN', -1)))
    multiallelic = is_mcnv(record)
    # Compute simple AF on the fly if AF is not present in record INFO
    if not multiallelic:
        if 'AF' in record.info.keys():
            af = record.info['AF']
        else:
            af = calc_af(record)
        svaf = tokenize_numeric(np.log10(af))
    nocalls = 0
    nonref = 0
    
    # Process each sample in serial
    for sample in record.samples.keys():

        # Add BS from boost
        if boost_r is not None:
            BS = boost_r.samples[sample].get('BS', None)
        else:
            BS = None
        record.samples[sample]['BS'] = BS

        # Rewrite GQ from gqrecal
        if gqrecal_r is not None:
            GQ = gqrecal_r.samples[sample].get('GQ', None)
            record.samples[sample]['GQ'] = GQ

        # Add SL from gqrecal
        if gqrecal_r is not None:
            SL = gqrecal_r.samples[sample].get('SL', None)
        else:
            SL = None
        record.samples[sample]['SL'] = SL

        # Do not modify GT for multiallelic variants
        if multiallelic:
            continue

        # Get list of filtering methods where this genotype passes
        pass_filts = set()
        for rec, tag in [(mingq_r, 'mingq'), (boost_r, 'boost'), (gqrecal_r, 'gqrecal')]:
            if rec is not None:
                GT = rec.samples[sample]['GT']

                # Don't consider Boost for homref GTs because Boost only scores non-ref GTs
                if tag == 'boost' and GT == (0, 0):
                    continue

                # If any GT other than ./. is observed, consider this GT passing
                if GT != (None, None):
                    pass_filts.add(tag)

        # Get list of filter combinations eligible for this GT
        # Don't consider Boost for homref GTs because Boost only scores non-ref GTs
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
    parser.add_argument('--mingq-vcf', help='VCF filtered by minGQ', required=True)
    parser.add_argument('--boost-vcf', help='VCF filtered by Boost', required=True)
    parser.add_argument('--gqrecalibrator-vcf', help='VCF filtered by GATK ' +
                        'GQRecalibrator', required=True)
    parser.add_argument('--rules', help='.tsv of rules for integrating genotypes. ' +
                        'If provided, must have the following five columns: SVTYPE, ' +
                        'min_log10_SVLEN, min_log10_AF, EV, passing_filter_combos. ' +
                        'If omitted, will keep all informative GTs over no-calls.')
    parser.add_argument('--strict', default=False, action='store_true',
                        help='Require intersection of all three filters for any SVs ' +
                        'not otherwise specified in --rules. [default: keep any ' +
                        'passing GT]')
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
    mingq_vcf = pysam.VariantFile(args.mingq_vcf)
    boost_vcf = pysam.VariantFile(args.boost_vcf)
    gqrecal_vcf = pysam.VariantFile(args.gqrecalibrator_vcf)

    # Update VCF header as needed
    raw_vcf.header.add_meta('FORMAT', 
                            items=[('ID', "BS"), ('Number', "1"), ('Type', "Float"), 
                                   ('Description', "lgBoost score")])
    raw_vcf.header.add_meta('FORMAT', 
                            items=[('ID', "SL"), ('Number', "1"), ('Type', "Integer"), 
                                   ('Description', "250 times the logits that " + \
                                                   "the genotype is correct")])
    raw_vcf.header.add_meta('INFO',
                            items=[('ID', "NCR"), ('Number', "1"), ('Type', "Float"),
                                   ('Description', "Proportion of no-call GTs")])

    # Load integration rules
    rules = load_rules(args.rules, args.strict)
    homref_rules = modify_homref_rules(rules)

    # Open connection to output VCF
    if args.vcf_out in '- stdout'.split():
        outvcf = pysam.VariantFile(sys.stdout, 'w', header=raw_vcf.header)
    else:
        outvcf = pysam.VariantFile(args.vcf_out, 'w', header=raw_vcf.header)

    # Load first records from all filtered VCFs into memory
    next_mingq_record = mingq_vcf.__next__()
    next_boost_record = boost_vcf.__next__()
    next_gqrecal_record = gqrecal_vcf.__next__()

    # Iterate over records
    k=0
    k_out=0
    if args.verbose:
        start_time = time()
    for record in raw_vcf:
        vid = record.id

        # Find matching records
        if next_mingq_record.id == vid:
            mingq_record = next_mingq_record.copy()
            next_mingq_record = mingq_vcf.__next__()
        else:
            mingq_record = None
        if next_boost_record.id == vid:
            boost_record = next_boost_record.copy()
            next_boost_record = boost_vcf.__next__()
        else:
            boost_record = None
        if next_gqrecal_record.id == vid:
            gqrecal_record = next_gqrecal_record.copy()
            next_gqrecal_record = gqrecal_vcf.__next__()
        else:
            gqrecal_record = None

        # Unify information from all versions of record
        record, n_nonref = unify_records(record, mingq_record, boost_record, 
                                         gqrecal_record, rules, homref_rules)

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

