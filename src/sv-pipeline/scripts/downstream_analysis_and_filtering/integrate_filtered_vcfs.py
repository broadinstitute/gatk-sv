#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2022 The Broad Institute of M.I.T. and Harvard
# Distributed under terms of the MIT license.
# Contact: Ryan Collins <rlcollins@g.harvard.edu>

"""
Integrate VCFs filtered by various different models
"""


from collections import defaultdict
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


def categorize_svlen(svlen):
    """
    Translates SVLEN values into string categories for simplicity
    """

    if svlen <= 250:
        return 'tiny'
    elif svlen > 250 and svlen <= 1000:
        return 'small'
    elif svlen > 1000 and svlen <= 5000:
        return 'medium'
    elif svlen > 5000:
        return 'large'


def tokenize_evidence(EV):
    """
    Convert EV tuple into a string token
    """

    return ','.join(sorted(EV))


def load_rules(rules_tsv):
    """
    Loads rules .tsv as a nested series of defaultdicts, keyed on:
    1. SVTYPE
    2. SVLEN
    3. EVIDENCE
    """

    # Set baseline: keep any model if no information is provided
    any_model = [set(['mingq']), set(['boost']), set(['gqrecal'])]
    rules = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: any_model)))

    if rules_tsv is not None:
        # TODO: update rules dict from .tsv
        import pdb; pdb.set_trace()

    return rules


def unify_records(record, mingq_r, boost_r, gqrecal_r, rules):
    """
    Unify information across all versions of the same variant record
    """

    # Get descriptive information about record
    svtype = record.info['SVTYPE']
    svsize = categorize_svlen(record.info.get('SVLEN', -1))
    multiallelic = is_mcnv(record)
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
        EV = tokenize_evidence(record.samples[sample]['EV'])
        elig_combos = rules[svtype][svsize][EV]

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
                        'If omitted, will keep all informative GTs over no-calls.')
    parser.add_argument('--keep-empty-records', default=False, action='store_true',
                        help='Retain records where 100\% of non-reference GTs are ' +
                        'removed during filtering. Default: drop records with AC=0.')
    parser.add_argument('-o','--vcf-out', help='Path to output VCF. Accepts "-" ' +
                        'and "stdout". Default: stdout', default='stdout')
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
    rules = load_rules(args.rules)

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
                                         gqrecal_record, rules)

        # Write updated record to output VCF if at least one non-reference GT was
        # retained (or unless overridden by --keep-empty-records)
        if n_nonref > 0 or args.keep_empty_records:
            outvcf.write(record)

    # Close connection to output VCF
    outvcf.close()

if __name__ == '__main__':
    main()

