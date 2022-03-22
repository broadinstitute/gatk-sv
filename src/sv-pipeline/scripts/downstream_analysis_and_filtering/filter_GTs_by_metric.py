#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Apply per-sample filters on GQ (or another FORMAT field) to an input VCF
"""

import argparse
import csv
import numpy as np
import pysam
import svtk.utils as svu
import sys
from collections import defaultdict


# Make dummy dict for all unique ranges of SVLEN
def _make_SVLEN_interval_dict(minMetricTable):
    # Prep lookup table
    SVLEN_table = {}
    i = 0

    with open(minMetricTable) as lt:
        reader = csv.reader(lt, delimiter='\t')

        for cond_id, minSVLEN, maxSVLEN, minAF, maxAF, includeSVTYPE, excludeSVTYPE, \
                includeFILTER, excludeFILTER, includeEV, excludeEV, minMetric, source in reader:

            # Skip header line
            if "#" in cond_id:
                continue

            # Add new range for minSVLEN to maxSVLEN if it doesn't already exist
            newrange = range(int(minSVLEN), int(maxSVLEN))
            if newrange not in SVLEN_table:
                SVLEN_table[newrange] = str(i)
                i += 1

    # Return SVLEN lookup table
    return SVLEN_table


# Helper function to index into SVLEN dummy table
def _lookup_SVLEN_key(SVLEN, SVLEN_table):
    out = None
    for key in SVLEN_table:
        if SVLEN in key:
            out = SVLEN_table[key]
            break
    return out


# Make dummy dict for all unique ranges of AF
def _make_AF_interval_dict(minMetricTable, scalar=10000):
    # Prep lookup table
    AF_table = {}
    i = 0

    with open(minMetricTable) as lt:
        reader = csv.reader(lt, delimiter='\t')

        for cond_id, minSVLEN, maxSVLEN, minAF, maxAF, includeSVTYPE, excludeSVTYPE, \
                includeFILTER, excludeFILTER, includeEV, excludeEV, minMetric, source in reader:

            # Skip header line
            if "#" in cond_id:
                continue

            # Add new range for minAF to maxAF if it doesn't already exist
            newrange = range(int(np.floor(scalar * float(minAF))), int(np.ceil(scalar * float(maxAF))))
            if newrange not in AF_table:
                AF_table[newrange] = str(i)
                i += 1

    # Return AF lookup table
    return AF_table


# Helper function to index into AF dummy table
def _lookup_AF_key(AF, AF_table, scalar=10000):
    val = int(np.round(scalar * float(AF)))
    if val >= scalar:
        val = scalar - 1
    out = None
    for key in AF_table:
        if val in key:
            out = AF_table[key]
            break
    return out


# Make dummy dict for all SVTYPEs
def _make_SVTYPE_dict(minMetricTable):
    # Prep lookup table
    SVTYPE_table = {}
    i = 0

    # Add all SVTYPE mappings from input minMetric table
    with open(minMetricTable) as lt:
        reader = csv.reader(lt, delimiter='\t')

        for cond_id, minSVLEN, maxSVLEN, minAF, maxAF, includeSVTYPE, excludeSVTYPE, \
                includeFILTER, excludeFILTER, includeEV, excludeEV, minMetric, source in reader:

            # Skip header line
            if "#" in cond_id:
                continue

            # Add new key for each SVTYPE to include specified that is not also excluded
            SVTYPES = [s for s in includeSVTYPE.split(',') if s not in excludeSVTYPE.split(',')]
            missing = len([s for s in SVTYPES if s not in SVTYPE_table])
            if missing > 0:
                for s in SVTYPES:
                    if s not in SVTYPE_table:
                        SVTYPE_table[s] = str(i)
                i += 1

                # Return SVTYPE lookup table
    return SVTYPE_table


# Helper function to index into SVTYPE table
def _lookup_SVTYPE_key(SVTYPE, SVTYPE_table):
    out = SVTYPE_table.get(SVTYPE, None)
    return out


# Make dummy dict for all FILTER combinations
def _make_FILTER_dict(minMetricTable, vcf):
    # Prep lookup table
    FILTER_table = {}
    i = 0

    # Get all FILTER statuses from VCF header
    vcf_filters = vcf.header.filters.keys()

    # Add all FILTER mappings from input minMetric table
    with open(minMetricTable) as lt:
        reader = csv.reader(lt, delimiter='\t')

        for cond_id, minSVLEN, maxSVLEN, minAF, maxAF, includeSVTYPE, excludeSVTYPE, \
                includeFILTER, excludeFILTER, includeEV, excludeEV, minMetric, source in reader:

            # Skip header line
            if "#" in cond_id:
                continue

            # Add new key for each FILTER to include specified that is not also excluded
            # If no FILTERs are explicitly included, default to all FILTERs in VCF header
            if includeFILTER == 'None':
                includeFILTER = (',').join(vcf_filters)
            FILTERS = [f for f in includeFILTER.split(',') if f not in excludeFILTER.split(',')]
            missing = len([f for f in FILTERS if f not in FILTER_table])
            if missing > 0:
                for f in FILTERS:
                    if f not in FILTER_table:
                        FILTER_table[f] = str(i)
                i += 1

    # Return FILTER lookup table
    return FILTER_table


# Helper function to index into FILTER table and default to lowest value
def _lookup_FILTER_key(FILTERs, FILTER_table):
    maxkey = str(np.max(list(map(int, set(FILTER_table.values())))))
    out = [FILTER_table.get(f, maxkey) for f in FILTERs.split(',')]
    if len(out) > 0:
        out = str(np.min(list(map(int, out))))
    return out


# Make dummy dict for all EV combinations
def _make_EV_dict(minMetricTable):
    # Prep lookup table
    EV_table = {}
    i = 0

    # Get list of universe of possible combinations of evidence
    ev_types = ['RD', 'PE', 'SR', 'RD,PE', 'RD,SR', 'PE,SR', 'RD,PE,SR']

    # Add all EV mappings from input minMetric table
    with open(minMetricTable) as lt:
        reader = csv.reader(lt, delimiter='\t')

        for cond_id, minSVLEN, maxSVLEN, minAF, maxAF, includeSVTYPE, excludeSVTYPE, \
                includeFILTER, excludeFILTER, includeEV, excludeEV, minMetric, source in reader:

            # Skip header line
            if "#" in cond_id:
                continue

            # Add new key for each EV to include specified that is not also excluded
            EVS = [e for e in includeEV.split(',') if e not in excludeEV.split(',')]
            missing = len([f for f in EVS if f not in EV_table and f != "None"])
            if missing > 0:
                for e in EVS:
                    if e not in EV_table:
                        EV_table[e] = str(i)
                i += 1

    # Add a final key for all evidence classes not already in dict
    for e in ev_types:
        if e not in EV_table:
            EV_table[e] = str(i)

    # Return EV lookup table
    return EV_table


# Helper function to index into EV table
def _lookup_EV_key(EV, EV_table):
    # Reorder EV string to match expected
    e = []
    for c in ['RD', 'PE', 'SR']:
        if c in EV:
            e.append(c)
    EV = ','.join(e)
    maxkey = str(np.max(list(map(int, set(EV_table.values())))))
    out = EV_table.get(EV, maxkey)
    return out


# Create master minMetric lookup table
def make_minMetric_dict(minMetricTable, SVLEN_table, AF_table, SVTYPE_table,
                        FILTER_table, EV_table):
    # Prep master minMetric lookup table
    def _nested_dict(n, type):
        if n == 1:
            return defaultdict(type)
        else:
            return defaultdict(lambda: _nested_dict(n - 1, type))

    minMetric_dict = _nested_dict(5, str)

    with open(minMetricTable) as lt:
        reader = csv.reader(lt, delimiter='\t')

        # Enter each line in the lookup table into the dictionary
        for cond_id, minSVLEN, maxSVLEN, minAF, maxAF, includeSVTYPE, excludeSVTYPE, \
                includeFILTER, excludeFILTER, includeEV, excludeEV, minMetric, source in reader:

            # Skip header line
            if "#" in cond_id:
                continue

            # Prep variables as needed
            SVLEN = int(np.mean([int(minSVLEN), int(maxSVLEN)]))
            SVLEN_idx = _lookup_SVLEN_key(SVLEN, SVLEN_table)
            AF = np.mean([float(minAF), float(maxAF)])
            AF_idx = _lookup_AF_key(AF, AF_table)
            FILTERs = (',').join([f for f in includeFILTER.split(',')
                                  if f not in excludeFILTER.split(',')])
            FILTER_idx = _lookup_FILTER_key(FILTERs, FILTER_table)
            EV_idx = _lookup_EV_key(includeEV, EV_table)
            minMetric = float(minMetric)

            # Update one line for each qualifying SVTYPE
            for SVTYPE in [s for s in includeSVTYPE.split(',')
                           if s not in excludeSVTYPE.split(',')]:
                SVTYPE_idx = _lookup_SVTYPE_key(SVTYPE, SVTYPE_table)

                # Add line to dict if all indexes are not None
                if SVLEN_idx is not None \
                        and AF_idx is not None \
                        and SVTYPE_idx is not None \
                        and FILTER_idx is not None:
                    minMetric_dict[SVLEN_idx][AF_idx][SVTYPE_idx][FILTER_idx][EV_idx] = int(minMetric)

    return minMetric_dict


# Query master minMetric lookup table
def _get_minMetric(record, minMetric_dict, SVLEN_table, AF_table, SVTYPE_table,
                   FILTER_table, EV_table, globalMin=0):
    # Get minMetric_dict indexes
    if 'SVLEN' in record.info.keys():
        SVLEN_idx = _lookup_SVLEN_key(record.info['SVLEN'], SVLEN_table)
    AF_idx = _lookup_AF_key(record.info['AF'][0], AF_table)
    if 'SVTYPE' in record.info.keys():
        SVTYPE_idx = _lookup_SVTYPE_key(record.info['SVTYPE'], SVTYPE_table)
    FILTER_query = ','.join(f for f in record.filter)
    FILTER_idx = _lookup_FILTER_key(FILTER_query, FILTER_table)

    # Get minMetric dict down to EV
    if SVLEN_idx is not None \
            and AF_idx is not None \
            and SVTYPE_idx is not None \
            and FILTER_idx is not None:
        minMetric = minMetric_dict[SVLEN_idx][AF_idx][SVTYPE_idx][FILTER_idx]
    # Otherwise, set all EV categories to globalMin
    else:
        minMetric = defaultdict(dict)
        for key in list(set(EV_table.values())):
            minMetric[key] = int(globalMin)

    return (minMetric)


def apply_minMetric_filter(record, metric, minMetric_dict, SVLEN_table, AF_table,
                           SVTYPE_table, FILTER_table, EV_table, globalMin=0,
                           maxNCR=0.05, highNCR_filter="COHORT"):
    # Get minMetric dict down to EV for this variant
    minMetric_by_ev = _get_minMetric(record, minMetric_dict, SVLEN_table,
                                     AF_table, SVTYPE_table, FILTER_table,
                                     EV_table, globalMin)

    n_samples = len(record.samples)

    bl = 0
    for s in record.samples:
        # Don't process homozygous genotypes
        if record.samples[s]['GT'] == (1, 1):
            continue

        # Get minMetric cutoff given that sample's evidence
        EV = record.samples[s]['EV']
        if isinstance(EV, tuple):
            EV = ','.join(list(EV))
        EV_key = _lookup_EV_key(EV, EV_table)
        minMetric = minMetric_by_ev.get(EV_key, globalMin)
        if record.samples[s][metric] is not None:
            if record.samples[s][metric] < minMetric:
                record.samples[s]['GT'] = (None, None)
                bl += 1

    if n_samples > 0:
        frac_bl = bl / float(n_samples)
        record.info['minGQ_NCR'] = round(frac_bl, 2)
        if frac_bl > maxNCR:
            record.filter.add(highNCR_filter)


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
    parser.add_argument('--metric', default='GQ', type=str, help='FORMAT field ' +
                                                                 'to which filtering should be applied [default: GQ]')
    parser.add_argument('--globalMin', help='Global minimum metric value',
                        type=int, default=0)
    parser.add_argument('--multiallelics', default=False, action='store_true',
                        help='Also apply filtering to multiallelic sites ' +
                             '(default: do not filter multiallelics).')
    parser.add_argument('--dropEmpties', default=False, action='store_true',
                        help='After GT reassignments, drop any SV with no remaining ' +
                             ' non-ref samples (default: keep all SV).')
    parser.add_argument('--simplify-INS-SVTYPEs', default=False, action='store_true',
                        help='Resets the SVTYPE of all INS variants, including MEIs, ' +
                             'to be SVTYPE=INS (default: keep original SVTYPEs).')
    parser.add_argument('--maxNCR', help='Max no-call rate among all ' +
                                         'samples before adding a flag to the record\'s FILTER field' +
                                         ' (default: 0.05)', type=float, default=0.05, dest='maxNCR')
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

    # Add HIGH_NOCALL_RATE filter and info field to vcf header
    NEW_FILTER = '##FILTER=<ID=HIGH_{0}_NOCALL_RATE,Description="More than '.format(args.prefix) + \
                 '{:.2%}'.format(args.maxNCR) + ' of {0} sample GTs were '.format(args.prefix) + \
                 'masked as no-call GTs due to low {0}. '.format(args.metric) + \
                 'Indicates a possibly noisy locus in {0} samples.>'.format(args.prefix)
    NEW_INFO = '##INFO=<ID=minGQ_NCR,Number=1,Type=Float,Description="Fraction of sample GTs that were ' \
               'masked as no-call GTs due to low GQ during minGQ filtering.>'
    vcf.header.add_line(NEW_FILTER)
    vcf.header.add_line(NEW_INFO)
    filter_text = 'HIGH_{0}_NOCALL_RATE'.format(args.prefix)

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

    # Make dummy lookup tables for SVLEN, AF, SVTYPE, FILTER, and EV
    SVLEN_table = _make_SVLEN_interval_dict(args.minMetricTable)
    AF_table = _make_AF_interval_dict(args.minMetricTable)
    SVTYPE_table = _make_SVTYPE_dict(args.minMetricTable)
    FILTER_table = _make_FILTER_dict(args.minMetricTable, vcf)
    EV_table = _make_EV_dict(args.minMetricTable)

    # Make filtering lookup table
    minMetric_dict = make_minMetric_dict(args.minMetricTable, SVLEN_table, AF_table,
                                         SVTYPE_table, FILTER_table, EV_table)

    # Iterate over records in vcf and apply filter
    for record in vcf.fetch():
        # Do not process multiallelic variants, unless optioned
        if args.multiallelics or (not args.multiallelics and
                                  not _is_multiallelic(record)):
            # Infer record's AF if AC & AN (but not AF) are provided
            if 'AF' in record.info.keys():
                pass
            else:
                AF = tuple([ac / record.info['AN'] for ac in record.info['AC']])
                record.info['AF'] = AF
            apply_minMetric_filter(record, args.metric, minMetric_dict,
                                   SVLEN_table, AF_table, SVTYPE_table,
                                   FILTER_table, EV_table, globalMin=args.globalMin,
                                   maxNCR=args.maxNCR, highNCR_filter=filter_text)

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
