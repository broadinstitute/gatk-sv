#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2017 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.

"""

"""

import argparse
import sys
import heapq
from collections import deque, defaultdict
import pysam
import svtk.utils as svu


# def merge_nested_depth_record(pesr_record, depth_record):
#     """Add samples from nested depth record to PE/SR record"""
#
#     pesr_record.info['MEMBERS'] = (pesr_record.info.get('MEMBERS', ()) +
#                                    (depth_record.id, ))
#     pesr_record.info['ALGORITHMS'] = pesr_record.info['ALGORITHMS'] + ('depth', )
#
#     def _quad(s):
#         return s.split('.')[0]
#
#     depth_samples = svu.get_called_samples(depth_record)
#     pesr_samples = svu.get_called_samples(pesr_record)
#
#     # If a sample is called in both the pe/sr record and the nested depth
#     # record, move any relatives called in the depth record to the pe/sr record
#     for quad, samples in itertools.groupby(depth_samples, _quad):
#         samples = list(samples)
#         if any([s in pesr_samples for s in samples]):
#             for sample in samples:
#                 svu.set_null(depth_record, sample)
#                 pesr_record.samples[sample]['GT'] = (0, 1)


def merge_pesr_depth(pesr_vcf, depth_vcf, frac=0.5, sample_overlap=0.5):
    # Memory inefficient but it's easier and shouldn't matter too much
    # now that the variants have been filtered down
    records = dict()
    records['pesr'] = {record.id: record for record in pesr_vcf}
    records['depth'] = {record.id: record for record in depth_vcf}

    # Wipe MEMBERS from prior clustering
    for source in 'pesr depth'.split():
        for ID, record in records[source].items():
            record.info['MEMBERS'] = [ID]

    # Reset for bedtool creation
    pesr_vcf.reset()
    base_record = next(pesr_vcf)

    # Reset for bedtool creation
    pesr_vcf.reset()
    depth_vcf.reset()
    pesr_bed = svu.vcf2bedtool(pesr_vcf, split_bnd=False,
                               include_samples=True,
                               include_strands=False,
                               report_alt=False)
    depth_bed = svu.vcf2bedtool(depth_vcf, split_bnd=False,
                                include_samples=True,
                                include_strands=False,
                                report_alt=False)

    # Remove records with no samples
    def _filter_allref(feature):
        "Returns False if feature has no called samples"
        exclude = False
        if len(feature.fields) == 6:
            samples = feature.fields[5]
            if samples not in ['.', '']:
                exclude = True
        return exclude

    pesr_bed = pesr_bed.filter(_filter_allref).saveas('filtered_pesr.bed')
    depth_bed = depth_bed.filter(_filter_allref).saveas('filtered_depth.bed')

    # Merge depth records with PE/SR records if they share 50% recip overlap
    sect = pesr_bed.intersect(depth_bed, wa=True, wb=True, r=True, f=frac)

    filtered_depth_IDs = deque()
    for pair in sect.intervals:
        # Check SV types match
        if pair.fields[4] != pair.fields[10]:
            continue

        # Get vcf records
        pesr_id, depth_id = pair.fields[3], pair.fields[9]
        pesr_record = records['pesr'][pesr_id]
        depth_record = records['depth'][depth_id]

        # Check for >=50% sample overlap
        samp_ovr = svu.samples_overlap(samplesA=pair.fields[5].split(','),
                                       samplesB=pair.fields[11].split(','))
        if not samp_ovr:
            continue

        # Note removal of depth ID
        filtered_depth_IDs.append(depth_id)

        # Update metadata and samples
        pesr_record.info['MEMBERS'] = (pesr_record.info.get('MEMBERS', ()) +
                                       (depth_record.id, ))
        pesr_record.info['ALGORITHMS'] = pesr_record.info['ALGORITHMS'] + \
            ('depth', )

        svu.update_best_genotypes(pesr_record,
                                  [pesr_record, depth_record],
                                  preserve_multiallelic=True)

        if 'varGQ' in pesr_record.info.keys() and 'varGQ' in depth_record.info.keys():
            pesr_record.info['varGQ'] = max(pesr_record.info['varGQ'],
                                            depth_record.info['varGQ'])

        for sample in pesr_record.samples:
            if 'EV' in pesr_record.samples[sample].keys() and 'EV' in depth_record.info.keys():
                pesr_ev = pesr_record.samples[sample]['EV']
                depth_ev = depth_record.samples[sample]['EV']
                pesr_record.samples[sample]['EV'] = tuple(
                    sorted(set(pesr_ev).union(depth_ev)))

    # Remove overlapping depth records (not performed in for loop to account
    # for double overlaps
    # TODO: handle double overlap of depth calls
    for ID in set(filtered_depth_IDs):
        records['depth'].pop(ID)

    # In remaining depth-only calls, add samples to PE/SR record if the
    # record covers 90% of the depth-only call.
    # SFARI ONLY - REMOVED FOR OTHER ANALYSES
#    sect = pesr_bed.intersect(depth_bed, wa=True, wb=True, F=0.9)
#
#    for pair in sect.intervals:
#        # Check SV types match
#        if pair.fields[4] != pair.fields[10]:
#            continue
#
#        pesr_id, depth_id = pair.fields[3], pair.fields[9]
#
#        # Skip depth records we already added with 50% reciprocal
#        if depth_id in filtered_depth_IDs:
#            continue
#
#        # If sample is in both depth record and pe/sr record, remove it from
#        # depth record
#        depth_record = records['depth'][depth_id]
#        pesr_record = records['pesr'][pesr_id]
#
#        merge_nested_depth_record(pesr_record, depth_record)

    # Merge records together
    def _sort_key(record):
        return (record.chrom, record.pos, record.info['CHR2'], record.stop)

    pesr_records = sorted(records['pesr'].values(), key=_sort_key)
    depth_records = sorted(records['depth'].values(), key=_sort_key)
    depth_records = [clean_depth_record(base_record, r) for r in depth_records]

    for record in heapq.merge(pesr_records, depth_records, key=_sort_key):
        # Clean out unwanted format keys
        # EDIT - this should be handled upstream by add_genotypes
        #  FORMATS = 'GT GQ RD_CN RD_GQ PE_GT PE_GQ SR_GT SR_GQ EV'.split()
        #  for key in record.format.keys():
        #  if key not in FORMATS:
        #  del record.format[key]

        record.info['ALGORITHMS'] = sorted(set(record.info['ALGORITHMS']))
        record.info['MEMBERS'] = sorted(set(record.info.get('MEMBERS', ())))

        # Skip emptied depth records
        if len(svu.get_called_samples(record)) == 0:
            continue

        yield record


def clean_depth_record(base_record, depth_record):
    base = base_record.copy()
    base.chrom = depth_record.chrom
    base.pos = depth_record.pos
    base.id = depth_record.id
    base.ref = depth_record.ref
    base.alts = depth_record.alts
    base.stop = depth_record.stop

    for key in base.info.keys():
        if key not in depth_record.info:
            base.info.pop(key)

    for key, val in depth_record.info.items():
        base.info[key] = val

    for sample in depth_record.samples:
        for key, val in depth_record.samples[sample].items():
            base.samples[sample][key] = val

    return base


def check_header(vcf):
    if 'MEMBERS' not in vcf.header.info:
        vcf.header.add_line(
            "##INFO=<ID=MEMBERS,Number=.,Type=String,Description=\"IDs of cluster's constituent records.\">")


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('pesr', help='VCF of PE/SR calls')
    parser.add_argument('depth', help='VCF of depth calls')
    parser.add_argument('fout')
    parser.add_argument('--prefix', default='pesr_rd_merged')
    args = parser.parse_args()

    pesr = pysam.VariantFile(args.pesr)
    depth = pysam.VariantFile(args.depth)
    check_header(pesr)
    check_header(depth)

    if args.fout in '- stdout'.split():
        fout = pysam.VariantFile(sys.stdout, 'w', header=pesr.header)
    else:
        fout = pysam.VariantFile(args.fout, 'w', header=pesr.header)

    counts = defaultdict(int)

    for record in merge_pesr_depth(pesr, depth):
        svtype = record.info['SVTYPE']
        counts[svtype] += 1
        record.id = '{0}_{1}_{2}'.format(args.prefix, svtype, counts[svtype])
        fout.write(record)


if __name__ == '__main__':
    main()
