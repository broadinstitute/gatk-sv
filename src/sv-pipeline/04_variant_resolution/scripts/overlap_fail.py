#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""

"""

import argparse
import sys
import pysam
import svtk.utils as svu
from svtk.vcfcluster import VCFCluster
from svtk.svfile import SVRecordCluster


def samples_overlap(samplesA, samplesB, upper_thresh=0.8, lower_thresh=0.5):
    # Get lists of called samples for each record
    samplesA = set(samplesA)
    samplesB = set(samplesB)

    # Compute fraction of each record's samples which are shared
    shared = samplesA & samplesB
    fracA = len(shared) / len(samplesA)
    fracB = len(shared) / len(samplesB)

    min_frac, max_frac = sorted([fracA, fracB])

    return min_frac >= lower_thresh and max_frac >= upper_thresh


def is_batch_only(cluster, batch):
    ids = [s.record.id for s in cluster.records]
    return all([batch in x for x in ids])


def get_sources(header):
    for record in header.records:
        rec = str(record)
        if rec.startswith('##source='):
            return rec.strip().split('=')[1].split(',')

    return []


def make_new_record(cluster, header):
    # Make new record and merge in pilot samples
    sources = get_sources(header)
    record = header.new_record()
    record = cluster.merge_record_data(record)
    record = cluster.merge_record_formats(record, sources, call_sources=True)
    record.info['MEMBERS'] = [r.record.id for r in cluster.records]
    return record


def is_pesr_only(svrecord):
    return 'depth' not in svrecord.record.info['SOURCES']


def is_depth_only(svrecord):
    return svrecord.record.info['SOURCES'] == ('depth',)


def is_pilot_denovo(svrecord):
    record = svrecord.record
    is_pilot = record.id.startswith('Pilot')

    called = svu.get_called_samples(record)
    is_private = len(called) == 1
    is_child = called[0].endswith('p1') or called[0].endswith('s1')

    return is_pilot and is_private and is_child


def overlap_fail(phase1, pilot, header, dist=300, frac=0.1):
    svc = VCFCluster([phase1, pilot], dist=dist, frac=frac, preserve_ids=True)

    for cluster in svc.cluster(merge=False):
        # Keep any pilot de novo variants
        other_records = []
        for record in cluster.records:
            if is_pilot_denovo(record):
                yield make_new_record(SVRecordCluster([record]), header)
            else:
                other_records.append(record)

        if len(other_records) > 0:
            cluster.records = other_records
        else:
            continue

        # skip phase1-only variants and pilot variants that overlap with
        # rejected phase1 variants
        if is_batch_only(cluster, 'Pilot'):
            records = cluster.records
            depth_only = all([is_depth_only(r) for r in records])

            if len(records) == 1:
                record = make_new_record(cluster, header)
                yield record

            # check that we're not overclustering Pilot variants
            elif len(records) == 2:
                samples = [svu.get_called_samples(r.record) for r in records]

                if samples_overlap(*samples) or depth_only:
                    record = make_new_record(cluster, header)
                    yield record
                else:
                    for svrecord in cluster.records:
                        newcluster = SVRecordCluster([svrecord])
                        record = make_new_record(newcluster, header)
                        yield record

            elif depth_only:
                record = make_new_record(cluster, header)
                yield record

            else:
                import ipdb
                ipdb.set_trace()
                raise Exception('Multiple Pilot variants clustered')

            #  record.id = cluster.records[0].record.id
            #  yield record


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('phase1', help='Failed Phase1 calls')
    parser.add_argument('pilot', help='Filtered pilot calls')
    parser.add_argument('fout')
    parser.add_argument('--prefix', default="SSC_pilotOnly")
    args = parser.parse_args()

    phase1 = pysam.VariantFile(args.phase1)
    pilot = pysam.VariantFile(args.pilot)

    header = phase1.header.copy()
    for sample in pilot.header.samples:
        header.add_sample(sample)

    # TEMPORARY; not all vcfs were clustered since adding MEMBERS
    if 'MEMBERS' not in header.info.keys():
        info = ('##INFO=<ID=MEMBERS,Number=.,Type=String,'
                'Description="IDs of cluster\'s constituent records.">')
        header.add_line(info)

    if args.fout in '- stdout'.split():
        fout = sys.stdout
    else:
        fout = open(args.fout, 'w')
    fout = pysam.VariantFile(fout, mode='w', header=header)

    for i, record in enumerate(overlap_fail(phase1, pilot, header)):
        record.id = args.prefix + '_' + str(i + 1)
        fout.write(record)


if __name__ == '__main__':
    main()
