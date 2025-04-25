#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""

"""

import argparse
import sys
import pysam
import numpy as np
import pandas as pd
from svtk.vcfcluster import VCFCluster
import svtk.utils as svu


def is_batch_only(cluster, batch):
    ids = [s.record.id for s in cluster.records]
    return all([batch in x for x in ids])


def get_sources(header):
    for record in header.records:
        rec = str(record)
        if rec.startswith('##source='):
            return rec.strip().split('=')[1].split(',')

    return []


def get_overcluster_stats(records):
    pilot = [r for r in records if 'Pilot' in r.id]
    phase1 = [r for r in records if 'Phase1' in r.id]

    dat = []
    for batch, records in zip('pilot phase1'.split(), [pilot, phase1]):
        names = ';'.join([r.id for r in records])
        sources = ';'.join([','.join(r.info['SOURCES']) for r in records])
        n = len(records)
        dat.append([batch, names, sources, n])

    df = pd.DataFrame(dat)
    df = df.replace('', np.nan)
    df.columns = 'batch records sources count'.split()
    return df


def to_bed(record, cluster):
    fmt = ('{chrom}\t{start}\t{end}\t{name}\t{svtype}\t{samples}\t{sources}\t'
           '{batch}\t{cluster}\n')

    batch = 'Pilot' if 'Pilot' in record.id else 'Phase1'

    return fmt.format(chrom=record.chrom, start=record.pos, end=record.stop,
                      name=record.id + '__' + record.chrom,
                      svtype=record.info['SVTYPE'],
                      samples=','.join(svu.get_called_samples(record)),
                      sources=','.join(record.info['SOURCES']),
                      batch=batch, cluster=cluster)


def overlap_pass(phase1, pilot, fout, dist=300, frac=0.1, prefix="SSC_merged"):
    svc = VCFCluster([phase1, pilot], dist=dist, frac=frac, preserve_ids=True)
    sources = get_sources(fout.header)

    # Helper for testing if SVRecord has pe/sr support
    pesr_sources = set('delly dragen lumpy manta wham'.split())

    def _has_pesr(record):
        sources = set(record.record.info['SOURCES'])
        return bool(sources & pesr_sources)

    for i, cluster in enumerate(svc.cluster(merge=False)):
        # skip pilot-only variants
        if is_batch_only(cluster, 'Pilot'):
            continue

        # Make new record and merge in pilot samples
        record = fout.header.new_record()
        record = cluster.merge_record_data(record)
        record = cluster.merge_record_formats(record, sources,
                                              call_sources=True)
        record.info['MEMBERS'] = [r.record.id for r in cluster.records]

        # Anecdotally, overclustering was observed to be small pe/sr-only
        # and a corresponding depth-only which didn't quite meet 80%
        records = cluster.records
        phase1_records = [r for r in records if 'Phase1' in r.record.id]
        if len(phase1_records) != 1:
            pesr_records = [r for r in phase1_records if _has_pesr(r)]
            depth_records = [r for r in phase1_records if not _has_pesr(r)]

            # Set record position to median of phase1, defaulting to pe/sr
            if len(pesr_records) > 0:
                record.pos = np.median([r.record.pos for r in pesr_records])
                record.stop = np.median([r.record.stop for r in pesr_records])
            else:
                record.pos = np.median([r.record.pos for r in depth_records])
                record.stop = np.median([r.record.stop for r in depth_records])
        else:
            phase1_record = phase1_records[0]
            record.pos = phase1_record.record.pos
            record.stop = phase1_record.record.stop

        record.id = prefix + '_' + str(i)
        fout.write(record)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('phase1')
    parser.add_argument('pilot')
    parser.add_argument('fout')
    parser.add_argument('--prefix')
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

    overlap_pass(phase1, pilot, fout, prefix=args.prefix)
    #  for record in overlap_pass(phase1, pilot, header):
    #      fout.write(record)


if __name__ == '__main__':
    main()
