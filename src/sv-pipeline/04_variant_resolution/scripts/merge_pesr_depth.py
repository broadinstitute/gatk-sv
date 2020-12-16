#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import sys
from collections import defaultdict
import pysam
import svtk.utils as svu


def merge_pesr_depth(vcf, fout, prefix, frac, sample_overlap, min_depth_only_size):

    def _get_shard_path(base_path, index):
        return "{}.shard_{}.vcf.gz".format(base_path, index)

    # Given one pesr record and one depth record, merge depth attributes into the pesr record
    def _merge_pair(record_a, record_b):
        is_depth_a = _record_is_depth(record_a)
        is_depth_b = _record_is_depth(record_b)
        if is_depth_a == is_depth_b:
            raise ValueError("Attempted to write pesr/pesr or depth/depth pair")
        if is_depth_a:
            depth_record = record_a
            pesr_record = record_b
        else:
            pesr_record = record_a
            depth_record = record_b

        pesr_record.info['ALGORITHMS'] = tuple(sorted(set(pesr_record.info['ALGORITHMS'] + ('depth',))))
        pesr_record.info['MEMBERS'] = tuple(sorted(set(pesr_record.info['MEMBERS'] + depth_record.info['MEMBERS'])))

        svu.update_best_genotypes(pesr_record, [pesr_record, depth_record], preserve_multiallelic=True)

        if 'varGQ' in pesr_record.info.keys() and 'varGQ' in depth_record.info.keys():
            pesr_record.info['varGQ'] = max(pesr_record.info['varGQ'], depth_record.info['varGQ'])

        for sample in pesr_record.samples:
            if 'EV' in pesr_record.samples[sample].keys() and 'EV' in depth_record.info.keys():
                pesr_ev = pesr_record.samples[sample]['EV']
                depth_ev = depth_record.samples[sample]['EV']
                pesr_record.samples[sample]['EV'] = tuple(
                    sorted(set(pesr_ev).union(depth_ev)))

        _cache_sample_overlap(pesr_record, force=True)

    def _write_record(record, salvaged):
        # We are done with this record
        _delete_from_sample_overlap_cache(record)
        # Only write depth record if not clustered with a pesr record
        if record.id in clustered_depth_ids:
            return
        svtype = record.info['SVTYPE']
        counts[svtype] += 1
        if _record_is_depth(record):
            new_record = clean_depth_record(base_record, record)
        else:
            new_record = record.copy()
        new_record.id = '{0}_{1}_{2}'.format(prefix, svtype, counts[svtype])
        if salvaged:
            new_record.id = new_record.id + "_salvaged"
        fout.write(new_record)

    def _flush_active_records():
        for r in active_records:
            _write_record(r, False)
        active_records.clear()
        _flush_sample_overlap_cache()

    def _record_is_depth(record):
        alg = record.info['ALGORITHMS']
        return len(alg) == 1 and alg[0] == 'depth'

    def _reciprocal_overlap(record_a, record_b):
        return svu.recip(record_a.pos, record_a.stop, record_b.pos, record_b.stop, frac=frac)

    def _cache_sample_overlap(record, force=False):
        if force or record.id not in sample_overlap_cache:
            _samples = svu.get_called_samples(record)
            sample_overlap_cache[record.id] = set(sample_id_to_index_dict[s] for s in _samples)
            return _samples
        else:
            return sample_overlap_cache[record.id]

    def _delete_from_sample_overlap_cache(record):
        if record.id in sample_overlap_cache:
            del sample_overlap_cache[record.id]

    def _flush_sample_overlap_cache():
        sample_overlap_cache.clear()

    def _sample_overlap(record_a, record_b):
        if sample_overlap == 0:
            return True
        _cache_sample_overlap(record_a)
        _cache_sample_overlap(record_b)
        return svu.samples_overlap(sample_overlap_cache[record_a.id], sample_overlap_cache[record_b.id],
                                   upper_thresh=sample_overlap, lower_thresh=sample_overlap)

    def _records_cluster_together(record_a, record_b):
        return _record_is_depth(record_a) != _record_is_depth(record_b) \
            and record_a.info['SVTYPE'] == record_b.info['SVTYPE'] \
            and _reciprocal_overlap(record_a, record_b) \
            and _sample_overlap(record_a, record_b)

    def _get_base_record(vcf):
        for record in vcf.fetch():
            if not _record_is_depth(record):
                vcf.reset()
                return record

    sample_overlap_cache = {}
    sample_id_to_index_dict = {s: i for i, s in enumerate(vcf.header.samples)}
    cnv_types = ['DEL', 'DUP']
    min_svlen = min_depth_only_size * frac

    base_record = _get_base_record(vcf)
    if base_record is None:
        raise ValueError("No PESR records were found")

    active_records = []
    counts = defaultdict(int)
    clustered_depth_ids = set()
    current_contig = None

    count = 0
    for record in vcf.fetch():

        if count > 0 and count % 1000 == 0:
            sys.stderr.write("Traversed {} records; {} active records; {} record sample sets cached\n"
                             .format(count, len(active_records), len(sample_overlap_cache)))
        count += 1

        # Seed MEMBERS info with original VID
        record.info['MEMBERS'] = (record.id,)

        if record.info['SVTYPE'] not in cnv_types \
                or record.info['SVLEN'] < min_svlen:
            _write_record(record, False)
            continue

        # Write all-ref sites as "salvaged"
        samples = _cache_sample_overlap(record)
        if len(samples) == 0:
            _write_record(record, True)
            continue

        finalized_record_ids = set()
        if current_contig is None or record.contig != current_contig:
            # Started a new contig
            _flush_active_records()
            current_contig = record.chrom
        else:
            # Check if this record belongs to an existing cluster
            for ar in active_records:
                # Upper-bound on current and future overlap (when > 0)
                ar_overlap_test = (ar.stop - record.start) / float(ar.stop - ar.start)
                if ar_overlap_test < frac:
                    # Since traversing in order, this cluster cannot have any more members
                    _write_record(ar, False)
                    finalized_record_ids.add(ar.id)
                elif _records_cluster_together(record, ar):
                    # Merges depth into the pesr record
                    _merge_pair(record, ar)
                    if _record_is_depth(record):
                        clustered_depth_ids.add(record.id)
                    if _record_is_depth(ar):
                        clustered_depth_ids.add(ar.id)
        active_records.append(record)
        active_records = [r for r in active_records if r.id not in finalized_record_ids]

    _flush_active_records()


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
    parser.add_argument('vcf', help='Combined but unmerged VCF of PE/SR calls')
    parser.add_argument('fout', help='Output VCF (unsorted!)')
    parser.add_argument('--interval-overlap', help='Interval reciprocal overlap fraction',
                        type=float, default=0.5)
    parser.add_argument('--sample-overlap', help='Sample overlap fraction',
                        type=float, default=0.5)
    parser.add_argument('--min-depth-only-size', help='Smallest depth only call SVLEN',
                        type=int, default=5000)
    parser.add_argument('--prefix', default='pesr_rd_merged')
    args = parser.parse_args()

    vcf = pysam.VariantFile(args.vcf)
    check_header(vcf)
    fout = pysam.VariantFile(args.fout, 'w', header=vcf.header)
    merge_pesr_depth(vcf, fout=fout, prefix=args.prefix,
                     frac=args.interval_overlap,
                     sample_overlap=args.sample_overlap,
                     min_depth_only_size=args.min_depth_only_size)
    fout.close()


if __name__ == '__main__':
    main()
