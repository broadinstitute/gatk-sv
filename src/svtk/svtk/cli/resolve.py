#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""
Resolve complex SV from inversion/translocation breakpoints and CNV intervals.
"""

import argparse
import sys
import subprocess
import numpy as np
import string
from collections import deque
from operator import attrgetter
import itertools
import pysam
import pybedtools as pbt
import svtk.utils as svu
from svtk.cxsv import link_cpx, ComplexSV, rescan_single_ender, link_cpx_V2
import datetime


CPX_INFO = [
    '##ALT=<ID=CTX,Description="Reciprocal chromosomal translocation">',
    '##ALT=<ID=CPX,Description="Complex SV">',
    '##ALT=<ID=INS,Description="Insertion">',
    '##ALT=<ID=UNR,Description="Unresolved breakend or complex SV">',
    '##INFO=<ID=SOURCE,Number=1,Type=String,Description="Source of inserted sequence.">',
    '##INFO=<ID=END2,Number=1,Type=Integer,Description="Position of breakpoint on CHR2">',
    '##INFO=<ID=CPX_TYPE,Number=1,Type=String,Description="Class of complex variant.">',
    '##INFO=<ID=CPX_INTERVALS,Number=.,Type=String,Description="Genomic intervals constituting complex variant.">',
    '##INFO=<ID=EVENT,Number=1,Type=String,Description="ID of event associated to breakend">',
    '##INFO=<ID=UNRESOLVED,Number=0,Type=Flag,Description="Variant is unresolved.">',
    '##INFO=<ID=UNRESOLVED_TYPE,Number=1,Type=String,Description="Class of unresolved variant.">'
]


def _merge_records(vcf, cpx_records, cpx_record_ids):
    """
    r1, r2 : iter of pysam.VariantRecord
    """
    def _next_record():
        # Skip VCF records that were included in complex event
        # get next record that's not already present in cpx_record_idss
        _rec = next(vcf, None)
        while _rec is not None and _rec.id in cpx_record_ids:
            _rec = next(vcf, None)
        return _rec

    def _next_cpx():
        try:
            return cpx_records.popleft()
        except IndexError:
            return None
    # Initialize merge
    curr_record = _next_record()
    curr_cpx = _next_cpx()
    while curr_record is not None and curr_cpx is not None:
        # Merge sort records not in complex event
        if curr_record.chrom == curr_cpx.chrom:
            if curr_record.pos <= curr_cpx.pos:
                yield curr_record
                curr_record = _next_record()
            else:
                yield curr_cpx
                curr_cpx = _next_cpx()
        elif svu.is_smaller_chrom(curr_record.chrom, curr_cpx.chrom):
            yield curr_record
            curr_record = _next_record()
        else:
            yield curr_cpx
            curr_cpx = _next_cpx()

    # At least one iterator is exhausted, return rest of other iterator, if any
    while curr_record is not None:
        yield curr_record
        curr_record = _next_record()
    while curr_cpx is not None:
        yield curr_cpx
        curr_cpx = _next_cpx()


def remove_CPX_from_INV(resolve_CPX, resolve_INV):
    """
    Return list of inversion calls not overlapped by list of complex calls
    """
    cpx_interval = [(i.chrom, i.pos, i.stop) for i in resolve_CPX]
    out = [
        inv for inv in resolve_INV
        if not any(cpx[0] == inv.chrom and cpx[1] <= i.stop and i.pos <= cpx[2] for cpx in cpx_interval)
    ]
    return out


def multisort(xs, specs):
    for key, reverse in reversed(specs):
        xs.sort(key=attrgetter(key), reverse=reverse)
    return xs


def cluster_INV(independent_INV):
    list_INV = [multisort(list(group), (('pos', False), ('stop', False)))
                for chrom, group in itertools.groupby(independent_INV, attrgetter('chrom'))]
    return [x for group in list_INV for x in _cluster_INV_list(group)]


def _cluster_INV_list(independent_INV):
    out = []
    rec = float('-inf')  # be sure to add first element
    for i in independent_INV:
        if i.pos > rec:
            out.append([i])
            rec = i.stop
        else:
            out[-1].append(i)
            if i.stop > rec:
                rec = i.stop
    return out


def cluster_single_cleanup(cluster):
    # for clusters including both FF and RR inversions, only keep inversions for CPX resolution
    cluster_svtype = [[i.info['SVTYPE'], i.info['STRANDS']] for i in cluster]
    if ['INV', '--'] in cluster_svtype and ['INV', '++'] in cluster_svtype:
        return deque(
            [
                cluster[cluster_svtype.index(['INV', '--'])],
                cluster[cluster_svtype.index(['INV', '++'])]
            ]
        )
    else:
        return cluster


def clusters_cleanup(clusters):
    return deque(cluster_single_cleanup(cluster) for cluster in clusters)


def get_random_string(random_string_len):
    """
    Produce string of random upper-case characters and digits, of requested length
    """
    return ''.join(np.random.choice(list(string.ascii_uppercase + string.digits))
                   for _ in range(random_string_len))


def resolve_complex_sv(vcf, cytobands, disc_pairs, mei_bed, variant_prefix='CPX_',
                       min_rescan_support=4, pe_blacklist=None, quiet=False,
                       SR_only_cutoff=1000, random_resolved_id_length=10):
    """
    Resolve complex SV from CNV intervals and BCA breakpoints.
    Yields all resolved events, simple or complex, in sorted order.
    Parameters
    ----------
    vcf : pysam.VariantFile
    cytobands : pysam.TabixFile
    disc_pairs : pysam.TabixFile
    mei_bed : pybedtools.BedTool
    variant_prefix : str
        Prefix to assign to resolved variants
    min_rescan_support : int
        Number of pairs required to count a sample as
        supported during PE rescan
    pe_blacklist : pysam.TabixFile, optional
        Blacklisted genomic regions. Anomalous pairs in these regions will be
        removed prior to clustering.
    quiet : boolean, optional
        Do not print status updates
    Yields
    ------
    sv : pysam.VariantRecord
    """

    clusters = link_cpx(vcf)
    clusters = clusters_cleanup(clusters)

    # Print number of candidate clusters identified
    if not quiet:
        now = datetime.datetime.now()
        print('svtk resolve @ ' + now.strftime("%H:%M:%S") + ': ' +
              'identified ' + str(len(clusters)) + ' candidate complex clusters ' +
              'during first pass', flush=True)

    # resolved_idx = unresolved_idx = 1

    if not variant_prefix.endswith('_'):
        variant_prefix += '_'

    cpx_records = deque()
    cpx_record_ids = set()
    np.random.seed(1)  # arbitrary fixed seed for reproducibility

    for cluster in clusters:
        # Print status for each cluster
        if not quiet:
            now = datetime.datetime.now()
            print('svtk resolve @ ' + now.strftime("%H:%M:%S") + ': ' +
                  'resolving candidate cluster containing the following records: ' +
                  ', '.join([e.id for e in cluster]), flush=True)

        # Try finding opposite strand support for single ender inversions
        if len(cluster) == 1 and cluster[0].info['SVTYPE'] == 'INV':
            rec, opp = rescan_single_ender(cluster[0], disc_pairs,
                                           min_rescan_support,
                                           pe_blacklist=pe_blacklist)
            if opp is not None:
                cluster = deque([rec, opp])

        # if cxsv overlap pulled in unrelated insertions, keep them separate
        if all(r.info['SVTYPE'] == 'INS' for r in cluster):
            for record in cluster:
                cpx = ComplexSV([record], cytobands, mei_bed, SR_only_cutoff)
                cpx_record_ids = cpx_record_ids.union(cpx.record_ids)

                # Assign random string as resolved ID to handle sharding
                cpx.vcf_record.id = variant_prefix + \
                    get_random_string(random_resolved_id_length)
                cpx_records.append(cpx.vcf_record)
                # resolved_idx += 1
            outcome = 'treated as separate unrelated insertions'
        else:
            cpx = ComplexSV(cluster, cytobands, mei_bed, SR_only_cutoff)
            cpx_record_ids = cpx_record_ids.union(cpx.record_ids)
            if cpx.svtype == 'UNR':
                # Assign random string as unresolved ID to handle sharding
                unresolved_vid = 'UNRESOLVED_' + \
                    get_random_string(random_resolved_id_length)
                for record in cpx.records:
                    record.info['EVENT'] = unresolved_vid
                    record.info['UNRESOLVED'] = True
                    cpx_records.append(record)
                # unresolved_idx += 1
                outcome = 'is unresolved'
            elif cpx.svtype == 'SPLIT':
                # check all CNVs for depth support and report the first
                # insertion record and all CNVs with depth support. CNVs without
                # depth support will have their IDs added to the MEMBERS field of
                # the INS record
                cnv_ids_to_append = []
                for cnv in cpx.cnvs:
                    if 'EVIDENCE' in cnv.info and 'RD' in cnv.info['EVIDENCE']:
                        cnv.info['MEMBERS'] = (cnv.id,)
                        cpx_records.append(cnv)
                    else:
                        cnv_ids_to_append.append(cnv.id)
                ins_rec = cpx.insertions[0]
                ins_rec.info['MEMBERS'] = (
                    ins_rec.id,) + tuple(cnv_ids_to_append)
                cpx_records.append(ins_rec)
                outcome = 'split into INS and CNV variants. ' + \
                          'The following records were merged into the INS record: ' + \
                          ', '.join(cnv_ids_to_append)
            else:
                cpx.vcf_record.id = variant_prefix + \
                    get_random_string(random_resolved_id_length)
                cpx_records.append(cpx.vcf_record)
                if 'CPX_TYPE' in cpx.vcf_record.info.keys():
                    outcome = 'resolved as ' + \
                        str(cpx.vcf_record.info['CPX_TYPE'])
                else:
                    outcome = 'resolved as ' + \
                        str(cpx.vcf_record.info['SVTYPE'])
                # resolved_idx += 1

        # Report outcome per cluster
        if not quiet:
            now = datetime.datetime.now()
            print('svtk resolve @ ' + now.strftime("%H:%M:%S") + ': ' +
                  'candidate cluster ' + outcome, flush=True)

    # Output all variants
    vcf.reset()

    for record in _merge_records(vcf, cpx_records, cpx_record_ids):
        # Clean all BNDs to ensure they are set to unresolved
        # Reason: some SR-only BNDs will escape being set as UNRESOLVED if they
        # are part of a multi-breakpoint complex cluster that remains UNRESOLVED
        # after excluding SR-only breakpoints
        if record.info['SVTYPE'] == 'BND':
            record.info['UNRESOLVED'] = True
            if 'UNRESOLVED_TYPE' not in record.info.keys():
                record.info['UNRESOLVED_TYPE'] = 'SINGLE_ENDER'
        if 'CPX_TYPE' in record.info.keys():
            if 'UNRESOLVED' in record.info.keys():
                record.info['UNRESOLVED_TYPE'] = record.info['CPX_TYPE']
                record.info.pop('CPX_TYPE')
            else:
                record.info.pop('STRANDS')
        if 'CIPOS' in record.info.keys():
            record.info.pop('CIPOS')
        if 'CIEND' in record.info.keys():
            record.info.pop('CIEND')
        if 'RMSSTD' in record.info.keys():
            record.info.pop('RMSSTD')
        yield record


def cluster_cleanup(clusters_v2):
    cluster_pos = []
    cluster_info = []
    for i in range(len(clusters_v2)):
        info = clusters_v2[i]
        info_pos = '_'.join(sorted(
            [','.join([str(k) for k in [j.pos, j.stop, j.info['SVTYPE']]]) for j in info]))
        if info_pos not in cluster_info:
            cluster_info.append(info_pos)
            cluster_pos.append(i)
    return [clusters_v2[i] for i in cluster_pos]


def resolve_complex_sv_v2(resolve_INV, cytobands, disc_pairs,
                          mei_bed, variant_prefix='CPX_', min_rescan_support=4,
                          pe_blacklist=None, quiet=False, SR_only_cutoff=1000,
                          random_resolved_id_length=10):
    linked_INV = cluster_INV(resolve_INV)
    clusters_v2 = link_cpx_V2(linked_INV, cpx_dist=2000)
    clusters_v2 = cluster_cleanup(clusters_v2)

    np.random.seed(0)  # arbitrary fixed seed for reproducibility

    # Print number of candidate clusters identified
    if not quiet:
        now = datetime.datetime.now()
        print('svtk resolve @ ' + now.strftime("%H:%M:%S") + ': ' +
              'identified ' + str(len(clusters_v2)) + ' candidate complex clusters ' +
              'during second pass', flush=True)

    cpx_records_v2 = deque()
    cpx_record_ids_v2 = set()
    for cluster in clusters_v2:
        # Print status for each cluster
        if not quiet:
            now = datetime.datetime.now()
            print('svtk resolve @ ' + now.strftime("%H:%M:%S") + ': ' +
                  'resolving candidate cluster containing the following records: ' +
                  ', '.join([e.id for e in cluster]), flush=True)

        # Try finding opposite strand support for single ender inversions
        if len(cluster) == 1 and cluster[0].info['SVTYPE'] == 'INV':
            rec, opp = rescan_single_ender(cluster[0], disc_pairs,
                                           min_rescan_support,
                                           pe_blacklist=pe_blacklist)
            if opp is not None:
                cluster = deque([rec, opp])

        # if cxsv overlap pulled in unrelated insertions, keep them separate
        if all(r.info['SVTYPE'] == 'INS' for r in cluster):
            for record in cluster:
                cpx = ComplexSV([record], cytobands, mei_bed, SR_only_cutoff)
                cpx_record_ids_v2.update(cpx.record_ids)

                # Assign random string as resolved ID to handle sharding
                cpx.vcf_record.id = variant_prefix + '_' + \
                    get_random_string(random_resolved_id_length)
                cpx_records_v2.append(cpx.vcf_record)
                # resolved_idx += 1
            outcome = 'treated as separate unrelated insertions'
        else:
            cpx = ComplexSV(cluster, cytobands, mei_bed, SR_only_cutoff)
            cpx_record_ids_v2.update(cpx.record_ids)
            if cpx.svtype == 'UNR':
                # Assign random string as unresolved ID to handle sharding
                unresolved_vid = 'UNRESOLVED_' + \
                    get_random_string(random_resolved_id_length)
                for record in cpx.records:
                    record.info['EVENT'] = unresolved_vid
                    record.info['UNRESOLVED'] = True
                    cpx_records_v2.append(record)
                # unresolved_idx += 1
                outcome = 'is unresolved'
            else:
                cpx.vcf_record.id = variant_prefix + '_' + \
                    get_random_string(random_resolved_id_length)
                cpx_records_v2.append(cpx.vcf_record)
                if 'CPX_TYPE' in cpx.vcf_record.info.keys():
                    outcome = 'resolved as ' + \
                        str(cpx.vcf_record.info['CPX_TYPE'])
                else:
                    outcome = 'resolved as ' + \
                        str(cpx.vcf_record.info['SVTYPE'])

        # Report outcome per cluster
        if not quiet:
            now = datetime.datetime.now()
            print('svtk resolve @ ' + now.strftime("%H:%M:%S") + ': ' +
                  'candidate cluster ' + outcome, flush=True)

    for i in cpx_records_v2:
        if i.info['SVTYPE'] == 'CPX':
            for info in 'UNRESOLVED EVENT UNRESOLVED_TYPE STRANDS'.split():
                if info in i.info.keys():
                    i.info.pop(info)
    return cpx_records_v2


def main(argv):
    parser = argparse.ArgumentParser(
        description=__doc__,
        prog='svtk resolve',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('raw', help='Filtered breakpoints and CNV intervals.')

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--discfile', default=None,
                       help='Scraped discordant pairs. Required '
                       'to attempt to resolve single-ender inversions.')
    group.add_argument('--discfile-list', default=None,
                       type=argparse.FileType('r'),
                       help='Tab-delimited list of discordant pair files '
                       'and indices')
    parser.add_argument('resolved', type=argparse.FileType('w'),
                        help='Resolved simple and complex variants.')
    parser.add_argument('--mei-bed', help='Mobile element insertion bed. '
                        'Required to classify inverted insertions.',
                        required=True)
    parser.add_argument('--cytobands', help='Cytoband file. Required to '
                        'correctly classify interchromosomal events.',
                        required=True)
    #  parser.add_argument('--bincov', help='Bincov file.', required=True)
    #  parser.add_argument('--medianfile', help='Medianfile', required=True)
    #  parser.add_argument('--famfile', help='Fam file', required=True)
    #  parser.add_argument('--cutoffs', help='Random forest cutoffs',
    #  required=True)
    parser.add_argument('--min-rescan-pe-support', type=int, default=4,
                        help='Minumum discordant pairs required during '
                        'single-ender rescan.')
    parser.add_argument('-x', '--pe-blacklist', metavar='BED.GZ',
                        default=None, help='Tabix indexed bed of blacklisted '
                        'regions. Any anomalous pair falling inside one '
                        'of these regions is excluded from PE rescanning.')
    parser.add_argument('-u', '--unresolved', type=argparse.FileType('w'),
                        help='Unresolved complex breakpoints and CNV.')
    parser.add_argument('-p', '--prefix', default='CPX_',
                        help='Variant prefix [CPX_]')
    parser.add_argument('-q', '--quiet', default=False,
                        help='Disable progress logging to stderr.')

    if len(argv) == 0:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args(argv)

    # Print status
    if not args.quiet:
        now = datetime.datetime.now()
        print('svtk resolve @ ' + now.strftime("%H:%M:%S") + ': ' +
              'starting variant resolution.', flush=True)

    vcf = pysam.VariantFile(args.raw)
    for line in CPX_INFO:
        vcf.header.add_line(line)

    resolved_pipe = subprocess.Popen(['vcf-sort', '-c'],
                                     stdin=subprocess.PIPE,
                                     stdout=args.resolved)

    resolved_f = pysam.VariantFile(resolved_pipe.stdin, 'w', header=vcf.header)
    unresolved_f = pysam.VariantFile(args.unresolved, 'w', header=vcf.header)

    cytobands = pysam.TabixFile(args.cytobands)

    mei_bed = pbt.BedTool(args.mei_bed)
    if args.pe_blacklist is not None:
        blacklist = pysam.TabixFile(args.pe_blacklist)
    else:
        blacklist = None
    #  cutoffs = pd.read_table(args.cutoffs)
    #  rdtest = svu.RdTest(args.bincov, args.medianfile, args.famfile,
        #  list(vcf.header.samples), cutoffs)

    if args.discfile is not None:
        disc_pairs = pysam.TabixFile(args.discfile)
    else:
        tabixfiles = []
        for line in args.discfile_list:
            fname, idx = line.strip().split()
            tabixfiles.append(pysam.TabixFile(fname, index=idx))
        disc_pairs = svu.MultiTabixFile(tabixfiles)

    resolved_records = []
    unresolved_records = []
    resolve_INV = []
    # cpx_dist = 20000

    for record in resolve_complex_sv(vcf, cytobands, disc_pairs, mei_bed, args.prefix,
                                     args.min_rescan_pe_support, blacklist, args.quiet):
        # Move members to existing variant IDs unless variant is complex
        if record.info['SVTYPE'] != 'CPX' and args.prefix not in record.id:
            # Don't alter MEMBERS if the prefix of record.id is already in MEMBERS
            if 'MEMBERS' in record.info.keys() and record.id not in record.info['MEMBERS']:
                record.info['MEMBERS'] = (record.id,)
        # Passes unresolved single-ender inversions to second-pass,
        # otherwise writes resolved records to output files
        if record.info['UNRESOLVED']:
            if record.info['SVTYPE'] == 'INV' and record.info['UNRESOLVED_TYPE'] != 'SR_ONLY_LARGE_INVERSION':
                resolve_INV.append(record)
            else:
                unresolved_records.append(record)
        else:
            resolved_records.append(record)

    # out_rec = resolve_complex_sv(vcf, cytobands, disc_pairs, mei_bed, args.prefix, args.min_rescan_pe_support, blacklist)
    # Print status
    if not args.quiet:
        now = datetime.datetime.now()
        print('svtk resolve @ ' + now.strftime("%H:%M:%S") + ': ' +
              'starting second pass through unresolved inversion single-enders ' +
              'for loose inversion linking', flush=True)

    # RLC: As of Sept 19, 2018, only considering inversion single-enders in second-pass
    # due to too many errors in second-pass linking and variant reporting
    cpx_records_v2 = resolve_complex_sv_v2(resolve_INV,
                                           cytobands, disc_pairs, mei_bed, args.prefix,
                                           args.min_rescan_pe_support, blacklist, args.quiet)

    for record in cpx_records_v2:
        # Move members to existing variant IDs unless variant is complex
        if record.info['SVTYPE'] != 'CPX' and 'CPX' not in record.id.split('_'):
            # Don't alter MEMBERS if record.id already in MEMBERS
            if 'MEMBERS' in record.info.keys() and record.id not in record.info['MEMBERS']:
                record.info['MEMBERS'] = (record.id,)
        if record.info['UNRESOLVED']:
            unresolved_records.append(record)
        else:
            resolved_records.append(record)

    # Add back all unresolved inversions from first pass that were not used
    # or discarded by second pass
    v2_used = {r.id for r in cpx_records_v2}
    for record in cpx_records_v2:
        if 'MEMBERS' in record.info.keys():
            v2_used.update(record.info['MEMBERS'])
    unresolved_records.extend(
        record for record in resolve_INV if record.id not in v2_used)

    # Final pass as sanity check: iterate over original VCF, find records that
    # do not appear in either resolved or unresolved VCFs
    k = 0
    used_vids = {r.id for r in resolved_records}
    used_vids.update(r.id for r in unresolved_records)
    for r in resolved_records:
        used_vids.update(r.info['MEMBERS'])
    for r in unresolved_records:
        used_vids.update(r.info['MEMBERS'])

    for r in vcf:
        if r.id not in used_vids:
            k += 1
            if r.info['SVTYPE'] in ('CNV', 'DEL', 'DUP', 'MCNV', 'INS'):
                # remove unresolved tags if present
                r.info.pop('UNRESOLVED', None)
                r.info.pop('UNRESOLVED_TYPE', None)
                resolved_records.append(r)
            else:
                r.info['UNRESOLVED'] = True
                r.info['UNRESOLVED_TYPE'] = 'POSTHOC_RESTORED'
                unresolved_records.append(r)
    # Print outcome
    if not args.quiet:
        now = datetime.datetime.now()
        print('svtk resolve @ ' + now.strftime("%H:%M:%S") + ': ' +
              'final sanity check of original VCF resulted in ' + str(k) +
              ' missing records being restored to final VCF.', flush=True)

    # Write records to file
    for r in resolved_records:
        resolved_f.write(r)
    resolved_f.close()
    for r in unresolved_records:
        r.info['UNRESOLVED'] = True
        if 'UNRESOLVED_TYPE' not in r.info.keys():
            r.info['UNRESOLVED_TYPE'] = 'UNKNOWN'
        if 'CPX_TYPE' in record.info.keys():
            record.info.pop('CPX_TYPE')
        if 'CPX_INTERVALS' in record.info.keys():
            record.info.pop('CPX_INTERVALS')
        unresolved_f.write(r)
    unresolved_f.close()

    stdout, stderr = resolved_pipe.communicate()


if __name__ == '__main__':
    main(sys.argv[1:])
