#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#

"""
Intersect SV called by PE/SR-based algorithms.

Paired-end and split-read callers provide a reasonably precise estimation of
an SV breakpoint. This program identifies variant calls that fall within
the expected margin of error made by these programs and clusters them together.
The cluster distance defaults to 500 bp but it is recommended to use the
maximum individual clustering distance across the libraries being analyzed.
(Generally median + 7 * MAD)
"""

import argparse
import os
import sys
from collections import deque
from pysam import VariantFile, TabixFile

from svtk.vcfcluster import VCFCluster


def flatten_pos(records, name, fout):
    for record in records:
        chrom = record.CHROM
        start = record.POS
        end = record.INFO['END']
        quad = record.samples[0].sample.split('.')[0]
        source = record.source
        ID = record.ID
        svtype = record.INFO['SVTYPE']

        entry = [str(x) for x in [chrom, start, end, quad, source, ID, name,
                                  svtype]]
        entry = '\t'.join(entry) + '\n'
        fout.write(entry)


def parse_filepaths(filepaths):
    """
    Parameters
    ----------
    filepaths : list of str
        List of paths to standardized VCFs

    Returns
    -------
    vcfs : list of pysam.VariantFile
    """

    vcfs = deque()
    for path in filepaths:
        if len(path.split()) != 1:
            raise ValueError('File list must be single column')
        if not os.path.isfile(path):
            raise FileNotFoundError('VCF {0} not found'.format(path))

        vcf = VariantFile(path)
        vcfs.append(vcf)

    return vcfs


def main(argv):
    parser = argparse.ArgumentParser(
        description=__doc__,
        prog='svtk vcfcluster',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('filelist', type=argparse.FileType('r'),
                        help='List of paths to standardized VCFS')
    parser.add_argument('fout', help='Clustered VCF.')
    parser.add_argument('-r', '--region', default=None,
                        help='Restrict clustering to genomic region.')
    parser.add_argument('-d', '--dist',
                        type=int, default=500,
                        help='Maximum clustering distance. Suggested to use '
                        'max of median + 7*MAD over samples. [500]')
    parser.add_argument('-f', '--frac',
                        type=float, default=0.1,
                        help='Minimum reciprocal overlap between variants. '
                        '[0.1]')
    parser.add_argument('-x', '--blacklist', metavar='BED.GZ',
                        type=TabixFile, default=None,
                        help='Tabix indexed bed of blacklisted regions. Any '
                        'SV with a breakpoint falling inside one of these '
                        'regions is filtered from output.')
    parser.add_argument('-z', '--svsize', type=int, default=500,
                        help='Minimum SV size to report for intrachromosomal '
                        'events. [0]')
    parser.add_argument('-p', '--prefix',
                        default='MERGED',
                        help='Prefix for merged variant IDs. [MERGED]')
    parser.add_argument('-t', '--svtypes', default='DEL,DUP,INV,BND',
                        help='Comma delimited list of svtypes to cluster '
                        '[DEL,DUP,INV,BND]')
    parser.add_argument('--ignore-svtypes', action='store_true', default=False,
                        help='Ignore svtypes when clustering.')
    parser.add_argument('-o', '--sample-overlap', type=float, default=0.0,
                        help='Minimum sample overlap for two variants to be '
                        'clustered together.')
    parser.add_argument('--preserve-ids', action='store_true', default=False,
                        help='Include list of IDs of constituent records in '
                        'each cluster.')
    parser.add_argument('--preserve-genotypes', action='store_true',
                        default=False,
                        help='In a set of clustered variants, report best '
                        '(highest GQ) non-reference genotype when available.')
    parser.add_argument('--preserve-header', action='store_true',
                        default=False,
                        help='Use header from clustering VCFs')
    parser.add_argument('--skip-merge', action='store_true',
                        default=False,
                        help='Do not merge clustered records. Adds CLUSTER info fields.')
    parser.add_argument('--merge-only', action='store_true',
                        default=False,
                        help='When run on a vcf generated with --skip-merge, only merges records '
                             'with identical CLUSTER fields.')
    parser.add_argument('--single-end', action='store_true',
                        default=False,
                        help='Require only one end to be within the minimum distance.')
    #  parser.add_argument('--cluster-bed', type=argparse.FileType('w'),
    #                      help='Bed of constituent calls in each cluster')

    # Print help if no arguments specified
    if len(argv) == 0:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args(argv)

    if args.skip_merge and args.merge_only:
        raise ValueError('Cannot use both --skip-merge and --merge-only')

    # Parse SV files and lists of samples and sources
    filepaths = [line.strip() for line in args.filelist.readlines()]
    vcfs = parse_filepaths(filepaths)

    svtypes = args.svtypes.split(',')
    match_svtypes = not args.ignore_svtypes

    do_merge = not args.skip_merge
    do_cluster = not args.merge_only
    svc = VCFCluster(vcfs, dist=args.dist, blacklist=args.blacklist,
                     frac=args.frac, svtypes=svtypes, region=args.region,
                     match_svtypes=match_svtypes,
                     preserve_ids=args.preserve_ids,
                     preserve_genotypes=args.preserve_genotypes,
                     sample_overlap=args.sample_overlap,
                     preserve_header=args.preserve_header,
                     do_cluster=do_cluster,
                     do_merge=do_merge,
                     single_end=args.single_end)

    # Open new file
    if args.fout in '- stdout'.split():
        fout = sys.stdout
    else:
        fout = open(args.fout, 'w')

    fout = VariantFile(fout, mode='w', header=svc.header)

    for i, cluster in enumerate(svc.cluster()):
        if args.prefix:
            cluster_id = [args.prefix]
        else:
            cluster_id = ['SV']
        if args.region:
            chrom = args.region.split(':')[0]
            cluster_id.append(chrom)
        if do_merge and do_cluster:
            cluster_index = i
        else:
            cluster_index = cluster[0].info['CLUSTER']
        cluster_id.append(str(cluster_index + 1))
        cluster_id = '_'.join(cluster_id)

        for record in cluster:
            # Name record
            if do_merge:
                name = cluster_id
            else:
                name = record.id

            record.id = name
            fout.write(record)

            # Size filter (CTX have size -1)
            if -1 < record.info['SVLEN'] < args.svsize:
                continue

            #  if args.cluster_bed is not None:
                #  flatten_pos(cluster, record.ID, args.cluster_bed)

    fout.close()


if __name__ == '__main__':
    main(sys.argv[1:])
