#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""
Convert a VCF to a BED.
"""

import argparse
import sys
import pysam
import pybedtools as pbt
import svtk.utils as svu


def vcf2bed(argv):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf', help='VCF to convert.')
    parser.add_argument('bed', help='Converted bed. Specify `-` or `stdout` '
                        'to write to stdout.')
    parser.add_argument('--no-samples', dest='include_samples',
                        action='store_false', default=True,
                        help='Don\'t include comma-delimited list of called '
                        'samples for each variant.')
    parser.add_argument('-i', '--info', action='append',
                        help='INFO field to include as column in output. '
                        'May be specified more than once. To include all INFO '
                        'fields, specify `--info ALL`. INFO fields are '
                        'reported in the order in which they are requested. '
                        'If ALL INFO fields are requested, they are reported '
                        'in the order in which they appear in the VCF header.')
    parser.add_argument('--include-filters', action='store_true', default=False,
                        help='Include FILTER status in output, with the same behavior an INFO field.')
    parser.add_argument('--split-bnd', action='store_true', default=False,
                        help='Report two entries in bed file for each BND.')
    parser.add_argument('--split-cpx', action='store_true', default=False,
                        help='Report entries for each CPX rearrangement interval.')
    parser.add_argument('--no-header', dest='header', action='store_false',
                        default=True, help='Suppress header.')
    parser.add_argument('--no-sort-coords', dest='no_sort_coords', action='store_true',
                        default=False, help='Do not sort start/end coordinates '
                        'per record before writing to bed.')
    parser.add_argument('--no-unresolved', dest='no_unresolved', action='store_true',
                        default=False, help='Do not output unresolved variants.')
    parser.add_argument('--simple-sinks', dest='simple_sinks', action='store_true',
                        default=False, help='Report all INS sinks as 1bp intervals.')

    # Print help if no arguments specified
    if len(argv) == 0:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args(argv)

    if args.vcf in '- stdin'.split():
        vcf = pysam.VariantFile(sys.stdin)
    else:
        vcf = pysam.VariantFile(args.vcf)

    header = '#chrom start end name svtype'.split()
    if args.include_samples:
        header.append('samples')
    if args.info:
        if 'ALL' in args.info:
            header = header + vcf.header.info.keys()
        else:
            header = header + args.info
    if args.include_filters:
        header = header + ['FILTER']
    header = '\t'.join(header)

    include_unresolved = not args.no_unresolved

    bt = svu.vcf2bedtool(vcf,
                         split_bnd=args.split_bnd,
                         include_samples=args.include_samples,
                         include_strands=False,
                         split_cpx=args.split_cpx,
                         include_infos=args.info,
                         annotate_ins=False,
                         report_alt=True,
                         no_sort_coords=args.no_sort_coords,
                         simple_sinks=args.simple_sinks,
                         include_unresolved=include_unresolved,
                         include_filters=args.include_filters)

    if args.bed in 'stdout -'.split():
        if args.header:
            sys.stdout.write(header + '\n')
        sys.stdout.write(str(bt))
    else:
        if args.header:
            bt.saveas(args.bed, trackline=header)
        else:
            bt.saveas(args.bed)


def remote_tabix(argv):
    parser = argparse.ArgumentParser(
        description="Tabix into a remotely hosted file",
        prog='svtk remote_tabix',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('url')
    parser.add_argument('index')
    parser.add_argument('region', nargs='?', default=None)
    parser.add_argument('-R', '--regions', default=None,
                        help='fetch all regions in bed file')
    parser.add_argument('--header', help='include header', action='store_true',
                        default=False)

    # Print help if no arguments specified
    if len(argv) == 0:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args(argv)

    # Create tabix file
    tbx = pysam.TabixFile(args.url, index=args.index)

    # Print header if requested
    if args.header:
        sys.stdout.write('\n'.join(tbx.header) + '\n')

    # Fetch and output region of interest
    if args.region is not None:
        if args.regions is not None:
            raise Exception("Must specify only one of region or regions file.")

        f = tbx.fetch(region=args.region)
        for line in f:
            sys.stdout.write(line + '\n')

    elif args.regions is not None:
        bed = pbt.BedTool(args.regions).merge()
        for i in bed.intervals:
            f = tbx.fetch(i.chrom, i.start, i.end)
            for line in f:
                sys.stdout.write(line + '\n')

    else:
        raise Exception('Must specify one of region or regions file.')
