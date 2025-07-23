#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""
Calculate enrichment of clipped reads or discordant pairs at SV breakpoints.
"""

import argparse
import sys
import pysam
import pandas as pd
from svtk.pesr import SRTestRunner, PETestRunner, PETest, SRTest


def sr_test(argv):
    parser = argparse.ArgumentParser(
        description="Calculate enrichment of clipped reads at SV breakpoints.",
        prog='svtk sr-test',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf',
                        help='VCF of variant calls. Standardized to include '
                        'CHR2, END, SVTYPE, STRANDS in INFO.')
    parser.add_argument('countfile', help='Tabix indexed file of split counts.'
                        ' Columns: chrom,pos,clip,count,sample')
    parser.add_argument('fout',
                        help='Output table of most significant start/end'
                        'positions.')
    parser.add_argument('-w', '--window', type=int, default=100,
                        help='Window around variant start/end to consider for '
                        'split read support. [100]')
    parser.add_argument('--insertion-window', type=int, default=50,
                        help='Maximum distance between left and right breakpoint for INS. Allows for crossing in '
                             'cases of microhomology. [50]')
    parser.add_argument('--common', default=False,
                        action='store_true', help='Ignore background for common AF')
    parser.add_argument('-b', '--background', type=int, default=160,
                        help='Number of background samples to choose for '
                        'comparison in t-test. [160]')
    parser.add_argument('-s', '--samples', type=argparse.FileType('r'),
                        default=None,
                        help='Whitelist of samples to restrict testing to.')
    parser.add_argument('--index', default=None,
                        help='Tabix index of discordant pair file. Required if '
                        'discordant pair file is hosted remotely.')
    # TODO: add normalization
    parser.add_argument('--medianfile', default=None,
                        help='Median coverage statistics for each library '
                        '(optional). If provided, each sample\'s split '
                        'counts will be normalized accordingly. '
                        'Same format as RdTest, one column per sample.')
    parser.add_argument('--log', action='store_true', default=False,
                        help='Print progress log to stderr.')
    parser.add_argument('--outlier-sample-ids', default=None,
                        help='Path to file containing outlier sample IDs.')
    parser.add_argument('--seed', type=int, default=0,
                        help='Random seed.')

    # Print help if no arguments specified
    if len(argv) == 0:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args(argv)

    vcf = pysam.VariantFile(args.vcf)

    if args.index is not None:
        countfile = pysam.TabixFile(args.countfile, index=args.index,
                                    parser=pysam.asTuple())
    else:
        if args.countfile.startswith('http'):
            raise Exception('Must provide tabix index with remote URL')
        countfile = pysam.TabixFile(args.countfile, parser=pysam.asTuple())

    if args.fout in '- stdout'.split():
        fout = sys.stdout
    else:
        fout = open(args.fout, 'w')

    header = 'name coord pos log_pval called_median bg_median bg_frac'.split()
    fout.write('\t'.join(header) + '\n')

    if args.samples is not None:
        whitelist = [s.strip() for s in args.samples.readlines()]
    else:
        whitelist = None

    if args.medianfile is not None:
        medians = pd.read_table(args.medianfile)
        medians = pd.melt(medians, var_name='sample', value_name='median_cov')
    else:
        medians = None

    outlier_sample_ids = None
    if args.outlier_sample_ids:
        outlier_sample_ids = args.outlier_sample_ids

    runner = SRTestRunner(vcf, countfile, fout, args.background, common=args.common,
                          window=args.window, ins_window=args.insertion_window,
                          whitelist=whitelist, medians=medians, log=args.log,
                          outlier_sample_ids=outlier_sample_ids, seed=args.seed)
    runner.run()


def pe_test(argv):
    parser = argparse.ArgumentParser(
        description="Calculate enrichment of discordant pairs at SV breakpoints.",
        prog='svtk pe-test',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf', help='Variants.')
    parser.add_argument('disc', help='Table of discordant pair coordinates.')
    parser.add_argument('fout', type=argparse.FileType('w'),
                        help='Output table of PE counts.')
    parser.add_argument('-o', '--window-out', type=int, default=500,
                        help='Window outside breakpoint to query for '
                        'discordant pairs. [500]')
    parser.add_argument('-i', '--window-in', type=int, default=50,
                        help='Window inside breakpoint to query for '
                        'discordant pairs. [50]')
    parser.add_argument('-b', '--background', type=int, default=160,
                        help='Number of background samples to sample for PE '
                        'evidence. [160]')
    parser.add_argument('--common', default=False,
                        action='store_true', help='Ignore background for common AF')
    parser.add_argument('-s', '--samples', type=argparse.FileType('r'),
                        default=None,
                        help='Whitelist of samples to restrict testing to.')
    parser.add_argument('--index', default=None,
                        help='Tabix index of discordant pair file. Required if '
                        'discordant pair file is hosted remotely.')
    parser.add_argument('--medianfile', default=None,
                        help='Median coverage statistics for each library '
                        '(optional). If provided, each sample\'s split '
                        'counts will be normalized accordingly. '
                        'Same format as RdTest, one column per sample.')
    parser.add_argument('--log', action='store_true', default=False,
                        help='Print progress log to stderr.')
    parser.add_argument('--outlier-sample-ids', default=None,
                        help='Path to file containing outlier sample IDs.')
    parser.add_argument('--seed', type=int, default=0,
                        help='Random seed.')

    if len(argv) == 0:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args(argv)

    if args.vcf in '- stdin'.split():
        vcf = pysam.VariantFile(sys.stdin)
    else:
        vcf = pysam.VariantFile(args.vcf)

    if args.fout in '- stdout'.split():
        fout = sys.stdout
    else:
        fout = args.fout

    header = 'name log_pval called_median bg_median bg_frac'.split()
    args.fout.write('\t'.join(header) + '\n')

    if args.samples is not None:
        whitelist = [s.strip() for s in args.samples.readlines()]
    else:
        whitelist = None

    if args.index is not None:
        discfile = pysam.TabixFile(args.disc, index=args.index)
    else:
        if args.disc.startswith('http'):
            raise Exception('Must provide tabix index with remote URL')
        discfile = pysam.TabixFile(args.disc)

    if args.medianfile is not None:
        medians = pd.read_table(args.medianfile)
        medians = pd.melt(medians, var_name='sample', value_name='median_cov')
    else:
        medians = None

    outlier_sample_ids = None
    if args.outlier_sample_ids:
        outlier_sample_ids = args.outlier_sample_ids

    runner = PETestRunner(vcf, discfile, fout, args.background, args.common, args.window_in, args.window_out,
                          whitelist, medians=medians, log=args.log, outlier_sample_ids=outlier_sample_ids,
                          seed=args.seed)

    runner.run()


def count_pe(argv):
    parser = argparse.ArgumentParser(
        description="Count discordant pairs supporting a SV breakpoints.",
        prog='svtk count-pe',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf', help='Variants.')
    parser.add_argument('disc', help='Table of discordant pair coordinates.')
    parser.add_argument('fout', type=argparse.FileType('w'),
                        help='Output table of PE counts.')
    parser.add_argument('-o', '--window-out', type=int, default=500,
                        help='Window outside breakpoint to query for '
                        'discordant pairs. [500]')
    parser.add_argument('-i', '--window-in', type=int, default=50,
                        help='Window inside breakpoint to query for '
                        'discordant pairs. [50]')
    parser.add_argument('--common', default=False,
                        action='store_true', help='Ignore background for common AF')
    parser.add_argument('-s', '--samples', type=argparse.FileType('r'),
                        default=None,
                        help='Whitelist of samples to restrict testing to.')
    parser.add_argument('--index', default=None,
                        help='Tabix index of discordant pair file. Required if '
                        'discordant pair file is hosted remotely.')
    parser.add_argument('--medianfile', default=None,
                        help='Median coverage statistics for each library '
                        '(optional). If provided, each sample\'s split '
                        'counts will be normalized accordingly. '
                        'Same format as RdTest, one column per sample.')

    if len(argv) == 0:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args(argv)

    if args.vcf in '- stdin'.split():
        vcf = pysam.VariantFile(sys.stdin)
    else:
        vcf = pysam.VariantFile(args.vcf)

    if args.fout in '- stdout'.split():
        fout = sys.stdout
    else:
        fout = args.fout

    header = 'name sample count'.split()
    args.fout.write('\t'.join(header) + '\n')

    if args.samples is not None:
        whitelist = [s.strip() for s in args.samples.readlines()]
    else:
        whitelist = [s for s in vcf.header.samples]

    if args.index is not None:
        discfile = pysam.TabixFile(args.disc, index=args.index)
    else:
        if args.disc.startswith('http'):
            raise Exception('Must provide tabix index with remote URL')
        discfile = pysam.TabixFile(args.disc)

    if args.medianfile is not None:
        medians = pd.read_table(args.medianfile)
        medians = pd.melt(medians, var_name='sample', value_name='median_cov')
    else:
        medians = None

    petest = PETest(discfile, args.common, args.window_in,
                    args.window_out, medians=medians)

    for record in vcf:
        counts = petest.load_counts(record, args.window_in, args.window_out)
        counts = petest.normalize_counts(counts)
        counts = counts.set_index('sample')
        counts = counts.reindex(whitelist).fillna(0).astype(int)
        counts = counts.reset_index()
        counts['name'] = record.id
        cols = 'name sample count'.split()

        for row in counts[cols].to_numpy():
            fout.write('\t'.join([str(x) for x in row]) + '\n')
        #  counts[cols].to_csv(fout, header=False, index=False, sep='\t', na_rep='NA')


def count_sr(argv):
    parser = argparse.ArgumentParser(
        description="Count clipped reads at SV breakpoints. Unwindowed.",
        prog='svtk count-sr',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf',
                        help='VCF of variant calls. Standardized to include '
                        'CHR2, END, SVTYPE, STRANDS in INFO.')
    parser.add_argument('countfile', help='Tabix indexed file of split counts.'
                        ' Columns: chrom,pos,clip,count,sample')
    parser.add_argument('fout',
                        help='Output table of split read counts.')
    parser.add_argument('--common', default=False,
                        action='store_true', help='Ignore background for common AF')
    parser.add_argument('-s', '--samples', type=argparse.FileType('r'),
                        default=None,
                        help='Whitelist of samples to restrict testing to.')
    parser.add_argument('--index', default=None,
                        help='Tabix index of discordant pair file. Required if '
                        'discordant pair file is hosted remotely.')
    # TODO: add normalization
    parser.add_argument('--medianfile', default=None,
                        help='Median coverage statistics for each library '
                        '(optional). If provided, each sample\'s split '
                        'counts will be normalized accordingly. '
                        'Same format as RdTest, one column per sample.')
    # Print help if no arguments specified
    if len(argv) == 0:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args(argv)

    vcf = pysam.VariantFile(args.vcf)

    if args.index is not None:
        countfile = pysam.TabixFile(args.countfile, index=args.index,
                                    parser=pysam.asTuple())
    else:
        if args.countfile.startswith('http'):
            raise Exception('Must provide tabix index with remote URL')
        countfile = pysam.TabixFile(args.countfile, parser=pysam.asTuple())

    if args.fout in '- stdout'.split():
        fout = sys.stdout
    else:
        fout = open(args.fout, 'w')

    header = 'name coord sample count'.split()
    fout.write('\t'.join(header) + '\n')

    if args.samples is not None:
        whitelist = [s.strip() for s in args.samples.readlines()]
    else:
        whitelist = [s for s in vcf.header.samples]

    if args.medianfile is not None:
        medians = pd.read_table(args.medianfile)
        medians = pd.melt(medians, var_name='sample', value_name='median_cov')
    else:
        medians = None
    srtest = SRTest(countfile, args.common, window=0, medians=medians)

    for record in vcf:
        for coord in 'start end'.split():
            if coord == 'start':
                pos, strand, chrom = record.pos, record.info['STRANDS'][0], record.chrom
            else:
                # TODO: With a properly formatted VCF, should be using END2 instead of END here
                pos, strand, chrom = record.stop, record.info['STRANDS'][1], record.info['CHR2']

            counts = srtest.load_counts(chrom, pos, strand)
            counts = srtest.normalize_counts(counts)
            counts = counts['sample count'.split()]
            counts = counts.set_index('sample')
            counts = counts.reindex(whitelist).fillna(0).astype(int)
            counts = counts.reset_index()
            counts['name'] = record.id
            counts['coord'] = coord

            for row in counts[header].values:
                fout.write('\t'.join([str(x) for x in row]) + '\n')
            #  counts[header].to_csv(fout, header=False, index=False, sep='\t', na_rep='NA')
