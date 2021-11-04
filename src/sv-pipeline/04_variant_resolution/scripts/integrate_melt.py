#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""

"""

import argparse
import sys
from collections import deque
import heapq
import pysam
import svtk.utils as svu


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


def integrate_melt(cxsv, melt, fout, window=100):
    cxsv_bed = svu.vcf2bedtool(cxsv, annotate_ins=False, include_samples=True)
    melt_bed = svu.vcf2bedtool(melt, annotate_ins=False, include_samples=True)

    sect = cxsv_bed.window(melt_bed, w=window)

    # Check breakpoints are within window
    def close_enough(interval):
        startA, endA = [int(x) for x in interval.fields[1:3]]
        startB, endB = [int(x) for x in interval.fields[8:10]]
        return abs(startA - startB) < window and abs(endA - endB) < window

    excluded_cxsv = deque()
    for interval in sect.intervals:
        samplesA = interval.fields[6].split(',')
        samplesB = interval.fields[13].split(',')

        if (samples_overlap(samplesA, samplesB) and
                close_enough(interval) and interval.fields[4] == 'INS'):
            excluded_cxsv.append(interval.fields[3])

    cxsv.reset()
    melt.reset()

    for record in heapq.merge(cxsv, melt, key=lambda record: record.pos):
        if record.id in excluded_cxsv:
            #  print(record.id)
            continue
        fout.write(record)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('cxsv', help='Resolved variant VCF')
    parser.add_argument('melt', help='Filtered MELT VCF')
    parser.add_argument('fout', help='Integrated VCF')
    parser.add_argument('-w', '--window', type=int, default=100,
                        help='Window around insertion point within which '
                        'variants will be considered the same.')
    args = parser.parse_args()

    cxsv = pysam.VariantFile(args.cxsv)
    melt = pysam.VariantFile(args.melt)

    if args.fout in '- stdout'.split():
        fout = pysam.VariantFile(sys.stdout, 'w', header=cxsv.header)
    else:
        fout = pysam.VariantFile(args.fout, 'w', header=cxsv.header)

    integrate_melt(cxsv, melt, fout, window=args.window)


if __name__ == '__main__':
    main()
