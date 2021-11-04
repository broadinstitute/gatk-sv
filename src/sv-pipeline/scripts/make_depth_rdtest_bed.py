#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""
Convert bedcluster output to RdTest format
"""

import argparse
import sys
import pandas as pd


def make_depth_rdtest_bed(svof):
    svof['#chrom'] = svof['#chrom'].astype(str)
    svof['start'] = svof.start.astype(int)
    svof['end'] = svof.end.astype(int)
    bed = svof['#chrom start end name svtype'.split()].drop_duplicates()

    # Add samples
    def agg_samples(samples):
        return ','.join(sorted(set(samples)))
    samples = svof.groupby('name')['sample'].agg(agg_samples)
    samples = samples.rename('samples').reset_index()
    bed = pd.merge(bed, samples, on='name', how='left')

    # Format
    bed['svtype'] = bed.svtype.str.upper()

    cols = '#chrom start end name samples svtype'.split()
    return bed[cols]


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('bed', help='Input BED')
    parser.add_argument('fout', help='Output BED', type=argparse.FileType('w'),
                        default=sys.stdout, nargs='?')
    args = parser.parse_args()

    clustered = pd.read_table(args.bed, dtype={'sample': str})

    bed = make_depth_rdtest_bed(clustered)

    bed.to_csv(args.fout, sep='\t', index=False)


if __name__ == '__main__':
    main()
