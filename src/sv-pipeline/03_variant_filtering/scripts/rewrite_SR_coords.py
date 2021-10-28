#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""

"""

import argparse
import sys
import pysam
import pandas as pd


def rewrite_SR_coords(record, metrics, pval_cutoff, bg_cutoff):
    row = metrics.loc[record.id]
    if row.SR_sum_log_pval >= pval_cutoff and row.SR_sum_bg_frac >= bg_cutoff:
        record.pos = int(row.SR_posA_pos)
        record.stop = int(row.SR_posB_pos)
        if record.info['SVTYPE'] == 'INV':
            record.pos, record.stop = sorted([record.pos, record.stop])
        if record.info['SVTYPE'] not in 'INS BND'.split():
            record.info['SVLEN'] = record.stop - record.pos


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf')
    parser.add_argument('metrics')
    parser.add_argument('cutoffs')
    parser.add_argument('fout')
    args = parser.parse_args()

    vcf = pysam.VariantFile(args.vcf)
    if args.fout in '- stdout'.split():
        fout = pysam.VariantFile(sys.stdout, 'w', header=vcf.header)
    else:
        fout = pysam.VariantFile(args.fout, 'w', header=vcf.header)

    # Load metrics
    metrics = pd.read_table(args.metrics).drop_duplicates()

    records = [r for r in vcf]
    IDs = [r.id for r in records]
    metrics = metrics.loc[metrics.name.isin(IDs)].copy()

    metrics = metrics.set_index('name')

    # Load cutoffs
    cutoffs = pd.read_table(args.cutoffs)
    pval_cutoff = cutoffs.loc[(cutoffs['test'] == 'SR1') &
                              (cutoffs['metric'] == 'SR_sum_log_pval'), 'cutoff'].iloc[0]
    bg_cutoff = cutoffs.loc[(cutoffs['test'] == 'SR1') &
                            (cutoffs['metric'] == 'SR_sum_bg_frac'), 'cutoff'].iloc[0]

    for record in records:
        rewrite_SR_coords(record, metrics, pval_cutoff, bg_cutoff)
        if record.info['SVTYPE'] in 'DEL DUP'.split():
            if record.info['SVLEN'] < 50:
                continue
        fout.write(record)


if __name__ == '__main__':
    main()
