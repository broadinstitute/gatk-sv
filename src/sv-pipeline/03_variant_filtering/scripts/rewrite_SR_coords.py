#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""

"""

import argparse
import sys
import pysam
import pandas as pd
import numpy as np


def rewrite_SR_coords(record, metrics, pval_cutoff, bg_cutoff):
    row = metrics.loc[record.id]
    if row.SRQ >= pval_cutoff and row.SRCS >= bg_cutoff:
        record.pos = int(row.SR1POS) if row.SR1POS and not np.isnan(row.SR1POS) else record.pos
        record.stop = int(row.SR2POS) if row.SR2POS and not np.isnan(row.SR2POS) else record.stop
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
                              (cutoffs['metric'] == 'SRQ'), 'cutoff'].iloc[0]
    bg_cutoff = cutoffs.loc[(cutoffs['test'] == 'SR1') &
                            (cutoffs['metric'] == 'SRCS'), 'cutoff'].iloc[0]

    for record in records:
        rewrite_SR_coords(record, metrics, pval_cutoff, bg_cutoff)
        if record.info['SVTYPE'] in 'DEL DUP'.split():
            if 'SVLEN' not in record.info or record.info['SVLEN'] < 50:
                continue
        fout.write(record)


if __name__ == '__main__':
    main()
