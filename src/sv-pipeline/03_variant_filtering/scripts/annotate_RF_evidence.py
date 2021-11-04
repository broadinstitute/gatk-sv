#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""

"""

import argparse
import pysam
import pandas as pd


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf')
    parser.add_argument('scores')
    parser.add_argument('fout')
    args = parser.parse_args()

    vcf = pysam.VariantFile(args.vcf)

    info = ('##INFO=<ID=EVIDENCE,Number=.,Type=String,'
            'Description="Classes of random forest support.">')
    vcf.header.add_line(info)

    fout = pysam.VariantFile(args.fout, 'w', header=vcf.header)

    scores = pd.read_table(args.scores)

    records = [r for r in vcf]
    IDs = [r.id for r in records]
    scores = scores.loc[scores.name.isin(IDs)].copy()

    scores = scores.set_index('name')

    for record in records:
        evidence = []
        score = scores.loc[record.id]

        if score.BAF1_prob >= 0.5:
            evidence.append('BAF')
        if score.PE_prob >= 0.5:
            evidence.append('PE')
        if score.RD_prob >= 0.5:
            evidence.append('RD')
        if score.SR1_prob >= 0.5:
            evidence.append('SR')

        record.info['EVIDENCE'] = evidence
        fout.write(record)


if __name__ == '__main__':
    main()
