#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""

"""

import argparse
import pandas as pd


def take_highest_effect(s):
    if 'LOF' in s.values:
        return 'LOF'
    if 'COPY_GAIN' in s.values:
        return 'COPY_GAIN'
    if 'DUP_PARTIAL' in s.values:
        return 'DUP_PARTIAL'
    if 'INV_SPAN' in s.values:
        return 'INV_SPAN'
    if 'INTRONIC' in s.values:
        return 'INTRONIC'
    if 'PROMOTER' in s.values:
        return 'PROMOTER'
    return 'OTHER'


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('protein_coding')
    parser.add_argument('promoter')
    parser.add_argument('fout')
    args = parser.parse_args()

    pc = pd.read_table(args.protein_coding)
    promoter = pd.read_table(args.promoter)

    genic = pd.concat([pc, promoter])
    genic = genic.loc[genic.effect != 'NEAREST_TSS'].copy()

    genic = genic.groupby('name gene'.split())['effect']\
                 .agg(take_highest_effect).reset_index()

    intergenic = pc.loc[~pc.name.isin(genic.name)].copy()
    intergenic['effect'] = 'INTERGENIC'
    intergenic['gene'] = 'NO_GENE'

    genic = pd.concat([genic, intergenic])

    cols = 'name effect gene'.split()
    genic[cols].sort_values('name')\
               .drop_duplicates()\
               .to_csv(args.fout, index=False, sep='\t')


if __name__ == '__main__':
    main()
