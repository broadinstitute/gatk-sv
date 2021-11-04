#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""

"""

import argparse
import sys
import pandas as pd
import pysam
import svtk.utils as svu


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf')
    parser.add_argument('dn_filter')
    parser.add_argument('fout')
    args = parser.parse_args()

    vcf = pysam.VariantFile(args.vcf)

    if args.fout in '- stdout'.split():
        fout = pysam.VariantFile(sys.stdout, 'w', header=vcf.header)
    else:
        fout = pysam.VariantFile(args.fout, 'w', header=vcf.header)

    dn_filter = pd.read_table(args.dn_filter)
    dn_filter = dn_filter.pivot_table(index='name',
                                      columns='operation',
                                      values='sample',
                                      aggfunc=lambda s: ','.join(s))
    for record in vcf:
        if record.id not in dn_filter.index:
            fout.write(record)
            continue

        # Add false negative parents
        parents = dn_filter.loc[record.id, 'add']
        if parents is not None:
            for sample in parents.split(','):
                record.samples[sample]['GT'] = (0, 1)

        # Remove false positive children
        children = dn_filter.loc[record.id, 'remove']
        if children is not None:
            for sample in children.split(','):
                svu.set_null(record, sample)

        # Skip variant if the child was only sample
        called = svu.get_called_samples(record)
        if len(called) == 0:
            continue

        fout.write(record)


if __name__ == '__main__':
    main()
