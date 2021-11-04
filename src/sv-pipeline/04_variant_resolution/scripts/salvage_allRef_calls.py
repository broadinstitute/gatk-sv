#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""
Reformat 100% reference calls salvaged during pesr + depth merging
"""

import argparse
import sys
import pysam
import random
import string


# Reformat salvaged records to be consistent with existing naming scheme
def reformat(record, prefix):
    record.info['MEMBERS'] = record.id
    record.id = prefix + '_' + \
        ''.join(random.choice(string.ascii_uppercase + string.digits)
                for _ in range(10))
    return record

# Main block


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf')
    parser.add_argument('fout')
    parser.add_argument('-p', '--prefix', default='salvaged',
                        help='Prefix to prepend to reformatted variant IDs.')
    args = parser.parse_args()

    if args.vcf in '- stdin'.split():
        vcf = pysam.VariantFile(sys.stdin)
    else:
        vcf = pysam.VariantFile(args.vcf)

    if args.fout in '- stdout'.split():
        fout = pysam.VariantFile(sys.stdout, 'w', header=vcf.header)
    else:
        fout = pysam.VariantFile(args.fout, 'w', header=vcf.header)

    for record in vcf:
        fout.write(reformat(record, args.prefix))


if __name__ == '__main__':
    main()
