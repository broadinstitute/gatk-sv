#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""
Strip CPX INFO fields from non-CPX variants
"""

import argparse
import sys
import pysam

# Cleanup function


def rmcpx(vcf, fout):

    # Iterate over records
    for record in vcf:

        # Get basic info about record
        svtype = record.info['SVTYPE']

        # Clear CPX info fields for non-CPX variants
        if svtype != "CPX" and svtype != "CTX":
            for info in 'CPX_TYPE CPX_INTERVALS'.split():
                if info in record.info.keys():
                    record.info.pop(info)

        # Write record to file
        fout.write(record)


# Main block
def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf')
    parser.add_argument('fout')

    args = parser.parse_args()

    if args.vcf in '- stdin'.split():
        vcf = pysam.VariantFile(sys.stdin)
    else:
        vcf = pysam.VariantFile(args.vcf)

    if args.fout in '- stdout'.split():
        fout = pysam.VariantFile(sys.stdout, 'w', header=vcf.header)
    else:
        fout = pysam.VariantFile(args.fout, 'w', header=vcf.header)

    rmcpx(vcf, fout)


if __name__ == '__main__':
    main()
