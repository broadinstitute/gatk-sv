#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""

"""

import argparse
import pysam


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('batch_vcf')
    parser.add_argument('cohort_vcf')
    parser.add_argument('fout', type=argparse.FileType('w'))
    args = parser.parse_args()

    # VCFs from other batches
    batch_vcf = pysam.VariantFile(args.batch_vcf)
    cohort_vcf = pysam.VariantFile(args.cohort_vcf)

    # Copy header of cohort VCF EXCEPT for samples line
    header_lines = str(cohort_vcf.header).strip().split('\n')
    for line in header_lines[:-1]:
        args.fout.write(line + '\n')

    # Copy samples line of batch vcf
    header_lines = str(batch_vcf.header).strip().split('\n')
    args.fout.write(header_lines[-1] + '\n')

    n_samples = len(batch_vcf.header.samples)

    # Write out null records for dedupped variants
    for record in cohort_vcf:
        base = '\t'.join(str(record).strip().split('\t')[:9] + ['GT'])
        null_gts = '\t'.join(['0/0' for i in range(n_samples - 1)])
        args.fout.write(base + '\t0/1\t' + null_gts + '\n')

    args.fout.close()


if __name__ == '__main__':
    main()
