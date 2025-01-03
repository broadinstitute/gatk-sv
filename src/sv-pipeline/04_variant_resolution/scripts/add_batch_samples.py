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
    parser.add_argument('fout')
    args = parser.parse_args()

    # VCFs from other batches
    with pysam.VariantFile(args.batch_vcf) as batch_vcf, pysam.VariantFile(args.cohort_vcf) as cohort_vcf:
        # Get samples from batch vcf and add to header
        samples = [s for s in batch_vcf.header.samples]
        header = cohort_vcf.header
        for s in samples:
            header.add_sample(s)
        gts = ['0/0' for _ in samples]
        gts[0] = '0/1'
        gt_string = '\tGT\t' + '\t'.join(gts) + '\n'
        # Write out null records for deduplicated variants
        with open(args.fout, mode='w') as fout:
            fout.write(str(header))
            for record in cohort_vcf:
                fout.write(str(record).strip() + gt_string)


if __name__ == '__main__':
    main()
