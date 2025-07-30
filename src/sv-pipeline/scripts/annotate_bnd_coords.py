#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""

"""

import pysam
import argparse


def main():
    parser = argparse.ArgumentParser(description="Update coordinates for BND variants in a VCF file.")
    parser.add_argument("input_vcf", help="Input VCF file path.")
    parser.add_argument("output_vcf", help="Output VCF file path.")
    args = parser.parse_args()

    input_vcf_path = args.input_vcf
    output_vcf_path = args.output_vcf

    with pysam.VariantFile(input_vcf_path, 'r') as vcf_in:
        header = vcf_in.header

        with pysam.VariantFile(output_vcf_path, 'w', header=header) as vcf_out:
            for record in vcf_in:
                if record.info['SVTYPE'] == 'BND' and 'END2' not in record.info:
                    record.info['END2'] = record.stop
                    record.stop = record.pos + 1
                vcf_out.write(record)


if __name__ == "__main__":
    main()
