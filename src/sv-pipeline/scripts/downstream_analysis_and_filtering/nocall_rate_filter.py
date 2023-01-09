#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2022 The Broad Institute of M.I.T. and Harvard
# Distributed under terms of the MIT license.
# Contact: Ryan Collins <rlcollins@g.harvard.edu>

"""
Update FILTERs in a VCF based on no-call GT rates
Expects VCF to have been pre-processed with annotate_nocall_rates.py
"""


import argparse
import pysam
import numpy as np
from time import time


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf_in', help='VCF to filter.')
    parser.add_argument('vcf_out', help='Path to output VCF. Accepts "-" ' +
                        'and "stdout". Default: stdout', default='stdout')
    parser.add_argument('--global-max-ncr', type=float, default=0.1,
                        help='Global maximum NCR to be permitted. [default %(default)s]')
    parser.add_argument('--global-filter-on', default='NCR', help='Key for INFO ' +
                        'field to use when applying --global-max-ncr. ' +
                        '[default %(default)s]')
    parser.add_argument('--chrx-filter-on', default=None, help='Key for INFO ' +
                        'field to use when applying --global-max-ncr on chrX. ' +
                        '[defaults to --global-filter-on]')
    parser.add_argument('--chry-filter-on', default=None, help='Key for INFO ' +
                        'field to use when applying --global-max-ncr on chrY. ' +
                        '[defaults to --global-filter-on]')
    parser.add_argument('--hard-filter', default=False, action='store_true',
                        help='Exclude failing records outright from vcf_out. ' +
                        '[default: retain failing records but tag them with ' +
                        'a corresponding FILTER value]')
    parser.add_argument('-v', '--verbose', default=False, action='store_true',
                        help='Print verbose logging [default: no logging]')
    args = parser.parse_args()

    # Open connection to input VCF
    vcf = pysam.VariantFile(args.vcf_in)

    # Check to ensure NCR is present in VCF header
    autosome_ncr_key = args.global_filter_on
    chrx_ncr_key = autosome_ncr_key if args.chrx_filter_on is None else args.chrx_filter_on
    chry_ncr_key = autosome_ncr_key if args.chry_filter_on is None else args.chry_filter_on
    for key in [autosome_ncr_key, chrx_ncr_key, chry_ncr_key]:
        if key not in vcf.header.records.header.info.keys():
            msg = 'ERROR: INFO "{}" not found in header of input VCF. Has this VCF ' + \
                  'been preprocessed by annotate_nocall_rates.py yet?'
            exit(msg.format(key))

    # Update VCF header as needed
    vcf.header.add_meta('FILTER',
                        items=[('ID', "HIGH_NCR"),
                               ('Description', "Unacceptably high rate of no-call GTs")])

    # Open connection to output VCF
    if args.vcf_out in '- stdout'.split():
        outvcf = pysam.VariantFile(sys.stdout, 'w', header=vcf.header)
    else:
        outvcf = pysam.VariantFile(args.vcf_out, 'w', header=vcf.header)

    # Iterate over records
    k = 0
    k_fail = 0
    k_missing = 0
    if args.verbose:
        start_time = time()
    for record in vcf:
        k += 1

        ncr_key = autosome_ncr_key
        if record.chrom == "chrX":
            ncr_key = chrx_ncr_key
        elif record.chrom == "chrY":
            ncr_key = chry_ncr_key

        # Check global NCR
        if ncr_key not in record.info.keys():
            k_missing += 1
        elif record.info[ncr_key] > args.global_max_ncr:
            k_fail += 1
            if args.hard_filter:
                continue
            else:
                record.filter.add('HIGH_NCR')

        # Unless caught by any exception above, write to output vcf
        outvcf.write(record)

        # Print log if optioned
        if args.verbose and k % 10 == 0:
            elapsed = time() - start_time
            msg = 'Progress: {} records parsed | {} seconds per record | {} ({:.1f}%) failed | {} ({:.1f}%) missing NCR values'
            print(msg.format(k, np.round(elapsed / k, 3), k_fail,
                             100 * k_fail / k, k_missing, 100 * k_missing / k))

    # Close connection to output VCF
    outvcf.close()


if __name__ == '__main__':
	main()

