#!/usr/bin/env python


import argparse
import sys
import pysam
import numpy as np
from time import time


_cnv_types = ('DEL', 'DUP')


def get_threshold(record, global_max_ncr, thresholds, med_size, large_size):
    svtype = record.info['SVTYPE']
    if svtype in _cnv_types:
        svlen = record.info['SVLEN']
        if svlen < med_size:
            size_index = 0
        elif svlen < large_size:
            size_index = 1
        else:
            size_index = 2
        threshold = thresholds[svtype][size_index]
    else:
        threshold = thresholds[svtype][0]
    if threshold is None:
        if svtype != "CNV" and global_max_ncr is not None:
            threshold = global_max_ncr  # global max NCR only used if class threshold not provided
        else:
            threshold = 1.0  # if no class or global threshold provided, do not filter on NCR
    return threshold


def get_ncr_key(record, autosome_ncr_key, chrx_ncr_key, chry_ncr_key, chr_x, chr_y):
    if record.chrom == chr_x:
        return chrx_ncr_key
    elif record.chrom == chr_y:
        return chry_ncr_key
    else:
        return autosome_ncr_key


def process(vcf, outvcf, global_max_ncr, thresholds, autosome_ncr_key, chrx_ncr_key, chry_ncr_key, 
            chr_x, chr_y, medium_size, large_size, hard_filter, verbose):
    k = 0
    k_fail = 0
    k_missing = 0
    start_time = time()
    for record in vcf:
        k += 1

        threshold = get_threshold(record, global_max_ncr, thresholds, medium_size, large_size)

        ncr_key = get_ncr_key(record, autosome_ncr_key, chrx_ncr_key, chry_ncr_key, chr_x, chr_y)

        # Check global NCR
        if ncr_key not in record.info.keys():
            k_missing += 1
        elif record.info[ncr_key] > threshold:
            k_fail += 1
            if hard_filter:
                continue
            else:
                record.filter.add('HIGH_NCR')

        # Unless caught by any exception above, write to output vcf
        outvcf.write(record)

        # Print log if optioned
        if verbose and k % 1000 == 0:
            elapsed = time() - start_time
            msg = 'Progress: {} records parsed | {} seconds per record | {} ({:.1f}%) failed | ' + \
                '{} ({:.1f}%) missing NCR values'
            print(msg.format(k, np.round(elapsed / k, 3), k_fail,
                             100 * k_fail / k, k_missing, 100 * k_missing / k))


def create_threshold_dict(args):
    return {
        'DEL': [args.small_del_threshold, args.medium_del_threshold, args.large_del_threshold],
        'DUP': [args.small_dup_threshold, args.medium_dup_threshold, args.large_dup_threshold],
        'INS': [args.ins_threshold],
        'INV': [args.inv_threshold],
        'BND': [args.bnd_threshold],
        'CPX': [args.cpx_threshold],
        'CTX': [args.ctx_threshold],
        'CNV': [None]
    }


def main():
    parser = argparse.ArgumentParser(
        description="Update FILTERs in a VCF based on no-call GT rates. " +
                    "Expects VCF to have been pre-processed with annotate_nocall_rates.py",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf_in', help='VCF to filter.')
    parser.add_argument('vcf_out', help='Path to output VCF. Accepts "-" ' +
                        'and "stdout". Default: stdout', default='stdout')
    parser.add_argument('--global-max-ncr', type=float,
                        help='Global maximum NCR to be permitted. ' +
                             'Overridden by subclass thresholds if provided.')
    parser.add_argument('--global-filter-on', default='NCR', help='Key for INFO ' +
                        'field to use when applying --global-max-ncr. ' +
                        '[default %(default)s]')
    parser.add_argument('--chrx-filter-on', default=None, help='Key for INFO ' +
                        'field to use when applying --global-max-ncr on chrX. ' +
                        '[defaults to --global-filter-on]')
    parser.add_argument('--chry-filter-on', default=None, help='Key for INFO ' +
                        'field to use when applying --global-max-ncr on chrY. ' +
                        '[defaults to --global-filter-on]')
    parser.add_argument('--chr-x', default="chrX", help="String for chrX")
    parser.add_argument('--chr-y', default="chrY", help="String for chrY")
    parser.add_argument('--hard-filter', default=False, action='store_true',
                        help='Exclude failing records outright from vcf_out. ' +
                        '[default: retain failing records but tag them with ' +
                        'a corresponding FILTER value]')
    parser.add_argument("--medium-size", type=int, default=500,
                        help="Min size for medium DEL/DUP")
    parser.add_argument("--large-size", type=int, default=10000,
                        help="Min size for large DEL/DUP")
    parser.add_argument("--small-del-threshold", type=float,
                        help="Threshold SL for small DELs")
    parser.add_argument("--medium-del-threshold", type=float,
                        help="Threshold SL for medium DELs")
    parser.add_argument("--large-del-threshold", type=float,
                        help="Threshold SL for large DELs")
    parser.add_argument("--small-dup-threshold", type=float,
                        help="Threshold SL for small DUPs")
    parser.add_argument("--medium-dup-threshold", type=float,
                        help="Threshold SL for medium DUPs")
    parser.add_argument("--large-dup-threshold", type=float,
                        help="Threshold SL for large DUPs")
    parser.add_argument("--ins-threshold", type=float,
                        help="Threshold SL for INS")
    parser.add_argument("--inv-threshold", type=float,
                        help="Threshold SL for INV")
    parser.add_argument("--bnd-threshold", type=float,
                        help="Threshold SL for BND")
    parser.add_argument("--cpx-threshold", type=float,
                        help="Threshold SL for CPX")
    parser.add_argument("--ctx-threshold", type=float,
                        help="Threshold SL for CTX")
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

    thresholds = create_threshold_dict(args)

    # Iterate over records
    process(vcf, outvcf, args.global_max_ncr, thresholds, autosome_ncr_key, chrx_ncr_key, chry_ncr_key,
            args.chr_x, args.chr_y, args.medium_size, args.large_size, args.hard_filter, args.verbose)

    # Close connection to output VCF
    outvcf.close()


if __name__ == '__main__':
    main()
