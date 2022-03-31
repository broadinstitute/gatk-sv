#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2022 The Broad Institute of M.I.T. and Harvard
# Distributed under terms of the MIT license.
# Contact: Ryan Collins <rlcollins@g.harvard.edu>

"""
Split an input VCF into multiple smaller shards based on input lists of variant IDs
"""


import argparse
import pysam
import sys


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf', help='Input vcf. Also accepts "stdin" and "-".')
    parser.add_argument('-i', '--id-lists', help='List containing the paths to ' +
                        'all variant ID lists for sharding.', required=True)
    parser.add_argument('-p','--prefix', help='Filename prefix for output VCF shards.',
                        type=str, default='vcf_shard_')
    args = parser.parse_args()

    # Load all variant ID lists
    vid_map = {}
    n_shards = 0
    with open(args.id_lists) as fin:
        for i, inpath in enumerate(fin.readlines()):
            with open(inpath.rstrip()) as flist:
                for vid in flist.readlines():
                    vid_map[vid.rstrip()] = i
                n_shards += 1

    # Open connections to input VCF
    if args.vcf in '- stdin'.split():
        vcf = pysam.VariantFile(sys.stdin) 
    else:
        vcf = pysam.VariantFile(args.vcf)

    # Open connections to all shard outfiles
    vout = args.prefix + '{}.vcf.gz'
    out_vcfs = {i : pysam.VariantFile(vout.format(i), 'w', header=vcf.header) \
                for i in range(n_shards)}

    # Iterate over input VCF and write records to file
    for record in vcf:
        shard_i = vid_map.get(record.id, None)
        if shard_i is not None:
            out_vcfs[shard_i].write(record)

    # Close all open file handles
    for ovcf in out_vcfs.values():
        ovcf.close()


if __name__ == '__main__':
    main()

