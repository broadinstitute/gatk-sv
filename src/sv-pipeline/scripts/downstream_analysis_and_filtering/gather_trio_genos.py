#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""
Extract trio allele counts & GQs for all variants in a vcf
"""

import argparse
import sys
import csv
import pysam


def read_ac_adj(infile, pro, fa, mo):
    """
    Read --ac-adj tsv into dictionary
    """

    ac_adj = {}

    with open(infile) as tsvin:
        reader = csv.reader(tsvin, delimiter='\t')
        for vid, ac_pro, ac_fa, ac_mo in reader:
            if vid not in ac_adj.keys():
                acs = [ac_pro, ac_fa, ac_mo]
                acs = [str(i) for i in acs]
                ac_adj[vid] = acs

    return ac_adj


def gather_info(vcf, fout, pro, fa, mo, ac_adj=None, no_header=False):
    GTs_to_skip = './. None/None 0/None None/0'.split()
    sex_chroms = 'X Y chrX chrY'.split()

    # Write header to output file
    if not no_header:
        header = '#VID\tSVLEN\tAF\tSVTYPE\tFILTER\tpro_EV\tpro_AC\tfa_AC\tmo_AC\tpro_GQ\tfa_GQ\tmo_GQ\n'
        fout.write(header)

    trio_samples = [pro, fa, mo]

    if ac_adj is not None:
        vids_to_correct = ac_adj.keys()
    else:
        vids_to_correct = []

    for record in vcf:
        # #Do not include UNRESOLVED variants
        # if 'UNRESOLVED' in record.info.keys() \
        # or 'UNRESOLVED_TYPE' in record.info.keys() \
        # or 'UNRESOLVED' in record.filter:
        #     continue

        # Do not include variants from sex chromosomes
        if record.chrom in sex_chroms:
            continue

        # Do not include multiallelic variants
        if 'MULTIALLELIC' in record.info.keys() \
                or 'MULTIALLELIC' in record.filter \
                or len(record.alts) > 1:
            continue

        # Get GTs for trio
        GTs = [get_GT(record, ID) for ID in trio_samples]

        # Skip sites that aren't het in proband
        if GTs[0] != '0/1':
            continue

        # Skip sites that are reference in all three samples
        if len([g for g in GTs if g == '0/0']) == 3:
            continue

        # Skip sites that are reference or missing in the proband, as these are
        # uninformative for the purposes of assessing de novo rate
        if GTs[0] in GTs_to_skip or GTs[0] == '0/0':
            continue

        # Skip sites that are missing in any sample
        if len([g for g in GTs if g in GTs_to_skip]) > 0:
            continue

        # Convert to ACs
        ACs = [get_AC(g) for g in GTs]

        # Overwrite ACs, if optioned
        if record.id in vids_to_correct:
            # oldACs = ACs
            newACs = ac_adj[record.id]
            for i in [0, 1, 2]:
                ACs[i] = str(max([int(ACs[i]), int(newACs[i])]))
            # oldACs_str = '(' + ', '.join(oldACs) + ')'
            # newACs_str = '(' + ', '.join(ACs) + ')'
            # print('Overwriting ACs for {0} from {1} to {2}\n'.format(record.id,
            #                                                          oldACs_str,
            #                                                          newACs_str))

        # Get genotype qualities for trio
        GQs = [record.samples[ID]['GQ'] for ID in trio_samples]

        # Skip sites that are missing integer GQs in any member of the trio
        # This shouldn't occur in theory, but somtimes does. Cause unclear.
        if len([g for g in GQs if g is None]) > 0:
            continue
        # Otherwise, convert GQs to string for writing to file
        else:
            GQs = [str(g) for g in GQs]

        # Get minimal variant info
        vid = record.id
        size = str(record.info['SVLEN'])
        if 'AF' in record.info.keys():
            freq = str(record.info['AF'][0])
        else:
            freq = 'NA'
        svtype = record.info['SVTYPE']
        filt = ','.join([f for f in record.filter])
        pro_ev = record.samples[trio_samples[0]]['EV']
        if isinstance(pro_ev, tuple):
            pro_ev = ','.join(list(pro_ev))
        vinfo = '{0}\t{1}\t{2}\t{3}\t{4}\t{5}'.format(vid, size, freq,
                                                      svtype, filt, pro_ev)

        # Write record to file
        newline = '{0}\t{1}\t{2}'.format(vinfo, '\t'.join(ACs), '\t'.join(GQs))
        fout.write(newline + '\n')


def get_GT(record, ID):
    GT = record.samples[ID]['GT']
    if GT is not None:
        GT_str = '/'.join([str(i) for i in GT])
    else:
        GT_str = './.'
    return GT_str


def get_AC(GT):
    AC = sum([int(a) for a in GT.split('/')])
    return str(AC)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf', help='Input vcf (supports "stdin").')
    parser.add_argument('fout', help='Output file (supports "stdout").')
    parser.add_argument('pro', help='Proband sample ID.')
    parser.add_argument('fa', help='Father sample ID.')
    parser.add_argument('mo', help='Mother sample ID.')
    parser.add_argument('--ac-adj', help='tsv with variant IDs and ' +
                        'pro/fa/mo AC to be manually overwritten.')
    parser.add_argument('--no-header', help='Do not write header line.',
                        action='store_true', default=False)

    args = parser.parse_args()

    if args.vcf in '- stdin'.split():
        vcf = pysam.VariantFile(sys.stdin)
    else:
        vcf = pysam.VariantFile(args.vcf)

    if args.fout in '- stdout'.split():
        fout = sys.stdout
    else:
        fout = open(args.fout, 'w')

    if args.ac_adj is not None:
        ac_adj = read_ac_adj(args.ac_adj, args.pro, args.fa, args.mo)
    else:
        ac_adj = None

    gather_info(vcf, fout, args.pro, args.fa, args.mo, ac_adj, args.no_header)

    fout.close()


if __name__ == '__main__':
    main()
