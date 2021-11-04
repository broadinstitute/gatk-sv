#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

# Retrieved from https://raw.githubusercontent.com/talkowski-lab/rCNV2/a8849b880e2c1e1d2c51f055d85fec07384886a2/data_curation/CNV/defragment_cnvs.py
# on May 14, 2020 and adapted for gatksv pipeline

"""
Perform local defragmentation of a CNV BED file
"""


import pybedtools as pbt
import numpy as np
import argparse


def load_cnvs(inbed):
    """
    Prepare CNVs for defragmentation
    """

    cnvs_orig = pbt.BedTool(inbed)

    def _hack_cnv(cnv):
        """
        Apply hack to chromosome field to allow for easy defragmentation
        """
        cnv.chrom = '_'.join([cnv.chrom, cnv.fields[4], cnv.fields[5]])
        return cnv

    cnvs = cnvs_orig.each(_hack_cnv).sort()

    return cnvs


def defragment_cnvs(cnvs, maxdist=0.25):
    """
    Defragment CNVs based on extending breakpoints by a fraction of CNV size
    """

    def _extend_cnv(cnv, maxdist):
        """
        Extend CNV as a fraction of size
        """
        ext = np.round(len(cnv) * maxdist)
        orig_start = cnv.start
        orig_end = cnv.end
        cnv.start = np.max([0, orig_start - ext])
        cnv.end = orig_end + ext
        cnv.append(orig_start)
        cnv.append(orig_end)
        return cnv

    cnvs_ext = cnvs.each(_extend_cnv, maxdist=maxdist).sort(). \
        merge(c=[4, 5, 6, 7, 8, 9], o='distinct', delim=',')

    def _clean_hit(cnv):
        ostarts = [int(x) for x in cnv[-2].split(',')]
        oends = [int(x) for x in cnv[-1].split(',')]

        cnv.start = np.min(ostarts)
        cnv.end = np.max(oends)

        # Uniquify caller list
        cnv[6] = ",".join(set(cnv[6].split(",")))

        return cnv[:-2]

    return cnvs_ext.each(_clean_hit)


def reformat_cnvs(cnvs):
    """
    Restructure sample & CNV type columns for final output
    """

    def _unhack_cnv(cnv):
        x = cnv.chrom.split('_')
        chrom = x[0]
        fields = [chrom, str(cnv.start), str(cnv.end)]
        fields.extend(cnv.fields[3:])
        newcnv = '\t'.join(fields)

        return newcnv

    return pbt.BedTool('\n'.join([_unhack_cnv(cnv) for cnv in cnvs]),
                       from_string=True)


def main():
    """
    Main block
    """

    # Parse command line arguments and options
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('inbed', help='BED file of CNVs to defragment. Must ' +
                                      'be BED7, where col4 = var ID, col5 = sample ID, '
                                      'col6 = cnv type and col7 = callers.')
    parser.add_argument(
        'outbed', help='Output BED file for defragmented CNVs.')
    parser.add_argument('--max-dist', dest='maxdist', type=float, default=0.25,
                        help='Maximum distance to extend each CNV during ' +
                             'defragmentation, specified as a fraction of total CNV ' +
                             'size.')

    args = parser.parse_args()

    # Load & reformat CNVs for defragmentation
    cnvs = load_cnvs(args.inbed)

    # Defragment CNVs
    defragged_cnvs = defragment_cnvs(cnvs, args.maxdist)

    # Reformat defragged CNVs and write to outfile
    reformat_cnvs(defragged_cnvs).saveas(args.outbed)


if __name__ == '__main__':
    main()
