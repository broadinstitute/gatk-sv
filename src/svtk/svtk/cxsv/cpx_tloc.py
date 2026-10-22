#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""
Classification of reciprocal translocations.
"""


def classify_insertion(plus, minus, mh_buffer=50):
    plus_A = plus.pos
    minus_A = minus.pos
    plus_B = plus.info['END2'] if plus.info['SVTYPE'] == 'BND' else plus.stop
    minus_B = minus.info['END2'] if minus.info['SVTYPE'] == 'BND' else minus.stop

    # Buffer comparisons
    def _greater(p1, p2):
        return p1 > p2 - mh_buffer

    if _greater(minus_A, plus_A) and _greater(minus_B, plus_B):
        return 'INS_B2A'
    elif _greater(plus_A, minus_A) and _greater(plus_B, minus_B):
        return 'INS_A2B'
    else:
        return 'INS_UNCLASSIFIED'


def classify_simple_translocation(plus, minus, mh_buffer=10):
    """
    Resolve a pair of interchromosomal breakends.

    Parameters
    ----------
    FF : pysam.VariantRecord
        FF inversion breakpoint.
    RR : pysam.VariantRecord
        RR inversion breakpoint.
    cnvs : list of pysam.VariantRecord
        List of overlapping CNVs.

    Returns
    -------
    svtype : str
        Complex SV class.
    """

    # plus refers to breakend whose strand begins with '+'
    if plus.chrom != minus.chrom or plus.info['CHR2'] != minus.info['CHR2']:
        return 'TLOC_MISMATCH_CHROM'

    # Reference chromosomes are labeled A and B
    # Breakpoints/Derivative chromosomes are labeled plus and minus, based on
    # ref chromosome A's strandedness on each breakpoint
    # plus_A = the breakend of ref chrom A on derivative chrom where A is
    # forward-stranded

    # get positions
    plus_A = plus.pos
    minus_A = minus.pos
    plus_B = plus.info['END2'] if plus.info['SVTYPE'] == 'BND' else plus.stop
    minus_B = minus.info['END2'] if minus.info['SVTYPE'] == 'BND' else minus.stop

    plus_strands = plus.info['STRANDS']

    # Buffer comparisons
    def _greater(p1, p2):
        return p1 > p2 - mh_buffer

    # Check for PE evidence
    def _hasPE(recA, recB):
        if 'EVIDENCE' in recA.info.keys() \
                and 'EVIDENCE' in recB.info.keys():
            if 'PE' in recA.info['EVIDENCE'] \
                    and 'PE' in recB.info['EVIDENCE']:
                return True
            else:
                return False
        else:
            return False

    if plus_strands == '+-':
        if _greater(minus_A, plus_A) and _greater(plus_B, minus_B):
            if _hasPE(plus, minus):
                return 'CTX_PP/QQ'
            else:
                return 'CTX_UNR'
        if _greater(minus_A, plus_A) and _greater(minus_B, plus_B):
            return 'CTX_INS_B2A'
        if _greater(plus_A, minus_A) and _greater(plus_B, minus_B):
            return 'CTX_INS_A2B'
    else:
        if _greater(minus_A, plus_A) and _greater(minus_B, plus_B):
            if _hasPE(plus, minus):
                return 'CTX_PQ/QP'
            else:
                return 'CTX_UNR'
        if _greater(minus_A, plus_A) and _greater(plus_B, minus_B):
            return 'CTX_INV_INS_B2A'
        if _greater(plus_A, minus_A) and _greater(minus_B, plus_B):
            return 'CTX_INV_INS_A2B'

    return 'CTX_UNR'
