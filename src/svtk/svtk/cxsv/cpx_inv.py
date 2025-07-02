#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""
Classification of complex inversion events.
"""

from collections import defaultdict
import numpy as np
import svtk.utils as svu


def breakpoint_ordering(FF, RR, mh_buffer=10):
    """
    Match paired breakpoints to known coordinate orderings.

    e.g. in a simple inversion, FF_start < RR_start < FF_end < RR_end.

    Parameters
    ----------
    FF : pysam.VariantRecord
    RR : pysam.VariantRecord
    mh_buffer : int, optional
        Microhomology buffer. Add window to coordinates to permit fuzzy match.

    Returns
    -------
    ordering : str
        SIMPLE/DEL, DUP5/INS3, DUP3/INS5, dupINVdup, or UNK.
    """

    ff_end = FF.info['END2'] if FF.info['SVTYPE'] == 'BND' else FF.stop
    rr_end = RR.info['END2'] if RR.info['SVTYPE'] == 'BND' else RR.stop
    # Check if breakpoints match simple/deletion ordering
    # (FF_start < RR_start < FF_end < RR_end)
    del_order = ((RR.pos > FF.pos - mh_buffer) and
                 (ff_end > RR.pos) and
                 (rr_end > ff_end - mh_buffer))

    # Check if breakpoints match 5' dup ordering
    # (RR_start < FF_start < FF_end < RR_end)
    dup5_order = ((RR.pos < FF.pos) and
                  (FF.pos < ff_end) and
                  (ff_end < rr_end + mh_buffer))

    # Check if breakpoints match 3' dup ordering
    # (FF_start < RR_start < RR_end < FF_end)
    dup3_order = ((FF.pos < RR.pos + mh_buffer) and
                  (RR.pos < rr_end) and
                  (rr_end < ff_end))

    # Check if breakpoints match dupINVdup ordering
    # (RR_start < FF_start < RR_end < FF_end)
    dupINVdup_order = (RR.pos < FF.pos < rr_end < ff_end)

    if del_order:
        return 'SIMPLE/DEL'
    elif dup5_order:
        return 'DUP5/INS3'
    elif dup3_order:
        return 'DUP3/INS5'
    elif dupINVdup_order:
        return 'dupINVdup'
    else:
        return 'UNK'


def breakpoints_match(FF, RR, svtype, mh_buffer=10):
    """
    Test if the coordinate ordering matches the class predicted by CNV overlap.

    Parameters
    ----------
    FF : pysam.VariantRecord
    RR : pysam.VariantRecord
    svtype : str
        Class of complex SV predicted by CNV overlap.
        (delINV, INVdel, delINVdel, dupINV, dupINVdel, INVdup, delINVdup,
         dupINVdup)
    mh_buffer : int, optional
        Microhomology buffer. Add window to coordinates to permit fuzzy match.

    Returns
    -------
    ordering : str
        SIMPLE/DEL, DUP5/INS3, DUP3/INS5, dupINVdup, or UNK
    """
    order = breakpoint_ordering(FF, RR, mh_buffer)

    if svtype in 'delINV INVdel delINVdel'.split():
        return order == 'SIMPLE/DEL'
    elif svtype in 'dupINV dupINVdel DUP5/INS3'.split():
        return order == 'DUP5/INS3'
    elif svtype in 'INVdup delINVdup DUP3/INS5'.split():
        return order == 'DUP3/INS5'
    else:
        return order == 'dupINVdup'


def classify_2_cnv(FF, RR, cnvs, min_frac=0.5):
    """
    Classify the cxSV class of a pair of inv bkpts and two associated CNVs.

    Matches each CNV to a 5' or 3' location, as constrained by the breakpoint
    coordinates.

    Parameters
    ----------
    FF : pysam.VariantRecord
    RR : pysam.VariantRecord
    cnvs : [pysam.VariantRecord, pysam.VariantRecord]
    min_frac : float, optional
        Minimum reciprocal overlap of each cnv with a candidate CNV interval
        defined by the breakpoint coordinates.

    Returns
    -------
    svtype : str
    """

    # Assign CNVs to 5' or 3' based on ordering
    cnv5, cnv3 = sorted(cnvs, key=lambda r: r.pos)

    # Check if 5' CNV matches breakpoints
    if cnv5.info['SVTYPE'] == 'DEL':
        interval5 = (FF.pos, RR.pos)
    else:
        interval5 = (RR.pos, FF.pos)
    cnv5_end = cnv5.info['END2'] if cnv5.info['SVTYPE'] == 'BND' else cnv5.stop
    frac5 = svu.reciprocal_overlap(cnv5.pos, cnv5_end, *interval5)

    # Check if 3' CNV matches breakpoints
    ff_end = FF.info['END2'] if FF.info['SVTYPE'] == 'BND' else FF.stop
    rr_end = RR.info['END2'] if RR.info['SVTYPE'] == 'BND' else RR.stop
    if cnv3.info['SVTYPE'] == 'DEL':
        interval3 = (ff_end, rr_end)
    else:
        interval3 = (rr_end, ff_end)
    cnv3_end = cnv3_end.info['END2'] if cnv3_end.info['SVTYPE'] == 'BND' else cnv3_end.stop
    frac3 = svu.reciprocal_overlap(cnv3.pos, cnv3_end, *interval3)

    # Report cxSV class based on whether CNVs matched intervals
    if frac5 >= min_frac and frac3 >= min_frac:
        svtype = (cnv5.info['SVTYPE'].lower() +
                  'INV' +
                  cnv3.info['SVTYPE'].lower())
    elif frac5 >= min_frac and frac3 < min_frac:
        return classify_1_cnv(FF, RR, cnv5)
    elif frac5 < min_frac and frac3 >= min_frac:
        return classify_1_cnv(FF, RR, cnv3)
    else:
        svtype = 'CNV_2_FAIL'

    return svtype, cnvs


def classify_1_cnv(FF, RR, cnv, min_frac=0.5,
                   min_bkpt_cnv_size=500, max_bkpt_cnv_size=4000):
    """
    Classify the cxSV class of a pair of inv bkpts and one associated CNV.

    Matches each CNV to a 5' or 3' location, as constrained by the breakpoint
    coordinates. After matching CNV, check if distance between breakpoints
    at other end is sufficient to call a second flanking CNV.

    Parameters
    ----------
    FF : pysam.VariantRecord
    RR : pysam.VariantRecord
    cnvs : [pysam.VariantRecord, pysam.VariantRecord]
    min_frac : float, optional
        Minimum reciprocal overlap of each cnv with a candidate CNV interval
        defined by the breakpoint coordinates.
    min_bkpt_cnv_size : int, optional
        Minimum distance between breakpoints to call flanking CNV.
    max_bkpt_cnv_size : int, optional
        Maximum distance between breakpoints to call flanking CNV.

    Returns
    -------
    svtype : str
    """

    # Make CNV class lowercase (for later concatenation with INV)
    cnv_type = cnv.info['SVTYPE'].lower()

    # Determine eligible 5'/3' CNV intervals defined by the breakpoints
    ff_end = FF.info['END2'] if FF.info['SVTYPE'] == 'BND' else FF.stop
    rr_end = RR.info['END2'] if RR.info['SVTYPE'] == 'BND' else RR.stop
    if cnv_type == 'del':
        interval5 = (FF.pos, RR.pos)
        interval3 = (ff_end, rr_end)
    else:
        interval5 = (RR.pos, FF.pos)
        interval3 = (rr_end, ff_end)

    # Check overlap of CNV against full inversion length
    start = min(FF.pos, RR.pos)
    end = max(ff_end, rr_end)
    cnv_end = cnv.info['END2'] if cnv.info['SVTYPE'] == 'BND' else cnv.stop
    total_frac = svu.reciprocal_overlap(cnv.pos, cnv_end, start, end)
    frac5 = svu.overlap_frac(*interval5, cnv.pos, cnv_end)
    frac3 = svu.overlap_frac(*interval3, cnv.pos, cnv_end)

    # If one CNV spans the entire event, it likely represents two CNV merged
    # during preprocessing or clustering
    if total_frac > 0.9 and frac5 > 0.95 and frac3 > 0.95:
        svtype = cnv_type + 'INV' + cnv_type  # + '_merged'
        return svtype, [cnv]

    # Otherwise, check whether it's 5' or 3'
    cnv_end = cnv.info['END2'] if cnv.info['SVTYPE'] == 'BND' else cnv.stop
    frac5 = svu.reciprocal_overlap(cnv.pos, cnv_end, *interval5)
    frac3 = svu.reciprocal_overlap(cnv.pos, cnv_end, *interval3)

    # 5' CNV; check 3' breakpoints for small flanking CNV
    if frac5 >= min_frac and frac3 < min_frac:
        svtype = cnv_type + 'INV'

        ff_end = FF.info['END2'] if FF.info['SVTYPE'] == 'BND' else FF.stop
        rr_end = RR.info['END2'] if RR.info['SVTYPE'] == 'BND' else RR.stop
        dist3 = rr_end - ff_end
        if min_bkpt_cnv_size <= dist3 < max_bkpt_cnv_size:
            svtype = svtype + 'del'
        elif min_bkpt_cnv_size <= -dist3 < max_bkpt_cnv_size:
            svtype = svtype + 'dup'

    # 3' CNV; check 5' breakpoints for small flanking CNV
    elif frac5 < min_frac and frac3 >= min_frac:
        svtype = 'INV' + cnv_type

        dist5 = RR.pos - FF.pos
        if min_bkpt_cnv_size <= dist5 < max_bkpt_cnv_size:
            svtype = 'del' + svtype
        elif min_bkpt_cnv_size <= -dist5 < max_bkpt_cnv_size:
            svtype = 'dup' + svtype

    # Couldn't match the CNV
    else:
        return classify_0_cnv(FF, RR)

    return svtype, [cnv]


def filter_multiple_cnvs(FF, RR, cnvs, min_frac=0.5):
    """
    For cases with 3 or more overlapping CNV, try to remove spurious hits
    by forcing 50% reciprocal with any of the possible CNV intervals. If
    multiple CNVs are present for a candidate interval (e.g. 5' deletion),
    their coordinates are merged by taking the median.

    Parameters
    ----------
    FF : pysam.VariantRecord
        FF inversion breakpoint.
    RR : pysam.VariantRecord.
        RR inversion breakpoint
    cnvs : list of pysam.VariantRecord
        List of CNVs overlapping breakpoints.

    Returns
    -------
    cnvs : list of pysam.VariantRecord
        Filtered and merged CNVs.
    """

    # Identify eligible intervals for flanking CNV, defined by inv breakpoints
    ff_end = FF.info['END2'] if FF.info['SVTYPE'] == 'BND' else FF.stop
    rr_end = RR.info['END2'] if RR.info['SVTYPE'] == 'BND' else RR.stop
    del5 = (FF.pos, RR.pos)
    del3 = (ff_end, rr_end)
    dup5 = (RR.pos, FF.pos)
    dup3 = (rr_end, ff_end)

    # Determine if CNV supports 5' CNV, 3' CNV, spans event, or fails overlap
    def _test_overlap(cnv):
        svtype = cnv.info['SVTYPE']
        cnv_end = cnv.info['END2'] if cnv.info['SVTYPE'] == 'BND' else cnv.stop
        if svtype == 'DEL':
            frac5 = svu.reciprocal_overlap(cnv.pos, cnv_end, *del5)
            frac3 = svu.reciprocal_overlap(cnv.pos, cnv_end, *del3)
        else:
            frac5 = svu.reciprocal_overlap(cnv.pos, cnv_end, *dup5)
            frac3 = svu.reciprocal_overlap(cnv.pos, cnv_end, *dup3)

        if frac5 >= min_frac and frac3 >= min_frac:
            return svtype + '_53'
        elif frac5 >= min_frac:
            return svtype + '_5'
        elif frac3 >= min_frac:
            return svtype + '_3'
        else:
            return 'no_hit'

    # Collect CNV of same overlap type (e.g., 5' deletion) for merging
    cnvlists = defaultdict(list)
    for cnv in cnvs:
        cnvtype = _test_overlap(cnv)
        if cnvtype == 'no_hit':
            continue
        cnvlists[cnvtype].append(cnv)

    # Keep original CNV if only one present,
    # else merge by taking median start/end
    cnvs = []
    for overlap in cnvlists.keys():
        if len(cnvlists[overlap]) == 1:
            cnvs.append(cnvlists[overlap][0])
        else:
            cnvlist = cnvlists[overlap]
            # Overwrite values in first VariantRecord
            # (can't add list of IDs yet)
            merged_cnv = cnvlist[0].copy()

            # get coordinates
            start = int(np.median([c.pos for c in cnvlist]))
            end = int(np.median([c.info['END2'] if c.info['SVTYPE'] == 'BND' else c.stop for c in cnvlist]))
            name = '__'.join([c.id for c in cnvlist])

            merged_cnv.pos = start
            merged_cnv.stop = start + 1 if merged_cnv.info['SVTYPE'] == 'BND' else end
            merged_cnv.id = name
            if merged_cnv.info['SVTYPE'] == 'BND':
                merged_cnv.info['END2'] = end

            cnvs.append(merged_cnv)

    return sorted(cnvs, key=lambda record: record.pos)


def classify_0_cnv(FF, RR, min_bkpt_cnv_size=300):
    """
    Classify the cxSV class of a pair of inv bkpts with no associated CNVs.

    Matches breakpoint ordering to a known class, then tests if breakpoint
    distances are sufficient to call a CNV that was missed by integration
    pipeline.

    Parameters
    ----------
    FF : pysam.VariantRecord
    RR : pysam.VariantRecord
    min_bkpt_cnv_size : int, optional
        Minimum distance between breakpoints to call flanking CNV.

    Returns
    -------
    svtype : str
    """

    # Identify breakpoint ordering
    order = breakpoint_ordering(FF, RR, mh_buffer=10)

    # Get end coordinates
    ff_end = FF.info['END2'] if FF.info['SVTYPE'] == 'BND' else FF.stop
    rr_end = RR.info['END2'] if RR.info['SVTYPE'] == 'BND' else RR.stop

    # Check for flanking deletions around a "simple" inversion
    if order == 'SIMPLE/DEL':
        start_dist = RR.pos - FF.pos
        end_dist = rr_end - ff_end
        if start_dist < min_bkpt_cnv_size and end_dist < min_bkpt_cnv_size:
            svtype = 'INV'
        elif start_dist >= min_bkpt_cnv_size and end_dist < min_bkpt_cnv_size:
            svtype = 'delINV'
        elif start_dist < min_bkpt_cnv_size and end_dist >= min_bkpt_cnv_size:
            svtype = 'INVdel'
        else:
            svtype = 'delINVdel'

    # Check for flanking dups
    elif order == 'dupINVdup':
        start_dist = FF.pos - RR.pos
        end_dist = ff_end - rr_end

        if start_dist >= min_bkpt_cnv_size and end_dist >= min_bkpt_cnv_size:
            svtype = 'dupINVdup'
        elif start_dist >= min_bkpt_cnv_size:
            svtype = 'DUP5/INS3'
        elif end_dist >= min_bkpt_cnv_size:
            svtype = 'DUP3/INS5'
        else:
            svtype = 'UNK'

    # DUP5/INS3, DUP3/INS5, and UNK don't require add'l check
    else:
        svtype = order

    return svtype, []


def classify_complex_inversion(FF, RR, cnvs):
    """
    Classify the complex class of an inversion and associated CNVs.

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
    cnvs : list of pysam.VariantRecord
        CNVs represented in resolved variant structure
    """

    raw_cnvs = cnvs

    if len(cnvs) > 2:
        cnvs = filter_multiple_cnvs(FF, RR, cnvs)

    # Get original CNV records after merging and filtering
    filtered_ids = [s for r in cnvs for s in r.id.split('__')]
    pass_raw_cnvs = [cnv for cnv in raw_cnvs if cnv.id in filtered_ids]

    if len(cnvs) == 0:
        return classify_0_cnv(FF, RR)
    elif len(cnvs) == 1:
        svtype, pass_raw_cnvs = classify_1_cnv(FF, RR, cnvs[0])
    elif len(cnvs) == 2:
        svtype, pass_raw_cnvs = classify_2_cnv(FF, RR, cnvs)
    else:
        return 'MULT_CNVS', pass_raw_cnvs

    if breakpoints_match(FF, RR, svtype, mh_buffer=10):
        return svtype, pass_raw_cnvs
    else:
        return 'COMPLEX_INS', pass_raw_cnvs
