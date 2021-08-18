#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2017 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.

"""
Link complex records
"""

import itertools
from collections import deque
import numpy as np
import scipy.sparse as sps
import natsort
import pysam
import svtk.utils as svu


def samples_overlap_records(recA, recB, upper_thresh=0.5, lower_thresh=0.5):
    samplesA = set(svu.get_called_samples(recA))
    samplesB = set(svu.get_called_samples(recB))
    return samples_overlap(samplesA, samplesB, upper_thresh=upper_thresh, lower_thresh=lower_thresh)


def samples_overlap(samplesA, samplesB, upper_thresh=0.5, lower_thresh=0.5):
    """
    Report if the samples called in two VCF records overlap sufficiently.
    The fraction of each record's samples which are shared with the other
    record is calculated. The record with a greater fraction of shared samples
    must exceed the upper threshold AND the record with a lesser fraction of
    shared samples must exceed the lower threshold. This is intended to
    maximize sensitivity in rare variants with a false negative in one
    breakpoint.
    Parameters
    ----------
    recA : pysam.VariantRecord
    recB : pysam.VariantRecord
    upper_thresh : float, optional
        Minimum sample overlap in record with greater overlap
    lower_thresh : float, optional
        Minimum sample overlap in record with lesser overlap
    Returns
    -------
    samples_overlap : bool
        Samples shared between records meet required thresholds.
    """
    # Compute fraction of each record's samples which are shared
    if len(samplesA) > 0 and len(samplesB) > 0:
        shared = samplesA & samplesB
        fracA = len(shared) / len(samplesA)
        fracB = len(shared) / len(samplesB)
        min_frac, max_frac = sorted([fracA, fracB])
    else:
        min_frac, max_frac = [0, 0]
    return min_frac >= lower_thresh and max_frac >= upper_thresh


def extract_breakpoints(vcf, bkpt_idxs):
    """
    Extract all VCF records in list of IDs
    (Assumes VCF is sorted by variant ID)
    Parameters
    ----------
    vcfpath : str
        Path to VCF
    bkpt_idxs : dict of {str : int}
        Mapping of variant IDs to array index
    Returns
    -------
    bkpts : list of pysam.VariantRecord
    """

    #  vcf = pysam.VariantFile(vcfpath)
    n_bkpts = len(bkpt_idxs)
    bkpts = np.empty(n_bkpts, dtype=object)

    for record in vcf:
        idx = bkpt_idxs.get(record.id)
        if idx is not None:
            bkpts[idx] = record

    return bkpts


def link_cpx(vcf, bkpt_window=300, cpx_dist=2000):
    """
    Parameters
    ----------
    vcfpath : str
        Path to breakpoint VCF
    """

    bt = svu.vcf2bedtool(vcf.filename, annotate_ins=False)

    # Identify breakpoints which overlap within specified window
    overlap = bt.window(bt, w=bkpt_window).saveas()

    # Exclude self-hits
    #  overlap = overlap.filter(lambda b: b.fields[3] != b.fields[9]).saveas()

    # Exclude intersections where two DELs or two DUPs cluster together
    # cnvtypes = 'DEL DUP'.split()
    overlap = overlap.filter(lambda b: not (
        b.fields[4] == "DEL" and b.fields[10] == "DEL")).saveas()
    overlap = overlap.filter(lambda b: not (
        b.fields[4] == "DUP" and b.fields[10] == "DUP")).saveas()

    # # Exclude intersections with annotated mobile elements (rather than BNDs)
    # overlap = overlap.filter(lambda b: b.fields[4] is not re.match(re.compile('INS\:ME\:*'), b.fields[4])).saveas()

    # Get linked variant IDs
    links = [(b[3], b[9]) for b in overlap.intervals]
    linked_IDs = natsort.natsorted(set(itertools.chain.from_iterable(links)))
    linked_IDs = np.array(linked_IDs)

    # Map variant IDs to indices
    bkpt_idxs = {ID: i for i, ID in enumerate(linked_IDs)}
    indexed_links = np.array([(bkpt_idxs[a], bkpt_idxs[b]) for a, b in links])

    # Extract VariantRecords corresponding to breakpoints
    n_bkpts = len(linked_IDs)
    bkpts = extract_breakpoints(vcf, bkpt_idxs)

    # Build called sample index
    # Get lists of called samples for each record
    sample_sets_dict = {idx: set(svu.get_called_samples(bkpts[idx])) for idx in set(indexed_links.flatten().tolist())}

    # Exclude wildly disparate overlaps
    # Build sparse graph from links
    G = sps.eye(n_bkpts, dtype=np.uint16, format='lil')
    for i, j in indexed_links:
        if (samples_overlap(sample_sets_dict[i], sample_sets_dict[j]) and
                close_enough(bkpts[i], bkpts[j])):
            G[i, j] = 1

    # Generate lists of clustered breakpoints
    n_comp, comp_list = sps.csgraph.connected_components(G)
    clusters = [deque() for x in range(n_comp)]
    for i, c_label in enumerate(comp_list):
        clusters[c_label].append(bkpts[i])

    # # Remove clusters of only CNV - leftover from shared sample filtering
    # def _ok_cluster(cluster):
    #     ok = any([record.info['SVTYPE'] not in cnvtypes for record in cluster])
    #     return ok

    # clusters = [c for c in clusters if _ok_cluster(c)]
    #  clusters = [c for c in clusters if len(c) > 1]

    return clusters


def unify_list(list):
    out = []
    for i in list:
        if i not in out:
            out.append(i)
    return out


def CNV_readin_from_resolved_vcf(resolved_name, inv_intervals):
    resolved_f = pysam.VariantFile(resolved_name, 'r')
    # rec_a = 0
    out = []
    for i in resolved_f:
        for j in inv_intervals:
            if i.chrom == j[0]:
                if (i.pos - j[1]) * (i.pos - j[2]) < 0 or (i.stop - j[1]) * (i.stop - j[2]) < 0:
                    if i.info['SVTYPE'] in ['DEL', 'DUP']:
                        out.append(i)
    resolved_f.close()
    return out


def link_cpx_V2(linked_INV, resolve_CNV, cpx_dist=2000):
    linked_INV_V2 = []
    for group in linked_INV:
        if len(group) > 1:
            for i in group:
                for j in group:
                    if ro_calu(i, j) > 0 and samples_overlap_records(i, j):
                        linked_INV_V2.append([i, j])
        else:
            linked_INV_V2.append([group[0]])
    inv_intervals = []
    for i in linked_INV_V2:
        if len(i) > 1:
            tmp = [i[0].chrom]
            for j in i:
                tmp += [j.pos, j.stop]
            inv_intervals.append(
                [tmp[0], min(unify_list(tmp[1:])), max(unify_list(tmp[1:]))])
        else:
            inv_intervals.append([i[0].chrom, i[0].pos, i[0].stop])
    inv_intervals = sorted(unify_list(inv_intervals))
    # out_rec = unify_list(CNV_readin_from_resolved_vcf(resolved_name,inv_intervals) + CNV_readin_from_resolved_vcf(unresolved_name,inv_intervals))
    out_rec = resolve_CNV
    cluster = []
    for i in linked_INV_V2:
        if len(i) > 1:
            if abs(i[1].pos - i[0].pos) > cpx_dist and abs(i[1].stop - i[0].stop) > cpx_dist:
                if 'STRANDS' in i[0].info.keys() and 'STRANDS' in i[1].info.keys():
                    if sorted(unify_list([i[0].info['STRANDS'], i[1].info['STRANDS']])) == ['++', '--']:
                        if i[0].pos < i[1].pos < i[0].stop < i[1].stop or i[1].pos < i[0].pos < i[1].stop < i[0].stop:
                            cpx_intervals = [[i[0].chrom, sorted([i[0].pos, i[0].stop, i[1].pos, i[1].stop])[0], sorted([i[0].pos, i[0].stop, i[1].pos, i[1].stop])[1]], [
                                i[0].chrom, sorted([i[0].pos, i[0].stop, i[1].pos, i[1].stop])[2], sorted([i[0].pos, i[0].stop, i[1].pos, i[1].stop])[3]]]
                            CNV_close = [j for j in out_rec if ro_calu_interval([j.chrom, j.pos, j.stop], cpx_intervals[0]) > .5 and abs(
                                j.pos - cpx_intervals[0][1]) < cpx_dist and abs(j.stop - cpx_intervals[0][2]) < cpx_dist]
                            CNV_close += [j for j in out_rec if ro_calu_interval([j.chrom, j.pos, j.stop], cpx_intervals[1]) > .5 and abs(
                                j.pos - cpx_intervals[1][1]) < cpx_dist and abs(j.stop - cpx_intervals[1][2]) < cpx_dist]
                            cluster.append(CNV_close + i)
        else:
            cluster.append(i)
    return cluster


def link_inv(vcf, bkpt_window=300, cpx_dist=2000):
    bt = svu.vcf2bedtool(vcf.filename, annotate_ins=False)
    overlap = bt.window(bt, w=bkpt_window).saveas()
    overlap = overlap.filter(lambda b: not (
        b.fields[4] == "DEL" and b.fields[10] == "DEL")).saveas()
    overlap = overlap.filter(lambda b: not (
        b.fields[4] == "DUP" and b.fields[10] == "DUP")).saveas()
    links = [(b[3], b[9]) for b in overlap.intervals]
    linked_IDs = natsort.natsorted(set(itertools.chain.from_iterable(links)))
    linked_IDs = np.array(linked_IDs)
    bkpt_idxs = {ID: i for i, ID in enumerate(linked_IDs)}
    indexed_links = np.array([(bkpt_idxs[a], bkpt_idxs[b]) for a, b in links])
    n_bkpts = len(linked_IDs)
    bkpts = extract_breakpoints(vcf, bkpt_idxs)
    # Exclude wildly disparate overlaps
    G = sps.eye(n_bkpts, dtype=np.uint16, format='lil')
    for i, j in indexed_links:
        if (ro_calu(bkpts[i], bkpts[j]) > 0 and samples_overlap_records(bkpts[i], bkpts[j])):
            G[i, j] = 1
    # Generate lists of clustered breakpoints
    n_comp, comp_list = sps.csgraph.connected_components(G)
    clusters = [deque() for x in range(n_comp)]
    for i, c_label in enumerate(comp_list):
        clusters[c_label].append(bkpts[i])
    return clusters


def close_enough(r1, r2, cpx_dist=2000):
    distA = np.abs(r1.pos - r2.pos)
    distB = np.abs(r1.stop - r2.stop)
    return distA < cpx_dist or distB < cpx_dist


def ro_calu(r1, r2):
    out = 0
    if not r1.chrom == r2.chrom:
        out = 0
    elif r1.pos > r2.stop or r1.stop < r2.pos:
        out = 0
    else:
        maxval = max([r1.stop - r1.pos, r2.stop - r2.pos])
        if maxval > 0:
            out = (sorted([r1.pos, r2.pos, r1.stop, r2.stop])[
                   2] - sorted([r1.pos, r2.pos, r1.stop, r2.stop])[1]) / maxval
        else:
            out = 0
    return out


def ro_calu_interval(r1, r2):
    out = 0
    if not r1[0] == r2[0]:
        out = 0
    elif r1[1] > r2[2] or r1[2] < r2[1]:
        out = 0
    else:
        out = (sorted(r1[1:] + r2[1:])[2] - sorted(r1[1:] +
                                                   r2[1:])[1]) / max([r1[2] - r1[1], r2[2] - r2[1]])
    return out
