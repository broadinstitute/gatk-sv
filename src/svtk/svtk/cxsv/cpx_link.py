#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""
Link complex records
"""

import itertools
from collections import deque
import numpy as np
import scipy.sparse as sps
import natsort
import svtk.utils as svu


def samples_overlap_records(recA, recB, called_samples_dict, upper_thresh=0.5, lower_thresh=0.5):
    if recA.id not in called_samples_dict:
        called_samples_dict[recA.id] = set(svu.get_called_samples(recA))
    if recB.id not in called_samples_dict:
        called_samples_dict[recB.id] = set(svu.get_called_samples(recB))
    return samples_overlap(called_samples_dict[recA.id], called_samples_dict[recB.id],
                           upper_thresh=upper_thresh, lower_thresh=lower_thresh)


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


def link_cpx(vcf, bkpt_window=300):
    """
    Parameters
    ----------
    vcfpath : str
        Path to breakpoint VCF
    """

    bt = svu.vcf2bedtool(vcf.filename, annotate_ins=False)

    # Identify breakpoints which overlap within specified window
    overlap = bt.window(bt, w=bkpt_window).saveas()

    # Exclude intersections where two DELs or two DUPs cluster together
    overlap = overlap.filter(lambda b: not (
        b.fields[4] == "DEL" and b.fields[10] == "DEL")).saveas()
    overlap = overlap.filter(lambda b: not (
        b.fields[4] == "DUP" and b.fields[10] == "DUP")).saveas()

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

    return clusters


def unify_list(list):
    out = []
    for i in list:
        if i not in out:
            out.append(i)
    return out


def link_cpx_V2(linked_INV, cpx_dist=2000):
    overlapping_inv = []
    called_samples_dict = {}
    for group in linked_INV:
        if len(group) > 1:
            for i, j in itertools.combinations(group, 2):
                if records_overlap(i, j) and samples_overlap_records(i, j, called_samples_dict):
                    overlapping_inv.append([i, j])
        else:
            overlapping_inv.append(group)
    cluster = []
    for inv in overlapping_inv:
        if len(inv) > 1:
            if abs(inv[1].pos - inv[0].pos) > cpx_dist and abs(inv[1].stop - inv[0].stop) > cpx_dist:
                if 'STRANDS' in inv[0].info.keys() and 'STRANDS' in inv[1].info.keys():
                    if inv[0].info['STRANDS'] != inv[1].info['STRANDS']:
                        if inv[0].pos < inv[1].pos < inv[0].stop < inv[1].stop \
                                or inv[1].pos < inv[0].pos < inv[1].stop < inv[0].stop:
                            cluster.append(inv)
        else:
            cluster.append(inv)
    return cluster


def close_enough(r1, r2, cpx_dist=2000):
    distA = np.abs(r1.pos - r2.pos)
    distB = np.abs(r1.stop - r2.stop)
    return distA < cpx_dist or distB < cpx_dist


def records_overlap(r1, r2):
    return r1.chrom == r2.chrom and not (r1.pos > r2.stop or r1.stop < r2.pos)
