#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""

"""

import argparse
from collections import defaultdict
import pysam
from svtk.genomeslink import GenomeSLINK, GSNode
import svtk.utils as svu
from statistics import median
import numpy as np
import datetime


class DiscPair(GSNode):
    def __init__(self, chrA, posA, strandA, chrB, posB, strandB, sample):
        self.strandA = strandA
        self.strandB = strandB
        self.sample = sample
        super().__init__(chrA, posA, chrB, posB)

    @property
    def is_inversion(self):
        return (self.chrA == self.chrB) and (self.strandA == self.strandB)

    def __str__(self):
        e = '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n'
        return e.format(self.chrA, self.posA, self.strandA,
                        self.chrB, self.posB, self.strandB, self.sample)


def match_cluster(record, cluster, dist=300):
    """
    Determine whether DiscPair cluster matches a VCF record of interest.

    Checks if pairs exist which match the record's strandedness, then checks
    whether the min/max coord of these pairs is within a specified distance
    of the record's coordinates.

    Arguments
    ---------
    record : pysam.VariantRecord
    cluster : list of DiscPair
    dist : int, optional

    Returns
    -------
    match : bool
    """
    r_start = record.pos
    r_end = record.info['END2'] if record.info['SVTYPE'] == 'BND' else record.stop

    if record.info['STRANDS'] == '++':
        # Choose max start/end of ++ pairs
        c_start = max((p.posA for p in cluster if (
            p.strandA == '+')), default=None)
        if c_start is None:
            # If no ++ pairs present, return False
            return False
        c_end = max(p.posB for p in cluster if (p.strandB == '+'))
    elif record.info['STRANDS'] == '--':
        # Choose min start/end of -- pairs
        c_start = min((p.posA for p in cluster if (
            p.strandA == '-')), default=None)
        if c_start is None:
            # If no -- pairs present, return False
            return False
        c_end = min(p.posB for p in cluster if (p.strandB == '-'))
    else:
        strands = record.info['STRANDS']
        raise Exception('Invalid inversion orientation: {0}'.format(strands))

    # Test if cluster start/end are sufficiently close to record start/end
    return abs(r_start - c_start) < dist and abs(r_end - c_end) < dist


def rescan_single_ender(record, pe, min_support=4, window=1000, dist=300,
                        min_frac_samples=0.5, pe_blacklist=None, max_samples=40,
                        quiet=False, min_span=50):
    """
    Test if a putative single-ender inversion has support from other strand.

    Selects discordant pairs in the neighborhood of the original record, then
    clusters them together. If enough samples have sufficient paired-end
    evidence supporting the opposite strand, we have found support for the
    other end of the record.

    Arguments
    ---------
    record : pysam.VariantRecord
    pe : pysam.TabixFile
        Scraped discordant pair metadata
    min_support : int, optional
        Number of pairs required to count a sample as supported
    window : int, optional
        Window around record start to search for pairs
    dist : int, optional
        Clustering distance for fetched pairs
    min_frac_samples : float, optional
        Fraction of called samples required to have opposite strand support in
        order to call the record as having both strands present. If 0, only one
        sample will be required.
    pe_blacklist : pysam.TabixFile, optional
        Blacklisted genomic regions. Anomalous pairs in these regions will be
        removed prior to clustering.
    quiet : boolean, optional
        Do not print status updates
    min_span : int, optional
        Minimum distance spanned between discordant read mapping positions in
        newly identified candidate breakpoints


    Returns
    -------
    opposite : pysam.VariantRecord
        Record corresponding to the pairs found supporting the other strand.
        None if no such pairs found
    """

    # Print statement that single ender rescan has been attempted
    if not quiet:
        now = datetime.datetime.now()
        print('svtk resolve @ ' + now.strftime("%H:%M:%S") + ': ' +
              'single-ender rescan procedure started for ' +
              record.id)

    # Select pairs nearby record
    search_start = max(0, record.pos - window)
    search_end = max(search_start + 1, record.pos + window)
    pairs = pe.fetch(
        '{0}:{1}-{2}'.format(record.chrom, search_start, search_end))
    pairs = [DiscPair(*p.split()) for p in pairs]

    # To protect against wasting time on particularly messy loci not captured
    # in the blacklist, automatically fail site if total number of discordant
    # pairs is > samples * min_support
    all_samples = record.samples.keys()
    print(len(pairs), flush=True)
    print(len(all_samples) * min_support, flush=True)
    if len(pairs) > len(all_samples) * min_support:
        return record, None

    # Subset to only inversion pairs
    pairs = [p for p in pairs if p.is_inversion]

    # If median number of pairs per sample not called in the original record
    # > min_support, fail record. Otherwise, keep going.
    called = svu.get_called_samples(record)
    if len(called) < len(all_samples):
        # Count number of pairs per sample for all samples
        sample_support_precluster = defaultdict(int)
        for pair in pairs:
            sample_support_precluster[pair.sample] += 1
        # compute median pairs per uncalled sample
        median_pairs_not_called = median(
            sample_support_precluster.get(s, 0) for s in all_samples if s not in called
        )
        if median_pairs_not_called > min_support:
            return record, None
        # Randomly subset pairs from all samples to max_samples
    if len(called) > max_samples:
        np.random.seed(2)  # arbitrary fixed seed for reproducibility
        called = np.random.choice(called, max_samples, replace=False).tolist()

    # Restrict to inversions in samples called in VCF record
    pairs = [p for p in pairs if p.sample in called and p.is_inversion]

    # Cluster pairs
    slink = GenomeSLINK(pairs, dist, blacklist=pe_blacklist)
    # choose the largest cluster
    cluster = max(
        (c for c in slink.cluster() if match_cluster(record, c, window)),
        key=len, default=None
    )
    if cluster is None:
        # if no clusters, fail site
        return record, None

    # Select clustered pairs which support the opposite strand as the record
    missing_strand = '+' if record.info['STRANDS'] == '--' else '-'
    supporting_pairs = [p for p in cluster if p.strandA == missing_strand]

    # Check span of supporting pairs from best cluster
    minA = round(np.percentile([p.posA for p in cluster], 10))
    maxA = round(np.percentile([p.posA for p in cluster], 90))
    spanA = maxA - minA

    minB = round(np.percentile([p.posB for p in cluster], 10))
    maxB = round(np.percentile([p.posB for p in cluster], 90))
    spanB = maxB - minB

    if min(spanA, spanB) < min_span:
        return record, None

    # Count number of supporting pairs in each called sample
    sample_support = defaultdict(int)
    for pair in supporting_pairs:
        sample_support[pair.sample] += 1

    # If enough samples were found to have support, make new variant record
    n_supported_samples = sum(sample_support[s] > min_support for s in called)
    if n_supported_samples >= min_frac_samples * len(called):
        opp_strand = make_new_record(supporting_pairs, record)

        same_strand_pairs = [p for p in cluster if p.strandA != missing_strand]
        same_strand = make_new_record(same_strand_pairs, record, True)
        same_strand.id = record.id

        # Print statement that single ender rescan has been successful
        if not quiet:
            now = datetime.datetime.now()
            print('svtk resolve @ ' + now.strftime("%H:%M:%S") + ': ' +
                  'single-ender rescan successful for ' +
                  record.id)

        return same_strand, opp_strand
    else:
        return record, None


def make_new_record(pairs, old_record, retain_algs=False):
    record = old_record.copy()

    record.id = record.id + '_OPPSTRAND'

    # Take third quartile of + read positions for +/+ breakpoints
    # Take first quartile of - read positions for -/- breakpoints
    if pairs[0].strandA == '+':
        record.pos = round(np.percentile([p.posA for p in pairs], 90), 0)
        end = round(np.percentile([p.posB for p in pairs], 90), 0)
        record.info['STRANDS'] = '++'
    else:
        record.pos = round(np.percentile([p.posA for p in pairs], 10), 0)
        end = round(np.percentile([p.posB for p in pairs], 10), 0)
        record.info['STRANDS'] = '--'

    if record.info.get('SVTYPE') == 'BND':
        record.info['END2'] = end
        record.stop = record.pos
    else:
        record.stop = end

    record.info['SVLEN'] = end - record.pos
    if retain_algs:
        old_algs = list(record.info['ALGORITHMS'])
        old_algs.append('rescan')
        record.info['ALGORITHMS'] = tuple(old_algs)

    return record


def rescan_single_enders(vcf, pe, min_support=4, window=500, pe_blacklist=None):
    for record in vcf:
        rescan_single_ender(record, pe, min_support, window,
                            pe_blacklist=pe_blacklist)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf', help='Single enders')
    parser.add_argument('pairs', help='Scraped discordant pair file.')
    parser.add_argument('--min-rescan-pe-support', type=int, default=4,
                        help='Minumum discordant pairs required during '
                        'single-ender rescan ')
    parser.add_argument('--window', type=int, default=500, help='Window around '
                        'single ender coordinates to search for pairs')
    parser.add_argument('-x', '--pe-blacklist', metavar='BED.GZ',
                        default=None, help='Tabix indexed bed of blacklisted '
                        'regions. Any anomalous pair falling inside one '
                        'of these regions is excluded from PE rescanning.')

    args = parser.parse_args()

    vcf = pysam.VariantFile(args.vcf)
    pe = pysam.TabixFile(args.pairs)
    blacklist = pysam.TabixFile(args.pe_blacklist)

    rescan_single_enders(vcf, pe, args.min_rescan_pe_support, args.window,
                         pe_blacklist=blacklist)


if __name__ == '__main__':
    main()
