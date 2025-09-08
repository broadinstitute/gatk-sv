#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""

"""

from collections import defaultdict
import numpy as np
import pandas as pd
import pysam
from .pesr_test import PESRTest, PESRTestRunner


class _DiscPair:
    def __init__(self, chrA, posA, strandA, chrB, posB, strandB, sample):
        self.chrA = chrA
        self.posA = int(posA)
        self.strandA = strandA
        self.chrB = chrB
        self.posB = int(posB)
        self.strandB = strandB
        self.sample = sample


class PETest(PESRTest):
    def __init__(self, discfile, common=False, window_in=50, window_out=500, medians=None):
        self.discfile = discfile
        self.window_in = window_in
        self.window_out = window_out
        self.common = common

        super().__init__(medians, common)

    def test_record(self, record, called, background):
        # Test SR support at all coordinates within window of start/end
        results = self.test(record, called, background)

        results = results.to_frame().transpose()

        # Clean up columns
        results['name'] = record.id
        results['bg_frac'] = results.called / \
            (results.background + results.called)
        results['bg_frac'] = results.bg_frac.fillna(0)
        cols = 'name log_pval called background bg_frac'.split()

        return results[cols]

    def test(self, record, called, background):
        """
        Test enrichment of discordant reads in a set of samples.

        Arguments
        ---------
        record : pysam.VariantFile
        window_in : int
        window_out : int
        called : list of str
            List of called samples to test
        background : list of str
            List of samples to use as background

        Returns
        -------
        called_median : float
        background_median : float
        log_pval : float
            Negative log10 p-value
        """

        # Load split counts.
        counts = self.load_counts(record, self.window_in, self.window_out)
        counts = self.normalize_counts(counts)
        return super().test(counts, called, background)

    def load_counts(self, record, window_in, window_out):
        """Load pandas DataFrame from tabixfile"""

        def _get_coords(pos, strand):
            if strand == '+':
                start, end = pos - window_out, pos + window_in
            else:
                start, end = pos - window_in, pos + window_out
            return start, end

        strandA, strandB = record.info['STRANDS']
        startA, endA = _get_coords(record.pos, strandA)
        end = record.info['END2'] if record.info['SVTYPE'] == 'BND' else record.stop
        startB, endB = _get_coords(end, strandB)

        # Add 1 because evidence is stored/indexed with 0-based coordinates
        region = '{0}:{1}-{2}'.format(record.chrom, startA + 1, endA + 1)

        try:
            pairs = self.discfile.fetch(region=region)
        except ValueError:
            pairs = []

        counts = defaultdict(int)
        i = 0
        for pair_record in pairs:
            pair = _DiscPair(*pair_record.split('\t'))
            if (i > 1000000):
                print(region)
                counts = defaultdict(int)
                break
            i += 1
            # Pairs were selected based on window around chrA;
            # just need to check chrB
            if pair.chrB != record.info['CHR2']:
                continue
            if not (startB <= pair.posB < endB):
                continue

            # Require pairs match breakpoint strand
            if pair.strandA != strandA or pair.strandB != strandB:
                continue

            counts[pair.sample] += 1

        counts = pd.DataFrame.from_dict({'count': counts})
        counts = counts.reset_index()
        counts = counts.rename(columns={'index': 'sample'})

        return counts


class PETestRunner(PESRTestRunner):
    def __init__(self, vcf, discfile, fout, n_background=160, common=False,
                 window_in=50, window_out=500, whitelist=None, blacklist=None,
                 medians=None, log=False, outlier_sample_ids=None, seed=42):
        """
        vcf : pysam.VariantFile
        discfile : pysam.TabixFile
        """
        self.petest = PETest(discfile, common=common, window_in=window_in,
                             window_out=window_out, medians=medians)
        self.fout = fout

        super().__init__(vcf, common, n_background, whitelist, blacklist, log, outlier_sample_ids, seed)

    def test_record(self, record):
        if not self._strand_check(record):
            counts = self.petest.null_score(null_val=np.nan)
        else:
            called, background = self.choose_background(record)
            counts = self.petest.test_record(record, called, background)

        counts = counts.rename(columns={'called': 'called_median',
                                        'background': 'bg_median'})
        counts.to_csv(self.fout, header=False, index=False,
                      sep='\t', na_rep='NA')

    @staticmethod
    def _strand_check(record):
        return ('STRANDS' in record.info.keys() and
                record.info['STRANDS'] in '++ +- -+ --'.split())
