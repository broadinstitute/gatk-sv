#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""

"""

import numpy as np
import scipy.stats as ss
import pandas as pd
import sys
from .pesr_test import PESRTest, PESRTestRunner


class SRTest(PESRTest):
    def __init__(self, countfile, common=False, window=50, ins_window=50, medians=None):
        self.countfile = countfile
        self.window = window
        self.ins_window = ins_window
        self.common = common

        super().__init__(medians, common)

    def test_record(self, record, called, background):
        # Test SR support at all coordinates within window of start/end
        resultA = self._test_coord(
            record.pos, record.info['STRANDS'][0], record.chrom,
            called, background, record.pos - self.window, record.pos + self.window)

        resultB = self._test_coord(
            record.stop, record.info['STRANDS'][1], record.info['CHR2'],
            called, background, record.stop - self.window, record.stop + self.window)

        posA = int(resultA.pos)
        posB = int(resultB.pos)
        log_pval_A = float(resultA.log_pval)
        log_pval_B = float(resultB.log_pval)
        if record.info['SVTYPE'] == 'INS':
            if posA > posB + self.ins_window or posA < posB - self.ins_window:
                # Invalid coordinates, need to re-optimize around the better coordinate within the insertion window size
                if log_pval_A >= log_pval_B:
                    # posA is better, so use posA as anchor and check for best posB in valid window
                    resultB = self._test_coord(
                        posA, record.info['STRANDS'][1], record.info['CHR2'],
                        called, background, posA - self.ins_window, posA + self.ins_window)
                else:
                    # vice versa
                    resultA = self._test_coord(
                        posB, record.info['STRANDS'][0], record.chrom,
                        called, background, posB - self.ins_window, posB + self.ins_window)
        elif record.chrom == record.info['CHR2']:
            # Check some corner cases for intrachromosomal variants
            if record.info['STRANDS'][0] == record.info['STRANDS'][1]:
                # If strands are the same, disallow case where posA = posB, which can happen for small variants
                # Note that the case posA > posB is allowed and corrected for in the SR coordinate rewriting step later
                if posA == posB and self.window > 0:
                    resultB = self._test_coord(
                        record.stop, record.info['STRANDS'][1], record.info['CHR2'],
                        called, background, record.stop - self.window, record.stop + self.window,
                        invalid_pos_list=[posA])
            elif posA >= posB:
                # Invalid coordinates, need to re-optimize around the better coordinate
                if log_pval_A >= log_pval_B:
                    # posA is better, so use posA as anchor and check for best posB in valid window
                    resultB = self._test_coord(
                        record.stop, record.info['STRANDS'][1], record.info['CHR2'],
                        called, background, posA, posA + self.window)
                else:
                    # vice versa
                    resultA = self._test_coord(
                        record.pos, record.info['STRANDS'][0], record.chrom,
                        called, background, posB - self.window, posB)

        resultA['coord'] = 'posA'
        resultB['coord'] = 'posB'
        results = pd.concat([resultA, resultB], ignore_index=True)
        # Add test for sum of posA and posB
        total = self._test_total(results)
        results = pd.concat([results, total], ignore_index=True)

        # Clean up columns
        results['name'] = record.id
        results['bg_frac'] = results.called / \
                             (results.background + results.called)
        results['bg_frac'] = results.bg_frac.fillna(0)
        cols = 'name coord pos log_pval called background bg_frac'.split()

        return results[cols]

    def test(self, chrom, pos, strand, called, background):
        """
        Test enrichment of clipped reads in a set of samples at a given coord.

        Arguments
        ---------
        chrom : str
        pos : int
        strand : str
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
        counts = self.load_counts(chrom, pos, strand)
        counts = self.normalize_counts(counts)

        return super().test(counts, called, background)

    def load_counts(self, chrom, pos, strand):
        """Load pandas DataFrame from tabixfile"""

        if pos > 0:
            # Add 1 because evidence is stored/indexed with 0-based coordinates
            region = '{0}:{1}-{1}'.format(chrom, pos + 1)
            try:
                lines = self.countfile.fetch(region)
            except ValueError:
                lines = []
        else:
            lines = []
        #  counts = io.StringIO('\n'.join([l for l in lines]))

        cols = 'chrom pos clip count sample'.split()
        #  dtypes = dict(chrom=str, pos=int, clip=str, count=int, sample=str)

        counts = pd.DataFrame.from_records(
            [l[:5] for l in lines], columns=cols)
        counts['count'] = counts['count'].astype(int)

        # Restrict to splits in orientation of interest
        clip = 'right' if strand == '+' else 'left'
        counts = counts.loc[counts['clip'] == clip].copy()

        return counts

    def _test_total(self, results):
        """Test enrichment of posA+posB"""
        total = results['called background'.split()].sum()
        pval = max(ss.poisson.cdf(total.background, total.called), sys.float_info.min)
        total['log_pval'] = np.abs(np.log10(pval))

        # format and add dummy metadata
        total = total.to_frame().transpose()
        total['coord'] = 'sum'
        total['pos'] = 0

        return total

    def _test_coord(self, coord, strand, chrom, samples, background, left_boundary, right_boundary,
                    invalid_pos_list=[]):
        """Test enrichment at all positions within window"""

        # Run SR test at each position
        results = []
        positions = [p for p in range(left_boundary, right_boundary + 1) if p not in invalid_pos_list]
        for pos in positions:
            result = self.test(chrom, pos, strand, samples, background)
            result = result.to_frame().transpose()
            result['pos'] = pos
            # make negative so it sorts correctly
            result['dist'] = -np.abs(pos - coord)
            results.append(result)

        results = pd.concat(results, ignore_index=True)

        # Choose most significant position, using distance to predicted
        # breakpoint as tiebreaker
        results = results.sort_values(['log_pval', 'dist'], ascending=False)
        best = results.iloc[0].to_frame().transpose()

        return best


class SRTestRunner(PESRTestRunner):
    def __init__(self, vcf, countfile, fout, n_background=160, common=False, window=100, ins_window=50,
                 whitelist=None, blacklist=None, medians=None, log=False, outlier_sample_ids=None, seed=42):
        """
        vcf : pysam.VariantFile
        countfile : pysam.TabixFile
        fout : writable file
        n_background : int
        window : int
        ins_window : int
        whitelist : list of str
        blacklist : list of str
        """
        self.srtest = SRTest(countfile, common=common, window=window, ins_window=ins_window, medians=medians)
        self.fout = fout

        super().__init__(vcf, common, n_background, whitelist, blacklist, log, outlier_sample_ids, seed)

    def test_record(self, record):
        called, background = self.choose_background(record)
        counts = self.srtest.test_record(record, called, background)
        counts = counts.rename(columns={'called': 'called_median',
                                        'background': 'bg_median'})
        counts.to_csv(self.fout, header=False, index=False, sep='\t')
