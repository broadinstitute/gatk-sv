#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""

"""

import numpy as np
import scipy.stats as ss
import pandas as pd
from .pesr_test import PESRTest, PESRTestRunner


class SRTest(PESRTest):
    def __init__(self, countfile, common=False, window=50, medians=None):
        self.countfile = countfile
        self.window = window
        self.common = common

        super().__init__(medians, common)

    def test_record_rec(self, record, called, background):
        # Test SR support at all coordinates within window of start/end
        results = []
        for coord in 'posA posB'.split():
            result = self._test_coord(record, coord, called, background)
            result['coord'] = coord
            results.append(result)
        results = pd.concat(results, ignore_index=True)

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

    def test_record(self, record, called, background):
        # Test SR support at all coordinates within window of start/end
        resultA = self._test_coord(
            record, 'posA', record.chrom, called, background)
        resultA['coord'] = 'posA'

        resultB = self._test_coord(
            record, 'posB', record.info['CHR2'], called, background)
        resultB['coord'] = 'posB'

        if (record.chrom == record.info['CHR2']) and (int(resultA.pos) >= int(resultB.pos)):
            rec = 1
            while True:
                result = self._test_coord_V2(
                    record, 'posB', record.chrom, called, background, rec)
                if int(result.pos) > int(resultA.pos.iloc[0]):
                    result['coord'] = 'posB'
                    resultB = result
                    break
                rec = rec + 1

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

    def _test_coord_V2(self, record, coord, chrom, samples, background, optimal=0):
        """Test enrichment at all positions within window"""
        if coord == 'posA':
            coord, strand = record.pos, record.info['STRANDS'][0]
        else:
            coord, strand = record.stop, record.info['STRANDS'][1]

        # Run SR test at each position
        results = []
        for pos in range(coord - self.window, coord + self.window + 1):
            result = self.test(chrom, pos, strand, samples, background)
            result = result.to_frame().transpose()
            result['pos'] = pos
            result['dist'] = np.abs(pos - coord)
            results.append(result)

        results = pd.concat(results, ignore_index=True)

        # Choose most significant position, using distance to predicted
        # breakpoint as tiebreaker
        results = results.sort_values(['log_pval', 'dist'], ascending=False)
        if optimal < len(results):
            best = results.iloc[optimal].to_frame().transpose()
        else:
            best = results.iloc[-1].to_frame().transpose()
            best['log_pval'] = 0
            best['pos'] = record.stop + self.window + 1
        return best

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
            region = '{0}:{1}-{1}'.format(chrom, pos)
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
        pval = ss.poisson.cdf(total.background, total.called)
        total['log_pval'] = np.abs(np.log10(pval))

        # format and add dummy metadata
        total = total.to_frame().transpose()
        total['coord'] = 'sum'
        total['pos'] = 0

        return total

    def _test_coord(self, record, coord, chrom, samples, background):
        """Test enrichment at all positions within window"""
        if coord == 'posA':
            coord, strand = record.pos, record.info['STRANDS'][0]
        else:
            coord, strand = record.stop, record.info['STRANDS'][1]

        # Run SR test at each position
        results = []
        for pos in range(coord - self.window, coord + self.window + 1):
            result = self.test(chrom, pos, strand, samples, background)
            result = result.to_frame().transpose()
            result['pos'] = pos
            result['dist'] = np.abs(pos - coord)
            results.append(result)

        results = pd.concat(results, ignore_index=True)

        # Choose most significant position, using distance to predicted
        # breakpoint as tiebreaker
        results = results.sort_values(['log_pval', 'dist'], ascending=False)
        best = results.iloc[0].to_frame().transpose()

        return best


class SRTestRunner(PESRTestRunner):
    def __init__(self, vcf, countfile, fout, n_background=160, common=False, window=100,
                 whitelist=None, blacklist=None, medians=None, log=False):
        """
        vcf : pysam.VariantFile
        countfile : pysam.TabixFile
        fout : writable file
        n_background : int
        window : int
        whitelist : list of str
        blacklist : list of str
        """
        self.srtest = SRTest(countfile, common, window, medians=medians)
        self.fout = fout

        super().__init__(vcf, common, n_background, whitelist, blacklist, log)

    def test_record(self, record):
        called, background = self.choose_background(record)
        counts = self.srtest.test_record(record, called, background)
        counts = counts.rename(columns={'called': 'called_median',
                                        'background': 'bg_median'})
        counts.to_csv(self.fout, header=False, index=False, sep='\t')
