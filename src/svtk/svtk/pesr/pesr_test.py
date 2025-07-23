#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""

"""

import sys
import datetime
import numpy as np
import scipy.stats as ss
import pandas as pd
import svtk.utils as svu


class PESRTest:
    def __init__(self, medians=None, common=False):
        self.medians = medians
        self.common = common

    def test(self, counts, samples, background):
        """
        Test enrichment of clipped reads in a set of samples at a given coord.

        Arguments
        ---------
        chrom : str
        pos : int
        strand : str
        samples : list of str
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

        # Restrict to called or background samples
        counts = counts.loc[counts['sample'].isin(samples + background)].copy()

        # Return null score if no eligible clipped reads present
        if counts.shape[0] == 0:
            return self.null_score()

        # Add called and background samples with no observed clipped reads
        counts = counts.set_index('sample')['count']\
                       .reindex(samples + background)\
                       .fillna(0).reset_index()

        # Label samples
        counts['is_called'] = counts['sample'].isin(samples)

        # Calculate enrichment
        result = counts.groupby('is_called')['count'].median()

        # Fill 0 if called in all samples
        result = result.reindex([True, False]).fillna(0)
        result.index = ['called', 'background']
        if self.common != "False":
            if len(samples) > len(background):
                result.background = 0.0
        pval = max(ss.poisson.cdf(result.background, result.called), sys.float_info.min)
        result['log_pval'] = np.abs(np.log10(pval))

        return result

    def normalize_counts(self, counts, target_cov=60):
        if self.medians is None:
            return counts

        counts = pd.merge(counts, self.medians, on='sample', how='left')
        counts['norm_count'] = counts['count'] * \
            target_cov / counts['median_cov']
        counts['count'] = counts['norm_count'].astype(float).round()
        counts.drop(['norm_count', 'median_cov'], axis=1, inplace=True)

        return counts

    def null_score(self, null_val=0.0):
        """Null score if no clipped reads observed"""
        score = pd.Series([null_val] * 3,
                          ['background', 'called', 'log_pval']).rename('count')
        score.index.name = 'status'

        return score


class PESRTestRunner:
    def __init__(self, vcf, common=False, n_background=160, whitelist=None, blacklist=None,
                 log=False, outlier_sample_ids=None, seed=0):
        self.vcf = vcf

        self.common = common
        self.samples = list(vcf.header.samples)
        self.n_background = n_background

        self.whitelist = whitelist if whitelist else self.samples
        self.blacklist = blacklist if blacklist else []

        self.log = log

        np.random.seed(seed)

        outlier_samples = set()
        if outlier_sample_ids:
            with open(outlier_sample_ids, 'r') as f:
                outlier_samples = set([line.strip() for line in f])
        self.outlier_sample_ids = outlier_samples

    def run(self):
        if self.log:
            start = datetime.datetime.now()

        for i, record in enumerate(self.vcf):
            t0 = datetime.datetime.now()
            self.test_record(record)
            t1 = datetime.datetime.now()

            if self.log:
                n_records = i + 1
                var_time = (t1 - t0).total_seconds()
                total_time = (t1 - start).total_seconds()
                hours, remainder = divmod(total_time, 3600)
                minutes, seconds = divmod(remainder, 60)

                msg = ('%d variants processed. '
                       'Time to process last variant: %0.2f seconds. '
                       'Total time elapsed: %d hours, %d minutes, %0.2f seconds.')
                msg = msg % (n_records, var_time, int(
                    hours), int(minutes), seconds)
                sys.stderr.write(msg + '\n')

    def test_record(self, record):
        called, background = self.choose_background(record)

    def choose_background(self, record, whitelist=None, blacklist=None):
        # Select called and background samples
        called = svu.get_called_samples(record)
        background = [s for s in self.samples if s not in called]

        # Create non-outlier sample lists
        non_outlier_called = [s for s in called if s not in self.outlier_sample_ids]
        non_outlier_background = [s for s in background if s not in self.outlier_sample_ids]

        # Exclude outlier samples only if non-outlier samples exist
        if len(non_outlier_called) > 0:
            called = non_outlier_called
        if len(non_outlier_background) > 0:
            background = non_outlier_background

        # Permit override of specified white/blacklists
        whitelist = whitelist if whitelist is not None else self.whitelist
        blacklist = blacklist if blacklist is not None else self.blacklist

        def _filter_whitelist(samples):
            return [s for s in samples if s in whitelist]

        def _filter_blacklist(samples):
            return [s for s in samples if s not in blacklist]

        called = _filter_whitelist(called)
        background = _filter_whitelist(background)

        called = _filter_blacklist(called)
        background = _filter_blacklist(background)

        if len(background) >= self.n_background:
            background = np.random.choice(background, self.n_background,
                                          replace=False).tolist()

        return called, background
