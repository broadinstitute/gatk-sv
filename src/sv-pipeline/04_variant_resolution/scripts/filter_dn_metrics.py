#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""
Filter putative de novo variants.

1) Evaluates the results of a single-sample test against metric cutoffs derived
from a random-forest classifier, removing false positives in children and
adding false negatives in parents.
2) For PE/SR BCA and small CNV, tests the residual between a child and each
parent to add additional false negatives in parents.
"""

import argparse
from collections import namedtuple
import itertools
import numpy as np
import scipy.stats as ss
import pandas as pd
from svtk.famfile import parse_famfile


class DenovoFilter:
    def __init__(self, child, father, mother, cutoffs):
        """
        Clean putative de novo variants with single sample and residual tests.

        Attributes
        ----------
        child : Metrics
        father : Metrics
        mother : Metrics
        cutoffs : pd.DataFrame
        samples : list of Metrics
        called_samples : list of Metrics
        cutoff_class : (str, str)
        """

        self.child = child
        self.father = father
        self.mother = mother

        # Choose cutoffs based on variant type and convert to {metric: cutoff}
        self.cutoff_class = get_cutoff_class(child)
        self.cutoffs = cutoffs

        self.samples = [child, father, mother]
        self.called_samples = []
        self.supports = []

    def filter(self):
        # Test child individually against cutoffs
        if self.test_child():
            self.called_samples.append(self.child)
            if hasattr(self, 'SR_support') and self.SR_support:
                self.supports.append('SR')
            if hasattr(self, 'PE_support') and self.PE_support:
                self.supports.append('PE')
            if hasattr(self, 'RD_support') and self.RD_support:
                self.supports.append('RD')
        # If child failed single-sample test, don't add parents
        else:
            return

        # Test each parent individually against cutoffs
        for parent in [self.father, self.mother]:
            if self.test_parent(parent):
                self.called_samples.append(parent)

        # Test if child stats remain significant after taking residual
        if self.cutoff_class == ('pesr', 'SV'):
            # Run residual test on both parents
            for parent in [self.father, self.mother]:
                if self.test_residual(parent):
                    self.called_samples.append(parent)

    def test_child(self):
        cutoffs = self.get_cutoffs(self.cutoff_class)

        PESR_cutoffs = self.get_cutoffs(('pesr', 'SV'))
        RD_cutoffs = self.get_cutoffs(('pesr', 'CNV'))

        # Check if variant is individually supported by PE, SR, or RD
        self.PE_support = self.child.PE_log_pval >= PESR_cutoffs['PE_log_pval']
        self.SR_support = (self.child.SR_sum_log_pval >=
                           PESR_cutoffs['SR_sum_log_pval'])
        self.RD_support = self.test_sample(self.child, RD_cutoffs, np.all)

        # BCA and sub 1kb CNV get the respective cutoffs,
        # plus additional filters
        if self.cutoff_class == ('pesr', 'SV'):
            # get default filter status
            filter_pass = self.test_sample(self.child, cutoffs, np.any)

            # Apply more stringent de novo filters for SR-supported variants
            SR_pass = (self.child.SR_posA_called_median >= 2 and
                       self.child.SR_posB_called_median >= 2)

            if self.child.svsize > 0:
                dist_pass = (self.child.SR_dist / self.child.svsize) >= 0.8
            else:
                dist_pass = False
            SR_pass = SR_pass and dist_pass

            # Apply more stringent de novo filters for PE-supported CNV
            if self.child.svtype in 'DEL DUP'.split():
                PE_pass = self.RD_support
            else:
                PE_pass = True

            # Add appropriate filters
            if self.PE_support and not self.SR_support:
                filter_pass = filter_pass and PE_pass
            elif not self.PE_support and self.SR_support:
                filter_pass = filter_pass and SR_pass
            elif self.PE_support and self.SR_support:
                filter_pass = filter_pass and (SR_pass or PE_pass)

            # Add more stringent background test
            BG_pass = (self.child.PE_bg_median < 2 and
                       self.child.SR_sum_bg_median < 2)
            filter_pass = filter_pass and BG_pass

        # >1 kb PE/SR CNV require RD and PE/SR cutoffs
        elif self.cutoff_class == ('pesr', 'CNV'):
            filter_pass = self.test_sample(self.child, cutoffs, np.all)

            # In addition to default depth, require pe/sr support
            PESR_cutoffs = self.get_cutoffs(('pesr', 'SV'))
            self.PESR_pass = self.test_sample(self.child, PESR_cutoffs, np.any)
            filter_pass = filter_pass and self.PESR_pass

            # If >5 kb, also permit more stringent depth-only RD support
            if self.child.svsize >= 5000:
                RD_cutoffs = self.get_cutoffs(('depth', self.child.svtype))
                RD_pass = self.test_sample(self.child, RD_cutoffs, np.all)
                RD_pass = RD_pass and (self.child.svsize >= 5000)
                filter_pass = filter_pass or RD_pass

        # everything else gets default cutoffs
        else:
            filter_pass = self.test_sample(self.child, cutoffs, np.all)

        return filter_pass

    def test_parent(self, parent):
        # BCA and sub 1kb CNV get the respective cutoffs,
        # plus additional filters
        if self.cutoff_class == ('pesr', 'SV'):
            cutoffs = self.get_cutoffs(('pesr', 'SV'))
            filter_pass = self.test_sample(parent, cutoffs, np.any)

            if self.SR_support:
                filter_pass = filter_pass or (parent.SR_sum_called_median >= 2)
            if self.PE_support:
                filter_pass = filter_pass or (parent.PE_called_median >= 2)

            cutoffs = self.get_cutoffs(('pesr', 'CNV'))
            filter_pass = (filter_pass or
                           parent.RD_log_pval >= cutoffs['RD_log_pval'])

        elif self.cutoff_class == ('pesr', 'CNV'):
            # First try checking against RD cutoffs
            cutoffs = self.get_cutoffs(('pesr', 'CNV'))
            filter_pass = parent.RD_log_pval >= cutoffs['RD_log_pval']

            # If >5 kb and child passed pe/sr, also check PE/SR pvals
            if parent.svsize >= 5000 and self.PESR_pass:
                PESR_cutoffs = self.get_cutoffs(('pesr', 'SV'))
                PESR_pass = self.test_sample(parent, PESR_cutoffs, np.any)
                filter_pass = filter_pass or PESR_pass
            # If <5 kb, check for any PE/SR evidence
            elif parent.svsize < 5000:
                if self.SR_support:
                    filter_pass = filter_pass or (
                        parent.SR_sum_called_median >= 2)
                if self.PE_support:
                    filter_pass = filter_pass or (parent.PE_called_median >= 2)

        # depth CNV only require pe/sr CNV >1kb cutoffs
        else:
            cutoffs = self.get_cutoffs(('pesr', 'CNV'))
            filter_pass = parent.RD_log_pval >= cutoffs['RD_log_pval']

        return filter_pass

    def get_cutoffs(self, cutoff_class):
        cutoffs = self.cutoffs.loc[cutoff_class]
        return cutoffs.squeeze().dropna().to_dict()

    def test_sample(self, sample, cutoffs, agg_func):
        """
        Test single sample's evidence against RF cutoffs

        Arguments
        ---------
        sample : namedtuple('Metrics')
        cutoffs : dict of {str: float}
        agg_func : function
            any or all

        Returns
        -------
        filter_pass : bool
            True if sample stats pass RF cutoffs
        """

        filters = []

        for metric, cutoff in cutoffs.items():
            if metric == 'svsize':
                continue
            filters.append(getattr(sample, metric) >= cutoff)

        return agg_func(filters)

    def test_residual(self, parent):
        """
        Test if residual of parent and child stats remains significant.

        If not, consider the variant inherited from the respective parent.

        Arguments
        ---------
        parent : namedtuple('Metrics')

        Returns
        -------
        resid_filter : bool
            True if variant appears inherited from given parent
        """

        # get names of variables to test
        signals = 'PE SR_sum PESR'.split()
        called_metrics = [signal + '_called_median' for signal in signals]
        bg_metrics = [signal + '_bg_median' for signal in signals]
        cutoff_cols = [signal + '_log_pval' for signal in signals]

        # extract counts to test
        parent = np.array([getattr(parent, m) for m in called_metrics])
        child = np.array([getattr(self.child, m) for m in called_metrics])
        bg = np.array([getattr(self.child, m) for m in bg_metrics])

        # get cutoffs
        cutoffs = self.get_cutoffs(self.cutoff_class)
        cutoffs = np.array([cutoffs.get(m) for m in cutoff_cols])

        # test if residual remains significant
        resid = np.abs((parent - child))

        resid_test = -np.log10(ss.poisson.cdf(bg, resid)) >= cutoffs

        # only test stats that were originally significant in child
        raw_test = -np.log10(ss.poisson.cdf(bg, child)) >= cutoffs
        resid_test = resid_test[np.where(raw_test)]

        # if any stat no longer appears significant, call it inherited
        return np.any(~resid_test)


def get_cutoff_class(sample):
    """
    Identify which set of cutoffs to use when classifying a given variant call

    Arguments
    ---------
    sample : namedtuple('Metrics')
        Row from variant metrics table

    Returns
    -------
    (source, svtype) : (str, str)
        Key to cutoffs table
    """

    if sample.sources == 'depth' and sample.svtype == 'DEL':
        return ('depth', 'DEL')
    elif sample.sources == 'depth' and sample.svtype == 'DUP':
        return ('depth', 'DUP')
    elif (sample.svtype in 'DEL DUP'.split() and sample.svsize >= 1000):
        return ('pesr', 'CNV')
    else:
        return ('pesr', 'SV')


def filter_denovo(metrics, fam, cutoffs):
    """

    Arguments
    ---------
    metrics : iterable of Metrics
    cutoffs : pd.DataFrame
        (source, svtype) as index

    Yields
    ------
    name, family, samples : str, str, str
    """

    # Get metrics of each family called in each variant
    grouped = itertools.groupby(metrics, key=lambda m: (m.name, m.family))

    for (name, family), samples in grouped:
        children = []

        for sample in samples:
            if fam.samples[sample.sample].has_parents:
                children.append(sample)
            else:
                child = fam.samples[sample.sample].children[0]
                if sample.sample == fam.samples[child].mother:
                    mother = sample
                elif sample.sample == fam.samples[child].father:
                    father = sample

        # test each child separately
        tests = [DenovoFilter(child, father, mother, cutoffs)
                 for child in children]

        # run de novo filtering to identify final set of called samples
        supports = []
        for test in tests:
            test.filter()
            supports += test.supports

        if len(supports) == 0:
            support = 'NA'
        else:
            support = ','.join(sorted(set(supports)))

        called_samples = [s for test in tests for s in test.called_samples]
        called_samples = sorted(set([s.sample for s in called_samples]))
        if len(called_samples) == 0:
            called_samples = ['NA']

        yield (name, family, ','.join(called_samples), support)


def metric_parser(metricfile, fam):
    """Convert table of variant metrics to iterator over Metrics tuples"""
    header = next(metricfile).strip().split()
    Metrics = namedtuple('Metrics', ['family'] + header)

    def _cast(data):
        for i, var in enumerate(header):
            if var == 'svsize':
                data[i] = np.int(data[i])
            elif var.startswith('PE') or var.startswith('SR'):
                if data[i] == 'NA':
                    data[i] = np.nan
                else:
                    data[i] = np.float(data[i])
            elif var.startswith('RD') and var != 'RD_Model':
                if data[i] == 'NA':
                    data[i] = np.nan
                else:
                    data[i] = np.float(data[i])
        return data

    for line in metricfile:
        data = line.strip().split()
        family = fam.samples[data[1]].family
        yield Metrics(family, *_cast(data))


def format_cutoffs(cutoffs):
    metrics = 'PESR_log_pval PE_log_pval SR_sum_log_pval RD_Median_Separation RD_log_2ndMaxP RD_log_pval'.split()
    cutoffs = cutoffs.loc[cutoffs.metric.isin(metrics)].copy()

    tests = 'PESR SR2 PE RD'.split()
    cutoffs = cutoffs.loc[cutoffs.test.isin(tests)].copy()

    cutoffs = cutoffs.loc[cutoffs.max_svsize.isnull()].copy()

    cutoffs.loc[cutoffs.test != 'RD', 'svtype'] = 'SV'
    cutoffs.algtype = cutoffs.algtype.str.lower()

    cutoffs = cutoffs.pivot_table(index='algtype svtype'.split(),
                                  values='cutoff', columns='metric')

    return cutoffs


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('metrics', type=argparse.FileType('r'))
    parser.add_argument('cutoffs')
    parser.add_argument('famfile', type=argparse.FileType('r'))
    parser.add_argument('fout', type=argparse.FileType('w'))
    args = parser.parse_args()

    fam = parse_famfile(args.famfile)
    metrics = metric_parser(args.metrics, fam)

    cutoffs = pd.read_table(args.cutoffs)
    cutoffs = format_cutoffs(cutoffs)

    fmt = '{0}\t{1}\t{2}\t{3}\n'
    for name, family, samples, support in filter_denovo(metrics, fam, cutoffs):
        args.fout.write(fmt.format(name, family, samples, support))


if __name__ == '__main__':
    main()
