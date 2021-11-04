#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""

"""

import os
import tempfile
import pkg_resources
import subprocess as sp
from collections import namedtuple
import numpy as np
import pandas as pd
from .utils import get_called_samples


def _record_to_bed(record):
    entry = '{chrom}\t{start}\t{end}\t{name}\t{samples}\t{svtype}\n'

    if record.info['SVTYPE'] in 'DEL DUP'.split():
        return entry.format(chrom=record.chrom,
                            start=record.pos,
                            end=record.stop,
                            name=record.id,
                            samples=','.join(get_called_samples(record)),
                            svtype=record.info['SVTYPE'])
    else:
        return ''


def _make_rdtest_bed(variants):
    """
    Make temporary bed file for RdTest

    Parameters
    ----------
    variants : list of pysam.VariantRecord

    Returns
    -------
    bed : tempfile.NamedTemporaryFile
    """

    bed = tempfile.NamedTemporaryFile(dir=os.getcwd())

    for variant in variants:
        bed_entry = _record_to_bed(variant)
        bed.write(bed_entry.encode('utf-8'))
        bed.flush()

    return bed


class RdTest:
    def __init__(self, bincov_file, medianfile, famfile, whitelist,
                 cutoffs=None):
        self.bincov_file = bincov_file
        self.medianfile = medianfile
        self.famfile = famfile
        self.whitelist = whitelist
        self.cutoffs = cutoffs

    def get_cutoffs(self, cutoff_type):
        if cutoff_type == 'pesr_lt1kb':
            cutoffs = self.cutoffs.loc[self.cutoffs.algtype == 'PESR']
            cutoffs = self.cutoffs.loc[self.cutoffs.max_svsize == 1000]
        elif cutoff_type == 'pesr_gt1kb':
            cutoffs = self.cutoffs.loc[self.cutoffs.algtype == 'PESR']
            cutoffs = self.cutoffs.loc[self.cutoffs.min_svsize == 1000]
        elif cutoff_type == 'depth':
            cutoffs = self.cutoffs.loc[self.cutoffs.algtype == 'Depth']
        else:
            msg = ('cutoff_type must be pesr_lt1kb, pesr_gt1kb, or depth, '
                   'not {0}')
            msg = msg.format(cutoff_type)
            raise Exception(msg)

        min_Median_Separation = cutoffs.loc[cutoffs.metric ==
                                            'RD_Median_Separation', 'cutoff'].iloc[0]
        min_log_pval = cutoffs.loc[cutoffs.metric ==
                                   'RD_log_pval', 'cutoff'].iloc[0]
        min_log_2ndMaxP = cutoffs.loc[cutoffs.metric ==
                                      'RD_log_2ndMaxP', 'cutoff'].iloc[0]

        Cutoffs = namedtuple('Cutoffs', ['min_Median_Separation',
                                         'min_log_pval', 'min_log_2ndMaxP'])
        return Cutoffs(min_Median_Separation, min_log_pval, min_log_2ndMaxP)

    def test_record(self, record, cutoff_type='pesr_gt1kb'):
        if self.cutoffs is None:
            raise Exception('Record testing not available without cutoffs')
        metrics = call_rdtest([record], self.bincov_file, self.medianfile,
                              self.famfile, self.whitelist, quiet=True)
        metrics = metrics.iloc[0]

        cutoffs = self.get_cutoffs(cutoff_type)

        # if separation is a string, coverage failure
        if isinstance(metrics.Median_Separation, str):
            return False
        else:
            return (metrics.Median_Separation >= cutoffs.min_Median_Separation and
                    -np.log10(metrics.P) >= cutoffs.min_log_pval and
                    -np.log10(metrics['2ndMaxP']) >= cutoffs.min_log_2ndMaxP)

    def test(self, records, quiet=True):
        metrics = call_rdtest(records, self.bincov_file, self.medianfile,
                              self.famfile, self.whitelist, quiet=True)

        return metrics


def call_rdtest(variants, bincov_file, medianfile, famfile, whitelist,
                quiet=False):
    """
    Utility wrapper around a basic RdTest call

    Parameters
    ----------
    variants : list of pysam.VariantRecord
    bincov_file : str
    medianfile : str
    famfile : str
    whitelist : str or list of str
        Filepath to sample whitelist or list of whitelisted sample IDs
    quiet : bool
        Suppress RdTest stdout/stderr

    Returns
    -------
    metrics : pd.DataFrame
        RdTest metrics for the provided variants
    """

    if not os.path.exists(bincov_file):
        raise Exception('Bincov file does not exist: {0}'.format(bincov_file))

    if not os.path.exists(medianfile):
        raise Exception('Medianfile does not exist: {0}'.format(medianfile))

    if not os.path.exists(famfile):
        raise Exception('Famfile does not exist: {0}'.format(famfile))

    if isinstance(whitelist, str):
        whitelist_filename = whitelist
    elif isinstance(whitelist, list):
        whitelist_file = tempfile.NamedTemporaryFile(dir=os.getcwd())
        for sample in whitelist:
            entry = sample + '\n'
            whitelist_file.write(entry.encode('utf-8'))
        whitelist_file.flush()
        whitelist_filename = whitelist_file.name
    else:
        msg = 'Invalid type for whitelist: {0}\n'.format(str(type(whitelist)))
        msg += 'Must be str or list of str'
        raise Exception(msg)

    for variant in variants:
        if variant.info['SVTYPE'] not in 'DEL DUP'.split():
            msg = 'Invalid svtype {0} for RdTest in record {1}'
            msg = msg.format(variant.info['SVTYPE'], variant.id)
            raise Exception(msg)

    output_dir = tempfile.TemporaryDirectory(dir=os.getcwd())

    bed = _make_rdtest_bed(variants)

    RdTest = pkg_resources.resource_filename('svtk', 'RdTest/RdTest.R')

    if quiet:
        FNULL = open(os.devnull, 'w')
        stdout = FNULL
        stderr = FNULL
    else:
        stdout = None
        stderr = None

    sp.run(['Rscript', RdTest,
            '-b', bed.name,
            '-o', output_dir.name + '/',
            '-n', 'tmp',
            '-c', bincov_file,
            '-m', medianfile,
            '-f', famfile,
            '-w', whitelist_filename],
           stdout=stdout,
           stderr=stderr)

    metrics = pd.read_table(os.path.join(output_dir.name, 'tmp.metrics'))

    return metrics


def filter_rdtest(variants, cutoffs):
    """
    Run RdTest on variants, return only those which met RF cutoffs
    """
