#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""
Wrapper for breakpoints stored in VCF records
"""

import numpy as np
import svtk.utils as svu


class Breakpoint:
    def __init__(self, chrA, posA, chrB, posB, name, samples, strands):
        """
        Attributes
        ----------
        chrA : str
        posA : int
        chrB : str
        posB : str
        name : str
        samples : list of str
        strands : str
            +-,-+,++,--
        background : list of str
        split_counts : pd.DataFrame
        """
        self.chrA = chrA
        self.posA = posA
        self.chrB = chrB
        self.posB = posB
        self.name = name
        self.samples = samples
        self.strands = strands

        # Constructed later
        self.background = []
        self.split_counts = None

    @classmethod
    def from_vcf(cls, record, whitelist=None):
        """
        Parameters
        ----------
        record : pysam.VariantRecord
        """

        chrA = record.chrom
        posA = record.pos
        chrB = record.info['CHR2']
        posB = record.info['END2'] if record.info['SVTYPE'] == 'BND' else record.stop

        name = record.id
        strands = record.info['STRANDS']

        samples = svu.get_called_samples(record)
        if whitelist is not None:
            samples = [s for s in samples if s in whitelist]

        return cls(chrA, posA, chrB, posB, name, samples, strands)

    @classmethod
    def from_bed(cls, chrom, start, end, name, samples, svtype):
        chrA = chrom
        posA = int(start)
        chrB = chrom
        posB = int(end)

        name = name
        samples = samples.split(',')

        if svtype.upper() not in 'DEL DUP'.split():
            msg = 'Invalid SV type: {0} [expected del,dup,DEL,DUP]'
            msg = msg.format(svtype)
            raise Exception(msg)

        if svtype.upper() == 'DEL':
            strands = '+-'
        else:
            strands = '-+'

        return cls(chrA, posA, chrB, posB, name, samples, strands)

    def choose_background(self, samples, n_background):
        """
        Choose background samples

        Parameters
        ----------
        samples : list of str
            All available samples
        n_background : int
            Max number of background samples
        """

        self.background = [s for s in samples if s not in self.samples]
        if len(self.background) >= n_background:
            self.background = np.random.choice(self.background, n_background,
                                               replace=False).tolist()
