# -*- coding: utf-8 -*-
#
"""
std_melt.py

Standardize a MELT record.
"""

from .standardize import VCFStandardizer


@VCFStandardizer.register('melt')
class MeltStandardizer(VCFStandardizer):
    def standardize_info(self, std_rec, raw_rec):
        """
        Standardize MELT record.

        1) Rename all SVTYPE fields to "INS" (keep subclass in ALT)
        2) Add END (POS + 1)
        3) Add CHR2 (CHROM)
        4) Add STRANDS (+-; treat as deletion since no start/end bkpt info)
        5) Add ALGORITHMS.
        """

        # Rename SVTYPE subclasses to INS
        std_rec.info['SVTYPE'] = 'INS'

        # Add END
        std_rec.stop = raw_rec.pos + 1

        # Add STRANDS
        std_rec.info['STRANDS'] = '+-'

        # Add CHR2
        std_rec.info['CHR2'] = raw_rec.chrom

        # Add SVLEN
        std_rec.info['SVLEN'] = raw_rec.info['SVLEN']

        # Add ALGORITHMS
        std_rec.info['ALGORITHMS'] = ['melt']

        return std_rec
