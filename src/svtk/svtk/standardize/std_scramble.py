# -*- coding: utf-8 -*-
#
"""
std_scramble.py

Standardize a Scramble record.
"""

from .standardize import VCFStandardizer


@VCFStandardizer.register('scramble')
class ScrambleStandardizer(VCFStandardizer):
    def standardize_info(self, std_rec, raw_rec):
        """
        Standardize Scramble record -- nothing to do.
        """
        return std_rec
