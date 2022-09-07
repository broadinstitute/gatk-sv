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
        Standardize Scramble record
        """
        std_rec.stop = raw_rec.pos + 1
        std_rec.info["SVTYPE"] = raw_rec.info["SVTYPE"]
        std_rec.info["SVLEN"] = raw_rec.info["SVLEN"]
        std_rec.info["STRANDS"] = "+-"
        std_rec.info["ALGORITHMS"] = ["scramble"]

        return std_rec
