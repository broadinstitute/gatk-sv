# -*- coding: utf-8 -*-
#
"""
std_dragen_sv_v3_7_6.py

Standardize a DRAGEN-SV record (DRAGEN-SV v3.7.6 VCF format).
"""

from .standardize import VCFStandardizer
from .std_dragen import DragenStandardizer


@VCFStandardizer.register('dragen_v3.7.6')
class DragenStandardizerV3_7_6(DragenStandardizer):
    def standardize_info(self, std_rec, raw_rec):
        """
        Standardize a DRAGEN-SV v3.7.6 record.

        Identical to DragenStandardizer except for DEL records: END is taken
        from the raw record. Without this, a symbolic <DEL> (single-base REF,
        END in INFO) inherits rlen=0 from the new record under pysam 0.15.x
        and is written with END=POS-1 and SVLEN=-1.
        """
        std_rec = super().standardize_info(std_rec, raw_rec)

        if std_rec.info['SVTYPE'] == 'DEL':
            std_rec.stop = raw_rec.stop
            std_rec.info['SVLEN'] = std_rec.stop - std_rec.pos

        return std_rec
