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
    # Currently identical to DragenStandardizer; override methods here if
    # DRAGEN-SV v3.7.6 ever needs version-specific standardization behavior.
    pass
