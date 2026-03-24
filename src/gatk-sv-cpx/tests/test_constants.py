"""
Tests for ``gatk_sv_cpx.constants``.
"""

from gatk_sv_cpx.constants import CpxType, CPX_TYPES_WITH_CNV, INSERTION_TYPE_CPX


class TestCpxTypes:
    def test_all_cpx_cnv_types_are_valid(self):
        all_values = {e.value for e in CpxType}
        for t in CPX_TYPES_WITH_CNV:
            assert t in all_values, f"{t} not in CpxType enum"

    def test_insertion_types_subset_of_cnv_types(self):
        assert INSERTION_TYPE_CPX <= CPX_TYPES_WITH_CNV
