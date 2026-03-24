"""
Tests for ``gatk_sv_cpx.evaluate_evidence``.
"""

from gatk_sv_cpx.evaluate_evidence import classify_pe_support, classify_depth_support


class TestPeClassification:
    def test_no_pe(self):
        assert classify_pe_support([[], []], min_pe=3) == "no_PE"

    def test_partial_pe_fwd_only(self):
        assert classify_pe_support([[2, "+", "100"], []], min_pe=3) == "partial_PE"

    def test_low_pe(self):
        assert classify_pe_support([[1, "+", "100"], [1, "-", "200"]], min_pe=3) == "low_PE"

    def test_high_pe(self):
        assert classify_pe_support([[5, "+", "100"], [5, "-", "200"]], min_pe=3) == "high_PE"


class TestDepthClassification:
    def test_na_empty(self):
        assert classify_depth_support([]) == "NA"

    def test_depth_above_threshold(self):
        assert classify_depth_support([0.8]) == "depth"

    def test_lack_depth(self):
        assert classify_depth_support([0.3]) == "lack_depth"

    def test_mixed(self):
        assert classify_depth_support([0.8, 0.2]) == "depth,lack_depth"
