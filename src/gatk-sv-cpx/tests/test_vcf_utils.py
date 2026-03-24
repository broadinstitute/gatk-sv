"""
Tests for ``gatk_sv_cpx.vcf_utils``.
"""

import pytest

from gatk_sv_cpx.vcf_utils import (
    get_overlap,
    has_reciprocal_overlap,
    parse_cpx_interval,
    parse_source,
)


class TestParseIntervals:
    def test_parse_cpx_interval_del(self):
        cnv, chrom, start, end = parse_cpx_interval("DEL_chr1:1000-2000")
        assert cnv == "DEL"
        assert chrom == "chr1"
        assert start == 1000
        assert end == 2000

    def test_parse_cpx_interval_dup(self):
        cnv, chrom, start, end = parse_cpx_interval("DUP_chrY:3125606-3125667")
        assert cnv == "DUP"
        assert chrom == "chrY"
        assert start == 3125606
        assert end == 3125667

    def test_parse_source_same_as_interval(self):
        result = parse_source("INV_chr2:5000-6000")
        assert result == ("INV", "chr2", 5000, 6000)


class TestOverlap:
    def test_no_overlap_diff_chrom(self):
        assert get_overlap("chr1", 100, 200, "chr2", 100, 200) == 0

    def test_no_overlap_same_chrom(self):
        assert get_overlap("chr1", 100, 200, "chr1", 300, 400) == 0

    def test_full_overlap(self):
        assert get_overlap("chr1", 100, 200, "chr1", 100, 200) == 100

    def test_partial_overlap(self):
        assert get_overlap("chr1", 100, 200, "chr1", 150, 250) == 50

    def test_reciprocal_overlap_true(self):
        assert has_reciprocal_overlap("chr1", 0, 100, "chr1", 10, 110, 0.5) is True

    def test_reciprocal_overlap_false(self):
        assert has_reciprocal_overlap("chr1", 0, 100, "chr1", 80, 200, 0.5) is False
