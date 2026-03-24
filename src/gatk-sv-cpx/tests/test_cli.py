"""
Tests for ``gatk_sv_cpx.cli`` argument parsing.
"""

import pytest

from gatk_sv_cpx.cli import main


class TestCliParsing:
    """Verify that the CLI parses subcommands and arguments correctly."""

    def test_resolve_requires_input_and_output(self):
        with pytest.raises(SystemExit):
            main(["resolve"])

    def test_evaluate_evidence_requires_input_and_output(self):
        with pytest.raises(SystemExit):
            main(["evaluate-evidence"])

    def test_genotype_and_refine_requires_evidence(self):
        with pytest.raises(SystemExit):
            main(["genotype-and-refine", "-i", "in.vcf", "-o", "out.vcf"])

    def test_no_subcommand_exits(self):
        with pytest.raises(SystemExit):
            main([])
