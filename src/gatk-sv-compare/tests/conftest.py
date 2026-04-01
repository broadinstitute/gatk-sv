from __future__ import annotations

import sys
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))


@pytest.fixture
def make_vcf(tmp_path):
    def _make_vcf(
        *,
        file_name: str,
        records: list[str],
        final_annotated: bool = False,
        include_members: bool = False,
        include_vargq: bool = False,
    ) -> Path:
        path = tmp_path / file_name
        header_lines = [
            "##fileformat=VCFv4.2",
            "##contig=<ID=chr1,length=1000000>",
            "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"SV type\">",
            "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"SV length\">",
            "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position\">",
            "##INFO=<ID=CHR2,Number=1,Type=String,Description=\"Partner chromosome\">",
            "##INFO=<ID=END2,Number=1,Type=Integer,Description=\"Partner end\">",
            "##INFO=<ID=CPX_TYPE,Number=1,Type=String,Description=\"Complex subtype\">",
            "##INFO=<ID=MULTIALLELIC,Number=0,Type=Flag,Description=\"Multiallelic marker\">",
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
            "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">",
            "##FORMAT=<ID=ECN,Number=1,Type=Integer,Description=\"Expected copy number\">",
        ]
        if include_members:
            header_lines.append("##INFO=<ID=MEMBERS,Number=.,Type=String,Description=\"Members\">")
        if include_vargq:
            header_lines.append("##INFO=<ID=varGQ,Number=1,Type=Integer,Description=\"Variant GQ\">")
        if final_annotated:
            header_lines.extend(
                [
                    "##INFO=<ID=N_BI_GENOS,Number=1,Type=Integer,Description=\"Biallelic genotype count\">",
                    "##INFO=<ID=N_HOMREF,Number=1,Type=Integer,Description=\"Hom ref\">",
                    "##INFO=<ID=N_HET,Number=1,Type=Integer,Description=\"Het\">",
                    "##INFO=<ID=N_HOMALT,Number=1,Type=Integer,Description=\"Hom alt\">",
                    "##INFO=<ID=gnomad_v4.1_sv_AF,Number=1,Type=Float,Description=\"gnomAD AF\">",
                    "##FORMAT=<ID=OGQ,Number=1,Type=Integer,Description=\"Original GQ\">",
                    "##FORMAT=<ID=SL,Number=1,Type=Float,Description=\"Log odds score\">",
                ]
            )
        header = "\n".join(header_lines)
        sample_header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n"
        path.write_text(header + "\n" + sample_header + "\n".join(records) + "\n")
        return path

    return _make_vcf
