from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional

import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from gatk_sv_compare.aggregate import aggregate
from gatk_sv_compare.config import AnalysisConfig


@dataclass
class ModuleTestContext:
    data: object
    config: AnalysisConfig


@pytest.fixture
def make_vcf(tmp_path):
    def _make_vcf(
        *,
        file_name: str,
        records: list[str],
        final_annotated: bool = False,
        include_members: bool = False,
        include_vargq: bool = False,
        extra_header_lines: Optional[List[str]] = None,
        sample_names: Optional[List[str]] = None,
    ) -> Path:
        path = tmp_path / file_name
        resolved_sample_names = sample_names or ["S1", "S2"]
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
        if extra_header_lines:
            header_lines.extend(extra_header_lines)
        header = "\n".join(header_lines)
        sample_header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(resolved_sample_names) + "\n"
        path.write_text(header + "\n" + sample_header + "\n".join(records) + "\n")
        return path

    return _make_vcf


@pytest.fixture
def module_test_context(tmp_path, make_vcf) -> ModuleTestContext:
    extra_headers = [
        "##INFO=<ID=AC,Number=1,Type=Integer,Description=\"Allele count\">",
        "##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele frequency\">",
        "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Allele number\">",
        "##INFO=<ID=ALGORITHMS,Number=.,Type=String,Description=\"Calling algorithms\">",
        "##INFO=<ID=EVIDENCE,Number=.,Type=String,Description=\"Evidence types\">",
        "##INFO=<ID=STATUS,Number=1,Type=String,Description=\"Match status\">",
        "##INFO=<ID=TRUTH_VID,Number=1,Type=String,Description=\"Truth variant id\">",
        "##INFO=<ID=segdup,Number=0,Type=Flag,Description=\"Segmental duplication overlap\">",
        "##INFO=<ID=simple_repeat,Number=0,Type=Flag,Description=\"Simple repeat overlap\">",
        "##INFO=<ID=repeatmasker,Number=0,Type=Flag,Description=\"RepeatMasker overlap\">",
        "##INFO=<ID=VAR_PPV,Number=1,Type=Float,Description=\"Variant PPV\">",
        "##INFO=<ID=CN_NONREF_FREQ,Number=1,Type=Float,Description=\"CNV nonref frequency\">",
        "##INFO=<ID=CN_NUMBER,Number=1,Type=Integer,Description=\"CNV sample count\">",
    ]

    sample_names = ["S1", "S2", "S3"]
    vcf_a = make_vcf(
        file_name="module_a.vcf",
        final_annotated=True,
        sample_names=sample_names,
        extra_header_lines=extra_headers,
        records=[
            "chr1\t100\ta_del\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=150;AC=3;AF=0.5;AN=6;ALGORITHMS=manta,wham;EVIDENCE=RD,PE;N_BI_GENOS=3;N_HOMREF=1;N_HET=1;N_HOMALT=1;STATUS=MATCHED;TRUTH_VID=b_del;segdup;VAR_PPV=0.9;gnomad_v4.1_sv_AF=0.1\tGT:GQ:ECN:OGQ:SL\t0/0:30:2:28:1.0\t0/1:50:2:45:1.0\t1/1:70:2:60:1.0",
            "chr1\t300\ta_ins\tN\t<INS:ME:ALU>\t.\tPASS\tSVTYPE=INS;SVLEN=305;AC=1;AF=0.1667;AN=6;ALGORITHMS=melt;EVIDENCE=SR;N_BI_GENOS=3;N_HOMREF=2;N_HET=1;N_HOMALT=0;STATUS=MATCHED;TRUTH_VID=b_ins;simple_repeat;VAR_PPV=0.95;gnomad_v4.1_sv_AF=0.02\tGT:GQ:ECN:OGQ:SL\t0/0:25:2:20:1.0\t0/1:65:2:55:1.0\t0/0:35:2:30:1.0",
            "chr1\t500\ta_dup\tN\t<DUP>\t.\tUNRESOLVED\tSVTYPE=DUP;SVLEN=1200;AC=2;AF=0.3333;AN=6;ALGORITHMS=depth;EVIDENCE=PE,SR;N_BI_GENOS=3;N_HOMREF=1;N_HET=2;N_HOMALT=0;STATUS=UNMATCHED;repeatmasker;VAR_PPV=0.4;gnomad_v4.1_sv_AF=0.15\tGT:GQ:ECN:OGQ:SL\t0/1:40:2:35:1.0\t0/1:45:2:40:1.0\t0/0:35:2:30:1.0",
            "chr1\t800\ta_cnv\tN\t<CNV>\t.\tMULTIALLELIC\tSVTYPE=CNV;SVLEN=5000;ALGORITHMS=cnmops,depth;EVIDENCE=RD,PE,SR,BAF;CN_NONREF_FREQ=0.3333;CN_NUMBER=3;N_BI_GENOS=3;N_HOMREF=0;N_HET=0;N_HOMALT=0;STATUS=UNMATCHED;gnomad_v4.1_sv_AF=0.0\tGT:GQ:ECN:OGQ:SL\t./.:20:2:18:1.0\t./.:22:2:19:1.0\t./.:24:2:20:1.0",
        ],
    )
    vcf_b = make_vcf(
        file_name="module_b.vcf",
        final_annotated=True,
        sample_names=sample_names,
        extra_header_lines=extra_headers,
        records=[
            "chr1\t120\tb_del\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=160;AC=2;AF=0.3333;AN=6;ALGORITHMS=manta;EVIDENCE=RD,PE;N_BI_GENOS=3;N_HOMREF=1;N_HET=2;N_HOMALT=0;STATUS=MATCHED;TRUTH_VID=a_del;segdup;VAR_PPV=0.88;gnomad_v4.1_sv_AF=0.12\tGT:GQ:ECN:OGQ:SL\t0/1:55:2:45:1.0\t0/1:52:2:44:1.0\t0/0:33:2:28:1.0",
            "chr1\t320\tb_ins\tN\t<INS:ME:ALU>\t.\tPASS\tSVTYPE=INS;SVLEN=295;AC=1;AF=0.1667;AN=6;ALGORITHMS=melt,wham;EVIDENCE=SR;N_BI_GENOS=3;N_HOMREF=2;N_HET=1;N_HOMALT=0;STATUS=MATCHED;TRUTH_VID=a_ins;simple_repeat;VAR_PPV=0.98;gnomad_v4.1_sv_AF=0.01\tGT:GQ:ECN:OGQ:SL\t0/0:20:2:18:1.0\t0/1:62:2:55:1.0\t0/0:40:2:32:1.0",
            "chr1\t650\tb_inv\tN\t<INV>\t.\tPASS\tSVTYPE=INV;SVLEN=2000;AC=1;AF=0.1667;AN=6;ALGORITHMS=depth;EVIDENCE=PE,SR;N_BI_GENOS=3;N_HOMREF=2;N_HET=1;N_HOMALT=0;STATUS=UNMATCHED;repeatmasker;VAR_PPV=0.7;gnomad_v4.1_sv_AF=0.05\tGT:GQ:ECN:OGQ:SL\t0/0:42:2:38:1.0\t0/0:39:2:32:1.0\t0/1:60:2:50:1.0",
        ],
    )

    config = AnalysisConfig(
        vcf_a_path=vcf_a,
        vcf_b_path=vcf_b,
        vcf_a_label="CallsetA",
        vcf_b_label="CallsetB",
        output_dir=tmp_path / "outputs",
        contigs=["chr1"],
        contig_lengths={"chr1": 1_000_000},
        n_workers=1,
    )
    data = aggregate(config)
    return ModuleTestContext(data=data, config=config)

