import importlib.util
from pathlib import Path

import pytest


def _load_pull_snps_module():
    root = Path(__file__).resolve().parents[1]
    module_path = root.joinpath("src", "gatk_sv_ploidy", "pull_snps.py")
    spec = importlib.util.spec_from_file_location("test_pull_snps_module", module_path)
    module = importlib.util.module_from_spec(spec)
    assert spec is not None
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


pull_snps_mod = _load_pull_snps_module()


def test_bgzip_and_index_vcf(tmp_path):
    pysam = pytest.importorskip("pysam")

    input_vcf = tmp_path / "sites.vcf"
    input_vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1,length=248956422>\n"
        "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t101\t.\tA\tG\t.\tPASS\tAF=0.6\n"
    )

    output_vcf = tmp_path / "sites.vcf.gz"
    pull_snps_mod._bgzip_and_index_vcf(str(input_vcf), str(output_vcf))

    assert output_vcf.exists()
    assert Path(f"{output_vcf}.tbi").exists()

    with pysam.BGZFile(str(output_vcf), "rb") as handle:
        contents = handle.read().decode()
    assert "#CHROM\tPOS" in contents

    with pysam.VariantFile(str(output_vcf)) as vcf:
        records = list(vcf.fetch("chr1", 100, 101))
    assert len(records) == 1
    assert records[0].contig == "chr1"
    assert records[0].pos == 101


def test_pull_snps_exports_bgzf_and_indexes(tmp_path, monkeypatch):
    pysam = pytest.importorskip("pysam")

    output_vcf = tmp_path / "gnomad_sites.vcf.gz"

    class DummyTable:
        freq = [type("Freq", (), {"AF": 0.7})()]
        alleles = ["A", "G"]
        filters = []
        locus = type(
            "Locus",
            (),
            {
                "in_autosome": staticmethod(lambda: True),
                "in_x_nonpar": staticmethod(lambda: False),
                "in_y_nonpar": staticmethod(lambda: False),
            },
        )()

        def index_globals(self):
            return type("Globals", (), {"freq_meta": None})()

        def filter(self, expr):
            return self

        def select(self, **kwargs):
            return self

    class DummyHl:
        @staticmethod
        def read_table(path):
            return DummyTable()

        @staticmethod
        def len(value):
            return 2 if isinstance(value, list) else 0

        @staticmethod
        def is_snp(ref, alt):
            return True

        @staticmethod
        def if_else(condition, when_true, when_false):
            return when_true if condition else when_false

        @staticmethod
        def all(fn, values):
            return all(values)

        @staticmethod
        def struct(**kwargs):
            return kwargs

        @staticmethod
        def float64(value):
            return float(value)

        @staticmethod
        def eval(value):
            return [{"group": "adj"}]

        @staticmethod
        def export_vcf(ht, path, metadata=None):
            Path(path).write_text(
                "##fileformat=VCFv4.2\n"
                "##contig=<ID=chr1,length=248956422>\n"
                "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                "chr1\t101\t.\tA\tG\t.\tPASS\tAF=0.7\n"
            )

    monkeypatch.setattr(pull_snps_mod, "_get_adj_freq_index", lambda hl, ht: 0)

    pull_snps_mod.pull_snps(DummyHl(), "dummy.ht", str(output_vcf))

    assert output_vcf.exists()
    assert Path(f"{output_vcf}.tbi").exists()

    with pysam.VariantFile(str(output_vcf)) as vcf:
        records = list(vcf.fetch("chr1", 100, 101))
    assert len(records) == 1
    assert records[0].info["AF"][0] == pytest.approx(0.7)
