import importlib.util
from pathlib import Path

import pysam
import pytest

MODULE_PATH = Path(__file__).with_name("annotate_moi.py")
MODULE_SPEC = importlib.util.spec_from_file_location("annotate_moi", MODULE_PATH)
assert MODULE_SPEC is not None
assert MODULE_SPEC.loader is not None
annotate_moi = importlib.util.module_from_spec(MODULE_SPEC)
MODULE_SPEC.loader.exec_module(annotate_moi)

SAMPLES = ["case", "mother", "father"]


def _write_vcf(path: Path, genotypes):
    """genotypes: list of (name, {sample: gt}) for the samples in SAMPLES order."""
    header = pysam.VariantHeader()
    header.contigs.add("chr1", length=1000)
    header.add_line('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Variant type">')
    header.add_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    for sample in SAMPLES:
        header.add_sample(sample)

    with pysam.VariantFile(str(path), "w", header=header) as vcf:
        for i, (name, gts) in enumerate(genotypes):
            record = vcf.new_record(contig="chr1", start=i * 100 + 1, alleles=("A", "C"))
            for sample in SAMPLES:
                record.samples[sample]["GT"] = gts[sample]
            record.info["SVTYPE"] = "DEL"
            vcf.write(record)


def _run_main(monkeypatch, tmp_path, vcf_path, case, mother=None, father=None):
    argv = [
        "annotate_moi.py",
        str(vcf_path),
        str(tmp_path / "out"),
        "--case", case,
    ]
    if mother is not None:
        argv += ["--mother", mother]
    if father is not None:
        argv += ["--father", father]
    monkeypatch.setattr("sys.argv", argv)
    annotate_moi.main()
    return tmp_path / "out.vcf.gz"


def _moi_map(out_path):
    with pysam.VariantFile(str(out_path)) as vcf:
        return {
            record.contig + ":" + str(record.start): (
                record.info.get("MOI"), record.info.get("MOI_CONFIDENCE")
            )
            for record in vcf
        }


def test_full_trio_moi_values(monkeypatch, tmp_path):
    genotypes = [
        ("de_novo", {"case": (0, 1), "mother": (0, 0), "father": (0, 0)}),
        ("mat", {"case": (0, 1), "mother": (0, 1), "father": (0, 0)}),
        ("pat", {"case": (0, 1), "mother": (0, 0), "father": (1, 1)}),
        ("both", {"case": (1, 1), "mother": (0, 1), "father": (0, 1)}),
        ("parent_only", {"case": (0, 0), "mother": (0, 0), "father": (1, 1)}),
        ("case_missing", {"case": (None, None), "mother": (0, 1), "father": (0, 0)}),
        ("all_ref", {"case": (0, 0), "mother": (0, 0), "father": (0, 0)}),
    ]
    vcf_path = tmp_path / "in.vcf.gz"
    _write_vcf(vcf_path, genotypes)
    out = _run_main(monkeypatch, tmp_path, vcf_path, "case", "mother", "father")

    moi = _moi_map(out)
    assert len(moi) == 7
    assert moi["chr1:1"] == ("DE_NOVO", "CONFIRMED")
    assert moi["chr1:101"] == ("INHERITED_FROM_MOTHER", "CONFIRMED")
    assert moi["chr1:201"] == ("INHERITED_FROM_FATHER", "CONFIRMED")
    assert moi["chr1:301"] == ("INHERITED_FROM_BOTH", "CONFIRMED")
    assert moi["chr1:401"] == ("PARENT_ONLY", "CONFIRMED")
    assert moi["chr1:501"] == ("UNASSESSABLE", "CONFIRMED")
    assert moi["chr1:601"] == ("UNASSESSABLE", "CONFIRMED")

    summary = (tmp_path / "out.moi_summary.tsv").read_text().splitlines()
    assert summary[0] == "MOI\tCOUNT"
    assert "DE_NOVO\t1" in summary
    assert "INHERITED_FROM_MOTHER\t1" in summary
    assert "UNASSESSABLE\t2" in summary
    assert (tmp_path / "out.vcf.gz.tbi").exists()


def test_half_trio_mother_only(monkeypatch, tmp_path):
    genotypes = [
        ("de_novo_unconfirmed", {"case": (0, 1), "mother": (0, 0), "father": (0, 1)}),
        ("mat_inherited", {"case": (0, 1), "mother": (0, 1), "father": (0, 0)}),
        ("parent_only", {"case": (0, 0), "mother": (0, 1), "father": (0, 0)}),
        ("parent_gt_missing", {"case": (0, 1), "mother": (None, None), "father": (0, 0)}),
    ]
    vcf_path = tmp_path / "in.vcf.gz"
    _write_vcf(vcf_path, genotypes)
    out = _run_main(monkeypatch, tmp_path, vcf_path, "case", "mother")

    moi = _moi_map(out)
    # father not assayed: de novo cannot be confirmed
    assert moi["chr1:1"] == ("DE_NOVO", "UNCONFIRMED")
    assert moi["chr1:101"] == ("INHERITED_FROM_MOTHER", "UNCONFIRMED")
    assert moi["chr1:201"] == ("PARENT_ONLY", "UNCONFIRMED")
    assert moi["chr1:301"] == ("DE_NOVO", "UNCONFIRMED")


def test_case_only_callset(monkeypatch, tmp_path):
    genotypes = [
        ("called", {"case": (0, 1), "mother": (0, 0), "father": (0, 0)}),
    ]
    vcf_path = tmp_path / "in.vcf.gz"
    _write_vcf(vcf_path, genotypes)
    out = _run_main(monkeypatch, tmp_path, vcf_path, "case")

    moi = _moi_map(out)
    assert moi["chr1:1"] == ("DE_NOVO", "UNCONFIRMED")


def test_missing_case_sample_errors(monkeypatch, tmp_path):
    genotypes = [("x", {"case": (0, 1), "mother": (0, 0), "father": (0, 0)})]
    vcf_path = tmp_path / "in.vcf.gz"
    _write_vcf(vcf_path, genotypes)
    monkeypatch.setattr("sys.argv", ["annotate_moi.py", str(vcf_path), str(tmp_path / "out"), "--case", "nobody"])
    with pytest.raises(SystemExit):
        annotate_moi.main()
