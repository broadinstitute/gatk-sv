import importlib.util
from pathlib import Path
import sys

import pysam


MODULE_PATH = Path(__file__).with_name("compute_AFs.py")
MODULE_SPEC = importlib.util.spec_from_file_location("compute_AFs", MODULE_PATH)
assert MODULE_SPEC is not None
assert MODULE_SPEC.loader is not None
compute_AFs = importlib.util.module_from_spec(MODULE_SPEC)
MODULE_SPEC.loader.exec_module(compute_AFs)


def _write_sex_chrom_vcf(path: Path, samples, gts):
    header = pysam.VariantHeader()
    header.contigs.add("chrX", length=156040895)
    header.add_line('##INFO=<ID=END,Number=1,Type=Integer,Description="End">')
    header.add_line('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">')
    header.add_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    for sample in samples:
        header.add_sample(sample)

    with pysam.VariantFile(path, "w", header=header) as vcf:
        record = vcf.new_record(
            contig="chrX",
            start=5000000,
            stop=5001000,
            alleles=("N", "<DEL>"),
            id="DEL_chrX_1",
            info={"SVTYPE": "DEL"},
        )
        for sample, gt in zip(samples, gts):
            record.samples[sample]["GT"] = gt
        vcf.write(record)


def _run_compute_afs(monkeypatch, vcf_path: Path, fam_path: Path, out_path: Path, pop_path: Path = None):
    argv = ["compute_AFs.py", "-f", str(fam_path)]
    if pop_path is not None:
        argv.extend(["-p", str(pop_path)])
    argv.extend([str(vcf_path), str(out_path)])
    monkeypatch.setattr(sys, "argv", argv)
    compute_AFs.main()
    with pysam.VariantFile(out_path) as out_vcf:
        return list(out_vcf.fetch())


def test_female_only_cohort_sex_chrom(tmp_path, monkeypatch):
    vcf_path = tmp_path / "in.vcf"
    fam_path = tmp_path / "female.fam"
    pop_path = tmp_path / "pop.tsv"
    out_path = tmp_path / "out.vcf"

    _write_sex_chrom_vcf(vcf_path, ["F1", "F2"], [(0, 1), (0, 0)])
    fam_path.write_text("FAM1\tF1\t0\t0\t2\t0\nFAM2\tF2\t0\t0\t2\t0\n")
    pop_path.write_text("F1\tEUR\nF2\tEUR\n")

    records = _run_compute_afs(monkeypatch, vcf_path, fam_path, out_path, pop_path)
    assert len(records) == 1
    info = records[0].info
    assert info["AN"] == 4
    assert info["AC"] == (1,)
    assert info["AF"] == (0.25,)
    assert info["AN_EUR"] == 4
    assert info["AC_EUR"] == (1,)
    assert info["AF_EUR"] == (0.25,)
    assert "AC_MALE" not in info


def test_male_only_cohort_sex_chrom(tmp_path, monkeypatch):
    vcf_path = tmp_path / "in.vcf"
    fam_path = tmp_path / "male.fam"
    pop_path = tmp_path / "pop.tsv"
    out_path = tmp_path / "out.vcf"

    _write_sex_chrom_vcf(vcf_path, ["M1", "M2"], [(0, 1), (0, 0)])
    fam_path.write_text("FAM1\tM1\t0\t0\t1\t0\nFAM2\tM2\t0\t0\t1\t0\n")
    pop_path.write_text("M1\tEUR\nM2\tEUR\n")

    records = _run_compute_afs(monkeypatch, vcf_path, fam_path, out_path, pop_path)
    assert len(records) == 1
    info = records[0].info
    assert info["AN"] == 2
    assert info["AC"] == (1,)
    assert info["AF"] == (0.5,)
    assert info["AN_EUR"] == 2
    assert info["AC_EUR"] == (1,)
    assert info["AF_EUR"] == (0.5,)
    assert "AC_FEMALE" not in info


def test_mixed_sex_cohort_sex_chrom(tmp_path, monkeypatch):
    vcf_path = tmp_path / "in.vcf"
    fam_path = tmp_path / "mixed.fam"
    out_path = tmp_path / "out.vcf"

    _write_sex_chrom_vcf(vcf_path, ["M1", "F1"], [(0, 1), (0, 1)])
    fam_path.write_text("FAM1\tM1\t0\t0\t1\t0\nFAM2\tF1\t0\t0\t2\t0\n")

    records = _run_compute_afs(monkeypatch, vcf_path, fam_path, out_path)
    assert len(records) == 1
    info = records[0].info
    assert info["AN"] == 3
    assert info["AC"] == (2,)
    assert abs(info["AF"][0] - (2 / 3)) < 1e-6
