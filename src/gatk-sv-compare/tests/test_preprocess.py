from __future__ import annotations

from pathlib import Path

import pysam

from gatk_sv_compare.config import AnalysisConfig
from gatk_sv_compare.preprocess import (
    build_svconcordance_command,
    build_svregionoverlap_command,
    parse_reference_dict,
    read_contig_list,
    run_preprocess,
)


def _write_bgzipped_vcf(path: Path, contig: str, variant_id: str) -> None:
    header = pysam.VariantHeader()
    header.add_meta("fileformat", "VCFv4.2")
    header.contigs.add(contig, length=1000)
    header.add_meta("INFO", items=[("ID", "SVTYPE"), ("Number", 1), ("Type", "String"), ("Description", "SV type")])
    header.add_meta("INFO", items=[("ID", "SVLEN"), ("Number", 1), ("Type", "Integer"), ("Description", "SV length")])
    header.add_meta("FORMAT", items=[("ID", "GT"), ("Number", 1), ("Type", "String"), ("Description", "Genotype")])
    with pysam.VariantFile(str(path), "wz", header=header) as out_vcf:
        record = out_vcf.new_record(contig=contig, start=99, stop=199, alleles=("N", "<DEL>"), id=variant_id, filter="PASS")
        record.info["SVTYPE"] = "DEL"
        record.info["SVLEN"] = 100
        out_vcf.write(record)
    pysam.tabix_index(str(path), preset="vcf", force=True)


def test_read_contig_list_and_parse_reference_dict(tmp_path) -> None:
    contig_list = tmp_path / "contigs.list"
    contig_list.write_text("# comment\nchr1\n\nchr2\n")
    reference_dict = tmp_path / "ref.dict"
    reference_dict.write_text("@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:2000\n")

    assert read_contig_list(contig_list) == ["chr1", "chr2"]
    assert parse_reference_dict(reference_dict) == {"chr1": 1000, "chr2": 2000}


def test_build_svconcordance_command_includes_tracks() -> None:
    command = build_svconcordance_command(
        eval_vcf=Path("eval.vcf.gz"),
        truth_vcf=Path("truth.vcf.gz"),
        contig="chr1",
        reference_dict=Path("ref.dict"),
        output_path=Path("out.vcf.gz"),
        gatk_path="gatk",
        java_options="-Xmx4g",
        track_names=["segdup"],
        track_intervals=[Path("segdup.bed")],
    )
    assert command[:4] == ["gatk", "--java-options", "-Xmx4g", "SVConcordance"]
    assert ["-L", "chr1"] == command[4:6]
    assert "--track-name" in command
    assert "segdup" in command
    assert "--track-intervals" in command


def test_build_svregionoverlap_command_includes_tracks() -> None:
    command = build_svregionoverlap_command(
        vcf=Path("in.vcf.gz"),
        reference_dict=Path("ref.dict"),
        output_path=Path("out.vcf.gz"),
        gatk_path="gatk",
        java_options="-Xmx4g",
        tracks={"segdup": Path("segdup.bed"), "simple_repeat": Path("repeat.bed")},
    )
    assert command[:4] == ["gatk", "--java-options", "-Xmx4g", "SVRegionOverlap"]
    assert command.count("--track-name") == 2
    assert command.count("--track-intervals") == 2


def test_run_preprocess_without_tracks(tmp_path, monkeypatch) -> None:
    calls = []

    def fake_get_gatk_version(_: str) -> str:
        return "4.6.0.0"

    def fake_run_command(command, timeout_seconds=3600):
        del timeout_seconds
        calls.append(list(command))
        output_path = Path(command[command.index("-O") + 1])
        contig = command[command.index("-L") + 1] if "-L" in command else "chr1"
        _write_bgzipped_vcf(output_path, contig=contig, variant_id=output_path.stem)

        class Result:
            returncode = 0
            stdout = ""
            stderr = ""

        return Result()

    monkeypatch.setattr("gatk_sv_compare.preprocess.get_gatk_version", fake_get_gatk_version)
    monkeypatch.setattr("gatk_sv_compare.preprocess.run_command", fake_run_command)

    reference_dict = tmp_path / "ref.dict"
    reference_dict.write_text("@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:2000\n")
    config = AnalysisConfig(
        vcf_a_path=tmp_path / "a.vcf.gz",
        vcf_b_path=tmp_path / "b.vcf.gz",
        output_dir=tmp_path / "out",
        reference_dict=reference_dict,
        contigs=["chr1", "chr2"],
        contig_lengths={"chr1": 1000, "chr2": 2000},
        n_workers=2,
    )

    annotated_a, annotated_b = run_preprocess(config)

    assert annotated_a.exists()
    assert annotated_b.exists()
    assert (tmp_path / "out" / "preprocess" / "concordance_a.vcf.gz").exists()
    assert len(calls) == 4
