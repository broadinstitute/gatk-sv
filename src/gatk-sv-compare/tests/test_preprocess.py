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


def _assert_vcf_has_no_str(path: Path) -> None:
    with pysam.VariantFile(str(path)) as vcf:
        for record in vcf:
            assert str(record.info.get("SVTYPE", "")) != "STR"
            assert record.alts != ("<STR>",)


def _write_mixed_str_sv_vcf(path: Path, *, str_id: str, del_id: str, repcn: str) -> None:
    header = pysam.VariantHeader()
    header.add_meta("fileformat", "VCFv4.2")
    header.contigs.add("chr1", length=1000)
    header.add_meta("INFO", items=[("ID", "SVTYPE"), ("Number", 1), ("Type", "String"), ("Description", "SV type")])
    header.add_meta("INFO", items=[("ID", "SVLEN"), ("Number", 1), ("Type", "Integer"), ("Description", "SV length")])
    header.add_meta("INFO", items=[("ID", "END"), ("Number", 1), ("Type", "Integer"), ("Description", "End")])
    header.add_meta("INFO", items=[("ID", "RU"), ("Number", 1), ("Type", "String"), ("Description", "Repeat unit")])
    header.add_meta("INFO", items=[("ID", "PERIOD"), ("Number", 1), ("Type", "Integer"), ("Description", "Repeat period")])
    header.add_meta("INFO", items=[("ID", "LOCUS"), ("Number", 1), ("Type", "String"), ("Description", "STR locus")])
    header.add_meta("FORMAT", items=[("ID", "GT"), ("Number", 1), ("Type", "String"), ("Description", "Genotype")])
    header.add_meta("FORMAT", items=[("ID", "REPCN"), ("Number", "."), ("Type", "Integer"), ("Description", "Repeat copy number")])
    header.add_sample("S1")
    with pysam.VariantFile(str(path), "wz", header=header) as out_vcf:
        del_record = out_vcf.new_record(contig="chr1", start=99, stop=199, alleles=("N", "<DEL>"), id=del_id, filter="PASS")
        del_record.info["SVTYPE"] = "DEL"
        del_record.info["SVLEN"] = 100
        del_record.samples["S1"]["GT"] = (0, 1)
        out_vcf.write(del_record)
        str_record = out_vcf.new_record(contig="chr1", start=299, stop=350, alleles=("N", "<STR>"), id=str_id, filter="PASS")
        str_record.info["SVTYPE"] = "STR"
        str_record.info["RU"] = "AC"
        str_record.info["PERIOD"] = 2
        str_record.info["LOCUS"] = "L1"
        str_record.samples["S1"]["GT"] = (0, 1)
        str_record.samples["S1"]["REPCN"] = tuple(int(value) for value in repcn.split(","))
        out_vcf.write(str_record)
    pysam.tabix_index(str(path), preset="vcf", force=True)


def test_read_contig_list_and_parse_reference_dict(tmp_path) -> None:
    contig_list = tmp_path / "contigs.list"
    contig_list.write_text("# comment\nchr1\n\nchr2\n")
    reference_dict = tmp_path / "ref.dict"
    reference_dict.write_text("@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:2000\n")

    assert read_contig_list(contig_list) == ["chr1", "chr2"]
    assert parse_reference_dict(reference_dict) == {"chr1": 1000, "chr2": 2000}


def test_build_svconcordance_command_omits_tracks() -> None:
    command = build_svconcordance_command(
        eval_vcf=Path("eval.vcf.gz"),
        truth_vcf=Path("truth.vcf.gz"),
        contig="chr1",
        reference_dict=Path("ref.dict"),
        output_path=Path("out.vcf.gz"),
        gatk_path="gatk",
        java_options="-Xmx4g",
    )
    assert command[:4] == ["gatk", "--java-options", "-Xmx4g", "SVConcordance"]
    assert ["-L", "chr1"] == command[4:6]
    assert "--track-name" not in command
    assert "--track-intervals" not in command


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
    _write_bgzipped_vcf(tmp_path / "a.vcf.gz", "chr1", "a_chr1")
    _write_bgzipped_vcf(tmp_path / "b.vcf.gz", "chr1", "b_chr1")
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


def test_run_preprocess_logs_stage_progress(tmp_path, monkeypatch, caplog) -> None:
    def fake_get_gatk_version(_: str) -> str:
        return "4.6.0.0"

    def fake_run_command(command, timeout_seconds=3600):
        del timeout_seconds
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
    _write_bgzipped_vcf(tmp_path / "a.vcf.gz", "chr1", "a_chr1")
    _write_bgzipped_vcf(tmp_path / "b.vcf.gz", "chr1", "b_chr1")
    config = AnalysisConfig(
        vcf_a_path=tmp_path / "a.vcf.gz",
        vcf_b_path=tmp_path / "b.vcf.gz",
        output_dir=tmp_path / "out",
        reference_dict=reference_dict,
        contigs=["chr1", "chr2"],
        contig_lengths={"chr1": 1000, "chr2": 2000},
        n_workers=2,
    )

    caplog.set_level("INFO")
    run_preprocess(config)

    assert "Starting preprocess into" in caplog.text
    assert "Starting SVConcordance scatter for 2 contigs" in caplog.text
    assert "Completed SVConcordance shard 2/2" in caplog.text
    assert "Copying concordance outputs to annotated outputs" in caplog.text
    assert "Preprocess complete:" in caplog.text


def test_run_preprocess_bypasses_and_merges_str_records(tmp_path, monkeypatch) -> None:
    calls = []

    def fake_get_gatk_version(_: str) -> str:
        return "4.6.0.0"

    def fake_run_command(command, timeout_seconds=3600):
        del timeout_seconds
        calls.append(list(command))
        eval_vcf = Path(command[command.index("--eval") + 1])
        truth_vcf = Path(command[command.index("--truth") + 1])
        _assert_vcf_has_no_str(eval_vcf)
        _assert_vcf_has_no_str(truth_vcf)
        output_path = Path(command[command.index("-O") + 1])
        with pysam.VariantFile(str(eval_vcf)) as in_vcf:
            header = in_vcf.header.copy()
            with pysam.VariantFile(str(output_path), "wz", header=header) as out_vcf:
                for record in in_vcf:
                    out_vcf.write(record)
        pysam.tabix_index(str(output_path), preset="vcf", force=True)

        class Result:
            returncode = 0
            stdout = ""
            stderr = ""

        return Result()

    monkeypatch.setattr("gatk_sv_compare.preprocess.get_gatk_version", fake_get_gatk_version)
    monkeypatch.setattr("gatk_sv_compare.preprocess.run_command", fake_run_command)

    reference_dict = tmp_path / "ref.dict"
    reference_dict.write_text("@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n")
    vcf_a = tmp_path / "a.vcf.gz"
    vcf_b = tmp_path / "b.vcf.gz"
    _write_mixed_str_sv_vcf(vcf_a, str_id="a_str", del_id="a_del", repcn="10,11")
    _write_mixed_str_sv_vcf(vcf_b, str_id="b_str", del_id="b_del", repcn="10,11")
    config = AnalysisConfig(
        vcf_a_path=vcf_a,
        vcf_b_path=vcf_b,
        output_dir=tmp_path / "out",
        reference_dict=reference_dict,
        contigs=["chr1"],
        contig_lengths={"chr1": 1000},
        n_workers=1,
    )

    annotated_a, annotated_b = run_preprocess(config)

    assert len(calls) == 2
    with pysam.VariantFile(str(annotated_a)) as vcf:
        records = {record.id: record for record in vcf}
        assert {"a_del", "a_str"}.issubset(records)
        assert str(records["a_str"].info["SVTYPE"]) == "STR"
        assert str(records["a_str"].info["STATUS"]) == "MATCHED"
        assert str(records["a_str"].info["TRUTH_VID"]) == "b_str"
    with pysam.VariantFile(str(annotated_b)) as vcf:
        records = {record.id: record for record in vcf}
        assert {"b_del", "b_str"}.issubset(records)
        assert str(records["b_str"].info["STATUS"]) == "MATCHED"
        assert str(records["b_str"].info["TRUTH_VID"]) == "a_str"
