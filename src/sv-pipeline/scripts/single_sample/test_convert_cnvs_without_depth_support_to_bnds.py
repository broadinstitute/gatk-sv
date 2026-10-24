"""Tests for the proband-union record rule in convert_cnvs_without_depth_support_to_bnds.py.

These cover what the trio de novo change to the record loop introduced: retention (keep / convert to
BND / drop) is decided over ALL proband samples rather than over the case alone, plus the check that a
requested proband name which is not a column in the VCF aborts the run.
"""

import importlib.util
import sys
import types
from pathlib import Path

import pysam
import pytest

MODULE_PATH = Path(__file__).with_name("convert_cnvs_without_depth_support_to_bnds.py")

SAMPLES = ["case", "mother", "father"]

# Pedigree sexes stand in for the parsed famfile. Every test record sits on chr1, which is absent from
# the allosome contigs file, so the script always takes the autosomal depth rule and never consults a
# sample's sex.
SAMPLE_SEXES = {"case": "1", "mother": "2", "father": "1"}

MIN_SIZE = 1000

HEADER_LINES = [
    '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Variant type">',
    '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of the event">',
    '##INFO=<ID=CHR2,Number=1,Type=String,Description="Chromosome for END2 coordinate">',
    '##INFO=<ID=END2,Number=1,Type=Integer,Description="End position of the event on chromosome CHR2">',
    '##INFO=<ID=STRANDS,Number=1,Type=String,Description="Strandness of the two breakends">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=RD_CN,Number=1,Type=Integer,Description="RD genotype (discrete)">',
    '##FORMAT=<ID=PE_GT,Number=1,Type=Integer,Description="PE genotype (discrete)">',
    '##FORMAT=<ID=SR_GT,Number=1,Type=Integer,Description="SR genotype (discrete)">',
]


def _load_script():
    """Load the converter by path with svtk.famfile stubbed out.

    The script imports svtk.famfile at module level, but svtk is not installed in this test
    environment (it only exists inside the sv_pipeline docker image), and nothing exercised here calls
    it: parse_famfile is used only to look up sample sexes for the allosome depth rule, which these
    autosome-only fixtures never reach. The stub stands in for that unexercised module and is removed
    from sys.modules again as soon as the script has been loaded.
    """
    saved = {name: sys.modules.get(name) for name in ("svtk", "svtk.famfile")}

    class _PedSample:
        def __init__(self, sex):
            self.sex = sex

    class _Ped:
        def __init__(self, samples):
            self.samples = samples

    svtk_module = types.ModuleType("svtk")
    famfile_module = types.ModuleType("svtk.famfile")
    famfile_module.parse_famfile = lambda _filehandle: _Ped(
        {name: _PedSample(sex) for name, sex in SAMPLE_SEXES.items()})
    svtk_module.famfile = famfile_module
    sys.modules["svtk"] = svtk_module
    sys.modules["svtk.famfile"] = famfile_module
    try:
        spec = importlib.util.spec_from_file_location(
            "convert_cnvs_without_depth_support_to_bnds", MODULE_PATH)
        assert spec is not None
        assert spec.loader is not None
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    finally:
        for name, previous in saved.items():
            if previous is None:
                del sys.modules[name]
            else:
                sys.modules[name] = previous
    return module


convert_cnvs = _load_script()


def _write_input_vcf(path: Path, records):
    """Write a genotyped PESR-style VCF from the per-record specs built by _del_spec.

    A spec is {"id", "pos", "svlen", "gts": {sample: GT}, "formats": {sample: {field: value}}}. pos is
    1-based. SVLEN is the only coordinate INFO field set: pysam/htslib fill in END = POS + SVLEN for the
    symbolic allele, so the record.stop/record.start that the loop reads are fixed by pos and svlen, and
    the coordinates are read back here to confirm that. A sample with no entry in "formats" carries "."
    for RD_CN/PE_GT/SR_GT, which the script reads as None.
    """
    header = pysam.VariantHeader()
    header.contigs.add("chr1", length=1000000)
    for line in HEADER_LINES:
        header.add_line(line)
    for sample in SAMPLES:
        header.add_sample(sample)

    with pysam.VariantFile(str(path), "w", header=header) as vcf:
        for spec in records:
            record = vcf.new_record(contig="chr1", start=spec["pos"] - 1, id=spec["id"],
                                    alleles=("N", "<DEL>"))
            record.info["SVTYPE"] = "DEL"
            record.info["SVLEN"] = spec["svlen"]
            for sample in SAMPLES:
                record.samples[sample]["GT"] = spec["gts"].get(sample, (0, 0))
                for field, value in spec.get("formats", {}).get(sample, {}).items():
                    record.samples[sample][field] = value
            vcf.write(record)

    with pysam.VariantFile(str(path)) as vcf:
        coords = {record.id: (record.pos, record.stop) for record in vcf}
    for spec in records:
        pos, stop = coords[spec["id"]]
        assert (pos, stop) == (spec["pos"], spec["pos"] + spec["svlen"])


def _run_converter(monkeypatch, tmp_path, records, probands):
    """Run the script's main() over the given records; returns {record id: snapshot} of the output."""
    vcf_path = tmp_path / "genotyped_pesr.vcf"
    _write_input_vcf(vcf_path, records)

    # Autosome-only allosome contigs file: it names a contig that none of the fixtures use, so every
    # record deterministically takes the autosomal depth rule.
    allosome_path = tmp_path / "allosome_contigs.txt"
    allosome_path.write_text("chrY\n")

    ped_path = tmp_path / "trio.ped"
    ped_path.write_text("".join(
        "TRIO\t{}\t0\t0\t{}\t.\n".format(name, SAMPLE_SEXES[name]) for name in SAMPLES))

    out_path = tmp_path / "out.vcf"
    monkeypatch.setattr("sys.argv", [
        str(MODULE_PATH), str(vcf_path), str(allosome_path), str(ped_path),
        *probands, str(MIN_SIZE), "-o", str(out_path),
    ])
    convert_cnvs.main()

    with pysam.VariantFile(str(out_path)) as out:
        return {
            record.id: {
                "info": dict(record.info),
                "alts": list(record.alts or []),
                "samples": {sample: dict(record.samples[sample]) for sample in SAMPLES},
            }
            for record in out
        }


def _del_spec(record_id, pos=1001, svlen=2000, **evidence):
    """One DEL record at/over min_size with per-sample RD_CN/PE_GT/SR_GT evidence.

    Evidence is given as <sample>_<field> keyword arguments, e.g. mother_rd_cn=1 or father_pe_gt=3; a
    field that is absent (or None) is left out of the record entirely, i.e. a missing genotype.
    """
    formats = {sample: {} for sample in SAMPLES}
    for key, value in evidence.items():
        sample, field = key.split("_", 1)
        assert sample in SAMPLES and field in ("rd_cn", "pe_gt", "sr_gt"), key
        if value is not None:
            formats[sample][field.upper()] = value
    return {"id": record_id, "pos": pos, "svlen": svlen,
            "gts": {sample: (0, 1) for sample in SAMPLES}, "formats": formats}


def test_parent_depth_support_keeps_case_record_as_del(monkeypatch, tmp_path):
    # Real discriminator for the union rule: the case has no RD_CN and no proband carries any PE/SR
    # evidence, so under the old case-only rule this record was DROPPED. Counting the mother's RD_CN=1
    # as depth support, it is retained, and unchanged as a DEL.
    out = _run_converter(monkeypatch, tmp_path,
                         [_del_spec("del_case_no_depth_mother_depth", mother_rd_cn=1)],
                         SAMPLES)

    assert list(out) == ["del_case_no_depth_mother_depth"]
    kept = out["del_case_no_depth_mother_depth"]
    assert kept["info"]["SVTYPE"] == "DEL"
    assert kept["alts"] == ["<DEL>"]
    # unchanged: nothing was added by a BND rewrite
    assert "CHR2" not in kept["info"]
    assert "END2" not in kept["info"]
    assert "STRANDS" not in kept["info"]
    assert kept["samples"]["mother"].get("RD_CN") == 1
    assert kept["samples"]["case"].get("RD_CN") is None


def test_no_depth_and_no_pesr_support_in_any_proband_is_dropped(monkeypatch, tmp_path):
    # Two records to show the decision is per record: the first has no depth genotype and no PE/SR
    # support in any proband and is dropped, the second is supported by the case's own depth.
    out = _run_converter(monkeypatch, tmp_path,
                         [_del_spec("del_no_evidence_anywhere"),
                          _del_spec("del_case_depth", pos=5001, case_rd_cn=1)],
                         SAMPLES)

    assert list(out) == ["del_case_depth"]


def test_parent_pesr_support_converts_record_to_bnd(monkeypatch, tmp_path):
    # No proband has depth support and only the father has PE evidence, so the conversion also
    # depends on a parent: case-only this record was neither supported nor converted, it was dropped.
    out = _run_converter(monkeypatch, tmp_path,
                         [_del_spec("del_father_pe_only", pos=9001, svlen=2000, father_pe_gt=3)],
                         SAMPLES)

    assert list(out) == ["del_father_pe_only"]
    converted = out["del_father_pe_only"]
    assert converted["info"]["SVTYPE"] == "BND"
    assert converted["alts"] == ["<BND>"]
    assert converted["info"]["STRANDS"] == "+-"
    # the DEL's span moves into the record-level CHR2/END2/SVLEN INFO fields
    assert converted["info"]["CHR2"] == "chr1"
    assert converted["info"]["END2"] == 11001          # pos + svlen, i.e. the DEL's END
    assert converted["info"]["SVLEN"] == 2001          # END2 - record.start, as the script computes it


def test_proband_sample_missing_from_vcf_errors(monkeypatch, tmp_path):
    # A mistyped parent name must not silently narrow the decision rule to the case alone.
    with pytest.raises(SystemExit) as excinfo:
        _run_converter(monkeypatch, tmp_path,
                       [_del_spec("del_case_no_depth_mother_depth", mother_rd_cn=1)],
                       ["case", "motherr"])

    message = str(excinfo.value)
    assert "motherr" in message
    assert "case" in message                        # the VCF's actual samples are named too
    assert (tmp_path / "out.vcf").exists() is False
