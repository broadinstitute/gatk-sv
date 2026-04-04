from __future__ import annotations

import pandas as pd
import pytest
from matplotlib.figure import Figure

from gatk_sv_compare.aggregate import aggregate
from gatk_sv_compare.config import AnalysisConfig
from gatk_sv_compare.modules.family_analysis import (
    FamilyAnalysisModule,
    _plot_dnr_heatmap,
    build_transmission_table,
    parse_ped_file,
    summarize_inheritance_stats,
)


def test_parse_ped_file_filters_to_present_members(tmp_path) -> None:
    ped = tmp_path / "test.fam"
    ped.write_text("F1\tPRO\tDAD\tMOM\t1\t1\nF2\tMISSING\t0\t0\t1\t1\n")

    families = parse_ped_file(ped, {"PRO", "DAD", "MOM"})

    assert len(families) == 1
    assert families[0].family_type == "trios"


def test_family_analysis_module_computes_expected_denovo_rate(tmp_path, make_vcf) -> None:
    extra_headers = [
        "##INFO=<ID=AC,Number=1,Type=Integer,Description=\"Allele count\">",
        "##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele frequency\">",
        "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Allele number\">",
    ]
    records = [
        "chr1\t100\tv1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=100;AC=2;AF=0.3333;AN=6\tGT:GQ:ECN\t0/1:60:2\t0/0:60:2\t0/1:60:2",
        "chr1\t200\tv2\tN\t<DUP>\t.\tPASS\tSVTYPE=DUP;SVLEN=300;AC=2;AF=0.3333;AN=6\tGT:GQ:ECN\t0/1:55:2\t0/1:55:2\t0/0:55:2",
        "chr1\t300\tv3\tN\t<INS>\t.\tPASS\tSVTYPE=INS;SVLEN=250;AC=1;AF=0.1667;AN=6\tGT:GQ:ECN\t0/1:50:2\t0/0:50:2\t0/0:50:2",
        "chr1\t400\tv4\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=150;AC=1;AF=0.1667;AN=6\tGT:GQ:ECN\t0/1:45:2\t0/0:45:2\t./.:.:2",
    ]
    vcf_a = make_vcf(
        file_name="family_a.vcf",
        sample_names=["PRO", "DAD", "MOM"],
        extra_header_lines=extra_headers,
        records=records,
    )
    vcf_b = make_vcf(
        file_name="family_b.vcf",
        sample_names=["PRO", "DAD", "MOM"],
        extra_header_lines=extra_headers,
        records=records,
    )
    ped = tmp_path / "family.fam"
    ped.write_text("F1\tPRO\tDAD\tMOM\t1\t1\n")

    config = AnalysisConfig(
        vcf_a_path=vcf_a,
        vcf_b_path=vcf_b,
        vcf_a_label="CallsetA",
        vcf_b_label="CallsetB",
        output_dir=tmp_path / "outputs",
        contigs=["chr1"],
        contig_lengths={"chr1": 1_000_000},
        n_workers=1,
        ped_file=ped,
    )
    data = aggregate(config)
    families = parse_ped_file(ped, set(data.sample_names_a))

    transmission = build_transmission_table(config.vcf_a_path, data.sites_a, config.vcf_a_label, families)
    inheritance = summarize_inheritance_stats(transmission)
    row = inheritance.loc[inheritance["label"] == "CallsetA"].iloc[0]
    assert "evidence_bucket" in transmission.columns
    assert row["site_de_novo_rate"] == pytest.approx(1.0 / 3.0)
    assert row["allele_de_novo_rate"] == pytest.approx(1.0 / 3.0)

    module = FamilyAnalysisModule()
    module.run(data, config)

    output_dir = config.output_dir / "family_analysis"
    inheritance_table = pd.read_csv(output_dir / "tables" / "inheritance_stats.trios.tsv.gz", sep="\t")
    by_class = pd.read_csv(output_dir / "tables" / "denovo_rate_by_class.tsv.gz", sep="\t")
    assert not inheritance_table.empty
    assert not by_class.empty
    assert (output_dir / "tables" / "denovo_rate_by_algorithm.tsv.gz").exists()
    assert (output_dir / "tables" / "denovo_rate_by_evidence.tsv.gz").exists()
    assert (output_dir / "inheritance.trios.all_sv.png").exists()
    assert (output_dir / "dnr_vs_algorithm.trios.variants.png").exists()
    assert (output_dir / "dnr_vs_evidence.trios.variants.png").exists()
    assert (output_dir / "dnr_vs_gq.trios.variants.png").exists()
    assert (output_dir / "dnr_heatmap.trios.variants.all_sv.png").exists()


def test_plot_dnr_heatmap_labels_colorbar(tmp_path, monkeypatch) -> None:
    table = pd.DataFrame(
        [
            {"size_bucket": "<100bp", "af_bucket": "AC=1", "site_denovo_rate": 0.25, "n_records": 4},
            {"size_bucket": "100bp-1kb", "af_bucket": "AC=2-5", "site_denovo_rate": 0.50, "n_records": 2},
        ]
    )
    captured = {}
    original_colorbar = Figure.colorbar

    def spy_colorbar(self, *args, **kwargs):
        colorbar = original_colorbar(self, *args, **kwargs)
        captured["colorbar"] = colorbar
        return colorbar

    monkeypatch.setattr(Figure, "colorbar", spy_colorbar)

    _plot_dnr_heatmap(table, tmp_path / "dnr_heatmap.png", "De novo heatmap", "site_denovo_rate")

    assert captured["colorbar"].ax.get_ylabel() == "De novo rate"


def test_plot_dnr_heatmap_includes_count_labels(tmp_path, monkeypatch) -> None:
    table = pd.DataFrame(
        [
            {"size_bucket": "<100bp", "af_bucket": "AC=1", "site_denovo_rate": 0.25, "n_records": 4},
        ]
    )
    captured = {}
    original_savefig = Figure.savefig

    def spy_savefig(self, *args, **kwargs):
        captured["texts"] = [text.get_text() for axis in self.axes for text in axis.texts]
        return original_savefig(self, *args, **kwargs)

    monkeypatch.setattr(Figure, "savefig", spy_savefig)

    _plot_dnr_heatmap(table, tmp_path / "dnr_heatmap_counts.png", "De novo heatmap", "site_denovo_rate")

    assert "n=4" in captured["texts"]
