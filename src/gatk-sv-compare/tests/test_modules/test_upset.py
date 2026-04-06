from __future__ import annotations

import pandas as pd
from matplotlib.figure import Figure

from gatk_sv_compare.modules.upset import UpSetModule, _MEMBERSHIP_FIELDS, _build_upset_series, _plot_upset, build_membership_count_table


def test_build_membership_count_table_aggregates_algorithm_combinations(module_test_context) -> None:
    table = build_membership_count_table(module_test_context.data.sites_a, "algorithms", pass_only=True)

    observed = dict(zip(table["combination"], table["n_variants"]))
    assert observed == {"manta,wham": 1, "melt": 1}


def test_build_membership_count_table_aggregates_evidence_combinations(module_test_context) -> None:
    table = build_membership_count_table(module_test_context.data.sites_a, "evidence", pass_only=True)

    observed = dict(zip(table["combination"], table["n_variants"]))
    assert observed == {"RD,PE": 1, "SR": 1}


def test_upset_module_writes_expected_plots_and_tables(module_test_context) -> None:
    module = UpSetModule()

    module.run(module_test_context.data, module_test_context.config)

    output_dir = module_test_context.config.output_dir / "upset"
    expected_files = [
        output_dir / "algorithms.upset.CallsetA.png",
        output_dir / "algorithms.upset.CallsetB.png",
        output_dir / "evidence.upset.CallsetA.png",
        output_dir / "evidence.upset.CallsetB.png",
        output_dir / "tables" / "algorithms_combinations.CallsetA.tsv.gz",
        output_dir / "tables" / "algorithms_combinations.CallsetB.tsv.gz",
        output_dir / "tables" / "evidence_combinations.CallsetA.tsv.gz",
        output_dir / "tables" / "evidence_combinations.CallsetB.tsv.gz",
    ]

    for path in expected_files:
        assert path.exists()

    algorithms_table = pd.read_csv(output_dir / "tables" / "algorithms_combinations.CallsetA.tsv.gz", sep="\t")
    evidence_table = pd.read_csv(output_dir / "tables" / "evidence_combinations.CallsetA.tsv.gz", sep="\t")
    assert not algorithms_table.empty
    assert not evidence_table.empty


def test_upset_module_titles_describe_observed_non_empty_subsets(module_test_context, monkeypatch) -> None:
    module = UpSetModule()
    captured = []
    original_savefig = Figure.savefig

    def spy_savefig(self, *args, **kwargs):
        suptitle = self._suptitle.get_text() if self._suptitle is not None else ""
        captured.append(suptitle)
        return original_savefig(self, *args, **kwargs)

    monkeypatch.setattr(Figure, "savefig", spy_savefig)

    module.run(module_test_context.data, module_test_context.config)

    assert any("observed combinations (all non-empty subsets)" in title for title in captured)


def test_build_upset_series_uses_multiindex_for_single_observed_category() -> None:
    sites = pd.DataFrame(
        {
            "algorithms": ["melt", "melt"],
            "evidence_bucket": ["SR", "SR"],
            "in_filtered_pass_view": [True, True],
        }
    )

    field = next(field for field in _MEMBERSHIP_FIELDS if field.name == "algorithms")
    series = _build_upset_series(sites, field, pass_only=True)

    assert isinstance(series.index, pd.MultiIndex)
    assert series.index.names == ["melt"]
    assert series.tolist() == [2]


def test_plot_upset_handles_single_observed_category(tmp_path) -> None:
    sites = pd.DataFrame(
        {
            "algorithms": ["melt", "melt"],
            "evidence_bucket": ["SR", "SR"],
            "in_filtered_pass_view": [True, True],
        }
    )

    field = next(field for field in _MEMBERSHIP_FIELDS if field.name == "algorithms")
    output_path = tmp_path / "single-category-upset.png"

    _plot_upset(sites, field, output_path, "SingleCategory", pass_only=True)

    assert output_path.exists()
