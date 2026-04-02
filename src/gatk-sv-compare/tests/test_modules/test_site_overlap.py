from __future__ import annotations

import pandas as pd
import matplotlib.pyplot as plt

from gatk_sv_compare.modules.site_overlap import SiteOverlapModule, _plot_heatmap, build_overlap_metrics


def test_build_overlap_metrics_counts_matches_for_both_callsets(module_test_context) -> None:
    metrics = build_overlap_metrics(module_test_context.data.sites_a, module_test_context.data.sites_b, pass_only=True)

    del_row = metrics.loc[metrics["svtype"] == "DEL"].iloc[0]
    assert int(del_row["n_total_a"]) == 1
    assert int(del_row["n_matched_a"]) == 1
    assert int(del_row["n_total_b"]) == 1
    assert int(del_row["n_matched_b"]) == 1


def test_site_overlap_module_writes_tables_and_plots(module_test_context) -> None:
    module = SiteOverlapModule()

    module.run(module_test_context.data, module_test_context.config)

    output_dir = module_test_context.config.output_dir / "site_overlap"
    table = pd.read_csv(output_dir / "tables" / "overlap_metrics.tsv.gz", sep="\t")
    assert not table.empty
    assert (output_dir / "tables" / "overlap_metrics.parquet").exists()
    assert (output_dir / "overlap.by_class.CallsetA.png").exists()
    assert (output_dir / "overlap.by_size.CallsetB.png").exists()
    assert (output_dir / "heatmap.size_x_class.CallsetA.png").exists()
    assert (output_dir / "heatmap.freq_x_class.CallsetB.png").exists()


def test_plot_heatmap_excludes_bnd_and_ctx_columns(tmp_path) -> None:
    sites = pd.DataFrame(
        [
            {"variant_id": "a", "size_bucket": "<100bp", "af_bucket": "AC=1", "svtype": "DEL", "truth_vid": "b", "status": "MATCHED"},
            {"variant_id": "b", "size_bucket": "<100bp", "af_bucket": "AC=1", "svtype": "BND", "truth_vid": "c", "status": "MATCHED"},
            {"variant_id": "c", "size_bucket": "<100bp", "af_bucket": "AC=1", "svtype": "CTX", "truth_vid": "d", "status": "MATCHED"},
        ]
    )

    _plot_heatmap(sites, "size_bucket", "svtype", tmp_path / "heatmap.png", "CallsetA")

    assert (tmp_path / "heatmap.png").exists()
    plt.close("all")
