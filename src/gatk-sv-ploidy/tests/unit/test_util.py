from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from gatk_sv_ploidy._util import (
    add_chromosome_labels,
    compute_contig_posterior_from_bin_posteriors,
    compute_plq_from_probabilities,
    format_count_summary,
    format_column_name,
    format_numeric_summary,
    get_chromosome_type,
    get_sample_columns,
    load_exclusion_ids,
    safe_path_label,
    save_and_close_plot,
)


def test_dataframe_and_format_helpers() -> None:
    df = pd.DataFrame(
        {
            "Chr": ["chr1"],
            "Start": [0],
            "End": [100],
            "source_file": ["x"],
            "Bin": ["chr1:0-100"],
            "S1": [1.0],
            "S2": [2.0],
        }
    )

    assert get_sample_columns(df) == ["S1", "S2"]
    assert get_chromosome_type("chrX") == "chrX"
    assert get_chromosome_type("Y") == "chrY"
    assert get_chromosome_type("chr8") == "Autosomal"
    assert format_column_name("mean_depth") == "Mean Depth"


def test_load_exclusion_ids(tmp_path) -> None:
    path = tmp_path / "items.txt"
    path.write_text("a\n\n b \n")

    assert load_exclusion_ids(str(path)) == ["a", "b"]


def test_safe_log_helpers_format_summaries_without_raw_paths() -> None:
    assert safe_path_label("/private/tmp/example/file.tsv.gz") == "file.tsv.gz"

    numeric = format_numeric_summary("Example", [1.0, 2.0, 3.0, 4.0])
    categorical = format_count_summary("States", [1, 3], ["A", "B"])

    assert numeric.startswith("Example: n=4")
    assert "median=" in numeric
    assert categorical == "States: total=4, A=1 (25.0%), B=3 (75.0%)"


def test_contig_posterior_and_plq_helpers() -> None:
    bin_posteriors = np.array(
        [
            [0.0, 0.0, 0.51, 0.49, 0.0, 0.0],
            [0.0, 0.0, 0.01, 0.99, 0.0, 0.0],
        ],
        dtype=np.float64,
    )

    contig_posterior = compute_contig_posterior_from_bin_posteriors(
        bin_posteriors,
    )
    plq = compute_plq_from_probabilities(contig_posterior)

    assert int(np.argmax(contig_posterior)) == 3
    assert contig_posterior[3] > 0.98
    assert int(plq) == 20


def test_plot_helpers_write_files_and_labels(tmp_path) -> None:
    fig, ax = plt.subplots()
    ax.plot([0, 1], [0, 1])
    add_chromosome_labels(ax, np.array(["chr1", "chr1", "chr2", "chr2"]))
    save_and_close_plot(str(tmp_path), "plot.png", subdir="plots")
    assert (tmp_path / "plots" / "plot.png").exists()

    fig, ax = plt.subplots()
    add_chromosome_labels(
        ax,
        np.array(["chr1", "chr1", "chr2", "chr2"]),
        x_transformed=np.array([0.0, 0.5, 1.0, 1.5]),
    )
    assert [tick.get_text() for tick in ax.get_xticklabels()] == ["1", "2"]
    plt.close(fig)
