from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import pytest

from gatk_sv_ploidy._util import (
    add_chromosome_labels,
    concatenate_tsvs,
    format_column_name,
    get_chromosome_type,
    get_sample_columns,
    load_exclusion_ids,
    read_file_list,
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


def test_read_file_list_and_load_exclusion_ids(tmp_path) -> None:
    path = tmp_path / "items.txt"
    path.write_text("a\n\n b \n")

    assert read_file_list(str(path)) == ["a", "b"]
    assert load_exclusion_ids(str(path)) == ["a", "b"]


def test_concatenate_tsvs_handles_missing_files(tmp_path, caplog) -> None:
    first = tmp_path / "a.tsv"
    second = tmp_path / "b.tsv"
    first.write_text("x\tvalue\n1\t10\n")
    second.write_text("x\tvalue\n2\t20\n")

    combined = concatenate_tsvs(
        [str(first), str(tmp_path / "missing.tsv"), str(second)],
    )

    assert combined["x"].tolist() == [1, 2]
    assert "File not found" in caplog.text


def test_concatenate_tsvs_raises_when_no_valid_files(tmp_path) -> None:
    with pytest.raises(ValueError, match="No valid TSV files"):
        concatenate_tsvs([str(tmp_path / "missing.tsv")])


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
