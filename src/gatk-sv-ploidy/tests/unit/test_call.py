from __future__ import annotations

import gzip

import pandas as pd

from gatk_sv_ploidy.call import (
    _classify_aneuploidy,
    _classify_sex,
    assign_sex_and_aneuploidy_types,
    export_aneuploid_data,
    save_sex_assignments,
)


def test_classify_sex_and_aneuploidy_tables() -> None:
    assert _classify_sex(1, 1) == "MALE"
    assert _classify_sex(2, 0) == "FEMALE"
    assert _classify_sex(3, 0) == "TRIPLE_X"
    assert _classify_sex(8, 1) == "OTHER"

    assert _classify_aneuploidy({}, {}, 2, 0) == "NORMAL"
    assert _classify_aneuploidy({"chr21": True}, {"chr21": 3}, 2, 0) == "TRISOMY_21"
    assert _classify_aneuploidy({"chrX": True}, {"chrX": 1, "chrY": 0}, 1, 0) == "TURNER"
    assert _classify_aneuploidy(
        {"chr13": True, "chr18": True},
        {"chr13": 3, "chr18": 3},
        2,
        0,
    ) == "MULTIPLE"


def test_assign_and_export_outputs(tmp_path) -> None:
    df = pd.DataFrame(
        [
            {
                "sample": "S1",
                "chromosome": "chrX",
                "copy_number": 2,
                "is_aneuploid": False,
                "mean_cn_probability": 0.95,
                "median_depth": 2.0,
                "sample_var_map": 0.04,
                "chr_type": 1,
            },
            {
                "sample": "S1",
                "chromosome": "chrY",
                "copy_number": 0,
                "is_aneuploid": False,
                "mean_cn_probability": 0.93,
                "median_depth": 0.0,
                "sample_var_map": 0.04,
                "chr_type": 2,
            },
            {
                "sample": "S1",
                "chromosome": "chr21",
                "copy_number": 2,
                "is_aneuploid": False,
                "mean_cn_probability": 0.97,
                "median_depth": 2.0,
                "sample_var_map": 0.04,
                "chr_type": 0,
            },
            {
                "sample": "S2",
                "chromosome": "chrX",
                "copy_number": 1,
                "is_aneuploid": True,
                "mean_cn_probability": 0.91,
                "median_depth": 1.0,
                "sample_var_map": 0.09,
                "chr_type": 1,
            },
            {
                "sample": "S2",
                "chromosome": "chrY",
                "copy_number": 0,
                "is_aneuploid": True,
                "mean_cn_probability": 0.92,
                "median_depth": 0.0,
                "sample_var_map": 0.09,
                "chr_type": 2,
            },
            {
                "sample": "S2",
                "chromosome": "chr21",
                "copy_number": 3,
                "is_aneuploid": True,
                "mean_cn_probability": 0.90,
                "median_depth": 3.0,
                "sample_var_map": 0.09,
                "chr_type": 0,
            },
        ]
    )

    pred_df = assign_sex_and_aneuploidy_types(df, {"S2": "MULTIPLE"})

    assert pred_df.set_index("sample").loc["S1", "predicted_aneuploidy_type"] == "NORMAL"
    assert pred_df.set_index("sample").loc["S2", "predicted_aneuploidy_type"] == "MULTIPLE"

    save_sex_assignments(pred_df, str(tmp_path))
    with gzip.open(tmp_path / "sex_assignments.txt.gz", "rt") as handle:
        written = handle.read()
    assert "sample_id" in written
    assert "Assignment" in written

    export_aneuploid_data(df, str(tmp_path / "aneuploid.tsv"))
    exported = pd.read_csv(tmp_path / "aneuploid.tsv", sep="\t")
    assert "inferred_sample_std" in exported.columns
    assert "chr_type" not in exported.columns
