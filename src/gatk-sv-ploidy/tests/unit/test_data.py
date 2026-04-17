from __future__ import annotations

import numpy as np
import pandas as pd
import torch

import pytest

from gatk_sv_ploidy.data import DepthData, load_site_data, read_depth_tsv


def test_read_depth_tsv_normalizes_chr_and_adds_metadata(tmp_path) -> None:
    path = tmp_path / "sample" / "depth.tsv"
    path.parent.mkdir(parents=True)
    path.write_text(
        "#Chr\tStart\tEnd\tS1\tS2\n"
        "chr21\t0\t100\t12\t15\n"
        "chrX\t100\t200\t7\t8\n"
    )

    df = read_depth_tsv(str(path))

    assert list(df.columns[:5]) == ["Chr", "Start", "End", "S1", "S2"]
    assert df.index.tolist() == ["chr21:0-100", "chrX:100-200"]
    assert df["source_file"].nunique() == 1


def test_depth_data_sorts_clamps_and_reorders_site_data(
    tiny_depth_df: pd.DataFrame,
    tiny_site_data: dict[str, np.ndarray],
) -> None:
    data = DepthData(
        tiny_depth_df,
        dtype=torch.float32,
        clamp_threshold=2.0,
        site_data=tiny_site_data,
    )

    assert data.chr.tolist() == ["chr18", "chr21", "chrX", "chrY"]
    assert data.sample_ids == ["SAMPLE_B", "SAMPLE_A"]
    assert float(data.depth.max()) == pytest.approx(2.0)
    assert data.site_alt.shape == (4, 3, 2)
    assert data.max_sites == 3


def test_depth_data_subsamples_deterministically(tiny_depth_df: pd.DataFrame) -> None:
    data_one = DepthData(
        tiny_depth_df,
        subsample_bins=2,
        subsample_samples=1,
        seed=11,
        clamp_threshold=None,
    )
    data_two = DepthData(
        tiny_depth_df,
        subsample_bins=2,
        subsample_samples=1,
        seed=11,
        clamp_threshold=None,
    )

    np.testing.assert_array_equal(data_one.chr, data_two.chr)
    assert data_one.sample_ids == data_two.sample_ids
    torch.testing.assert_close(data_one.depth, data_two.depth)


def test_depth_data_rejects_misaligned_site_data(
    tiny_depth_df: pd.DataFrame,
    tiny_site_data: dict[str, np.ndarray],
) -> None:
    bad_site_data = dict(tiny_site_data)
    bad_site_data["bin_start"] = bad_site_data["bin_start"].copy()
    bad_site_data["bin_start"][0] += 1

    with pytest.raises(ValueError, match="bin_start"):
        DepthData(tiny_depth_df, site_data=bad_site_data)


def test_load_site_data_round_trips_metadata(tmp_path) -> None:
    np.savez_compressed(
        tmp_path / "site_data.npz",
        site_alt=np.ones((2, 3, 1), dtype=np.int32),
        site_total=np.ones((2, 3, 1), dtype=np.int32) * 2,
        site_pop_af=np.full((2, 3), 0.25, dtype=np.float32),
        site_mask=np.ones((2, 3, 1), dtype=bool),
        bin_chr=np.array(["chr18", "chr21"], dtype=object),
        bin_start=np.array([0, 100], dtype=np.int64),
        bin_end=np.array([100, 200], dtype=np.int64),
    )

    loaded = load_site_data(str(tmp_path / "site_data.npz"))

    assert set(loaded) == {
        "site_alt",
        "site_total",
        "site_pop_af",
        "site_mask",
        "bin_chr",
        "bin_start",
        "bin_end",
    }
    assert loaded["site_alt"].shape == (2, 3, 1)
