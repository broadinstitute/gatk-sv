from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

PACKAGE_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = PACKAGE_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

FIXTURE_ROOT = Path(__file__).resolve().parent / "fixtures"
MEDIUM_FIXTURE_ROOT = FIXTURE_ROOT / "medium"
LARGE_FIXTURE_ROOT = FIXTURE_ROOT / "large"
TINY_FIXTURE_ROOT = FIXTURE_ROOT / "tiny"


def get_sample_columns(df: pd.DataFrame) -> list[str]:
    metadata_cols = {"Chr", "Start", "End", "source_file", "Bin"}
    return [column for column in df.columns if column not in metadata_cols]


@pytest.fixture
def fixture_root() -> Path:
    return FIXTURE_ROOT


@pytest.fixture
def medium_fixture_root() -> Path:
    return MEDIUM_FIXTURE_ROOT


@pytest.fixture
def large_fixture_root() -> Path:
    return LARGE_FIXTURE_ROOT


@pytest.fixture
def tiny_depth_df() -> pd.DataFrame:
    df = pd.DataFrame(
        {
            "Chr": ["chr21", "chr18", "chrX", "chrY"],
            "Start": [200, 100, 50, 25],
            "End": [300, 200, 150, 125],
            "SAMPLE_B": [2.1, 1.9, 2.0, 0.0],
            "SAMPLE_A": [1.9, 2.2, 1.0, 1.0],
        }
    )
    df["Bin"] = df["Chr"].astype(str) + ":" + df["Start"].astype(str) + "-" + df["End"].astype(str)
    return df.set_index("Bin")


@pytest.fixture
def tiny_site_data() -> dict[str, np.ndarray]:
    site_alt = np.array(
        [
            [[7, 6], [4, 0], [0, 0]],
            [[8, 7], [5, 0], [0, 0]],
            [[3, 12], [0, 0], [0, 0]],
            [[0, 4], [0, 0], [0, 0]],
        ],
        dtype=np.int32,
    )
    site_total = np.array(
        [
            [[14, 12], [8, 0], [0, 0]],
            [[16, 14], [10, 0], [0, 0]],
            [[6, 24], [0, 0], [0, 0]],
            [[0, 8], [0, 0], [0, 0]],
        ],
        dtype=np.int32,
    )
    site_pop_af = np.array(
        [
            [0.40, 0.25, 0.0],
            [0.45, 0.35, 0.0],
            [0.50, 0.0, 0.0],
            [0.50, 0.0, 0.0],
        ],
        dtype=np.float32,
    )
    site_mask = site_total > 0
    return {
        "site_alt": site_alt,
        "site_total": site_total,
        "site_pop_af": site_pop_af,
        "site_mask": site_mask,
        "bin_chr": np.array(["chr21", "chr18", "chrX", "chrY"], dtype=object),
        "bin_start": np.array([200, 100, 50, 25], dtype=np.int64),
        "bin_end": np.array([300, 200, 150, 125], dtype=np.int64),
    }


@pytest.fixture
def tiny_truth_dict() -> dict[str, str]:
    return {
        "SAMPLE_A": "TURNER",
        "SAMPLE_B": "MALE",
    }


@pytest.fixture
def tiny_truth_json(tmp_path: Path, tiny_truth_dict: dict[str, str]) -> Path:
    out = tmp_path / "truth.json"
    out.write_text(json.dumps(tiny_truth_dict, indent=2, sort_keys=True))
    return out


@pytest.fixture
def medium_expected_truth() -> dict[str, str]:
    path = MEDIUM_FIXTURE_ROOT / "truth.json"
    return json.loads(path.read_text())


@pytest.fixture
def medium_expected_sex() -> dict[str, str]:
    path = MEDIUM_FIXTURE_ROOT / "sex_truth.json"
    return json.loads(path.read_text())


@pytest.fixture
def large_expected_truth() -> dict[str, str]:
    path = LARGE_FIXTURE_ROOT / "truth.json"
    return json.loads(path.read_text())


@pytest.fixture
def large_expected_sex() -> dict[str, str]:
    path = LARGE_FIXTURE_ROOT / "sex_truth.json"
    return json.loads(path.read_text())


@pytest.fixture
def medium_depth_df(medium_fixture_root: Path) -> pd.DataFrame:
    return pd.read_csv(medium_fixture_root / "preprocessed_depth.tsv", sep="\t", index_col=0)


@pytest.fixture
def large_raw_depth_df(large_fixture_root: Path) -> pd.DataFrame:
    return pd.read_csv(large_fixture_root / "raw_depth.tsv.gz", sep="\t", compression="gzip")


@pytest.fixture
def medium_sample_ids(medium_depth_df: pd.DataFrame) -> list[str]:
    return get_sample_columns(medium_depth_df)
