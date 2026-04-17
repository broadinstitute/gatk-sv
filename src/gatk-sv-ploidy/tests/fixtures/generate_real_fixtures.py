from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

METADATA_COLS = {"Chr", "Start", "End", "Bin", "source_file"}
TARGET_CHROMS = ["chr13", "chr18", "chr21", "chrX", "chrY"]

MEDIUM_SPEC = {
    "female_count": 3,
    "male_count": 3,
    "max_sites": 8,
    "bins_per_chr": {
        "chr13": 15,
        "chr18": 15,
        "chr21": 15,
        "chrX": 15,
        "chrY": 9,
    },
}

LARGE_SPEC = {
    "female_count": 6,
    "male_count": 6,
    "max_sites": 16,
    "bins_per_chr": {
        "chr13": 40,
        "chr18": 40,
        "chr21": 31,
        "chrX": 40,
        "chrY": 9,
    },
}


def get_sample_columns(df: pd.DataFrame) -> list[str]:
    return [column for column in df.columns if column not in METADATA_COLS]


def choose_sample_ids(call_df: pd.DataFrame, female_count: int, male_count: int) -> list[str]:
    females = sorted(call_df.loc[call_df["sex"] == "FEMALE", "sample"].tolist())
    males = sorted(call_df.loc[call_df["sex"] == "MALE", "sample"].tolist())
    selected = females[:female_count] + males[:male_count]
    if len(selected) != female_count + male_count:
        raise ValueError("Unable to select the requested number of male/female samples")
    return selected


def evenly_spaced_positions(length: int, target_count: int) -> np.ndarray:
    if target_count >= length:
        return np.arange(length, dtype=np.int64)
    return np.unique(np.round(np.linspace(0, length - 1, num=target_count)).astype(np.int64))


def choose_bin_indices(df: pd.DataFrame, bins_per_chr: dict[str, int]) -> np.ndarray:
    selected: list[int] = []
    for chrom in TARGET_CHROMS:
        chrom_idx = np.flatnonzero(df["Chr"].to_numpy() == chrom)
        if chrom_idx.size == 0:
            continue
        rel = evenly_spaced_positions(chrom_idx.size, bins_per_chr[chrom])
        selected.extend(chrom_idx[rel].tolist())
    return np.array(sorted(selected), dtype=np.int64)


def subset_preprocessed(pre_df: pd.DataFrame, sample_ids: list[str], bin_idx: np.ndarray) -> pd.DataFrame:
    keep_cols = ["Chr", "Start", "End", *sample_ids]
    if "source_file" in pre_df.columns:
        keep_cols.append("source_file")
    subset = pre_df.iloc[bin_idx][keep_cols].copy()
    subset["Bin"] = subset["Chr"].astype(str) + ":" + subset["Start"].astype(str) + "-" + subset["End"].astype(str)
    return subset.set_index("Bin")


def subset_site_data(site_data: np.lib.npyio.NpzFile, sample_ids: list[str], bin_idx: np.ndarray, max_sites: int) -> dict[str, np.ndarray]:
    all_sample_ids = site_data["sample_ids"].tolist()
    sample_lookup = {sample_id: idx for idx, sample_id in enumerate(all_sample_ids)}
    sample_idx = np.array([sample_lookup[sample_id] for sample_id in sample_ids], dtype=np.int64)
    return {
        "site_alt": site_data["site_alt"][bin_idx, :max_sites][:, :, sample_idx],
        "site_total": site_data["site_total"][bin_idx, :max_sites][:, :, sample_idx],
        "site_pop_af": site_data["site_pop_af"][bin_idx, :max_sites],
        "site_mask": site_data["site_mask"][bin_idx, :max_sites][:, :, sample_idx],
        "sample_ids": np.array(sample_ids, dtype=object),
        "bin_chr": site_data["bin_chr"][bin_idx],
        "bin_start": site_data["bin_start"][bin_idx],
        "bin_end": site_data["bin_end"][bin_idx],
    }


def subset_raw_depth(raw_df: pd.DataFrame, sample_ids: list[str]) -> pd.DataFrame:
    keep_cols = ["#Chr", "Start", "End", *sample_ids]
    return raw_df.loc[:, keep_cols].copy()


def write_truth_jsons(
    output_dir: Path,
    call_df: pd.DataFrame,
    sample_ids: list[str],
) -> None:
    truth = {sample_id: "NORMAL" for sample_id in sample_ids}
    sex_truth = {
        sample_id: str(
            call_df.loc[call_df["sample"] == sample_id, "sex"].iloc[0]
        )
        for sample_id in sample_ids
    }
    (output_dir / "truth.json").write_text(json.dumps(truth, indent=2, sort_keys=True))
    (output_dir / "sex_truth.json").write_text(
        json.dumps(sex_truth, indent=2, sort_keys=True)
    )


def write_manifest(output_dir: Path, source_run: Path, spec_name: str, sample_ids: list[str], bin_idx: np.ndarray, max_sites: int) -> None:
    manifest = {
        "source_run": str(source_run),
        "spec": spec_name,
        "target_chromosomes": TARGET_CHROMS,
        "sample_ids": sample_ids,
        "n_samples": len(sample_ids),
        "n_bins": int(len(bin_idx)),
        "max_sites_per_bin": max_sites,
    }
    (output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True))


def build_medium_fixture(source_run: Path, output_dir: Path) -> None:
    pre_df = pd.read_csv(source_run / "preprocess" / "preprocessed_depth.tsv", sep="\t", index_col=0)
    site_data = np.load(source_run / "preprocess" / "site_data.npz", allow_pickle=True)
    call_df = pd.read_csv(source_run / "call" / "aneuploidy_type_predictions.tsv", sep="\t")

    sample_ids = choose_sample_ids(call_df, MEDIUM_SPEC["female_count"], MEDIUM_SPEC["male_count"])
    bin_idx = choose_bin_indices(pre_df.reset_index(drop=True), MEDIUM_SPEC["bins_per_chr"])
    subset_df = subset_preprocessed(pre_df.reset_index(drop=True), sample_ids, bin_idx)
    subset_df.to_csv(output_dir / "preprocessed_depth.tsv", sep="\t")

    subset_npz = subset_site_data(site_data, sample_ids, bin_idx, MEDIUM_SPEC["max_sites"])
    np.savez_compressed(output_dir / "site_data.npz", **subset_npz)

    write_truth_jsons(output_dir, call_df, sample_ids)
    write_manifest(output_dir, source_run, "medium", sample_ids, bin_idx, MEDIUM_SPEC["max_sites"])


def build_large_fixture(source_run: Path, output_dir: Path) -> None:
    raw_df = pd.read_csv(source_run / "data" / "all_samples_condensed_depth.rd.txt.gz", sep="\t", compression="gzip")
    call_df = pd.read_csv(source_run / "call" / "aneuploidy_type_predictions.tsv", sep="\t")
    sample_ids = choose_sample_ids(call_df, LARGE_SPEC["female_count"], LARGE_SPEC["male_count"])
    subset_df = subset_raw_depth(raw_df, sample_ids)
    subset_df.to_csv(output_dir / "raw_depth.tsv.gz", sep="\t", index=False, compression="gzip")
    write_truth_jsons(output_dir, call_df, sample_ids)
    write_manifest(output_dir, source_run, "large", sample_ids, np.arange(len(subset_df)), LARGE_SPEC["max_sites"])


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate gatk-sv-ploidy real-data test fixtures")
    parser.add_argument(
        "--source-run",
        default="/Users/markw/Work/talkowski/sv-pipe-testing/mw_ploidy/gatk-sv-ploidy/site_depth",
        help="Path to the example site_depth run",
    )
    parser.add_argument(
        "--output-root",
        default=str(Path(__file__).resolve().parent),
        help="Fixture output root",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    source_run = Path(args.source_run)
    output_root = Path(args.output_root)
    medium_dir = output_root / "medium"
    large_dir = output_root / "large"
    medium_dir.mkdir(parents=True, exist_ok=True)
    large_dir.mkdir(parents=True, exist_ok=True)
    build_medium_fixture(source_run, medium_dir)
    build_large_fixture(source_run, large_dir)


if __name__ == "__main__":
    main()
