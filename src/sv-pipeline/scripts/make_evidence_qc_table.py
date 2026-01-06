"""
This script creates a single table of sample QC & batching metrics from EvidenceQC.
"""

import argparse
import pandas as pd
from collections import Counter
from functools import reduce
from pathlib import Path
import numpy as np

ID_COL = "sample_id"
EMPTY_OUTLIERS = "EMPTY_ROWS_DROP"


def read_ploidy(filename: str) -> pd.DataFrame:
    """
    Args:
        filename: A tab-delimited file containing estimated copy numbers.
    Returns:
        A pandas DataFrame containing the following columns:
        [id, chr1_CopyNumber, ..., chr22_CopyNumber, chrX_CopyNumber, chrY_CopyNumber, chrX_CopyNumber_rounded].
    """
    df_ploidy = pd.read_csv(filename, sep="\t")
    df_ploidy.loc[round(df_ploidy["chrX_CopyNumber"]) < 2, "chrX_CopyNumber_rounded"] = 1
    df_ploidy.loc[round(df_ploidy["chrX_CopyNumber"]) >= 2, "chrX_CopyNumber_rounded"] = 2
    return df_ploidy


def read_sex_assignments(filename: str) -> pd.DataFrame:
    """
    Args:
        filename: A tab-delimited file containing estimated sex assignments.
    Returns:
        A pandas DataFrame containing the estimated sex assignment for each sample.
    """
    df_assignments = pd.read_csv(filename, sep="\t")
    df_assignments = df_assignments[[ID_COL, "Assignment"]]
    df_assignments.rename(columns={'Assignment': 'sex_assignment'}, inplace=True)
    return df_assignments


def read_bincov_median(filename: str) -> pd.DataFrame:
    """
    Median coverage (bincov_median)
    Args:
        filename: bincov_median output from EvidenceQC (medianCov.transposed) in a tab separated format
    Returns:
        a pandas DataFrame containing the transposed table given in the input file.
    """
    df_median = pd.read_csv(filename, sep="\t").T
    df_median = df_median.rename_axis(ID_COL).reset_index()
    df_median.columns = [ID_COL, "median_coverage"]
    df_median = df_median.reset_index(drop=True)
    return df_median


def read_wgd_scores(filename: str) -> pd.DataFrame:
    """
    Args:
        filename: a tab-delimited file containing wgd scores.
    Returns:
        A pandas DataFrame containing wgd scores for each sample.
    """
    df_wgd_scores = pd.read_csv(filename, sep="\t")
    df_wgd_scores.rename(columns={'score': 'wgd_score'}, inplace=True)
    return df_wgd_scores


def read_non_diploid(filename: str) -> pd.DataFrame:
    """
    EvidenceQC → untar → /ploidy_est/binwise_CNV_qValues.bed.gz →
    count number of bins with q-value < 0.05 for each sample
    Args:
        filename: a tab-delimited file containing binwise CNV qValues.
    Returns:
        A pandas DataFrame containing the number of bins with q-value < 0.05 for each sample.
    """
    df_non_diploid = pd.read_csv(filename, sep="\t").T
    df_non_diploid = df_non_diploid.T.iloc[0:, 3:]
    nondiploid_counts = df_non_diploid[df_non_diploid < 0.05].count()
    nondiploid_counts_df = pd.DataFrame(nondiploid_counts, columns=["nondiploid_bins"]).rename_axis(
        ID_COL).reset_index()
    return nondiploid_counts_df


def read_melt_insert_size(filename: str, col_name="mean_insert_size") -> pd.DataFrame:
    """
    Reads MELT insert size into a Pandas DataFrame.
    Args:
        filename: melt_insert_size in sample data table
        col_name: Set title to use for the column containing MELT mean insert size.
    Returns:
        MELT insert size for each sample given in a Pandas DataFrame with two columns: sample ID and mean insert size.
    """
    columns_names = [ID_COL, col_name]
    if filename is None:
        return pd.DataFrame(columns=columns_names)
    df = pd.read_csv(filename, names=columns_names, header=0, sep="\t")
    df[col_name] = df[col_name].astype(float)
    return df


def get_col_name(caller: str, outlier_type: str) -> str:
    return f"{caller}_{outlier_type}_outlier"


def read_outlier(filename: str, outlier_col_label: str) -> pd.DataFrame:
    """
    Args:
        filename: A tab-delimited file containing the output of
        Evidence QC on a caller (e.g., Manta) with low/high.
        outlier_col_label: Column header for the outliers.
    Returns:
        A pandas DataFrame containing the number of times a sample
        appears in the QC of a caller output.
    """
    if filename is None:
        return pd.DataFrame(columns=[ID_COL, outlier_col_label])
    df = pd.read_csv(filename, sep="\t")
    if df["Outlier_Sample"].empty:
        df[ID_COL] = df.apply(lambda _: EMPTY_OUTLIERS, axis=0)
        outlier_sample = df.pivot_table(columns=[ID_COL], aggfunc="size").astype(int)
    else:
        df['Outlier_Sample'] = df['Outlier_Sample'].apply(Path)
        df[ID_COL] = df['Outlier_Sample'].apply(lambda x: x.stem if '.' in x.suffix else np.nan)
        outlier_sample = df.pivot_table(columns=[ID_COL], aggfunc="size").astype(int)
    outlier_df = outlier_sample.reset_index()
    outlier_df.columns = [ID_COL, outlier_col_label]
    return outlier_df


def read_all_outlier(outlier_manta_df: pd.DataFrame, outlier_melt_df: pd.DataFrame, outlier_wham_df: pd.DataFrame,
                     outlier_scramble_df: pd.DataFrame, outlier_dragen_df: pd.DataFrame, outlier_type: str) -> pd.DataFrame:
    """
    Args:
        outlier_manta_df: Outliers determined in EvidenceQC for Manta.
        outlier_melt_df: Outliers determined in EvidenceQC for MELT.
        outlier_wham_df: Outliers determined in EvidenceQC for Wham.
        outlier_scramble_df: Outliers determined in EvidenceQC for Scramble.
        outlier_dragen_df: Outliers determined in EvidenceQC for Dragen.
        outlier_type: high or low. Determined in EvidenceQC for each of the three callers.
    Returns:
        The total number of times that a sample appears as an outlier
        across all three algorithms for each sample.
    """
    # Manta:
    col_name = get_col_name("manta", outlier_type)
    dict_manta = dict(zip(outlier_manta_df[ID_COL], outlier_manta_df[col_name]))

    # Melt:
    col_name = get_col_name("melt", outlier_type)
    dict_melt = dict(zip(outlier_melt_df[ID_COL], outlier_melt_df[col_name]))

    # Wham:
    col_name = get_col_name("wham", outlier_type)
    dict_wham = dict(zip(outlier_wham_df[ID_COL], outlier_wham_df[col_name]))

    # Scramble:
    col_name = get_col_name("scramble", outlier_type)
    dict_scramble = dict(zip(outlier_scramble_df[ID_COL], outlier_scramble_df[col_name]))

    # Dragen:
    col_name = get_col_name("dragen", outlier_type)
    dict_dragen = dict(zip(outlier_dragen_df[ID_COL], outlier_dragen_df[col_name]))

    # merging all the dictionaries
    outlier_dicts = [dict_manta, dict_melt, dict_wham, dict_scramble, dict_dragen]
    merged_dicts = Counter()
    for counted in outlier_dicts:
        merged_dicts.update(counted)
    all_outliers = dict(merged_dicts)
    if len(all_outliers) == 0:
        all_outliers_df = pd.DataFrame(columns=[ID_COL, "overall_" + outlier_type + "_outlier"])
    else:
        all_outliers_df = pd.DataFrame.from_dict(all_outliers, orient="index").reset_index()
        all_outliers_df.columns = [ID_COL, "overall_" + outlier_type + "_outlier"]
    return all_outliers_df


def read_variant_counts(filename: str) -> pd.DataFrame:
    """
    Reads variant counts from a TSV file.
    Args:
        filename: TSV file containing variant counts by sample, SV type, and chromosome.
    Returns:
        A pandas DataFrame containing variant counts by sample.
    """
    if filename is None:
        return pd.DataFrame(columns=[ID_COL])

    # Read the file and ensure sample_id is used as a regular column
    df = pd.read_csv(filename, sep="\t")

    # If there are no rows or no matching column, return empty DataFrame
    if df.empty or ID_COL not in df.columns:
        return pd.DataFrame(columns=[ID_COL])

    # Ensure all columns have integer values (replace NaN with 0)
    value_columns = [col for col in df.columns if col != ID_COL]
    if value_columns:
        df[value_columns] = df[value_columns].fillna(0).astype(int)

    return df


def merge_evidence_qc_table(
        filename_estimated_cn: str,
        filename_sex_assignments: str,
        filename_mediancov: str,
        filename_wgd: str,
        filename_cnv_qvalues: str,
        filename_high_manta: str,
        filename_high_melt: str,
        filename_high_wham: str,
        filename_high_scramble: str,
        filename_high_dragen: str,
        filename_low_manta: str,
        filename_low_melt: str,
        filename_low_wham: str,
        filename_low_scramble: str,
        filename_low_dragen: str,
        filename_manta_variant_counts: str,
        filename_melt_variant_counts: str,
        filename_wham_variant_counts: str,
        filename_scramble_variant_counts: str,
        filename_dragen_variant_counts: str,
        filename_melt_insert_size: str,
        output_prefix: str) -> None:
    """
    Reads the provided TSV files (tab-delimited) and merges all the information in one table
    serialized to the given output filename.
    """
    df_ploidy = read_ploidy(filename_estimated_cn)
    df_sex_assignments = read_sex_assignments(filename_sex_assignments)
    df_bincov_median = read_bincov_median(filename_mediancov)
    df_wgd_scores = read_wgd_scores(filename_wgd)
    df_non_diploid = read_non_diploid(filename_cnv_qvalues)
    df_manta_high_outlier = read_outlier(filename_high_manta, get_col_name("manta", "high"))
    df_melt_high_outlier = read_outlier(filename_high_melt, get_col_name("melt", "high"))
    df_wham_high_outlier = read_outlier(filename_high_wham, get_col_name("wham", "high"))
    df_scramble_high_outlier = read_outlier(filename_high_scramble, get_col_name("scramble", "high"))
    df_dragen_high_outlier = read_outlier(filename_high_dragen, get_col_name("dragen", "high"))
    df_total_high_outliers = read_all_outlier(df_manta_high_outlier, df_melt_high_outlier, df_wham_high_outlier,
                                              df_scramble_high_outlier, df_dragen_high_outlier, "high")
    df_manta_low_outlier = read_outlier(filename_low_manta, get_col_name("manta", "low"))
    df_melt_low_outlier = read_outlier(filename_low_melt, get_col_name("melt", "low"))
    df_wham_low_outlier = read_outlier(filename_low_wham, get_col_name("wham", "low"))
    df_scramble_low_outlier = read_outlier(filename_low_scramble, get_col_name("scramble", "low"))
    df_dragen_low_outlier = read_outlier(filename_low_dragen, get_col_name("dragen", "low"))
    df_total_low_outliers = read_all_outlier(df_manta_low_outlier, df_melt_low_outlier, df_wham_low_outlier,
                                             df_scramble_low_outlier, df_dragen_low_outlier, "low")
    df_manta_variant_counts = read_variant_counts(filename_manta_variant_counts)
    df_melt_variant_counts = read_variant_counts(filename_melt_variant_counts)
    df_wham_variant_counts = read_variant_counts(filename_wham_variant_counts)
    df_scramble_variant_counts = read_variant_counts(filename_scramble_variant_counts)
    df_dragen_variant_counts = read_variant_counts(filename_dragen_variant_counts)
    df_melt_insert_size = read_melt_insert_size(filename_melt_insert_size)

    # outlier column names
    callers = ["wham", "melt", "manta", "scramble", "dragen", "overall"]
    types = ["high", "low"]
    outlier_cols = [get_col_name(caller, type) for caller in callers for type in types]

    # all data frames
    dfs = [df_ploidy, df_sex_assignments, df_bincov_median, df_wgd_scores, df_non_diploid, df_melt_insert_size,
           df_manta_high_outlier, df_melt_high_outlier, df_wham_high_outlier, df_scramble_high_outlier, df_dragen_high_outlier, df_total_high_outliers,
           df_manta_low_outlier, df_melt_low_outlier, df_wham_low_outlier, df_scramble_low_outlier, df_dragen_low_outlier, df_total_low_outliers,
           df_manta_variant_counts, df_melt_variant_counts, df_wham_variant_counts, df_scramble_variant_counts, df_dragen_variant_counts]
    for df in dfs:
        df[ID_COL] = df[ID_COL].astype(object)
    output_df = reduce(lambda left, right: pd.merge(left, right, on=ID_COL, how="outer"), dfs)
    output_df = output_df[output_df[ID_COL] != EMPTY_OUTLIERS]
    output_df[outlier_cols] = output_df[outlier_cols].replace([None, np.nan], 0.0)

    # save the file
    output_df.to_csv(f"{output_prefix}.evidence_qc_table.tsv", sep="\t", header=True, index=False, na_rep=np.nan)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-y", "--estimated-copy-number-filename",
        help="Sets the filename containing estimated copy numbers per contig.")

    parser.add_argument(
        "-x", "--sex-assignments-filename",
        help="Sets the filename containing copy number-based sex assignments.")

    parser.add_argument(
        "-d", "--median-cov-filename",
        help="Sets the filename containing median coverage.")

    parser.add_argument(
        "-g", "--wgd-scores-filename",
        help="Sets the filename containing WGD scores.")

    parser.add_argument(
        "-b", "--binwise-cnv-qvalues-filename",
        help="Sets the filename containing bin-wise CNV q-values.")

    parser.add_argument(
        "-z", "--manta-qc-outlier-high-filename",
        help="Sets the filename containing Manta QC outlier high.")

    parser.add_argument(
        "-e", "--melt-qc-outlier-high-filename",
        help="Sets the filename containing Melt QC outlier high.")

    parser.add_argument(
        "-w", "--wham-qc-outlier-high-filename",
        help="Sets the filename containing Wham QC outlier high.")

    parser.add_argument(
        "-t", "--scramble-qc-outlier-high-filename",
        help="Sets the filename containing Scramble QC outlier high.")

    parser.add_argument(
        "-i", "--dragen-qc-outlier-high-filename",
        help="Sets the filename containing Dragen QC outlier high.")

    parser.add_argument(
        "-a", "--manta-qc-outlier-low-filename",
        help="Sets the filename containing Manta QC outlier low.")

    parser.add_argument(
        "-s", "--melt-qc-outlier-low-filename",
        help="Sets the filename containing Melt QC outlier low.")

    parser.add_argument(
        "-r", "--wham-qc-outlier-low-filename",
        help="Sets the filename containing Wham QC outlier low.")

    parser.add_argument(
        "-c", "--scramble-qc-outlier-low-filename",
        help="Sets the filename containing Scramble QC outlier low.")

    parser.add_argument(
        "-j", "--dragen-qc-outlier-low-filename",
        help="Sets the filename containing Dragen QC outlier low.")

    parser.add_argument(
        "-v", "--manta-variant-counts-filename",
        help="Sets the filename containing Manta variant counts per sample.")

    parser.add_argument(
        "-k", "--melt-variant-counts-filename",
        help="Sets the filename containing Melt variant counts per sample.")

    parser.add_argument(
        "-l", "--wham-variant-counts-filename",
        help="Sets the filename containing Wham variant counts per sample.")

    parser.add_argument(
        "-p", "--scramble-variant-counts-filename",
        help="Sets the filename containing Scramble variant counts per sample.")

    parser.add_argument(
        "-q", "--dragen-variant-counts-filename",
        help="Sets the filename containing DRAGEN-SV variant counts per sample.")

    parser.add_argument(
        "-m", "--melt-insert-size", dest="melt_insert_size_filename",
        help="Sets the filename containing Melt insert size. "
             "This file is expected to have two columns, containing "
             "sample ID and mean insert size in the first and second "
             "columns respectively. Also, the first line of the file "
             "is considered as header and is ignored.")

    parser.add_argument(
        "-o", "--output-prefix", required=True,
        help="Sets a prefix to be added to the output file(s).")

    args = parser.parse_args()

    merge_evidence_qc_table(
        args.estimated_copy_number_filename,
        args.sex_assignments_filename,
        args.median_cov_filename,
        args.wgd_scores_filename,
        args.binwise_cnv_qvalues_filename,
        args.manta_qc_outlier_high_filename,
        args.melt_qc_outlier_high_filename,
        args.wham_qc_outlier_high_filename,
        args.scramble_qc_outlier_high_filename,
        args.dragen_qc_outlier_high_filename,
        args.manta_qc_outlier_low_filename,
        args.melt_qc_outlier_low_filename,
        args.wham_qc_outlier_low_filename,
        args.scramble_qc_outlier_low_filename,
        args.dragen_qc_outlier_low_filename,
        args.manta_variant_counts_filename,
        args.melt_variant_counts_filename,
        args.wham_variant_counts_filename,
        args.scramble_variant_counts_filename,
        args.dragen_variant_counts_filename,
        args.melt_insert_size_filename,
        args.output_prefix)


if __name__ == "__main__":
    main()
