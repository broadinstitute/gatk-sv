"""
This script creates a single table of sample QC & batching metrics from EvidenceQC.
"""

import argparse
import pandas as pd
from collections import Counter
from numpy import array
from functools import reduce

def read_ploidy(filename: str) -> pd.DataFrame:
    """

    Args:
        filename: A tab-delimited file containing estimated copy numbers.

    Returns:
        A pandas DataFrame containing the following columns:
        [id, chr1_CopyNumber, ..., chrX_CopyNumber, chrX_CopyNumber_rounded].
    """
    df_ploidy = pd.read_csv(filename, sep="\t")
    df_ploidy.loc[df_ploidy["chrX_CopyNumber"] <= 2, "chrX_CopyNumber_rounded"] = 1
    df_ploidy.loc[df_ploidy["chrX_CopyNumber"] > 2, "chrX_CopyNumber_rounded"] = 2
    return df_ploidy


def read_bincov_median(filename: str) -> pd.DataFrame:
    """
    Median coverage (bincov_median)
    Args:
        filename: bincov_median output from EvidenceQC (medianCov.transposed) in a tab separated format

    Returns:
        a pandas DataFrame containing the transposed table given in the input file.
    """
    df_median = pd.read_csv(filename, sep="\t").T
    df_median = df_median.rename_axis("#ID").reset_index()
    df_median.columns = ["#ID", "median_coverage"]
    df_median = df_median.reset_index(drop=True)
    return df_median


def read_wgd_scores(filename: str) -> pd.DataFrame:
    return pd.read_csv(filename, sep="\t")


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
    df_non_diploid = df_non_diploid.T.iloc[0:, 3:]  # .rename_axis("#ID").reset_index()
    new_c = df_non_diploid.iloc[0:, 0:]
    nondiploid_counts = new_c[new_c < 0.05].count()
    nondiploid_counts_df = pd.DataFrame(nondiploid_counts, columns=["nondiploid_counts"]).rename_axis(
        "#ID").reset_index()
    return nondiploid_counts_df


def read_melt_insert_size(filename: str) -> array:
    """
    Reads MELT insert size into an array.
    Args:
        filename: melt_insert_size in sample data table

    Returns:
        Array[Float] of melt_insert_size for each sample
    """
    insert_size = pd.read_csv(filename, sep="\t")
    melt_insert_size_array = pd.DataFrame(insert_size, columns=["melt_insert_size"]).astype(float)
    return melt_insert_size_array


def read_manta_outlier(filename: str) -> pd.DataFrame:
    """
    Outliers from Manta.
    Args:
        filename: A tab-delimited file containing Manta QC Low/High.

    Returns:
        A pandas DataFrame containing the number of times a sample appears
        in Manta QC low and/or high.
    """
    outlier_manta = filename.split(".")[-4] + "_" + filename.split(".")[-1]  # + "_" + filename.split(".")[-3]
    df_manta = pd.read_csv(filename, sep="\t")
    if df_manta["Outlier_Sample"].empty:
        df_manta["#ID"] = df_manta.apply(lambda _: " ", axis=0)
        outlier_manta_sample = df_manta.pivot_table(columns=["#ID"], aggfunc="size").astype(int)
    else:
        df_manta["#ID"] = df_manta["Outlier_Sample"].str.split("/", expand=True)[10].str.split(".", expand=True)[0]
        outlier_manta_sample = df_manta.pivot_table(columns=["#ID"], aggfunc="size").astype(int)
    outlier_manta_df = outlier_manta_sample.reset_index()
    outlier_manta_df.columns = ["#ID", str(outlier_manta) + "_outlier"]
    return outlier_manta_df


def read_melt_outlier(filename: str) -> pd.DataFrame:
    """
    Outliers from Melt.
    Args:
        filename: A tab-delimited file containing the output of EvidenceQC on MELT QC Low/High.

    Returns:
        A pandas DataFrame containing the number of times a sample appears
        in MELT QC low and/or high.
    """
    outlier_melt = filename.split(".")[-4] + "_" + filename.split(".")[-1]
    df_melt = pd.read_csv(filename, sep="\t")
    if df_melt["Outlier_Sample"].empty:
        df_melt["#ID"] = df_melt.apply(lambda _: " ", axis=0)
        outlier_melt_sample = df_melt.pivot_table(columns=["#ID"], aggfunc="size").astype(int)
    else:
        df_melt["#ID"] = df_melt["Outlier_Sample"].str.split("/", expand=True)[10].str.split(".", expand=True)[0]
        outlier_melt_sample = df_melt.pivot_table(columns=["#ID"], aggfunc="size").astype(int)
    outlier_melt_df = outlier_melt_sample.reset_index()
    outlier_melt_df.columns = ["#ID", str(outlier_melt) + "_outlier"]
    return outlier_melt_df


def read_wham_outlier(filename: str) -> pd.DataFrame:
    """
    Outliers from Wham.
    Args:
        filename: A tab-delimited file containing the output of EvidenceQC on Wham QC Low/High.

    Returns:
        A pandas DataFrame containing the number of times a sample appears
        in MELT QC low and/or high.
    """
    outlier_wham = filename.split(".")[-4] + "_" + filename.split(".")[-1]  # + "_" + filename.split(".")[-3]
    df_wham = pd.read_csv(filename, sep="\t")
    if df_wham["Outlier_Sample"].empty:
        df_wham["#ID"] = df_wham.apply(lambda _: " ", axis=0)
        outlier_wham_sample = df_wham.pivot_table(columns=["#ID"], aggfunc="size").astype(int)
    else:
        df_wham["#ID"] = df_wham["Outlier_Sample"].str.split("/", expand=True)[11].str.split(".", expand=True)[0]
        outlier_wham_sample = df_wham.pivot_table(columns=["#ID"], aggfunc="size").astype(int)
    outlier_wham_df = outlier_wham_sample.reset_index()
    outlier_wham_df.columns = ["#ID", str(outlier_wham) + "_outlier"]
    return outlier_wham_df


def read_all_outlier(filename_manta: str, filename_melt: str, filename_wham: str, outlier_type: str) -> pd.DataFrame:
    """

    Args:
        filename_manta: Outliers determined in EvidenceQC for Manta.
        filename_melt: Outliers determined in EvidenceQC for MELT.
        filename_wham: Outliers determined in EvidenceQC for Wham.
        outlier_type:

    Returns:
        The total number of times that a sample appears as an outlier
        across all four algorithms for each sample.
    """
    # Manta:
    outlier_manta = filename_manta.split(".")[-4] + "_" + filename_manta.split(".")[
        -1]
    outlier_manta_df = read_manta_outlier(filename_manta)
    dict_manta = dict(list(zip(outlier_manta_df["#ID"], outlier_manta_df[str(outlier_manta) + "_outlier"])))
    # Melt:
    outlier_melt = filename_melt.split(".")[-4] + "_" + filename_melt.split(".")[
        -1]
    outlier_melt_df = read_melt_outlier(filename_melt)
    dict_melt = dict(list(zip(outlier_melt_df["#ID"], outlier_melt_df[str(outlier_melt) + "_outlier"])))
    # Wham:
    outlier_wham = filename_wham.split(".")[-4] + "_" + filename_wham.split(".")[
        -1]
    outlier_wham_df = read_wham_outlier(filename_wham)
    dict_wham = dict(list(zip(outlier_wham_df["#ID"], outlier_wham_df[str(outlier_wham) + "_outlier"])))
    # merging all the dictionaries
    outlier_dicts = [dict_manta, dict_melt, dict_wham]
    merged_dicts = Counter()
    for counted in outlier_dicts:
        merged_dicts.update(counted)
    all_outliers = dict(merged_dicts)
    all_outliers_df = pd.DataFrame.from_dict(all_outliers, orient="index").reset_index()  # index_col=None
    all_outliers_df.columns = ["#ID", outlier_type + "_overall_outliers"]
    return all_outliers_df


###########################################################################
# merging all the dataframes:
def merge_evidence_qc_table(
        filename_estimated_cn: str,
        filename_medianCov: str,
        filename_wgd: str,
        filename_cnv_qvalues: str,
        filename_high_manta: str,
        filename_high_melt: str,
        filename_high_wham: str,
        filename_low_manta: str,
        filename_low_melt: str,
        filename_low_wham: str,
        output_prefix: str) -> str:
    """
    Reads the provided TSV files (tab-delimited) and merges all the information in one table
    serialized to the given output filename.
    """
    df_ploidy = read_ploidy(filename_estimated_cn)
    df_bincov_median = read_bincov_median(filename_medianCov)
    df_wgd_scores = read_wgd_scores(filename_wgd)
    df_non_diploid = read_non_diploid(filename_cnv_qvalues)
    df_manta_high_outlier = read_manta_outlier(filename_high_manta)
    df_melt_high_outlier = read_melt_outlier(filename_high_melt)
    df_wham_high_outlier = read_wham_outlier(filename_high_wham)
    df_total_high_outliers = read_all_outlier(filename_high_manta, filename_high_melt, filename_high_wham, "high")
    df_manta_low_outlier = read_manta_outlier(filename_low_manta)
    df_melt_low_outlier = read_melt_outlier(filename_low_melt)
    df_wham_low_outlier = read_wham_outlier(filename_low_wham)
    df_total_low_outliers = read_all_outlier(filename_low_manta, filename_low_melt, filename_low_wham, "low")

    # all data frames
    dfs = [df_ploidy, df_bincov_median, df_wgd_scores, df_non_diploid, df_manta_high_outlier,
           df_melt_high_outlier, df_wham_high_outlier, df_total_high_outliers,
           df_manta_low_outlier, df_melt_low_outlier, df_wham_low_outlier, df_total_low_outliers]
    output_df = reduce(lambda left, right: pd.merge(left, right, on="#ID",
                                                    how="outer"), dfs).fillna(0)  # "none"

    # save the file
    return output_df.to_csv(f"{output_prefix}evidence_qc_table.tsv", sep="\t", header=True, index=False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("ploidy_test_filename", help="...")
    parser.add_argument("median_cov_transposed_filename", help="...")
    parser.add_argument("wgd_scores_filename", help="...")
    parser.add_argument("binwise_cnv_qValues_filename", help="...")
    parser.add_argument("manta_qc_outlier_high_filename", help="...")
    parser.add_argument("melt_qc_outlier_high_filename", help="...")
    parser.add_argument("wham_qc_outlier_high_filename", help="...")
    parser.add_argument("manta_qc_outlier_low_filename", help="...")
    parser.add_argument("melt_qc_outlier_low_filename", help="...")
    parser.add_argument("wham_qc_outlier_low_filename", help="...")
    parser.add_argument("melt_insert_size_file", help="...")
    parser.add_argument("--output-prefix", default="", help="...")

    args = parser.parse_args()

    merge_evidence_qc_table(
        args.ploidy_test_filename,
        args.median_cov_transposed_filename,
        args.wgd_scores_filename,
        args.binwise_cnv_qValues_filename,
        args.manta_qc_outlier_high_filename,
        args.melt_qc_outlier_high_filename,
        args.wham_qc_outlier_high_filename,
        args.manta_qc_outlier_low_filename,
        args.melt_qc_outlier_low_filename,
        args.wham_qc_outlier_low_filename,
        args.output_prefix)


if __name__ == "__main__":
    main()
