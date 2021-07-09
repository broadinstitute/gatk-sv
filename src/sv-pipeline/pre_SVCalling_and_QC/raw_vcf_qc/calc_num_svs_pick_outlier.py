#!/usr/bin/env python

import sys
from typing import Sequence, Set
import argparse
import numpy
import pandas


_zero_svs_are_outliers = True
_outlier_std_threshold = 5.0
_column_order = ["CHROM", "SVTYPE", "Mean", "Median", "STD",
                 "Outlier_Sample", "Outlier_Number", "Outlier_Cate"]


def read_statfile(statfile: str) -> pandas.DataFrame:
    """
    Special function needed to read in stats data table because
    a) pandas doesn't understand that the '#' means header
    b) there are multiple stats files concatenated together, resulting in headers being randomly mixed in
    Args:
        statfile: str
            File name with concatenated tab-separated tables of variant stats
    Returns:
        stats_data: pandas.DataFrame
            Table of variant stats
    """
    with open(statfile, 'r') as f_in:
        # get column header from first line, stripping '#'
        columns = f_in.readline().lstrip('#').split()
        # read rest of tsv file, using these columns as header and ignoring any future lines starting with '#'
        return pandas.read_csv(statfile, sep='\t', comment='#', names=columns)


def pick_outliers_by_group(
        chrom: str,
        sv_type: str,
        check_stats: pandas.DataFrame,
        all_samples: Set[str],
        zero_svs_are_outliers: bool = _zero_svs_are_outliers,
        outlier_std_threshold: float = _outlier_std_threshold
) -> pandas.DataFrame:
    """
    For given combination of contig and SV type, find samples that have outlier number of SVs. Return table of outliers
    along with statistics about SV count.
    Args:
        chrom: str
            Contig for checking SV counts
        sv_type: str
            SV type for checking SV counts
        check_stats: pandas.DataFrame
            Table with SV counts on this contig with this sv_type
        all_samples: Set[str]
            Set of all sample IDs in cohort
        zero_svs_are_outliers: bool
            Whether to treat samples with no counts as automatic outliers, or explicitly code as zero counts
        outlier_std_threshold: float
            Threshold for outlier status as multiple of standard deviation of SV counts
    Returns:
        outliers: pandas.DataFrame
            Table of outliers
    """
    # find samples that are missing: they have 0 SVs of this type on this contig
    missing_samples = pandas.DataFrame(
        tuple(
            {"CHROM": chrom, "SVTYPE": sv_type, "SAMPLE": sample_id, "NUM": 0}
            for sample_id in all_samples.difference(check_stats["SAMPLE"])
        )
    )

    if zero_svs_are_outliers:
        # THIS IS THE ORIGINAL PIPELINE BEHAVIOR
        # compute basic stats about observed nonzero SV counts
        count_mean = check_stats["NUM"].mean()
        count_median = check_stats["NUM"].median()
        count_std = check_stats["NUM"].std()
        # Amongst samples that have SVs, find counts deviating by more than set multiple of std from the median
        is_outlier = numpy.abs(
            check_stats["NUM"] - count_median) > outlier_std_threshold * count_std
        # Treat missing samples as outliers.
        outliers = pandas.concat(
            (missing_samples, check_stats.loc[is_outlier]), axis=0)
    else:
        # THIS FINDS FEWER, MORE MEANINGFUL OUTLIERS
        # Which samples are missing / included but have zero counts is unpredictable.
        # 1) concatenate all samples together
        check_stats = pandas.concat((check_stats, missing_samples), axis=0)
        # 2) compute stats from non-zero SV counts
        nonzero = check_stats["NUM"] > 0
        count_mean = check_stats.loc[nonzero, "NUM"].mean()
        count_median = check_stats.loc[nonzero, "NUM"].median()
        count_std = check_stats.loc[nonzero, "NUM"].std()
        # 3) check outliers by usual means from those stats
        # Set threshold to be set multiple of greater of: std of counts, sqrt(median of counts)
        #  (i.e. greater of std or expected Poisson std)
        # Find counts those deviating by more than threshold from the median (including zeros)
        is_outlier = (
            numpy.abs(check_stats["NUM"] - count_median) >
            outlier_std_threshold * numpy.maximum(count_std, numpy.sqrt(count_median))
        )
        outliers = check_stats.loc[is_outlier].copy()

    if outliers.empty:
        return pandas.DataFrame([], columns=_column_order)
    # augment outlier table with some statistics
    outliers["Mean"] = count_mean
    outliers["Median"] = count_median
    outliers["STD"] = count_std
    outliers["Outlier_Cate"] = numpy.where(
        outliers["NUM"] > count_median, "high", "low")
    # rename and re-order columns
    return outliers.rename({"NUM": "Outlier_Number", "SAMPLE": "Outlier_Sample"}, axis=1).reindex(_column_order, axis=1)


def pick_outliers(
        stats_data: pandas.DataFrame,
        zero_svs_are_outliers: bool = _zero_svs_are_outliers,
        outlier_std_threshold: float = _outlier_std_threshold
) -> pandas.DataFrame:
    """
    Find samples that have outlier number of SVs when broken down by contig and SV type. Return table of outliers
    along with statistics about SV count.
    Args:
        stats_data: pandas.DataFrame
            Table with SV counts
        zero_svs_are_outliers: bool
            Whether to treat samples with no counts as automatic outliers, or explicitly code as zero counts
        outlier_std_threshold: float
            Threshold for outlier status as multiple of standard deviation of SV counts
    Returns:
        outliers: pandas.DataFrame
            Table of outliers
    """
    # get set of all samples in stats data
    all_samples = set(stats_data["SAMPLE"])

    # loop over unique combinations of contig and sv type
    #    find outliers from each unique combination
    #    and concatenate those outliers into one table
    outliers = pandas.concat(
        tuple(
            pick_outliers_by_group(
                chrom=chrom, sv_type=sv_type, check_stats=check_stats, all_samples=all_samples,
                zero_svs_are_outliers=zero_svs_are_outliers, outlier_std_threshold=outlier_std_threshold
            )
            for (chrom, sv_type), check_stats in stats_data.groupby(
                ["CHROM", "SVTYPE"], sort=False, as_index=False, group_keys=False
            )
        ),
        axis=0
    )
    return outliers


def write_outliers_file(
        outliers: pandas.DataFrame,
        outname: str,
        outlier_type: str
):
    """
    Write outliers of the appropriate type ("low" or "high") to TSV file.
    Args:
        outliers: pandas.DataFrame
            Table of outlier data
        outname: str
            Base name of outlier TSV file. Final file name will have ".low" or ".high" appended to it.
        outlier_type: str
            "low" or "high".
    """
    # write outliers to tsv. Add "#" in front of header
    with open(outname + "." + outlier_type, 'w') as f_out:
        f_out.write("#")  # add '#' in front of header
        outlier_wanted = outliers["Outlier_Cate"] == outlier_type
        outliers.loc[outlier_wanted].to_csv(f_out, sep='\t', index=False)


def calc_num_svs_pick_outlier(
        statfile: str,
        outname: str,
        zero_svs_are_outliers: bool = _zero_svs_are_outliers,
        outlier_std_threshold: float = _outlier_std_threshold
):
    """
    Find samples that have outlier number of SVs when broken down by contig and SV type.
    Write two tables of outliers, along with statistics about SV count: one for those with above-median counts ("high")
    and one for those at median or below ("low").
    Args:
        statfile: str
            TSV file with table with SV counts
        outname: str
            Base name for saving outlier files. Low file will have ".low" appended to the name, and high file will have
            ".high"
        zero_svs_are_outliers: bool
            Whether to treat samples with no counts as automatic outliers, or explicitly code as zero counts
        outlier_std_threshold: float
            Threshold for outlier status as multiple of standard deviation of SV counts
    """
    stats_data = read_statfile(statfile)
    outliers = pick_outliers(stats_data, zero_svs_are_outliers=zero_svs_are_outliers,
                             outlier_std_threshold=outlier_std_threshold)
    write_outliers_file(outliers, outname, "low")
    write_outliers_file(outliers, outname, "high")


def _parse_arguments(argv: Sequence[str]) -> argparse.Namespace:
    # noinspection PyTypeChecker
    parser = argparse.ArgumentParser(
        description="Find outliers in SV counts broken down by contig and SV type",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("statfile", type=str,
                        help="name of stats concatinated from all samples")
    parser.add_argument("outname", type=str, help="name of output file")
    parser.add_argument("-z", "--zero-counts-are-not-outliers", action="store_true",
                        help="don't make zero SV counts an automatic outlier, check deviation from median as usual")
    parser.add_argument("-t", "--outlier-std-threshold", type=float, default=_outlier_std_threshold,
                        help="threshold multiple of std of counts for outliers")
    return parser.parse_args(argv[1:])


if __name__ == "__main__":
    args = _parse_arguments(sys.argv)
    calc_num_svs_pick_outlier(statfile=args.statfile, outname=args.outname,
                              zero_svs_are_outliers=not args.zero_counts_are_not_outliers,
                              outlier_std_threshold=args.outlier_std_threshold)
