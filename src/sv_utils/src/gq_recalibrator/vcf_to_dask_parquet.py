import sys
from pathlib import Path
import tempfile
import argparse

import numpy
import pandas
import pandas.core.dtypes.dtypes
import pandas.core.dtypes.cast
import dask
import logging
from dask.distributed import Client, LocalCluster
import dask.dataframe
from sv_utils import common, genomics_io
from gq_recalibrator import dask_utils
from gq_recalibrator import tarred_properties_to_parquet
from gq_recalibrator.tarred_properties_to_parquet import CompressionAlgorithms
from typing import Optional


class Default:
    num_threads = common.num_physical_cpus
    temp_dir = Path(tempfile.gettempdir())
    compression_algorithm = tarred_properties_to_parquet.Default.compression_algorithm
    column_levels = genomics_io.Default.column_levels
    drop_rows_with_missing_values = False


class Keys:
    contig = genomics_io.Keys.contig
    gt = genomics_io.Keys.gt
    allele_count = genomics_io.Keys.allele_count
    allele_frequency = genomics_io.Keys.allele_frequency
    is_autosome = tarred_properties_to_parquet.Keys.is_autosome
    id = genomics_io.Keys.id


def vcf_to_dask_parquet(
        input_vcf: Path,
        output_parquet: Path,
        wanted_properties: list[str],
        compression_algorithm: str = Default.compression_algorithm,
        drop_rows_with_missing_values: bool = Default.drop_rows_with_missing_values,
        reconcile_dtypes: bool = True,
        num_threads: int = Default.num_threads
) -> None:
    """
    Converts a VCF file to a Dask Parquet file.

    Args:
        input_vcf: Path to the input VCF file.
        output_parquet: Path to the output Parquet file.
        wanted_properties: List of properties to extract from the VCF.
        compression_algorithm: Compression algorithm to use for parquet.
        drop_rows_with_missing_values: If True, drop rows/variants with missing values
                                       If False, they will be included as missing values in the df.
        reconcile_dtypes: If True, reconcile dtypes of integer columns across dask partitions.
                          I think this is only necessary because of a bug in either dask or pyarrow
                          but I'm not sure which one and don't have time for them to fix it.
        num_threads: number of threads to use for decompression of the VCF
    """
    raw_vcf_properties = _get_wanted_vcf_properties(wanted_properties)
    raw_properties_df = genomics_io.vcf_to_dask(
        input_vcf, wanted_properties=raw_vcf_properties, num_threads=num_threads
    )

    if raw_vcf_properties == wanted_properties:
        # we have the properties we want, just use them as is
        if reconcile_dtypes:
            dask_properties_df = dask.dataframe.from_delayed(
                dask_utils.reconcile_dtypes(list(dask.partitions)),
                divisions=raw_properties_df.divisions,
                verify_meta=False,
            )
        else:
            dask_properties_df = raw_properties_df
    else:
        # we need to compute / drop some properties
        partitions = [
            dask.delayed(_compute_wanted_partition_properties)(
                raw_vcf_partition=_partition,
                raw_vcf_properties=raw_vcf_properties,
                wanted_properties=wanted_properties,
                drop_rows_with_missing_values=drop_rows_with_missing_values,
            )
            for _partition in raw_properties_df.partitions
        ]
        dask_properties_df = dask.dataframe.from_delayed(
            dask_utils.reconcile_dtypes(partitions) if reconcile_dtypes else partitions,
            divisions=raw_properties_df.divisions,
            verify_meta=False,
        )

    tarred_properties_to_parquet.df_to_parquet(
        dask_properties_df, output_path=output_parquet, compression_algorithm=compression_algorithm
    )


def _get_wanted_vcf_properties(wanted_dask_properties: list[str]) -> list[str]:
    """
    Get the list of properties to load from the VCF, given the list of final wanted properties.
    """
    return list({
        (
            Keys.contig
            if prop == Keys.is_autosome
            else Keys.gt
            if prop == Keys.allele_count or prop == Keys.allele_frequency
            else prop
        )
        for prop in wanted_dask_properties
    })


def _compute_wanted_partition_properties(
        raw_vcf_partition: pandas.DataFrame,
        raw_vcf_properties: list[str],
        wanted_properties: list[str],
        drop_rows_with_missing_values: bool,
) -> pandas.DataFrame:
    """
    Compute the final wanted properties for this partition of the dask DataFrame.
    """
    drop_properties = set(raw_vcf_properties).difference(wanted_properties)
    concat_dfs: list[pandas.DataFrame] = []
    if Keys.allele_count in wanted_properties or Keys.allele_frequency in wanted_properties:
        allele_counts_df = _get_partition_allele_counts(raw_vcf_partition)
        if Keys.allele_count in wanted_properties:
            concat_dfs.append(allele_counts_df)
        if Keys.allele_frequency in wanted_properties:
            allele_frequencies_df = _get_partition_allele_frequencies(allele_counts_df)
            concat_dfs.append(allele_frequencies_df)
    if Keys.is_autosome in wanted_properties:
        concat_dfs.append(
            _get_partition_is_autosome(raw_vcf_partition)
        )
    if len(concat_dfs) == 0:
        # don't need to add new properties
        properties_partition = _drop_properties(raw_vcf_partition, drop_properties)
    else:
        concat_dfs.insert(0, _drop_properties(raw_vcf_partition, drop_properties))
        properties_partition = genomics_io.standard_sort_properties(
            pandas.concat(concat_dfs, axis=1)
        )
    if drop_rows_with_missing_values:
        has_missing_values = properties_partition.isna().any(axis=1)
        if has_missing_values.any():
            properties_partition = properties_partition.loc[~has_missing_values]
            for col in properties_partition.columns:
                col_property = properties_partition[col]
                if pandas.api.types.is_extension_array_dtype(col_property):
                    properties_partition[col] = col_property.astype(
                        _make_non_nullable_dtype(col_property.dtype)
                    )
        if properties_partition.isna().any().any():
            raise ValueError("Still have missing values")
    else:
        if properties_partition.isna().any().any():
            print("allowing missing values")

    return properties_partition


def _make_non_nullable_dtype(
        dtype: pandas.core.dtypes.dtypes.Dtype
) -> pandas.core.dtypes.dtypes.Dtype:
    """Convert a nullable dtype to a non-nullable dtype (typically after dropping null values)."""
    return pandas.core.dtypes.cast.pandas_dtype(dtype.name.replace("I", "i").replace("U", "u"))


def _get_partition_allele_counts(raw_vcf_partition: pandas.DataFrame) -> pandas.DataFrame:
    """Get the allele counts """
    allele_counts_df: pandas.DataFrame = genomics_io.get_allele_counts(
        raw_vcf_partition, no_call_ac=-1
    )
    allele_counts_df.columns = pandas.MultiIndex.from_tuples(
        [(sample_id, Keys.allele_count) for sample_id in allele_counts_df.columns],
        name=Default.column_levels
    )
    return allele_counts_df


def _get_partition_allele_frequencies(allele_counts_df: pandas.DataFrame) -> pandas.DataFrame:
    num_alleles = 2 * (allele_counts_df >= 0).sum(axis=1)
    num_non_ref_alleles = allele_counts_df.mask(allele_counts_df < 0, 0).sum(axis=1)

    allele_frequencies_df = (
        (num_non_ref_alleles / num_alleles).fillna(0).astype(numpy.float32).to_frame()
    )
    allele_frequencies_df.columns = pandas.MultiIndex.from_tuples(
        [(None, Keys.allele_frequency)], name=Default.column_levels
    )
    return allele_frequencies_df


def _get_partition_is_autosome(raw_vcf_partition: pandas.DataFrame) -> pandas.DataFrame:
    is_autosome: pandas.DataFrame = genomics_io.get_is_autosome(raw_vcf_partition).to_frame()
    is_autosome.columns = pandas.MultiIndex.from_tuples(
        [(None, Keys.is_autosome)], name=Default.column_levels
    )
    return is_autosome


def _drop_properties(
    raw_vcf_partition: pandas.DataFrame, drop_properties: set[str]
) -> pandas.DataFrame:
    # don't need to add new properties
    if len(drop_properties) == 0:
        # don't need to drop any properties
        return raw_vcf_partition
    else:
        # need to drop some properties
        keep_cols = [col for col in raw_vcf_partition.columns if col[1] not in drop_properties]
        return raw_vcf_partition[keep_cols]


def __parse_arguments(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract properties from a VCF file and save them in a parquet file",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        prog=argv[0]
    )
    parser.add_argument("--input_vcf", "-i", type=Path, help="Path to the input VCF file"),
    parser.add_argument(
        "--output-parquet", "-o",
        type=Path,
        help=(
            "path to output SV data in parquet format. If path ends in .tar, the parquet folder "
            "will be archived into a single tar file"
        )
    )
    parser.add_argument("--wanted-properties", type=str, nargs="+",
                        help="Comma-separated list of properties to extract from the VCF.")
    parser.add_argument("--temp-dir", "-t", type=Path, default=Default.temp_dir,
                        help="Path to preferred temporary directory")
    parser.add_argument(
        "--compression-algorithm", type=str, default=Default.compression_algorithm,
        choices=CompressionAlgorithms.list(), help="compression algorithm for parquet data"
    )
    parser.add_argument(
        "--drop-rows-with-missing-values", action="store_true",
        help="If true, drop rows/variants with missing values, if false, they will be included as"
             " missing values in the df."
    )
    parser.add_argument("--threads", type=int, default=Default.num_threads,
                        help="Number of parallel processes to use")
    return parser.parse_args(argv[1:] if len(argv) > 1 else ["--help"])


def main(argv: Optional[list[str]] = None):
    args = __parse_arguments(sys.argv if argv is None else argv)
    cluster = LocalCluster(n_workers=args.threads, silence_logs=logging.ERROR)
    client = Client(cluster)
    wanted_properties = (
        [args.wanted_properties]
        if isinstance(args.wanted_properties, str)
        else args.wanted_properties
    )
    wanted_properties = [p.strip() for arg in wanted_properties for p in arg.split(",")]
    with dask.config.set(temporary_directory=args.temp_dir):
        vcf_to_dask_parquet(
            input_vcf=args.input_vcf,
            output_parquet=args.output_parquet,
            wanted_properties=wanted_properties,
            compression_algorithm=args.compression_algorithm,
            drop_rows_with_missing_values=args.drop_rows_with_missing_values,
            num_threads=args.threads
        )


if __name__ == "__main__":
    main()
