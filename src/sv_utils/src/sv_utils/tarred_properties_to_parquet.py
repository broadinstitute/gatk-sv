#!/usr/bin/env python

import sys
import os
import ast
import glob
import warnings
import argparse
import gzip
import shutil
import tarfile
import tempfile
import json
import numpy
import pandas
import dask
import pyarrow

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import dask.dataframe
    import dask.array
from sv_utils import common, genomics_io
from typing import Optional, Union
from collections.abc import Mapping, Iterator
from types import MappingProxyType

PropertySummary = dict[str, Union[int, str, list[str]]]
DataFrame = dask.dataframe.DataFrame


class Keys:
    num_columns = "columns"
    num_rows = "rows"
    codes = "codes"
    type_name = "type"
    sample_id = genomics_io.Keys.sample_id
    property = genomics_io.Keys.property
    float_type = "float"
    double_type = "double"
    byte_type = "byte"
    short_type = "short"
    int_type = "int"
    row = "row"


class CompressionAlgorithms:
    uncompressed = "uncompressed"
    brotli = "brotli"
    gzip = "gzip"
    lz4 = "lz4"
    lz4_raw = "lz4_raw"
    snappy = "snappy"
    zstd = "zstd"

    @classmethod
    def list(cls) -> list[str]:
        # noinspection PyTypeChecker
        return [name for name, val in cls.__dict__.items() if isinstance(val, str) and not name.startswith('_')]

    @staticmethod
    def assert_valid(algorithm: str):
        allowed_vals = set(CompressionAlgorithms.list())
        if algorithm not in allowed_vals:
            raise ValueError(f"Invalid compression algorithm: {algorithm}. Allowed values are {allowed_vals}")


class Default:
    temp_dir = tempfile.gettempdir()
    remove_input_tar = False
    column_levels = genomics_io.Default.column_levels
    dask_partition_size_mb = 75.0
    compression_algorithm = CompressionAlgorithms.zstd,
    error_on_missing_property = True
    error_on_missing_sample = True


_dtype_map = MappingProxyType({
    Keys.float_type: numpy.float32,
    Keys.double_type: numpy.float64,
    Keys.byte_type: numpy.int8,
    Keys.short_type: numpy.int16,
    Keys.int_type: numpy.int32
})


_extracted_tar_files = {}


def __parse_arguments(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert properties produced by ExtractSV properties from tarred-gzipped TSV to a parquet data set",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        prog=argv[0]
    )
    parser.add_argument("input_tar", type=str, help="full path to .tar file to convert")
    parser.add_argument("output_parquet", type=str,
                        help="path to output parquet data. If path ends in .tar, the parquet folder will be archived "
                             "into a single tar file")
    parser.add_argument("--temp-dir", "-t", type=str, default=Default.temp_dir,
                        help="full path to preferred temporary directory")
    parser.add_argument("--remove-input-tar", type=bool, default=Default.remove_input_tar,
                        help="if true, remove input tar to save space, if false, leave it in place")
    parser.add_argument("--dask-partition-size-mb", type=float, default=Default.dask_partition_size_mb,
                        help="approximate in-memory size of dask partitions, in MiB")
    parser.add_argument("--compression-algorithm", type=str, default=Default.compression_algorithm,
                        choices=CompressionAlgorithms.list(), help="compression algorithm for parquet data")
    return parser.parse_args(argv[1:] if len(argv) > 1 else ["--help"])


def main(argv: Optional[list[str]] = None):
    args = __parse_arguments(sys.argv if argv is None else argv)
    os.environ["PYTHONWARNINGS"] = "ignore"
    tarred_properties_to_parquet(
        input_tar=args.input_tar,
        output_path=args.output_parquet,
        temp_dir=args.temp_dir,
        remove_input_tar=args.remove_input_tar,
        dask_partition_size_mb=args.dask_partition_size_mb,
        compression_algorithm=args.compression_algorithm
    )


def tarred_properties_to_parquet(
        input_tar: str,
        output_path: str,
        temp_dir: str = Default.temp_dir,
        remove_input_tar: bool = Default.remove_input_tar,
        dask_partition_size_mb: float = Default.dask_partition_size_mb,
        compression_algorithm: str = Default.compression_algorithm
):
    CompressionAlgorithms.assert_valid(compression_algorithm)
    input_folder = _extract_tar_to_folder(input_tar=input_tar, base_dir=temp_dir, remove_input_tar=remove_input_tar)
    print(f"input_folder={input_folder}")
    try:
        property_summaries = _get_property_summaries(input_folder)
    except Exception as exception:
        common.add_exception_context(exception, f"Getting property summaries JSON from {input_tar}")
        raise
    try:
        sample_ids = _get_sample_ids(input_folder)
    except Exception as exception:
        common.add_exception_context(exception, f"Getting sample IDs from {input_tar}")
        raise

    with dask.config.set(temporary_directory=temp_dir, scheduler='processes'):
        properties = load_properties_from_map(
            input_folder=input_folder, property_summaries=property_summaries, sample_ids=sample_ids,
            dask_partition_size_mb=dask_partition_size_mb
        )
        with warnings.catch_warnings():
            warnings.filterwarnings(action="ignore", category=FutureWarning)
            print(properties.head())
        df_to_parquet(df=properties, output_path=output_path, compression_algorithm=compression_algorithm)


def flatten_columns(df: DataFrame) -> pandas.Index:
    old_columns = df.columns
    if df.columns.nlevels > 1:
        df.columns = pandas.Index(
            [str(column) for column in old_columns.to_flat_index()]
        )
    return old_columns


def unflatten_columns(df: DataFrame):
    old_columns = df.columns
    try:
        multicolumns = pandas.MultiIndex.from_tuples(
            [ast.literal_eval(x.replace("(nan,", "(None,")) for x in old_columns]
        )
        if multicolumns.nlevels > 1:
            if multicolumns.nlevels == 2:
                multicolumns.names = Default.column_levels
            df.columns = multicolumns
    except (SyntaxError, ValueError):
        # couldn't unflatten, so probably not a flattened DataFrame. just return df unchanged
        pass
    return old_columns


def df_to_parquet(
        df: DataFrame,
        output_path: str,
        compression_algorithm: str = Default.compression_algorithm
):
    CompressionAlgorithms.assert_valid(compression_algorithm)
    # flatten columns to be compatible with parquet
    old_columns = flatten_columns(df)
    if output_path.endswith(".tar"):
        tar_file = output_path
        output_path = output_path.rsplit(".tar")[0] + ".parquet"
    else:
        tar_file = None
    with warnings.catch_warnings():
        warnings.filterwarnings(action="ignore", category=FutureWarning)
        df.to_parquet(output_path, engine="pyarrow", compression=compression_algorithm)

    # restore the old columns
    df.columns = old_columns

    if tar_file is not None:
        archive_to_tar(output_path, tar_file, remove_files=True)


def archive_to_tar(input_folder: str, tar_file: Optional[str] = None, remove_files: bool = False):
    if tar_file is None:
        tar_file = input_folder + ".tar"
    tar_root_folder = os.path.dirname(input_folder)
    with tarfile.open(tar_file, "w") as tar_out:
        for file_to_archive in recursive_ls_dir(input_folder):
            tar_out.add(file_to_archive, arcname=os.path.relpath(file_to_archive, tar_root_folder))
            if remove_files:
                os.remove(file_to_archive)
    if remove_files:
        os.rmdir(input_folder)


def recursive_ls_dir(folder_path: str) -> Iterator[str]:
    for root, dirs, files in os.walk(folder_path):
        for f in files:
            yield os.path.join(root, f)


def parquet_to_df(
        input_path: str,
        base_dir: Optional[str] = None,
        remove_input_tar: bool = Default.remove_input_tar,
        wanted_properties: Optional[list[str]] = None,
        wanted_samples: Optional[list[str]] = None,
        error_on_missing_property: bool = Default.error_on_missing_property,
        error_on_missing_sample: bool = Default.error_on_missing_sample
) -> DataFrame:
    if input_path.endswith(".tar"):
        input_path = _extract_tar_to_folder(input_path, base_dir=base_dir, remove_input_tar=remove_input_tar)
        print(f"extracted tarfile to {input_path}")
    with warnings.catch_warnings():
        df = dask.dataframe.read_parquet(
            input_path, engine="pyarrow", index=Keys.row,
            columns=_get_wanted_columns(input_path, wanted_properties, wanted_samples,
                                        error_on_missing_property=error_on_missing_property,
                                        error_on_missing_sample=error_on_missing_sample)
        )
    unflatten_columns(df)
    return df


def _extract_tar_to_folder(
        input_tar: str,
        base_dir: Optional[str] = None,
        remove_input_tar: bool = Default.remove_input_tar
) -> str:
    if input_tar in _extracted_tar_files:
        if remove_input_tar and os.path.isfile(input_tar):
            os.remove(input_tar)
        return _extracted_tar_files[input_tar]
    if base_dir is None:
        base_dir = os.path.dirname(input_tar)
    with tarfile.open(input_tar, 'r') as tar_in:
        tar_folder = os.path.join(base_dir, os.path.basename(os.path.commonpath(tar_in.getnames())))
        _extracted_tar_files[input_tar] = tar_folder
        tar_in.extractall(path=base_dir)
    if remove_input_tar:
        os.remove(input_tar)
    return tar_folder


def _get_wanted_columns(
        input_path: str,
        wanted_properties: Optional[list[str]] = None,
        wanted_samples: Optional[list[str]] = None,
        error_on_missing_property: bool = Default.error_on_missing_property,
        error_on_missing_sample: bool = Default.error_on_missing_sample,
        file_columns: Optional[list[str]] = None
) -> Optional[list[str]]:
    if wanted_samples is None:
        if wanted_properties is None:
            return None  # get all columns
        # get columns with desired properties
        if file_columns is None:
            file_columns = get_parquet_file_columns(input_path)
        wanted_columns = [
            col for col in file_columns
            if any(f"'{prop}'" in col for prop in wanted_properties)
        ]
        if error_on_missing_property:
            for prop in wanted_properties:
                if not any(f"'{prop}'" in col for col in wanted_columns):
                    raise ValueError(f"Wanted property '{prop}' not present in {input_path}")
    elif wanted_properties is None:
        # get columns with desired samples (or no sample / per-variant data)
        if file_columns is None:
            file_columns = get_parquet_file_columns(input_path)
        wanted_columns = [
            col for col in file_columns
            if any(f"'{sample}'" in col for sample in wanted_samples)
        ]
        if error_on_missing_sample:
            for sample in wanted_samples:
                if not any(f"'{sample}'" in col for col in wanted_columns):
                    raise ValueError(f"Wanted sample {sample} not present in {input_path}")
    else:
        # get columns with specified samples and properties
        file_columns = get_parquet_file_columns(input_path)
        sample_wanted_columns = set(
            _get_wanted_columns(
                input_path, wanted_properties=None, wanted_samples=wanted_samples,
                error_on_missing_sample=error_on_missing_sample, file_columns=file_columns
            )
        )
        wanted_columns = [
            col for col in _get_wanted_columns(
                input_path, wanted_properties=wanted_properties, wanted_samples=None,
                error_on_missing_property=error_on_missing_property, file_columns=file_columns
            )
            if col in sample_wanted_columns
        ]
    return wanted_columns


def get_parquet_file_columns(
        input_path: str,
        base_dir: Optional[str] = None,
        remove_input_tar: bool = Default.remove_input_tar
) -> list[str]:
    if input_path.endswith(".tar"):
        input_path = _extract_tar_to_folder(input_path, base_dir=base_dir, remove_input_tar=remove_input_tar)
        print(f"extracted tarfile to {input_path}")
    pq_file = glob.glob(f"{input_path}/*.parquet")[0]
    schema = pyarrow.parquet.read_schema(pq_file)
    return schema.names


def _get_property_summaries(input_folder: str) -> dict[str, PropertySummary]:
    summary_json = next(
        (os.path.join(input_folder, name) for name in os.listdir(input_folder) if name.endswith(".json")),
        None
    )
    if summary_json is None:
        raise ValueError("No summaries .json present")
    with open(summary_json, 'r') as f_in:
        return json.load(f_in)


def _get_sample_ids(input_folder: str) -> list[str]:
    samples_list_file = next(
        (os.path.join(input_folder, name) for name in os.listdir(input_folder) if name.endswith(".list")),
        None
    )
    if samples_list_file is None:
        raise ValueError("No sample IDs list is present")
    with open(samples_list_file, 'r') as f_in:
        return [line.strip() for line in f_in]


def load_properties_from_map(
        input_folder: str,
        property_summaries: Mapping[str, PropertySummary],
        sample_ids: list[str],
        dask_partition_size_mb: float = Default.dask_partition_size_mb
) -> DataFrame:
    # unzip TSVs in parallel, wait until all complete
    dask.compute(
        *(
            dask.delayed(_decompress_gzipped_file)(os.path.join(input_folder, f"{property_name}.tsv.gz"))
            for property_name in property_summaries.keys()
        )
    )

    divisions = _get_dask_divisions(property_summaries, dask_partition_size_mb)
    division_ranges = [(divisions[_idx - 1], divisions[_idx]) for _idx in range(1, len(divisions))]

    # sort columns so that all properties for given sample ID are adjacent?
    return dask.dataframe.from_map(
        _tsvs_to_df_partition,
        division_ranges,
        divisions=divisions,
        input_folder=input_folder,
        property_summaries=property_summaries,
        sample_ids=sample_ids
    )


def _get_dask_divisions(
        property_summaries: Mapping[str, PropertySummary],
        dask_partition_size_mb: float = Default.dask_partition_size_mb
) -> list[int]:
    num_rows = next(iter(property_summaries.values()))[Keys.num_rows]
    bytes_per_row = sum(_get_bytes_per_row(property_name, property_summary)
                        for property_name, property_summary in property_summaries.items())
    rows_per_division = int(numpy.ceil(dask_partition_size_mb * (2 ** 20) / bytes_per_row))
    return list(range(0, num_rows, rows_per_division)) + [num_rows]


def _get_bytes_per_row(property_name: str, property_summary: PropertySummary) -> int:
    if property_name == "ID":
        return 20  # wild guess as to string length, won't be a dominant factor
    else:
        _dtype = _dtype_map.get(property_summary[Keys.type_name], numpy.uint8)
        bits_per_column = numpy.finfo(_dtype).bits if numpy.issubdtype(_dtype, numpy.floating) \
            else numpy.iinfo(_dtype).bits
        return (bits_per_column // 8) * property_summary[Keys.num_columns]


def _tsvs_to_df_partition(
        division_range: tuple[int, int],
        input_folder: str,
        property_summaries: Mapping[str, PropertySummary],
        sample_ids: list[str]
) -> pandas.DataFrame:
    def _summary_sort_key(item):
        _name, _summary = item
        return _summary[Keys.num_columns], _name

    row_start, row_end = division_range

    df = pandas.concat(
        tuple(
            _read_csv_chunk(
                input_folder=input_folder, property_name=property_name, property_summary=property_summary,
                sample_ids=sample_ids, row_start=row_start, row_end=row_end
            )
            for property_name, property_summary in sorted(property_summaries.items(), key=_summary_sort_key)
        ),
        axis=1
    )
    num_rows = next(iter(property_summaries.values()))[Keys.num_rows]
    df.index = pandas.RangeIndex(row_start, row_end, dtype=numpy.dtype(numpy.min_scalar_type(-num_rows)), name=Keys.row)
    return df


def _read_csv_chunk(
        input_folder: str,
        property_name: str,
        property_summary: PropertySummary,
        sample_ids: list[str],
        row_start: int,
        row_end: int
) -> pandas.DataFrame:
    tsv = os.path.join(input_folder, f"{property_name}.tsv")
    if property_summary[Keys.num_columns] > 1:
        # this is a format / per-sample property
        column_names = sample_ids
        columns = pandas.MultiIndex.from_product([sample_ids, [property_name]], names=Default.column_levels)
    else:
        # this is an INFO / per-variant property
        column_names = [property_name]
        columns = pandas.MultiIndex.from_product([[None], [property_name]], names=Default.column_levels)
    type_name = property_summary[Keys.type_name]
    _dtype = "string[pyarrow]" if property_name == "ID" else _dtype_map.get(type_name, numpy.uint8)
    dtypes = {column_name: _dtype for column_name in column_names}
    df = pandas.read_csv(
        tsv, sep='\t', low_memory=True, engine='c', memory_map=False, names=column_names,
        dtype=dtypes, nrows=row_end - row_start, skiprows=row_start
    )
    df.columns = columns
    return df


def _decompress_gzipped_file(
        compressed_file: str,
        remove_compressed: bool = True
) -> str:
    if not compressed_file.endswith(".gz"):
        raise ValueError("compressed file does not end with .gz extension")
    decompressed_file = compressed_file.rsplit(".", 1)[0]
    with gzip.open(compressed_file, "rb") as f_in, open(decompressed_file, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
    if remove_compressed:
        os.remove(compressed_file)
    return decompressed_file


def get_value_counts(property_df: DataFrame) -> dict[int, int]:
    def _red(_vcs: list[dict[int, int]]) -> dict[int, int]:
        _out = _vcs[0]
        for _vc in _vcs[1:]:
            for _x, _v in _vc.items():
                if _x in _out:
                    _out[_x] += _v
                else:
                    _out[_x] = _v
        return _out

    return property_df.reduction(
        lambda p: pandas.Series(p.values.ravel()).value_counts().sort_index().to_dict(), _red, meta={}
    ).compute()[0]


if __name__ == "__main__":
    main()
