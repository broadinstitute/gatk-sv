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
import pyarrow.parquet

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import dask.dataframe
from sv_utils import common, genomics_io
from typing import Optional, Union
from collections.abc import Mapping, Iterator
from types import MappingProxyType

# Mapping holds data for summary of an individual property
PropertySummary = dict[str, Union[int, str, list[str]]]
# Mapping holds data for summaries of all properties
PropertiesSummary = dict[str, PropertySummary]
DataFrame = dask.dataframe.DataFrame


class Keys:
    num_samples = "num_samples"
    num_rows = "num_rows"
    codes = "codes"
    type_name = "type"
    start_row = "start_row"
    sample_id = genomics_io.Keys.sample_id
    property = genomics_io.Keys.property
    float_type = "float"
    double_type = "double"
    byte_type = "byte"
    short_type = "short"
    int_type = "int"
    string_type = "String"
    string_set_type = "StringSet"
    row = "row"
    category = "category"
    one_hot = "one-hot"


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
    compression_algorithm = CompressionAlgorithms.snappy
    error_on_missing_property = True
    error_on_missing_sample = True
    category_encoding = Keys.one_hot
    property_names_to_lower = True
    num_jobs = common.num_physical_cpus


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
    parser.add_argument("--input_tar", "-i", type=str, help="full path to .tar file to convert", action="extend",
                        required=True, nargs='+')
    parser.add_argument("--output_parquet", "-o", type=str, required=True,
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
    parser.add_argument("--category-encoding", type=str, default=Default.category_encoding,
                        choices=[Keys.category, Keys.one_hot],
                        help=f"Method to encode categorical data: {Keys.category}: as pandas Categorical, or "
                             f"{Keys.one_hot}: one-hot encode. Note: contig is always kept as Categorical.")
    parser.add_argument("--property-names-to-lower", type=bool, default=Default.property_names_to_lower,
                        help="Convert all property names to lower case, to make interactions simpler.")
    parser.add_argument("--threads", type=int, default=Default.num_jobs, help="Number of parallel processes to use")
    return parser.parse_args(argv[1:] if len(argv) > 1 else ["--help"])


def main(argv: Optional[list[str]] = None):
    args = __parse_arguments(sys.argv if argv is None else argv)
    tarred_properties_to_parquet(
        input_tar=args.input_tar,
        output_path=args.output_parquet,
        temp_dir=args.temp_dir,
        remove_input_tar=args.remove_input_tar,
        dask_partition_size_mb=args.dask_partition_size_mb,
        compression_algorithm=args.compression_algorithm,
        category_encoding=args.category_encoding,
        property_names_to_lower=args.property_names_to_lower,
        num_jobs=args.threads
    )


def tarred_properties_to_parquet(
        input_tar: Union[str, list[str]],
        output_path: str,
        temp_dir: str = Default.temp_dir,
        remove_input_tar: bool = Default.remove_input_tar,
        dask_partition_size_mb: float = Default.dask_partition_size_mb,
        compression_algorithm: str = Default.compression_algorithm,
        category_encoding: str = Default.category_encoding,
        property_names_to_lower: bool = Default.property_names_to_lower,
        num_jobs: int = Default.num_jobs
):
    CompressionAlgorithms.assert_valid(compression_algorithm)
    properties_summaries, sample_ids, input_folders = _get_gzipped_metadata(
        input_tar=input_tar, temp_dir=temp_dir, remove_input_tar=remove_input_tar
    )

    with dask.config.set(temporary_directory=temp_dir, scheduler='processes', num_workers=num_jobs):
        properties = load_properties_from_map(
            input_folder=input_folders, properties_summary=properties_summaries, sample_ids=sample_ids,
            dask_partition_size_mb=dask_partition_size_mb, category_encoding=category_encoding,
            property_names_to_lower=property_names_to_lower
        )
        with warnings.catch_warnings():
            warnings.filterwarnings(action="ignore", category=FutureWarning)
            print(properties.head())
        df_to_parquet(df=properties, output_path=output_path, compression_algorithm=compression_algorithm)


def _get_gzipped_metadata(
        input_tar: Union[str, list[str]],
        temp_dir: str = Default.temp_dir,
        remove_input_tar: bool = Default.remove_input_tar,
) -> (list[PropertiesSummary], list[str], list[str]):
    print("Getting metadata")
    # explicitly convert to list
    input_tars = [input_tar] if isinstance(input_tar, str) else input_tar
    properties_summaries = []
    sample_ids = None
    input_folders = []
    global_codes_map = {}
    for input_tar in input_tars:
        try:
            print(f"\tdecompress {input_tar}")
            input_folder = extract_tar_to_folder(input_tar=input_tar, base_dir=temp_dir,
                                                 remove_input_tar=remove_input_tar)
            input_folders.append(input_folder)
            _sample_ids = _get_sample_ids(input_folder)
            if sample_ids is None:
                sample_ids = _sample_ids
            else:
                if sample_ids != _sample_ids:
                    raise ValueError(f"sample IDs in {input_tar} differ from those in {input_tars[0]}")
            properties_summary = _get_property_summaries(input_folder)
            _validate_properties_summary(properties_summary, properties_summaries, len(sample_ids))
            # if there are multiple tar files, some codes may not be present in some of them: if a particular variant
            # type was not present in a chunk of the VCF. The order will bary too, but we can ensure that every codes
            # dict at least has all the same codes.
            for prop_name, property_summary in properties_summary.items():
                codes = property_summary.get(Keys.codes, [])
                if codes:
                    global_codes_map[prop_name] = global_codes_map.get(prop_name, set()).union(codes)
            properties_summaries.append(properties_summary)
        except (Exception, IsADirectoryError) as exception:
            common.add_exception_context(exception, f"Getting metadata from {input_tar}")
            raise

    # add any missing codes
    for prop_name, all_codes in global_codes_map.items():
        for properties_summary in properties_summaries:
            existing_codes = properties_summary[prop_name].get(Keys.codes, [])
            missing_codes = all_codes.difference(existing_codes)
            if missing_codes:
                properties_summary[prop_name][Keys.codes] = existing_codes + sorted(missing_codes)

    print("\tgetting metadata complete.")
    return properties_summaries, sample_ids, input_folders


def extract_tar_to_folder(
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


def _get_property_summaries(input_folder: str) -> PropertiesSummary:
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


def _validate_properties_summary(
        properties_summary: PropertiesSummary,
        previous_properties_summaries: list[PropertiesSummary],
        num_samples: int
):
    num_rows = None
    first_row_prop = None
    for prop_name, prop_summary in properties_summary.items():
        if num_rows is None:
            num_rows = prop_summary[Keys.num_rows]
            first_row_prop = prop_name
        elif prop_summary[Keys.num_rows] != num_rows:
            raise ValueError(f"Inconsistent {Keys.num_rows}: {first_row_prop}={num_rows}, "
                             f"{prop_name}={prop_summary[Keys.num_rows]}")
        num_prop_samples = prop_summary[Keys.num_samples]
        if num_prop_samples > 0 and num_prop_samples != num_samples:
            raise ValueError(
                f"Property {prop_name} has {num_prop_samples} columns but there are {num_samples} samples"
            )
        codes = prop_summary.get(Keys.codes, [])
        if codes:
            if len(set(codes)) != len(codes):
                raise ValueError(f"Property {prop_name} has non-unique codes {codes}")

    if previous_properties_summaries:
        # ensure that property names, columns, types, and codes are consistent
        compare_summary = previous_properties_summaries[0]
        if set(compare_summary.keys()) != set(properties_summary.keys()):
            raise ValueError("Property names are not consistent")
        for prop_name, prop_summary in compare_summary.items():
            num_prop_samples = properties_summary[prop_name][Keys.num_samples]
            num_compare_samples = prop_summary[Keys.num_samples]
            if num_prop_samples != num_compare_samples:
                raise ValueError(
                    f"Inconsistent {Keys.num_samples} between data sets for {prop_name}: {num_compare_samples} -> "
                    f"{num_prop_samples}"
                )
            type_name = properties_summary[prop_name][Keys.type_name]
            compare_type = prop_summary[Keys.type_name]
            if type_name != compare_type:
                raise ValueError(
                    f"Inconsistent {Keys.type_name} between data sets for {prop_name}: {type_name} -> {compare_type}"
                )


def load_properties_from_map(
        input_folder: Union[str, list[str]],
        properties_summary: Union[PropertiesSummary, list[PropertiesSummary]],
        sample_ids: list[str],
        dask_partition_size_mb: float = Default.dask_partition_size_mb,
        category_encoding: str = Default.category_encoding,
        property_names_to_lower: bool = Default.property_names_to_lower
) -> DataFrame:
    if isinstance(input_folder, str):
        if not isinstance(properties_summary, Mapping):
            raise ValueError("Must specify the same numbers of input_folder and property_summaries")
        # explicitly convert to list
        input_folders = [input_folder]
        properties_summaries = [properties_summary]
    else:
        if not (isinstance(input_folder, list) and isinstance(properties_summary, list) and
                len(input_folder) == len(properties_summary)):
            raise ValueError("Must specify the same numbers of input_folder and property_summaries")
        input_folders = input_folder
        properties_summaries = properties_summary

    # decompress files, and form structure of the map
    divisions = [0]
    map_parameters = []
    start_row = 0
    for input_folder, properties_summary in zip(input_folders, properties_summaries):
        # unzip TSVs in parallel, wait until all complete
        dask.compute(
            *(
                dask.delayed(_decompress_gzipped_file)(os.path.join(input_folder, f"{property_name}.tsv.gz"))
                for property_name in properties_summary.keys()
            )
        )
        input_divisions = _get_dask_divisions(properties_summary, dask_partition_size_mb)
        divisions.extend(division + start_row for division in input_divisions[1:])
        map_parameters.extend(
            (input_divisions[_idx - 1], input_divisions[_idx], start_row, input_folder, properties_summary)
            for _idx in range(1, len(input_divisions))
        )
        start_row = start_row + input_divisions[-1]

    # form dask DataFrame from this map
    return dask.dataframe.from_map(
        _tsvs_to_df_partition,
        map_parameters,
        divisions=divisions,
        sample_ids=sample_ids,
        total_num_rows=start_row,
        category_encoding=category_encoding,
        property_names_to_lower=property_names_to_lower
    )


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
        return (bits_per_column // 8) * max(1, property_summary[Keys.num_samples])


def _tsvs_to_df_partition(
        division_parameters: tuple[int, int, int, str, Mapping[str, PropertySummary]],
        sample_ids: list[str],
        total_num_rows: int,
        category_encoding: str,
        property_names_to_lower: bool
) -> pandas.DataFrame:
    row_start, row_end, data_set_row_start, input_folder, property_summaries = division_parameters

    def _summary_sort_key(item):
        _name, _summary = item
        return _summary[Keys.num_samples], _name
    df = pandas.concat(
        tuple(
            _read_csv_chunk(
                input_folder=input_folder, property_name=property_name, property_summary=property_summary,
                sample_ids=sample_ids, row_start=row_start, row_end=row_end, category_encoding=category_encoding,
                property_names_to_lower=property_names_to_lower
            )
            for property_name, property_summary in sorted(property_summaries.items(), key=_summary_sort_key)
        ),
        axis=1
    )
    df.index = pandas.RangeIndex(
        row_start + data_set_row_start, row_end + data_set_row_start,
        dtype=numpy.dtype(numpy.min_scalar_type(-total_num_rows)), name=Keys.row
    )
    # return DataFrame sorted so that:
    #  1) per-variant columns are first, followed by per-sample/GT properties in blocks
    #  2) all the properties for a given sample are sequential, sorted in alphabetical order
    #  3) the samples are in the original order listed in the VCF
    # e.g. if the VCF had two samples: NA19088 and NA19001, listed in that order, the properties might look like:
    # None             NA19088         190001
    # SVTYPE   SVLEN   AC       GQ     AC       GQ
    return df[
        sorted(df.columns.to_list(),
               key=lambda c: (sample_ids.index(c[0]) if isinstance(c[0], str) else -1, c[1]))
    ]


def _read_csv_chunk(
        input_folder: str,
        property_name: str,
        property_summary: PropertySummary,
        sample_ids: list[str],
        row_start: int,
        row_end: int,
        category_encoding: str,
        property_names_to_lower: bool
) -> pandas.DataFrame:
    tsv = os.path.join(input_folder, f"{property_name}.tsv")
    if property_names_to_lower:
        property_name = property_name.lower()
    is_format_property = property_summary[Keys.num_samples] > 0
    if is_format_property:
        # this is a format / per-sample property
        column_names = sample_ids
        columns = pandas.MultiIndex.from_product([sample_ids, [property_name]], names=Default.column_levels)
    else:
        # this is an INFO / per-variant property
        column_names = [property_name]
        columns = pandas.MultiIndex.from_product([[None], [property_name]], names=Default.column_levels)
    type_name = property_summary[Keys.type_name]
    codes = property_summary.get(Keys.codes, [])
    _dtype = "string[pyarrow]" if property_name == "id" else \
        _dtype_map.get(type_name, numpy.min_scalar_type(len(codes)))
    dtypes = {column_name: _dtype for column_name in column_names}
    df = pandas.read_csv(
        tsv, sep='\t', low_memory=True, engine='c', memory_map=False, names=column_names,
        dtype=dtypes, nrows=row_end - row_start, skiprows=row_start
    )
    if codes:
        if category_encoding == Keys.category:
            # need to do convert ordinal ints to pandas Categorical
            df = pandas.DataFrame(
                pandas.Categorical.from_codes(df, categories=codes),
                index=df.index, columns=columns
            )
        else:
            # need to one-hot this DataFrame
            codewords_map = _get_codewords_map(codes=codes, type_name=type_name)
            df = codes_to_one_hot(df, codewords_map=codewords_map, property_name=property_name,
                                  is_format_property=is_format_property)
    else:
        df.columns = columns
    return df


def _get_codewords_map(codes: list[str], type_name: str) -> dict[str, numpy.array]:
    if type_name == Keys.string_type:
        # easy: the nth code is all false except for the nth value in the map
        num_codewords = len(codes)

        def _basis(_ind: int) -> numpy.ndarray:
            _arr = numpy.zeros(num_codewords, dtype=bool)
            _arr[_ind] = True
            return _arr
        return {code: _basis(ind) for ind, code in enumerate(codes)}
    else:
        code_sets = [set(code.split(',')) for code in codes]
        all_codewords = sorted(set.union(*code_sets))

        def _contains(_codeword: str) -> numpy.ndarray:
            return numpy.array([_codeword in code_set for code_set in code_sets], dtype=bool)
        return {codeword: _contains(codeword) for codeword in all_codewords}


def codes_to_one_hot(
        df: DataFrame,
        codewords_map: dict[str, numpy.array],
        property_name: str,
        is_format_property: bool
) -> DataFrame:
    mapped_dfs = {codeword: df.applymap(lambda code: code_included[code])
                  for codeword, code_included in codewords_map.items()}
    for codeword, mapped_df in mapped_dfs.items():
        one_hot_property_name = f"{property_name}={codeword}"
        mapped_df.columns = pandas.MultiIndex.from_product(
            [df.columns.tolist(), [one_hot_property_name]] if is_format_property else [[None], [one_hot_property_name]],
            names=Default.column_levels
        )

    return pandas.concat(
        [mapped_df for __, mapped_df in sorted(mapped_dfs.items(), key=lambda _item: _item[0])],
        axis=1
    )


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


def flatten_columns(df: DataFrame) -> pandas.Index:
    old_columns = df.columns
    if df.columns.nlevels > 1:
        df.columns = pandas.Index(
            [str(column) for column in old_columns.to_flat_index()]
        )
    return old_columns


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
        input_path = extract_tar_to_folder(input_path, base_dir=base_dir, remove_input_tar=remove_input_tar)
        print(f"extracted tarfile to {input_path}")
    with warnings.catch_warnings():
        df = dask.dataframe.read_parquet(
            input_path, engine="pyarrow", index=Keys.row, calculate_divisions=True,
            columns=_get_wanted_columns(input_path, wanted_properties, wanted_samples,
                                        error_on_missing_property=error_on_missing_property,
                                        error_on_missing_sample=error_on_missing_sample)
        )
    unflatten_columns(df)
    return df


def unflatten_column_name(flat_column_name: str) -> tuple[Optional[str], str]:
    try:
        return ast.literal_eval(flat_column_name.replace("(nan,", "(None,"))
    except ValueError as value_error:
        common.add_exception_context(value_error, f"unflattening '{flat_column_name}'")
        raise


def unflatten_columns(df: DataFrame):
    old_columns = df.columns
    try:
        multicolumns = pandas.MultiIndex.from_tuples(
            [unflatten_column_name(flat_column_name) for flat_column_name in old_columns]
        )
        if multicolumns.nlevels > 1:
            if multicolumns.nlevels == 2:
                multicolumns.names = Default.column_levels
            df.columns = multicolumns
    except (SyntaxError, ValueError):
        # couldn't unflatten, so probably not a flattened DataFrame. just return df unchanged
        pass
    return old_columns


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
        input_path = extract_tar_to_folder(input_path, base_dir=base_dir, remove_input_tar=remove_input_tar)
        print(f"extracted tarfile to {input_path}")
    pq_file = glob.glob(f"{input_path}/*.parquet")[0]
    schema = pyarrow.parquet.read_schema(pq_file)
    return schema.names


def get_parquet_file_num_rows(
        input_path: str,
        base_dir: Optional[str] = None,
        remove_input_tar: bool = Default.remove_input_tar
) -> int:
    if input_path.endswith(".tar"):
        input_path = extract_tar_to_folder(input_path, base_dir=base_dir, remove_input_tar=remove_input_tar)
        print(f"extracted tarfile to {input_path}")

    @dask.delayed
    def _get_nrows(_pq_file: str) -> int:
        return pyarrow.parquet.ParquetFile(_pq_file).metadata.num_rows

    return dask.compute(sum(_get_nrows(_pq_file) for _pq_file in glob.glob(f"{input_path}/*.parquet")))[0]


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
