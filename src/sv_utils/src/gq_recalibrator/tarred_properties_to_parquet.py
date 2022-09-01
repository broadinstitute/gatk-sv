#!/usr/bin/env python

import sys
import os
from pathlib import Path
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
import attr

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import dask.dataframe
from sv_utils import common, genomics_io
from gq_recalibrator import dask_utils
from typing import Optional, Union, Any
from collections.abc import Mapping, Iterable, Iterator
from types import MappingProxyType
import enum


class Keys:
    num_samples = "num_samples"
    num_rows = "num_rows"
    codes = "codes"
    property_type = "type"
    start_row = "start_row"
    sample_id = genomics_io.Keys.sample_id
    property = genomics_io.Keys.property
    row = genomics_io.Keys.row
    category = "category"
    one_hot = "one-hot"
    positive_value = "positive_value"
    negative_value = "negative_value"
    baseline = "baseline"
    scale = "scale"
    id = genomics_io.Keys.id
    begin = genomics_io.Keys.begin
    end = genomics_io.Keys.end
    contig = genomics_io.Keys.contig
    svtype = genomics_io.Keys.svtype
    svlen = genomics_io.Keys.svlen
    allele_frequency = genomics_io.Keys.allele_frequency
    variant_weights = "variant_weights"
    is_multiallelic = "is_multiallelic"
    is_autosome = "is_autosome"
    parquet_suffix = ".parquet"
    properties_scaling_json = "properties_scaling.json"
    properties_summary_json = "properties_summary.json"


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
        return [
            name for name, val in cls.__dict__.items()
            if isinstance(val, str) and not name.startswith('_')
        ]

    @staticmethod
    def assert_valid(algorithm: str):
        allowed_vals = set(CompressionAlgorithms.list())
        if algorithm not in allowed_vals:
            raise ValueError(
                f"Invalid compression algorithm: {algorithm}. Allowed values are {allowed_vals}"
            )


class Default:
    temp_dir = Path(tempfile.gettempdir())
    compute_weights = True
    remove_input_tar = False
    column_levels = genomics_io.Default.column_levels
    dask_partition_size_mb = 250.0
    compression_algorithm = CompressionAlgorithms.snappy
    error_on_missing_property = True
    error_on_missing_sample = True
    category_encoding = Keys.category
    property_names_to_lower = True
    num_jobs = common.num_physical_cpus
    training_weight_properties = ()
    omit_properties_stats = False


non_ml_properties = frozenset((
    Keys.id, Keys.begin, Keys.end, Keys.contig, Keys.variant_weights, Keys.is_multiallelic,
    Keys.is_autosome
))
training_weight_bins = MappingProxyType({
    Keys.svlen: numpy.array([50, 500, 5000]),
    Keys.allele_frequency: numpy.array([0.1, 0.9])
})


@attr.frozen(slots=True, weakref_slot=False)
class PropertyScaling:
    baseline: float
    scale: float
    positive_value: Optional[float] = None
    negative_value: Optional[float] = None

    def to_json(self) -> dict[str, float]:
        return {
            key: value for key, value in (
                (key, getattr(self, key)) for key in self.__slots__
            ) if value is not None
        }

    @classmethod
    def from_json(cls, json_dict: dict[str, float]) -> "PropertyScaling":
        return PropertyScaling(**json_dict)

    @property
    def is_bool_property(self) -> bool:
        return self.positive_value is not None


class PropertiesScaling(dict[str, PropertyScaling]):
    @classmethod
    def from_json(
            cls, json_dict: dict[str, float] | dict[str, PropertyScaling]
    ) -> Union[PropertyScaling, "PropertiesScaling"]:
        if isinstance(next(iter(json_dict.values())), PropertyScaling):
            return PropertiesScaling(json_dict)
        else:
            return PropertyScaling.from_json(json_dict)


# Mapping holds data for summary of an individual property
PropertySummaryJsonDict = dict[str, int | str | list[str]]
# Mapping holds baseline and scale for an individual property
PropertyStatsDict = dict[str | int, int | list[float]]
# Mapping holds statistics for all properties
AllPropertiesStatsDict = dict[str, PropertyStatsDict]
DataFrame = dask.dataframe.DataFrame
DType = numpy.dtype | pandas.api.extensions.ExtensionDtype


class PropertyType(enum.Enum):
    String = "String"
    StringSet = "StringSet"
    boolean = "boolean"
    byte = "byte"
    short = "short"
    int = "int"
    float = "float"
    double = "double"

    @property
    def dtype(self) -> DType:
        # noinspection PyProtectedMember
        return PropertyType._dtype_map.get(self)

    @property
    def has_codes(self) -> bool:
        return self == PropertyType.String or self == PropertyType.StringSet

    @property
    def is_floating_type(self) -> bool:
        return not self.has_codes and numpy.issubdtype(self.dtype, numpy.floating)

    @property
    def is_integer_type(self) -> bool:
        return not self.has_codes and numpy.issubdtype(self.dtype, numpy.integer)


# MappingProxyType[PropertyType, numpy.dtype | pandas.api.extensions.ExtensionDtype]
PropertyType._dtype_map = MappingProxyType({
    PropertyType.String: pandas.api.types.CategoricalDtype,
    PropertyType.float: numpy.dtype(numpy.float32),
    PropertyType.double: numpy.dtype(numpy.float64),
    PropertyType.boolean: numpy.dtype(bool),
    PropertyType.byte: numpy.dtype(numpy.int8),
    PropertyType.short: numpy.dtype(numpy.int16),
    PropertyType.int: numpy.dtype(numpy.int32)
})


@attr.frozen(slots=True, weakref_slot=False)
class PropertySummary:
    property_type: PropertyType
    num_rows: int
    num_samples: int = 0
    codes: list[str] = []

    @classmethod
    def from_json(
            cls, json_dict: PropertySummaryJsonDict) -> "PropertySummary":
        return PropertySummary(
            property_type=PropertyType(json_dict[Keys.property_type]),
            num_rows=json_dict[Keys.num_rows],
            num_samples=json_dict[Keys.num_samples],
            codes=json_dict.get(Keys.codes, [])
        )

    def to_json(self) -> dict[str, float]:
        return {
            Keys.property_type: self.property_type.value,
            Keys.num_rows: self.num_rows,
            Keys.num_samples: self.num_samples,
            Keys.codes: self.codes
        }

    @property
    def encoding_dtype(self) -> numpy.dtype:
        """get the pure numpy dtype that would be used for ordinal encoding of String/StringSet,
        for storing to disk; hence booleans should be represented as uint8 0 or 1
        """
        if self.has_codes:
            # use signed integers so that -1 can be reserved for missing data
            return numpy.dtype(numpy.min_scalar_type(-len(self.codes)))
        elif self.is_bool_type:
            return numpy.dtype(numpy.uint8)
        else:
            return self.property_type.dtype

    @property
    def decoding_dtype(self) -> numpy.dtype:
        """get the pure numpy dtype that would be used for ordinal encoding of String/StringSet,
        for loading from disk; hence booleans should be represented as bool.
        """
        if self.has_codes:
            return numpy.dtype(numpy.min_scalar_type(len(self.codes)))
        else:
            return self.property_type.dtype

    @property
    def dtype(self) -> DType:
        if self.has_codes:
            return pandas.CategoricalDtype(categories=self.codes, ordered=False)
        else:
            return self.property_type.dtype

    @property
    def has_codes(self) -> bool:
        return self.property_type.has_codes

    @property
    def is_floating_type(self) -> bool:
        return self.property_type.is_floating_type

    @property
    def is_integer_type(self) -> bool:
        return self.property_type.is_integer_type

    @property
    def is_bool_type(self) -> bool:
        return self.property_type == PropertyType.boolean

    def with_added_codes(self, added_codes: list[str]) -> "PropertySummary":
        return PropertySummary(
            property_type=self.property_type, num_rows=self.num_rows, num_samples=self.num_samples,
            codes=self.codes + added_codes
        )

    @staticmethod
    def combine(property_summaries: list["PropertySummary"]) -> "PropertySummary":
        property_types = {
            property_summary.property_type for property_summary in property_summaries
        }
        if len(property_types) != 1:
            raise ValueError(f"Inconsistent .property_type: {property_types}")
        num_samples = {
            property_summary.num_samples for property_summary in property_summaries
        }
        if len(num_samples) != 1:
            raise ValueError(f"Inconsistent num_samples: {num_samples}")
        codes = sorted(
            {code for property_summary in property_summaries for code in property_summary.codes}
        )
        return PropertySummary(
            property_type=property_types.pop(),
            num_rows=sum(summary.num_rows for summary in property_summaries),
            num_samples=num_samples.pop(),
            codes=codes
        )


class PropertiesSummary(dict[str, PropertySummary]):
    @classmethod
    def from_json(
            cls, json_dict: PropertySummaryJsonDict | dict[str, PropertySummary]
    ) -> Union[PropertySummary, "PropertiesSummary"]:
        if isinstance(next(iter(json_dict.values())), PropertySummary):
            return PropertiesSummary(json_dict)
        else:
            return PropertySummary.from_json(json_dict)

    @staticmethod
    def combine(properties_summaries: list["PropertiesSummary"]) -> "PropertiesSummary":
        """Combine the properties_summaries from each shard to make them consistent. Also convert
        from VCF names (generally upper-case) to table/DataFrame names (snake case) since this will
        be used to process Dataframes, not VCFs
        """
        first_summary = next(iter(properties_summaries))
        return PropertiesSummary(
            {
                key.lower(): PropertySummary.combine(
                    [properties_summary[key] for properties_summary in properties_summaries]
                ) for key in first_summary.keys()
            }
        )


# keep track of tar files that have already been extracted
_extracted_tar_files = {}


def __parse_arguments(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert properties produced by ExtractSV properties from tarred-gzipped TSV"
                    " to a parquet data set",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        prog=argv[0]
    )
    parser.add_argument("--input-tar", "-i", type=Path, help="full path to .tar file to convert",
                        action="extend", required=True, nargs='+')
    parser.add_argument(
        "--output-parquet", "-o", type=Path, required=True,
        help="path to output SV data in parquet format. If path ends in .tar, the parquet folder "
             "will be archived into a single tar file")
    parser.add_argument(
        "--compute-weights", type=common.argparse_bool, default=Default.compute_weights,
        help="If true, compute the weight that each variant should have for training."
    )
    parser.add_argument(
        "--omit-properties-stats", type=bool, default=Default.omit_properties_stats,
        help="The properties stats JSONs are only needed for training (filtering will use the "
             "values stored from the training set). If these data will ONLY be used for filtering "
             "then this argument can safely be set to True and time can be saved not calculating "
             "these values."
    )
    parser.add_argument("--temp-dir", "-t", type=Path, default=Default.temp_dir,
                        help="Path to preferred temporary directory")
    parser.add_argument(
        "--remove-input-tar", type=common.argparse_bool, default=Default.remove_input_tar,
        help="if true, remove input tar to save space, if false, leave it in place"
    )
    parser.add_argument(
        "--wanted-properties", type=str, default=None,
        help="comma-separated list of properties to load. If omitted, take all properties"
    )
    parser.add_argument(
        "--dask-partition-size-mb", type=float, default=Default.dask_partition_size_mb,
        help="approximate in-memory size of dask partitions, in MiB"
    )
    parser.add_argument(
        "--compression-algorithm", type=str, default=Default.compression_algorithm,
        choices=CompressionAlgorithms.list(), help="compression algorithm for parquet data"
    )
    parser.add_argument(
        "--category-encoding", type=str, default=Default.category_encoding,
        choices=[Keys.category, Keys.one_hot],
        help=f"Method to encode categorical data: {Keys.category}: as pandas Categorical, or "
             f"{Keys.one_hot}: one-hot encode. Note: contig is always kept as Categorical."
    )
    parser.add_argument(
        "--property-names-to-lower", type=common.argparse_bool,
        default=Default.property_names_to_lower,
        help="Convert all property names to lower case, to make interactions simpler."
    )
    parser.add_argument("--threads", type=int, default=Default.num_jobs,
                        help="Number of parallel processes to use")
    return parser.parse_args(argv[1:] if len(argv) > 1 else ["--help"])


def main(argv: Optional[list[str]] = None):
    args = __parse_arguments(sys.argv if argv is None else argv)
    tarred_properties_to_parquet(
        input_tar=args.input_tar,
        output_path=args.output_parquet,
        compute_weights=args.compute_weights,
        omit_properties_stats=args.omit_properties_stats,
        wanted_properties=(
            None if args.wanted_properties is None else set(args.wanted_properies.split(","))
        ),
        temp_dir=args.temp_dir,
        remove_input_tar=args.remove_input_tar,
        dask_partition_size_mb=args.dask_partition_size_mb,
        compression_algorithm=args.compression_algorithm,
        category_encoding=args.category_encoding,
        property_names_to_lower=args.property_names_to_lower,
        num_jobs=args.threads
    )


def tarred_properties_to_parquet(
        input_tar: Path | list[Path],
        output_path: Path,
        compute_weights: bool = Default.compute_weights,
        omit_properties_stats: bool = Default.omit_properties_stats,
        wanted_properties: Optional[set[str]] = None,
        temp_dir: Path = Default.temp_dir,
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
    with dask.config.set(temporary_directory=temp_dir, scheduler='processes',
                         num_workers=num_jobs):
        # load all the properties extracted by
        properties = load_properties_from_map(
            input_folder=input_folders, properties_summary=properties_summaries,
            sample_ids=sample_ids, wanted_properties=wanted_properties,
            dask_partition_size_mb=dask_partition_size_mb, category_encoding=category_encoding,
            property_names_to_lower=property_names_to_lower
        )
        if compute_weights:
            # add variant_weights property
            properties = add_variant_weights(properties)
        with warnings.catch_warnings():
            warnings.filterwarnings(action="ignore", category=FutureWarning)

        # we want to add metadata into the parquet folder, so if it will be tarred, delay that
        # action
        if output_path.suffix == ".tar":
            tar_path = output_path
            parquet_folder = output_path.with_suffix("")
        else:
            tar_path = None
            parquet_folder = output_path
        df_to_parquet(
            df=properties, output_path=parquet_folder, compression_algorithm=compression_algorithm
        )

        if not omit_properties_stats:
            properties_scaling = get_properties_scaling(properties)
            with open(parquet_folder / Keys.properties_scaling_json, 'w') as f_out:
                json.dump(properties_scaling, f_out, indent="  ", default=PropertyScaling.to_json)
            with open(parquet_folder / Keys.properties_summary_json, 'w') as f_out:
                json.dump(PropertiesSummary.combine(properties_summaries), f_out, indent="  ",
                          default=PropertySummary.to_json)

        if tar_path is not None:
            archive_to_tar(folder_to_archive=parquet_folder, tar_file=tar_path, remove_files=True)


def _get_gzipped_metadata(
        input_tar: Path | list[Path],
        temp_dir: Path = Default.temp_dir,
        remove_input_tar: bool = Default.remove_input_tar,
) -> (list[PropertiesSummary], list[str], list[Path]):
    print("Getting metadata")
    # explicitly convert to list
    input_tars = [input_tar] if isinstance(input_tar, Path) else input_tar
    properties_summaries = []
    sample_ids = None
    input_folders = []
    global_codes_map = {}
    for input_tar in input_tars:
        try:
            print(f"\tdecompress {input_tar}")
            input_folder = extract_tar_to_folder(
                input_tar=input_tar, base_dir=temp_dir, remove_input_tar=remove_input_tar
            )
            input_folders.append(input_folder)
            _sample_ids = _get_sample_ids(input_folder)
            if sample_ids is None:
                sample_ids = _sample_ids
            else:
                if sample_ids != _sample_ids:
                    raise ValueError(
                        f"sample IDs in {input_tar} differ from those in {input_tars[0]}"
                    )
            properties_summary = _get_property_summaries(input_folder)
            _validate_properties_summary(properties_summary, properties_summaries, len(sample_ids))
            # if there are multiple tar files, some codes may not be present in some of them: if a
            # particular variant type was not present in a chunk of the VCF. The order will vary
            # too, but we can ensure that every codes dict at least has all the same codes.
            for prop_name, property_summary in properties_summary.items():
                codes = property_summary.codes
                if codes:
                    global_codes_map[prop_name] = \
                        global_codes_map.get(prop_name, set()).union(codes)
            properties_summaries.append(properties_summary)
        except (Exception, IsADirectoryError) as exception:
            common.add_exception_context(exception, f"Getting metadata from {input_tar}")
            raise

    # add any missing codes
    for prop_name, all_codes in global_codes_map.items():
        for properties_summary in properties_summaries:
            existing_codes = properties_summary[prop_name].codes
            missing_codes = all_codes.difference(existing_codes)
            if missing_codes:
                properties_summary[prop_name] = properties_summary[prop_name].with_added_codes(
                    sorted(missing_codes)
                )

    print("\tgetting metadata complete.")
    return properties_summaries, sample_ids, input_folders


def _get_tar_top_level_folder(tar_in: tarfile.TarFile) -> str:
    root_folders = {
        archive_parents[-2]
        for archive_parents in (
            Path(archive_file).parents for archive_file in tar_in.getnames()
        )
        if len(archive_parents) > 1
    }
    # if there is exactly one base folder in the archive, return that, otherwise return '.'
    return root_folders.pop() if len(root_folders) == 1 else '.'


def extract_tar_to_folder(
        input_tar: Path,
        base_dir: Optional[Path] = None,
        remove_input_tar: bool = Default.remove_input_tar
) -> Path:
    """Extract contents of tarfile to a folder and return the path to the extracted folder. The
    folder name will be set by the root path of the tarred archive.

    Args:
        input_tar: path to input tar file
        base_dir: If specifie, use as the parent dir of extracted tarred folder. If None, then
                  use the parent dir of input_tar
        remove_input_tar: if True, remove the input tar file after extracting.
    """
    input_tar = input_tar.expanduser()
    if input_tar in _extracted_tar_files:
        if remove_input_tar:
            input_tar.unlink(missing_ok=True)
        return _extracted_tar_files[input_tar]
    if base_dir is None:
        base_dir = input_tar.parent
    base_dir.mkdir(parents=True, exist_ok=True)
    with tarfile.open(input_tar, 'r') as tar_in:
        tar_folder = base_dir / _get_tar_top_level_folder(tar_in)
        shutil.rmtree(tar_folder, ignore_errors=True)
        _extracted_tar_files[input_tar] = tar_folder
        tar_in.extractall(path=base_dir)
    if remove_input_tar:
        input_tar.unlink()
    print(f"Extracted tar file to {tar_folder}")
    return tar_folder


def _get_property_summaries(input_folder: Path) -> PropertiesSummary:
    summary_json = next(input_folder.glob("*.json"), None)
    if summary_json is None:
        raise ValueError(f"No summaries .json present in {input_folder}")
    with open(summary_json, 'r') as f_in:
        return json.load(f_in, object_hook=PropertiesSummary.from_json)


def _get_sample_ids(input_folder: Path) -> list[str]:
    samples_list_file = next(input_folder.glob("*.list"), None)
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
            num_rows = prop_summary.num_rows
            first_row_prop = prop_name
        elif prop_summary.num_rows != num_rows:
            raise ValueError(f"Inconsistent {Keys.num_rows}: {first_row_prop}={num_rows}, "
                             f"{prop_name}={prop_summary.num_rows}")
        num_prop_samples = prop_summary.num_samples
        if num_prop_samples > 0 and num_prop_samples != num_samples:
            raise ValueError(
                f"Property {prop_name} has {num_prop_samples} columns but there are {num_samples}"
                " samples"
            )
        codes = prop_summary.codes
        if codes:
            if len(set(codes)) != len(codes):
                raise ValueError(f"Property {prop_name} has non-unique codes {codes}")

    if previous_properties_summaries:
        # ensure that property names, columns, types, and codes are consistent
        compare_summary = previous_properties_summaries[0]
        if set(compare_summary.keys()) != set(properties_summary.keys()):
            raise ValueError("Property names are not consistent")
        for prop_name, prop_summary in compare_summary.items():
            num_prop_samples = properties_summary[prop_name].num_samples
            num_compare_samples = prop_summary.num_samples
            if num_prop_samples != num_compare_samples:
                raise ValueError(
                    f"Inconsistent {Keys.num_samples} between data sets for {prop_name}: "
                    f"{num_compare_samples} -> {num_prop_samples}"
                )
            type_name = properties_summary[prop_name].property_type
            compare_type = prop_summary.property_type
            if type_name != compare_type:
                raise ValueError(
                    f"Inconsistent {Keys.property_type} between data sets for {prop_name}: "
                    f"{type_name.name} -> {compare_type.name}"
                )


def load_properties_from_map(
        input_folder: Path | list[Path],
        properties_summary: PropertiesSummary | list[PropertiesSummary],
        sample_ids: list[str],
        wanted_properties: Optional[set[str]] = None,
        dask_partition_size_mb: float = Default.dask_partition_size_mb,
        category_encoding: str = Default.category_encoding,
        property_names_to_lower: bool = Default.property_names_to_lower
) -> DataFrame:
    if isinstance(input_folder, Path):
        if not isinstance(properties_summary, Mapping):
            raise ValueError(
                "Must specify the same numbers of input_folder and property_summaries"
            )
        # explicitly convert to list
        input_folders = [input_folder]
        properties_summaries = [properties_summary]
    else:
        if not (isinstance(input_folder, list) and isinstance(properties_summary, list) and
                len(input_folder) == len(properties_summary)):
            raise ValueError(
                "Must specify the same numbers of input_folder and property_summaries"
            )
        input_folders = input_folder
        properties_summaries = properties_summary

    if wanted_properties is not None:
        properties_summaries = [
            {
                property_name: property_summary
                for property_name, property_summary in properties_summary.items()
                if property_name in wanted_properties
            }
            for properties_summary in properties_summaries
        ]

    # decompress files, and form structure of the map
    divisions = [0]
    map_parameters = []
    start_row = 0
    for input_folder, properties_summary in zip(input_folders, properties_summaries):
        # unzip TSVs in parallel, wait until all complete
        dask.compute(
            *(
                dask.delayed(_decompress_gzipped_file)(input_folder / f"{property_name}.tsv.gz")
                for property_name in properties_summary.keys()
            )
        )
        input_divisions = _get_dask_divisions(properties_summary, dask_partition_size_mb)
        divisions.extend(division + start_row for division in input_divisions[1:])
        map_parameters.extend(
            (
                input_divisions[_idx - 1], input_divisions[_idx], start_row, input_folder,
                properties_summary
            )
            for _idx in range(1, len(input_divisions))
        )
        start_row = start_row + input_divisions[-1]

    # form dask DataFrame from this map
    properties = dask.dataframe.from_map(
        _tsvs_to_df_partition,
        map_parameters,
        divisions=divisions,
        sample_ids=sample_ids,
        total_num_rows=start_row,
        category_encoding=category_encoding,
        property_names_to_lower=property_names_to_lower
    )
    properties.columns.names = Default.column_levels
    return properties


def _decompress_gzipped_file(
        compressed_file: Path,
        remove_compressed: bool = True
) -> Path:
    if not compressed_file.suffix == ".gz":
        raise ValueError("compressed file does not end with .gz extension")
    decompressed_file = compressed_file.with_suffix("")
    if decompressed_file.is_file() and not compressed_file.is_file():
        # decompression was already performed, or file was never compressed:
        return decompressed_file
    with gzip.open(compressed_file, "rb") as f_in, open(decompressed_file, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
    if remove_compressed:
        compressed_file.unlink()
    return decompressed_file


def _get_dask_divisions(
        property_summaries: Mapping[str, PropertySummary],
        dask_partition_size_mb: float = Default.dask_partition_size_mb
) -> list[int]:
    num_rows = next(iter(property_summaries.values())).num_rows
    bytes_per_row = sum(_get_bytes_per_row(property_name, property_summary)
                        for property_name, property_summary in property_summaries.items())
    # compute the number of rows per division that should fit within our memory budget, but set to
    # a minimum of 5 so that basic operations can work (e.g. .head())
    rows_per_division = max(
        5,
        int(numpy.ceil(dask_partition_size_mb * (2 ** 20) / bytes_per_row))
    )
    return list(range(0, num_rows, rows_per_division)) + [num_rows]


def _get_bytes_per_row(property_name: str, property_summary: PropertySummary) -> int:
    if property_name == "ID":
        return 20  # wild guess as to string length, won't be a dominant factor
    else:
        _dtype = property_summary.encoding_dtype
        bits_per_column = numpy.finfo(_dtype).bits if numpy.issubdtype(_dtype, numpy.floating) \
            else numpy.iinfo(_dtype).bits
        return (bits_per_column // 8) * max(1, property_summary.num_samples)


def _tsvs_to_df_partition(
        division_parameters: tuple[int, int, int, Path, PropertiesSummary],
        sample_ids: list[str],
        total_num_rows: int,
        category_encoding: str,
        property_names_to_lower: bool
) -> pandas.DataFrame:
    row_start, row_end, data_set_row_start, input_folder, properties_summary = division_parameters

    def _summary_sort_key(_item: tuple[str, PropertySummary]) -> tuple[int, str]:
        _name, _summary = _item
        return _summary.num_samples, _name

    # noinspection PyTypeChecker
    df = pandas.concat(
        tuple(
            _read_csv_chunk(
                input_folder=input_folder, property_name=property_name,
                property_summary=property_summary, sample_ids=sample_ids, row_start=row_start,
                row_end=row_end, category_encoding=category_encoding,
                property_names_to_lower=property_names_to_lower
            )
            for property_name, property_summary in sorted(properties_summary.items(),
                                                          key=_summary_sort_key)
        ),
        axis=1
    )
    df.index = pandas.Index(
        numpy.arange(row_start + data_set_row_start, row_end + data_set_row_start),
        dtype=numpy.dtype(numpy.min_scalar_type(-total_num_rows)), name=Keys.row
    )
    return genomics_io.standard_sort_properties(df, sample_ids=sample_ids)


def _read_csv_chunk(
        input_folder: Path,
        property_name: str,
        property_summary: PropertySummary,
        sample_ids: list[str],
        row_start: int,
        row_end: int,
        category_encoding: str,
        property_names_to_lower: bool
) -> pandas.DataFrame:
    tsv = input_folder / f"{property_name}.tsv"
    if property_names_to_lower:
        property_name = property_name.lower()
    is_format_property = property_summary.num_samples > 0
    if is_format_property:
        # this is a format / per-sample property
        column_names = sample_ids
        columns = pandas.MultiIndex.from_product([sample_ids, [property_name]],
                                                 names=Default.column_levels)
    else:
        # this is an INFO / per-variant property
        column_names = [property_name]
        columns = pandas.MultiIndex.from_product([[None], [property_name]],
                                                 names=Default.column_levels)
    # noinspection PyTypeChecker
    _dtype = "string[pyarrow]" if property_name == "id" else \
        property_summary.decoding_dtype
    dtypes = {column_name: _dtype for column_name in column_names}
    df = pandas.read_csv(
        tsv, sep='\t', low_memory=True, engine='c', memory_map=False, names=column_names,
        dtype=dtypes, nrows=row_end - row_start, skiprows=row_start
    )
    if property_summary.codes:
        if category_encoding == Keys.category:
            # need to do convert ordinal ints to pandas Categorical
            df = pandas.DataFrame(
                pandas.Categorical.from_codes(df, categories=property_summary.codes),
                index=df.index, columns=columns
            )
        else:
            # need to one-hot this DataFrame
            codewords_map = _get_codewords_map(property_summary)
            df = codes_to_one_hot(df, codewords_map=codewords_map, property_name=property_name,
                                  is_format_property=is_format_property)
    else:
        df.columns = columns
    return df


def categorical_properties_to_one_hot(
        properties_dataframe: pandas.DataFrame,
        properties_summary: PropertiesSummary,
        property_names_to_lower: bool = Default.property_names_to_lower,
        do_not_1_hot_properties: set[str] = frozenset(
            {Keys.variant_weights, Keys.id, Keys.contig}
        ),
        excluded_1_hot_properties: set[str] = frozenset({})
) -> pandas.DataFrame:
    """
    Convert a properties DataFrame with MultiIndex columns, where encoding is categorical to one
    where encoding is 1-hot.
    Utility function mainly useful for testing, but here in case it's useful for debugging
    """
    property_names = set(properties_dataframe.columns.get_level_values(level=Keys.property))

    def _iter_properties() -> tuple[str, pandas.DataFrame, PropertySummary]:
        for property_name in property_names:
            property_values = properties_dataframe.loc[:, (slice(None), property_name)]
            property_summary = None if property_name in do_not_1_hot_properties \
                else properties_summary[property_name]
            if property_summary is not None and property_summary.codes:
                inverse_codes = {code_word: ind for ind, code_word in
                                 enumerate(property_summary.codes)}
                # expect that property_values is probably CategoricalDtype. If so, need to convert
                # to integer dtypes corresponding to the appropriate codes
                property_values = property_values.applymap(inverse_codes.get).astype(
                    property_summary.encoding_dtype
                )
                property_values.columns = (
                    [property_name.lower()] if property_summary.num_samples == 0
                    else property_values.columns.get_level_values(level=Keys.sample_id)
                )

            yield (
                property_name.lower() if property_names_to_lower else property_name,
                property_values,
                property_summary
            )

    return genomics_io.standard_sort_properties(
        pandas.concat(
            (
                codes_to_one_hot(
                    property_values,
                    codewords_map=_get_codewords_map(property_summary),
                    property_name=property_name,
                    is_format_property=property_summary.num_samples > 0,
                    excluded_properties=excluded_1_hot_properties
                ) if property_summary is not None and property_summary.codes else property_values
                for property_name, property_values, property_summary in _iter_properties()
            ),
            axis=1
        )
    )


def _get_codewords_map(property_summary: PropertySummary) -> dict[str, numpy.array]:
    """
    Create a map from value of String/StringSet property to 1-hot encoding for that value, for
    fast lookup
    Args:
        property_summary: PropertySummary
            Summary info for this property
    Returns:
        codewords_map: dict[str, numpy.array]
            Map from code to 1-hot encoding for that code
    """
    if property_summary.property_type == PropertyType.String:
        # this property is a single string
        # easy: the nth code is all false except for the nth value in the map
        num_codewords = len(property_summary.codes)

        def _basis(_ind: int) -> numpy.ndarray:
            # add last code as always False for "missing" code (-1), corresponding to unwanted
            # properties
            _arr = numpy.zeros(num_codewords + 1, dtype=bool)
            _arr[_ind] = True
            return _arr
        return {code: _basis(ind) for ind, code in enumerate(property_summary.codes)}
    else:
        # this property is a set of strings
        code_sets = [set(code.split(',')) for code in property_summary.codes]

        all_codewords = sorted(set.union(*code_sets))

        def _contains(_codeword: str) -> numpy.ndarray:
            # add last code as always False for "missing" code (-1), corresponding to unwanted
            # properties
            return numpy.array([_codeword in code_set for code_set in code_sets] + [False],
                               dtype=bool)
        return {codeword: _contains(codeword) for codeword in all_codewords}


def codeword_to_one_hot_property_name(codeword: str, property_name: str) -> str:
    return (
        codeword if codeword.lower().startswith(f"{property_name}=".lower())
        else f"{property_name}={codeword}"
    )


def codes_to_one_hot(
        df: DataFrame,
        codewords_map: dict[str, numpy.array],
        property_name: str,
        is_format_property: bool,
        excluded_properties: set[str] = frozenset({})
) -> DataFrame:
    """
    Convert DataFrame of codes (words for String property, comma-separated words for StringSet
    property) to 1-hot encoded DataFrame
    Args:
        df: pandas.DataFrame
            DataFrame with values as string codes
        codewords_map: dict[str, numpy.array]
            Map from code to boolean numpy array of the 1-hot encoding corresponding to that code
        property_name: str
            Name of the property
        is_format_property: bool
            True if a per-sample / Format property, False otherwise
        excluded_properties: set[str]
            List of any properties to omit from final DataFrame
    """
    # make sensible but non-redundant property names
    one_hot_property_names = [
        codeword_to_one_hot_property_name(codeword=codeword, property_name=property_name)
        for codeword in codewords_map.keys()
    ]
    try:
        # make one dataframe for each code, indicating if that code is present for that row
        mapped_dfs = {
            one_hot_property_name: df.applymap(code_is_included.take)
            for one_hot_property_name, (codeword, code_is_included) in zip(
                one_hot_property_names, codewords_map.items()
            ) if one_hot_property_name not in excluded_properties
        }
    except IndexError as index_error:
        common.add_exception_context(
            index_error,
            f"codewords_map: {codewords_map}, df_codes: {numpy.unique(df.values[:])}"
        )
        raise

    for one_hot_property_name, mapped_df in mapped_dfs.items():
        mapped_df.columns = pandas.MultiIndex.from_product(
            [df.columns.tolist(), [one_hot_property_name]] if is_format_property
            else [[None], [one_hot_property_name]],
            names=Default.column_levels
        )

    return pandas.concat(
        [mapped_df for __, mapped_df in sorted(mapped_dfs.items(), key=lambda _item: _item[0])],
        axis=1
    )


def df_to_parquet(
        df: DataFrame,
        output_path: Path,
        compression_algorithm: str = Default.compression_algorithm,
        write_metadata: Optional[bool] = None,
        write_index: bool = True
):
    CompressionAlgorithms.assert_valid(compression_algorithm)
    # flatten columns to be compatible with parquet
    old_columns = df.columns
    dask_utils.flatten_columns(df)
    output_path = output_path.expanduser()
    if output_path.suffix == ".tar":
        # request final output as a tarred parquet folder
        tar_file = output_path
        # remove .tar from output folder
        parquet_folder = output_path.with_suffix("")
    else:
        # don't produce a tar_file
        tar_file = None
        parquet_folder = output_path
    with warnings.catch_warnings():
        warnings.filterwarnings(action="ignore", category=FutureWarning)
        print(f"writing parquet files to {parquet_folder}")
        df.to_parquet(parquet_folder, engine="pyarrow", compression=compression_algorithm,
                      write_metadata_file=write_metadata, write_index=write_index)
        os.system(f"ls -alh {parquet_folder}")

    # restore the old columns
    df.columns = old_columns

    if tar_file is not None:
        archive_to_tar(folder_to_archive=parquet_folder, tar_file=tar_file, remove_files=True)


def archive_to_tar(
        folder_to_archive: Path,
        tar_file: Optional[Path] = None,
        remove_files: bool = False
) -> Path:
    folder_to_archive = folder_to_archive.expanduser()
    if not folder_to_archive.is_dir():
        raise ValueError(f"{folder_to_archive} is not a directory")
    tar_file = folder_to_archive.with_suffix(".tar") if tar_file is None else tar_file.expanduser()
    tar_root_folder = folder_to_archive.parent
    with tarfile.open(tar_file, "w") as tar_out:
        for file_to_archive in folder_to_archive.rglob("*"):
            if file_to_archive.is_dir():
                continue  # only need to archive files
            tar_out.add(file_to_archive, arcname=os.path.relpath(file_to_archive, tar_root_folder))
            if remove_files:
                file_to_archive.unlink()
    if remove_files:
        shutil.rmtree(folder_to_archive)
    return tar_file


def parquet_to_df(
        input_path: Path | list[Path],
        base_dir: Optional[Path] = None,
        remove_input_tar: bool = Default.remove_input_tar,
        index_column: Optional[str] = Keys.row,
        wanted_properties: Optional[list[str]] = None,
        wanted_samples: Optional[list[str]] = None,
        error_on_missing_property: bool = Default.error_on_missing_property,
        error_on_missing_sample: bool = Default.error_on_missing_sample,
        filters: Union[list[tuple[str, str, Any]], list[list[tuple[str, str, Any]]]] = None,
        engine: str = "pyarrow",
        unflatten_columns: bool = True,
) -> DataFrame:
    """Load parquet folder into dask DataFrame with multi-index columns. If the parquet folder is
    tarred, manage untarring.

    Args:
        input_path: path to parquet folder (or tarfile of parquet folder). Can pass a list of
                    such parquet folders, which will then be combined into a single dask dataframe.
        base_dir: Optional path to extract tarred parquet folder into
        remove_input_tar: if True and the input_path is a parquet folder, remove tarfile after
                          extracting
        index_column: which column to use for index if available
        wanted_properties: list of properties to extract. If None, take all
        wanted_samples: list of samples to take genotype data for. If None, take all
        error_on_missing_property: if True and a requested property is not present, throw error. If
                                 False, continue anyway
        error_on_missing_sample: if True and a requested sample is not present, throw error. If
                                 False, continue anyway
        filters: Row filters to specify which variants to take. See documentation for
                 pyarrow.parquet.read_table
                 https://arrow.apache.org/docs/python/generated/pyarrow.parquet.read_table.html
        engine: which engine to use. In general this code expects pyarrow.
        unflatten_columns: if True, return columns in multi-index form
    """
    if isinstance(input_path, Path):
        input_path = [input_path]

    input_path = [
        extract_tar_to_folder(
            parquet_path.expanduser(), base_dir=base_dir, remove_input_tar=remove_input_tar
        ) if parquet_path.suffix == ".tar" else parquet_path
        for parquet_path in input_path
    ]

    if index_column is None:
        index = None
    elif index_column in get_parquet_folder_columns(input_path[0]):
        index = index_column
    else:
        index = None
        warnings.warn(
            f"{index_column} is not a column in {input_path[0]}, "
            f"available columns: {get_parquet_folder_columns(input_path[0])}"
        )

    # convert input argument to form dask can definitely read: input files listed in order if there
    # is more than one input path, otherwise just the folder if there is only one.
    input_file_arg = [
        f"{input_file}"
        for parquet_path in input_path
        for input_file in get_parquet_files(parquet_path)
    ] if len(input_path) > 1 else f"{input_path[0]}"

    with warnings.catch_warnings():
        df = dask.dataframe.read_parquet(
            input_file_arg, engine=engine, index=index, calculate_divisions=True,
            columns=_get_wanted_columns(input_path[0], wanted_properties, wanted_samples,
                                        error_on_missing_property=error_on_missing_property,
                                        error_on_missing_sample=error_on_missing_sample),
            filters=filters
        )

    return dask_utils.unflatten_columns(df) if unflatten_columns else df


def _get_wanted_columns(
        input_path: Path,
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
            file_columns = get_parquet_folder_columns(input_path)
        wanted_columns = [
            col for col in file_columns
            if any(f"'{prop.lower()}'" in col.lower() for prop in wanted_properties)
        ]
        if error_on_missing_property:
            for prop in wanted_properties:
                if not any(f"'{prop}'" in col for col in wanted_columns):
                    raise ValueError(f"Wanted property '{prop}' not present in {input_path}")
    elif wanted_properties is None:
        # get columns with desired samples (or no sample / per-variant data)
        if file_columns is None:
            file_columns = get_parquet_folder_columns(input_path)
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
        file_columns = get_parquet_folder_columns(input_path)
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


def get_arbritrary_parquet_file(
        input_path: Path,
        base_dir: Optional[Path] = None,
        remove_input_tar: bool = Default.remove_input_tar
) -> Path:
    if input_path.suffix == ".tar":
        input_path = extract_tar_to_folder(input_path, base_dir=base_dir,
                                           remove_input_tar=remove_input_tar)
    return next(input_path.glob("*.parquet"))


def get_parquet_file_schema(parquet_file: Path) -> pyarrow.Schema:
    try:
        return pyarrow.parquet.read_schema(parquet_file)
    except Exception as exception:
        common.add_exception_context(exception, f"Reading {parquet_file}")
        raise


def get_parquet_folder_sample_ids(
        input_path: Path,
        base_dir: Optional[Path] = None,
        remove_input_tar: bool = Default.remove_input_tar
) -> list[str]:
    return get_parquet_file_sample_ids(
        get_arbritrary_parquet_file(
            input_path=input_path, base_dir=base_dir, remove_input_tar=remove_input_tar
        )
    )


def get_parquet_folder_columns(
        input_path: Path,
        base_dir: Optional[Path] = None,
        remove_input_tar: bool = Default.remove_input_tar,
) -> list[str]:
    parquet_file = get_arbritrary_parquet_file(
        input_path=input_path, base_dir=base_dir, remove_input_tar=remove_input_tar
    )
    return get_parquet_file_columns(parquet_file)


def get_parquet_file_columns(parquet_file: Path) -> list[str]:
    return get_parquet_file_schema(parquet_file).names


def get_parquet_folder_unflat_columns(
    input_path: Path,
    base_dir: Optional[Path] = None,
    remove_input_tar: bool = Default.remove_input_tar,
) -> list[tuple[Optional[str], str]]:
    parquet_file = get_arbritrary_parquet_file(
        input_path=input_path, base_dir=base_dir, remove_input_tar=remove_input_tar
    )
    return get_parquet_file_unflat_columns(parquet_file)


def get_parquet_file_unflat_columns(
    parquet_file: Path
) -> list[tuple[Optional[str], str]]:
    return [
        dask_utils.unflatten_column_name(flat_column_name)
        for flat_column_name in get_parquet_file_columns(parquet_file)
        if flat_column_name != Keys.row
    ]


def get_parquet_file_sample_ids(parquet_file: Path) -> list[str]:
    seen_samples = set()

    def _first_time_seen(_sample_id: str) -> bool:
        if _sample_id in seen_samples:
            return False
        else:
            seen_samples.add(_sample_id)
            return True

    return [
        sample_id for (sample_id, property_name) in get_parquet_file_unflat_columns(parquet_file)
        if sample_id is not None and _first_time_seen(sample_id)
    ]


def get_parquet_folder_num_rows(
        input_path: Path,
        base_dir: Optional[Path] = None,
        remove_input_tar: bool = Default.remove_input_tar
) -> int:
    if input_path.suffix == ".tar":
        input_path = extract_tar_to_folder(input_path, base_dir=base_dir,
                                           remove_input_tar=remove_input_tar)
    return get_parquet_files_num_rows(parquet_files=input_path.glob("*.parquet"))


def get_parquet_files(parquet_folder: Path) -> list[Path]:
    return sorted(
        parquet_folder.glob(f"*{Keys.parquet_suffix}"),
        key=lambda p: int(p.name.split('.')[1])
    )


def get_parquet_files_num_rows0(parquet_files: Iterable[Path]) -> int:
    return dask.compute(
        sum(
            dask.delayed(get_parquet_file_num_rows)(parquet_file)
            for parquet_file in parquet_files
        )
    )[0]


def get_parquet_files_num_rows(parquet_files: Iterable[Path]) -> int:
    return sum(
        get_parquet_file_num_rows(parquet_file)
        for parquet_file in parquet_files
    )


def get_parquet_file_num_rows(parquet_file: Path) -> int:
    return pyarrow.parquet.ParquetFile(parquet_file).metadata.num_rows


def get_parquet_file_dtypes(
        parquet_file: Path,
        properties_summary: PropertiesSummary
) -> dict[str, DType] | dict[tuple[Optional[str], str], DType]:
    schema = get_parquet_file_schema(parquet_file)

    def _to_pandas_dtype(
            _name: tuple[Optional[str], str], _pyarrow_dtype: pyarrow.DataType
    ) -> DType:
        try:
            return _pyarrow_dtype.to_pandas_dtype()
        except NotImplementedError:  # hand-override for category dtypes
            return pandas.CategoricalDtype(properties_summary[_name[1]].codes)

    unflattened_names = [name if name == Keys.row else dask_utils.unflatten_column_name(name)
                         for name in schema.names]

    return {
        name: _to_pandas_dtype(name, pyarrow_dtype)
        for name, pyarrow_dtype in zip(unflattened_names, schema.types)
        if name != Keys.row
    }


def get_properties_scaling(properties: DataFrame) -> PropertiesScaling:
    def _get_bool_counts_dict(_partition: pandas.DataFrame) -> dict[str, int]:
        n_elements = _partition.size
        n_true = _partition.values.ravel().sum()
        return {
            Keys.positive_value: n_true,
            Keys.negative_value: n_elements - n_true
        }

    def _get_value_counts_dict(_partition: pandas.DataFrame) -> dict[int, int]:
        return pandas.Series(_partition.values.ravel()).value_counts().sort_index().to_dict()

    def _get_percentiles_dict(_partition: pandas.DataFrame) -> dict[str, list[float]]:
        quantiles = numpy.nanquantile(_partition.values.ravel(), [0.1, 0.5, 0.9])
        return {"low": [quantiles[0]], Keys.baseline: [quantiles[1]], "high": [quantiles[2]]}

    def _get_property_stats_dict(
            _property_partition: pandas.DataFrame,
            property_name: str
    ) -> PropertyStatsDict:
        dtypes = set(_property_partition.dtypes.values)
        if dtypes == {numpy.dtype(bool)}:
            return _get_bool_counts_dict(_property_partition)
        elif any(
            pandas.api.types.is_float_dtype(_dtype)
            for _dtype in set(_property_partition.dtypes.values)
        ):
            return _get_percentiles_dict(_property_partition)
        else:
            try:
                return _get_value_counts_dict(_property_partition)
            except TypeError as type_error:
                common.add_exception_context(type_error, f"getting stats for {property_name}")
                raise

    def _get_all_stats_dict(_partition: pandas.DataFrame) -> AllPropertiesStatsDict:
        return {
            property_name: _get_property_stats_dict(
                _partition.loc[:, (slice(None), property_name)], property_name
            )
            for property_name in set(_partition.columns.get_level_values(Keys.property))
            if property_name not in non_ml_properties
        }

    def _reduce_all_stats_dict(
            _all_properties_stats_dicts: list[AllPropertiesStatsDict]
    ) -> AllPropertiesStatsDict:
        _reduced_properties_dict = _all_properties_stats_dicts[0]
        for _all_properties_stats_dict in _all_properties_stats_dicts[1:]:
            for _property_name, _stats_dict in _all_properties_stats_dict.items():
                if Keys.baseline in _stats_dict:
                    # float property, concatenate values (by adding lists)
                    _reduced_dict = _reduced_properties_dict[_property_name]
                    for _value, _counts_list in _stats_dict.items():
                        _reduced_dict[_value] += _counts_list
                else:
                    # property with value counts: add values
                    _reduced_dict = _reduced_properties_dict[_property_name]
                    for _value, _counts_list in _stats_dict.items():
                        _reduced_dict[_value] = _reduced_dict.get(_value, 0) + _counts_list
        return _reduced_properties_dict

    stats_properties = [
        property_name for property_name in set(properties.columns.get_level_values(Keys.property))
        if property_name not in non_ml_properties
    ]
    print(f"Getting stats for {stats_properties}")

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        all_properties_stats = properties.reduction(
            _get_all_stats_dict, _reduce_all_stats_dict, meta={}
        ).compute()[0]

    # now need to do final computation of baseline and scale

    def _iter_property_scaling(
        _all_properties_stats: dict[str, Any]
    ) -> Iterator[tuple[str, PropertyScaling]]:
        for _property_name, _property_stats in _all_properties_stats.items():
            if Keys.positive_value in _property_stats:
                # boolean property
                if (
                        _property_stats[Keys.positive_value] == 0 or
                        _property_stats[Keys.negative_value] == 0
                ):
                    # there's no real data in this property
                    positive_value = 0.0
                    negative_value = 0.0
                    mean_value = 0.0
                    scale = 1.0
                else:
                    mean_value = (
                            _property_stats[Keys.positive_value] / (
                                _property_stats[Keys.positive_value] +
                                _property_stats[Keys.negative_value]
                            )
                    )
                    scale = (mean_value * (1.0 - mean_value)) ** 0.5
                    positive_value = (1.0 - mean_value) / scale
                    negative_value = -mean_value / scale
                yield _property_name, PropertyScaling(
                    baseline=mean_value, scale=scale,
                    positive_value=positive_value, negative_value=negative_value
                )
            elif Keys.baseline in _property_stats:
                # float property
                # noinspection PyTypeChecker
                baseline: float = numpy.median(_property_stats[Keys.baseline])
                scale = numpy.median(_property_stats["high"]) \
                    - numpy.median(_property_stats["low"])
                yield _property_name, PropertyScaling(
                    baseline=baseline, scale=1.0 if scale == 0.0 else scale
                )
            elif all(isinstance(key, int) for key in _property_stats.keys()):
                # int property
                num_counts = sum(_property_stats.values())

                def _get_value_at_quantile(_quantile: float) -> int:
                    _target_counts = _quantile * num_counts
                    for _value, _counts in sorted(_property_stats.items(), key=lambda t: t[0]):
                        _target_counts -= _counts
                        if _target_counts <= 0:
                            return _value

                baseline = _get_value_at_quantile(0.5)  # baseline is median
                # scale is the range of the central 80% of data
                scale = _get_value_at_quantile(0.9) - _get_value_at_quantile(0.1)
                if scale == 0:
                    # if the data is mostly constant, set scale to 1.0
                    scale = 1.0

                yield _property_name, PropertyScaling(baseline=baseline, scale=scale)
            else:
                # object / categorical properties: make stats for 1-hot encoding
                num_entries = sum(_property_stats.values())
                sub_category_counts_dict = {
                    sub_category: 0 for category in _property_stats.keys()
                    for sub_category in category.split(',')
                }
                for category, category_counts in _property_stats.items():
                    for sub_category in category.split(','):
                        sub_category_counts_dict[sub_category] += category_counts
                one_hot_stats = {
                    codeword_to_one_hot_property_name(sub_category, property_name=_property_name):
                        {Keys.positive_value: sub_category_counts,
                         Keys.negative_value: num_entries - sub_category_counts}
                    for sub_category, sub_category_counts in sub_category_counts_dict.items()
                }
                yield from _iter_property_scaling(one_hot_stats)

    return PropertiesScaling({
        property_name: property_scaling
        for property_name, property_scaling in _iter_property_scaling(all_properties_stats)
    })


def add_variant_weights(properties: DataFrame, temp_dir: Path = Default.temp_dir) -> DataFrame:
    properties = dask.dataframe.concat(
        [get_variant_weights(properties, temp_dir=temp_dir), properties],
        axis=1, interleave_partitions=True, ignore_order=True
    )
    properties.columns.names = Default.column_levels
    properties.index = properties.index.rename(Keys.row)
    return properties


def get_variant_weights(
        properties: DataFrame, temp_dir: Path = Default.temp_dir
) -> dask.dataframe.Series:
    variant_category_columns = {
        prop for prop in properties.columns.get_level_values(level="property")
        if prop in training_weight_bins or prop == Keys.svtype
    }

    def _get_categories_dataframe(_df: DataFrame) -> DataFrame:
        _raw_categories_df = _df[[(None, prop) for prop in variant_category_columns]]
        _raw_categories_df.columns = _raw_categories_df.columns.droplevel(Keys.sample_id)
        dtypes = _raw_categories_df.dtypes
        for prop_name in training_weight_bins.keys():
            dtypes[prop_name] = numpy.int64
        _categories_df = _raw_categories_df.map_partitions(
            _bin_categories,
            meta=pandas.DataFrame(
                {
                    _name: pandas.Series(dtype=_dtype, name=_name)
                    for _name, _dtype in dtypes.items()
                }
            )
        )
        _categories_df.index = _df.index
        return _categories_df

    def _bin_categories(_partition: pandas.DataFrame):
        for prop_name, prop_bins in training_weight_bins.items():
            _partition[prop_name] = numpy.searchsorted(
                prop_bins, _partition.loc[:, prop_name], side="right"
            )
        return _partition.copy()

    def _get_category_counts(_partition: pandas.DataFrame) -> CategoryCounts:
        # explicitly convert to dict so that it can be edited (otherwise get Mapping)
        return dict(_partition.value_counts(sort=False).to_dict())

    def _reduce_category_counts(_category_counts_list: list[CategoryCounts]) -> CategoryCounts:
        _reduced_counts = _category_counts_list[0]
        for _category_counts in _category_counts_list[1:]:
            for _category, _counts in _category_counts.items():
                _reduced_counts[_category] = _reduced_counts.get(_category, 0) + _counts
        # sort in category order
        return {
            _category: _counts
            for _category, _counts in sorted(_reduced_counts.items(), key=lambda t: t[0])
        }

    # get dataframe of binned categories
    categories_dataframe = _get_categories_dataframe(properties)

    variant_category_counts = categories_dataframe.reduction(
        _get_category_counts, _reduce_category_counts, meta={}
    ).compute()[0]

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        variant_category_weights = get_variant_category_weights(variant_category_counts)
    # we want to compress weights as a categorical Series, but it can't handle non-unique values,
    # so we will need to merge the weights that are equal

    unique_weights = []
    variant_bin_map = dict()
    next_bin = 0
    for bin_number, category in enumerate(sorted(variant_category_weights.keys())):
        weight = variant_category_weights[category]
        for _first_bin, _first_category in enumerate(sorted(variant_category_weights.keys())):
            if _first_bin >= bin_number:
                unique_weights.append(weight)
                variant_bin_map[category] = next_bin
                next_bin += 1
                break
            elif variant_category_weights[_first_category] == weight:
                variant_bin_map[category] = variant_bin_map[_first_category]
                break

    _dtype = pandas.CategoricalDtype(unique_weights)

    def _get_partition_variant_weights_categorical(_partition: pandas.DataFrame) -> pandas.Series:
        return pandas.Series(
            pandas.Categorical.from_codes(
                [
                    variant_bin_map[tuple(_category)]
                    for _category in _partition.itertuples(index=False)
                ],
                categories=unique_weights
            ),
            index=_partition.index, name=Keys.variant_weights, dtype=_dtype
        )

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        variant_weights = categories_dataframe.map_partitions(
            _get_partition_variant_weights_categorical,
            meta=pandas.Series(
                [], dtype=pandas.CategoricalDtype(categories=unique_weights),
                name=Keys.variant_weights,
                index=pandas.Index([], name=Keys.row)
            )
        )

    variant_weights.name = (None, Keys.variant_weights)
    variant_weights.index = variant_weights.index.rename(Keys.row)
    return variant_weights


Category = tuple[bool | int, ...]
CategoryCounts = dict[Category, int]
CategoryWeights = dict[Category, float]


def get_variant_category_weights(
        variant_category_counts: CategoryCounts,
        max_weight: float = 100.0
) -> CategoryWeights:
    """
    Given counts of number of variants in each variant category, compute weights such that:
     - weight is inversely proportional to the counts in each category
     - the mean weight (over variants) is 1.0
    After computing initial weights, restrain values so that no weight is greater than max_weight,
    or less than 1.0 / max_weight

    Args:
        variant_category_counts:
        max_weight:
    Returns:
    sum_c N_c * w_c = N_v
    w_c = A / N_c
    sum_c A = N_v
    """
    num_categories = len(variant_category_counts)
    num_variants = sum(variant_category_counts.values())
    min_weight = 1.0 / max_weight
    return {
        category: max(min(num_variants / float(num_categories * counts), max_weight), min_weight)
        for category, counts in variant_category_counts.items()
    }


if __name__ == "__main__":
    main()
