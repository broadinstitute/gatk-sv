#!/usr/bin/env python
import logging
import sys
import tarfile
from pathlib import Path
import warnings
import argparse
import shutil
import tempfile
import json
import numpy
import pandas
import attr
from tqdm.auto import tqdm as tqdm
import dask
from dask.distributed import Client, LocalCluster
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import dask.dataframe
from sv_utils import common, genomics_io
from gq_recalibrator import tarred_properties_to_parquet
from gq_recalibrator.tarred_properties_to_parquet import (
    PropertiesSummary, PropertySummary, PropertiesScaling, PropertyScaling, CompressionAlgorithms
)
from gq_recalibrator import dask_utils
from typing import Optional

DaskDataFrame = dask.dataframe.DataFrame
DaskSeries = dask.dataframe.Series


class Keys:
    row = genomics_io.Keys.row
    svtype = genomics_io.Keys.svtype
    variant_weights = tarred_properties_to_parquet.Keys.variant_weights
    parquet_suffix = tarred_properties_to_parquet.Keys.parquet_suffix
    properties_scaling_json = tarred_properties_to_parquet.Keys.properties_scaling_json
    properties_summary_json = tarred_properties_to_parquet.Keys.properties_summary_json


class Default:
    temp_dir = Path(tempfile.gettempdir())
    compute_weights = True
    remove_input_parquets = False
    compression_algorithm = CompressionAlgorithms.snappy
    num_jobs = common.num_physical_cpus


def __parse_arguments(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Concatenate tarred parquet folders from different shards into a single valid"
                    "parquet folder",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        prog=argv[0]
    )
    parser.add_argument("--input-parquet", "-i", type=str,
                        help="comma-separated path to parquet folder/tarfile to convert",
                        action="extend", required=True, nargs='+')
    parser.add_argument(
        "--output-parquet", "-o", type=Path, required=True,
        help="path to output SV data in parquet format. If path ends in .tar, the parquet folder "
             "will be archived into a single tar file")
    parser.add_argument("--temp-dir", "-t", type=Path, default=Default.temp_dir,
                        help="Path to preferred temporary directory")
    parser.add_argument(
        "--compute-weights", type=common.argparse_bool, default=Default.compute_weights,
        help="If true, compute the weight that each variant should have for training."
    )
    parser.add_argument(
        "--remove-input-parquets", type=common.argparse_bool,
        default=Default.remove_input_parquets,
        help="if true, remove each input tar to save space, if false, leave it in place"
    )
    parser.add_argument(
        "--compression-algorithm", type=str, default=Default.compression_algorithm,
        choices=CompressionAlgorithms.list(), help="compression algorithm for parquet data"
    )
    parser.add_argument("--threads", type=int, default=Default.num_jobs,
                        help="Number of parallel processes to use")
    return parser.parse_args(argv[1:] if len(argv) > 1 else ["--help"])


def main(argv: Optional[list[str]] = None):
    args = __parse_arguments(sys.argv if argv is None else argv)
    cluster = LocalCluster(n_workers=args.threads, silence_logs=logging.ERROR)
    # noinspection PyUnusedLocal
    client = Client(cluster)
    with dask.config.set(temporary_directory=args.temp_dir):
        concat_parquet_shards(
            input_parquets=[
                Path(separated_tar) for input_arg in args.input_parquet
                for separated_tar in input_arg.split(',')
            ],
            output_path=args.output_parquet,
            temp_dir=args.temp_dir,
            compute_weights=args.compute_weights,
            remove_input_parquets=args.remove_input_parquets,
            compression_algorithm=args.compression_algorithm
        )


@attr.define(slots=True, weakref_slot=False)
class RunningParquetStats:
    num_parquet_files: int = 0
    num_parquet_rows: int = 0
    shard_properties_summary: list[Optional[PropertiesSummary]] = []
    shard_properties_scaling: list[Optional[PropertiesScaling]] = []

    @property
    def has_properties_scaling(self) -> bool:
        if any(
            properties_scaling is not None for properties_scaling in self.shard_properties_scaling
        ):
            if any(
                properties_scaling is None for properties_scaling in self.shard_properties_scaling
            ):
                raise ValueError("Some PropertiesScalings are present, some are missing")
            return True
        else:
            return False

    @property
    def has_properties_summary(self) -> bool:
        if any(
            properties_summary is not None for properties_summary in self.shard_properties_summary
        ):
            if any(
                properties_summary is None for properties_summary in self.shard_properties_summary
            ):
                raise ValueError("Some PropertiesSummaries are present, some are missing")
            return True
        else:
            return False

    @property
    def combined_properties_summary(self) -> PropertiesSummary:
        return PropertiesSummary.combine(self.shard_properties_summary)

    def update_from_parquet_folder(
        self, input_parquet: Path, temp_dir: Path, remove_input_tar: bool
    ) -> None:
        if input_parquet.suffix == ".tar":
            input_folder = tarred_properties_to_parquet.extract_tar_to_folder(
                input_tar=input_parquet, base_dir=temp_dir, remove_input_tar=remove_input_tar
            )
        else:
            input_folder = input_parquet
        self.num_parquet_rows += tarred_properties_to_parquet.get_parquet_folder_num_rows(
            input_folder
        )
        self.num_parquet_files += len(tarred_properties_to_parquet.get_parquet_files(input_folder))

        properties_summary_file = input_folder / Keys.properties_summary_json
        if properties_summary_file.is_file():
            with open(properties_summary_file, 'r') as f_in:
                self.shard_properties_summary.append(
                    json.load(f_in, object_hook=PropertiesSummary.from_json)
                )

        properties_scaling_file = input_folder / Keys.properties_scaling_json
        if properties_scaling_file.is_file():
            with open(properties_scaling_file, 'r') as f_in:
                self.shard_properties_scaling.append(
                    json.load(f_in, object_hook=PropertiesScaling.from_json)
                )


def concat_parquet_shards(
        input_parquets: list[Path],
        output_path: Path,
        temp_dir: Path = Default.temp_dir,
        compute_weights: bool = Default.compute_weights,
        remove_input_parquets: bool = Default.remove_input_parquets,
        compression_algorithm: str = Default.compression_algorithm,
):
    f"""Concatenate input tarred parquet folders into a single output parquet folder.
    If the input parquet folders contain {Keys.properties_scaling_json} a combined
        {Keys.properties_scaling_json} will be calculated.
    If the input parquet folders contain {Keys.properties_summary_json} a combined
        {Keys.properties_summary_json} will be formed. 
    
    Args:
        input_parquets: list of paths to input parquet folders (tarred or not), in dataframe order
        output_path: path to combined parquet dataset. If it ends with ".tar", the parquet folder
                     will be tarred.
        temp_dir: path to use for temp files and folders
        compute_weights: if True, compute the weight that each variant should have for training
        remove_input_parquets: if True, remove the input parquets after concatenating
        compression_algorithm: name of algorithm to use for parquet files.
                               default={Default.compression_algorithm} 
    """
    weights = _compute_combined_weights(
        input_parquets=input_parquets, remove_input_tar=remove_input_parquets
    ) if compute_weights else None
    output_folder = (
        output_path.with_suffix("") if output_path.suffix == ".tar"
        else output_path
    )
    # need a temp output folder because we'll need to read the combined parquet file in and set the
    # index so that dask knows its divisions
    temp_output_folder = Path(tempfile.mkdtemp(suffix=Keys.parquet_suffix, dir=temp_dir))
    shutil.rmtree(output_folder, ignore_errors=True)
    output_folder.mkdir(exist_ok=False, parents=True)
    running_parquet_stats = RunningParquetStats()
    num_total_parquet_files = sum(
        get_num_parquet_files(input_parquet) for input_parquet in input_parquets
    )
    progress_bar = tqdm(desc="transferring parquet files", total=num_total_parquet_files)
    for input_parquet in input_parquets:
        if input_parquet.suffix == ".tar":
            input_tar = input_parquet
            input_folder = tarred_properties_to_parquet.extract_tar_to_folder(
                input_tar=input_tar, base_dir=temp_dir, remove_input_tar=remove_input_parquets
            )
        else:
            input_tar = None
            input_folder = input_parquet

        _transfer_parquet_folder(
            progress_bar=progress_bar,
            input_folder=input_folder,
            output_folder=output_folder,
            running_parquet_stats=running_parquet_stats,
            weights=weights,
            compression_algorithm=compression_algorithm
        )
        if remove_input_parquets or input_tar is not None:
            # remove unarchived folder
            shutil.rmtree(input_folder)

    dask_df = tarred_properties_to_parquet.parquet_to_df(
        input_path=output_folder, index_column=Keys.row
    )
    if not dask_df.known_divisions:
        raise RuntimeError(
            "Dataframe has unknown divisions. This is a bug, because at this point the divisions "
            "should be known."
        )

    if running_parquet_stats.has_properties_summary:
        print("Writing properties summary")
        output_folder.mkdir(parents=True, exist_ok=True)
        with open(output_folder / Keys.properties_summary_json, 'w') as f_out:
            json.dump(
                running_parquet_stats.combined_properties_summary,
                f_out,
                indent="  ",
                default=PropertySummary.to_json
            )

    if running_parquet_stats.has_properties_scaling:
        print("Getting properties scaling")
        output_folder.mkdir(parents=True, exist_ok=True)
        combined_properties_scaling = tarred_properties_to_parquet.get_properties_scaling(
            properties=dask_df
        )
        print("Writing properties scaling")
        with open(output_folder / Keys.properties_scaling_json, 'w') as f_out:
            json.dump(
                combined_properties_scaling,
                f_out,
                indent="  ",
                default=PropertyScaling.to_json
            )
    print("Removing temp folder")
    shutil.rmtree(temp_output_folder)
    print("Archiving to tar")
    if output_path.suffix == ".tar":
        tarred_properties_to_parquet.archive_to_tar(
            folder_to_archive=output_folder, tar_file=output_path, remove_files=True
        )


def _compute_combined_weights(input_parquets: list[Path], remove_input_tar: bool) -> DaskSeries:
    """Compute variant weights for final combined data set by concatenating the shards in a dask
    DataFrame, taking the necessary properties to compute weights, and then calling the weight
    calculation function
    """
    print("BEGIN Computing weights")
    wanted_properties = set(
        tarred_properties_to_parquet.training_weight_bins.keys()
    ).update(Keys.svtype)
    properties = tarred_properties_to_parquet.parquet_to_df(
        input_path=input_parquets,
        remove_input_tar=remove_input_tar,
        wanted_properties=wanted_properties,
        error_on_missing_property=True
    )
    # had .persist() here, but I don't think it helps:
    weights = tarred_properties_to_parquet.get_variant_weights(properties=properties)
    print("FINISHED Computing weights")
    return weights


def get_num_parquet_files(input_parquet: Path) -> int:
    if input_parquet.suffix == ".tar":
        with tarfile.open(input_parquet, "r") as tar_in:
            return sum(1 for member in tar_in.getmembers()
                       if member.name.endswith(Keys.parquet_suffix))
    else:
        return sum(1 for __ in input_parquet.glob(f"*{Keys.parquet_suffix}"))


def _transfer_parquet_folder(
        progress_bar: tqdm,
        input_folder: Path,
        output_folder: Path,
        running_parquet_stats: RunningParquetStats,
        weights: Optional[DaskSeries],
        compression_algorithm: str = Default.compression_algorithm
) -> None:
    """Orchestrate file transfer"""
    parquet_files = sorted(
        input_folder.glob(f"*{Keys.parquet_suffix}"),
        key=lambda p: int(p.name.split('.')[1])
    )

    for parquet_file in parquet_files:
        new_parquet_file = \
            output_folder / f"part.{running_parquet_stats.num_parquet_files}{Keys.parquet_suffix}"
        num_file_rows = tarred_properties_to_parquet.get_parquet_file_num_rows(parquet_file)
        row_end = running_parquet_stats.num_parquet_rows + num_file_rows
        _transfer_parquet_file(
            old_parquet_file=parquet_file,
            new_parquet_file=new_parquet_file,
            row_start=running_parquet_stats.num_parquet_rows,
            row_end=row_end,
            weights=(
                None if weights is None
                else weights.partitions[running_parquet_stats.num_parquet_files].compute()
            ),
            compression_algorithm=compression_algorithm
        )
        running_parquet_stats.num_parquet_files += 1
        running_parquet_stats.num_parquet_rows += num_file_rows

        progress_bar.update(n=1)

    properties_summary_file = input_folder / Keys.properties_summary_json
    if properties_summary_file.is_file():
        with open(properties_summary_file, 'r') as f_in:
            running_parquet_stats.shard_properties_summary.append(
                json.load(f_in, object_hook=PropertiesSummary.from_json)
            )

    properties_scaling_file = input_folder / Keys.properties_scaling_json
    if properties_scaling_file.is_file():
        with open(properties_scaling_file, 'r') as f_in:
            running_parquet_stats.shard_properties_scaling.append(
                json.load(f_in, object_hook=PropertiesScaling.from_json)
            )


def _transfer_parquet_file(
        old_parquet_file: Path,
        new_parquet_file: Path,
        row_start: int,
        row_end: int,
        weights: Optional[pandas.Series],
        compression_algorithm: str
) -> None:
    """Move parquet file by loading, adjusting starting row, adding weights, and re-saving"""
    # load old parquet file
    df = pandas.read_parquet(path=old_parquet_file, engine="pyarrow")

    # update row if needed
    save_index = False
    if Keys.row in df.columns:
        # need to update row
        df[Keys.row] = numpy.arange(row_start, row_end, dtype=numpy.int64)
    elif df.index.name == Keys.row:
        # update index, it's the row
        df.index = pandas.Index(numpy.arange(row_start, row_end, dtype=numpy.int64), name=Keys.row)
        save_index = True

    # add weights if needed
    if weights is not None:
        df = dask_utils.unflatten_columns(df)
        # noinspection PyTypeChecker
        col = pandas.MultiIndex.from_tuples(
            [(None, Keys.variant_weights)],
            names=tarred_properties_to_parquet.Default.column_levels
        )[0]
        df[col] = weights.values
        df = dask_utils.flatten_columns(
            genomics_io.standard_sort_properties(properties_dataframe=df)
        )

    # save new parquet file
    df.to_parquet(
        path=new_parquet_file,
        engine="pyarrow",
        compression=compression_algorithm,
        index=save_index
    )


if __name__ == "__main__":
    main()
