import os
from abc import ABC
from pathlib import Path
import atexit
import shutil
import tempfile
import contextlib
import time
import json
import concurrent.futures
import numpy
import pandas
import pandas.core.arrays
import torch
from typing import Optional, Any, IO, TypeVar
from collections.abc import Iterator, Collection
from gq_recalibrator import tarred_properties_to_parquet, training_utils
from sv_utils import common, get_truth_overlap, benchmark_variant_filter

ArrowStringArray = pandas.core.arrays.string_arrow.ArrowStringArray
ConfidentVariants = get_truth_overlap.ConfidentVariants
_debug_file_locks = False


class Keys:
    row = tarred_properties_to_parquet.Keys.row
    id = tarred_properties_to_parquet.Keys.id
    variant_weights = tarred_properties_to_parquet.Keys.variant_weights
    positive_value = tarred_properties_to_parquet.Keys.positive_value
    negative_value = tarred_properties_to_parquet.Keys.negative_value
    baseline = tarred_properties_to_parquet.Keys.baseline
    scale = tarred_properties_to_parquet.Keys.scale
    gq = benchmark_variant_filter.Keys.gq
    training = "training"
    validation = "validation"
    pickle_suffix = ".pickle"
    lock_suffix = ".lock"
    remainder = "remainder"
    remainder_sizes = f"{remainder}_sizes{pickle_suffix}"
    training_remainder = f"{training}_{remainder}"
    validation_remainder = f"{validation}_{remainder}"
    manifest_file = "manifest.json"


class Default:
    shuffle = True
    random_seed = 0
    non_ml_properties = tarred_properties_to_parquet.non_ml_properties
    # don't exclude variant_weights: they're not fed into the neural net, but we use them in
    # training:
    non_ml_needed_for_training = frozenset({Keys.id, Keys.variant_weights})
    excluded_properties = non_ml_properties.difference(non_ml_needed_for_training)
    scale_bool_properties = False,
    validation_proportion = 0.2
    temp_dir = Path(tempfile.gettempdir())
    sleep_time_s = 0.1
    delete_after_read = True
    max_look_ahead_batches = 3


def _get_lock_path(file_path: Path) -> Path:
    """ Unified function for getting lock-file Path that corresponds to a target path """
    return file_path.with_suffix(file_path.suffix + Keys.lock_suffix)


if _debug_file_locks:
    @contextlib.contextmanager
    def locked_write(file_name: Path, mode: str) -> IO:
        """
        Get context manager for writing to a file, generating lock file during write to indicate
        that reading should wait until the write operation is finished
        """
        lock_file = _get_lock_path(file_name)
        with open(lock_file, 'w'):
            # just create lock file
            print(f"created lock {lock_file} for {file_name}")
        with open(file_name, mode) as f_out:
            yield f_out
        print(f"wrote file {file_name}")
        lock_file.unlink()
        print(f"unlocked file {file_name}")

    @contextlib.contextmanager
    def locked_read(
            file_name: Path,
            mode: str,
            delete_after_read: bool = Default.delete_after_read,
            sleep_time_s: float = Default.sleep_time_s
    ) -> IO:
        # wait until file exists
        n_wait = 0
        while not file_name.is_file():
            n_wait += 1
            if n_wait % 100 == 0:
                print(f"waiting for {file_name} to exist")
            time.sleep(sleep_time_s)
        # wait until lock file is deleted
        lock_file = _get_lock_path(file_name)
        n_wait = 0
        while lock_file.is_file():
            n_wait += 1
            if n_wait % 100 == 0:
                print(f"waiting for {file_name} to be unlocked")
            time.sleep(sleep_time_s)
        # read
        with open(file_name, mode) as f_in:
            print(f"opened {file_name}")
            yield f_in
        print(f"read {file_name}")
        if delete_after_read:
            file_name.unlink()

else:
    @contextlib.contextmanager
    def locked_write(file_name: Path, mode: str) -> IO:
        """
        Get context manager for writing to a file, generating lock file during write to indicate
        that reading should wait until the write operation is finished
        """
        lock_file = _get_lock_path(file_name)
        with open(lock_file, 'w'):
            pass  # just create lock file
        with open(file_name, mode) as f_out:
            yield f_out
        lock_file.unlink()


    @contextlib.contextmanager
    def locked_read(
            file_name: Path,
            mode: str,
            delete_after_read: bool = Default.delete_after_read,
            sleep_time_s: float = Default.sleep_time_s
    ) -> IO:
        # wait until file exists
        while not os.path.isfile(file_name):
            time.sleep(sleep_time_s)
        # wait until lock file is deleted
        lock_file = _get_lock_path(file_name)
        while os.path.isfile(lock_file):
            time.sleep(sleep_time_s)
        # read
        with open(file_name, mode) as f_in:
            yield f_in
        if delete_after_read:
            file_name.unlink()


BatchIterator = TypeVar("BatchIterator", bound="BatchIteratorBase")


class BatchIteratorBase(Iterator, ABC):
    __slots__ = ("pickle_folder", "variants_per_batch", "parquet_file")
    __state_keys__ = ("pickle_folder", "variants_per_batch", "parquet_file")

    def __init__(
            self: BatchIterator,
            pickle_folder: Path | str,
            variants_per_batch: int,
            parquet_file: Path | str
    ):
        self.pickle_folder = Path(pickle_folder)
        self.variants_per_batch = variants_per_batch
        self.parquet_file = Path(parquet_file)

    def __iter__(self: BatchIterator) -> BatchIterator:
        return self

    def set_pickle_folder(self: BatchIterator, pickle_folder: Path) -> BatchIterator:
        self.pickle_folder = pickle_folder
        return self

    @property
    def state_dict(self: BatchIterator) -> dict[str, Any]:
        return {key: getattr(self, key) for key in self.__class__.__state_keys__}

    @property
    def json_dict(self: BatchIterator) -> dict[str, Any]:
        return {
            key: (
                f"{value}" if isinstance(value, Path)
                else value.tolist() if hasattr(value, "tolist")
                else value
            )
            for key, value in self.state_dict.items()
        }

    def to_json(self: BatchIterator) -> None:
        json_file = self.pickle_folder / Keys.manifest_file
        with locked_write(json_file, 'w') as f_out:
            json.dump(self.json_dict, f_out)

    @classmethod
    def from_json(
            cls: BatchIterator.__class__,
            pickle_folder: Path,
            delete_after_read: bool = Default.delete_after_read,
            sleep_time_s: float = Default.sleep_time_s
    ) -> BatchIterator:
        json_file = pickle_folder / Keys.manifest_file
        with locked_read(
                json_file, mode='r', delete_after_read=delete_after_read, sleep_time_s=sleep_time_s
        ) as f_in:
            json_dict = json.load(f_in)
        return cls(**json_dict)


class VcfTensorDataLoaderBase:
    """
    A DataLoader-like object for a set of tensors that can be much faster than
    TensorDataset + DataLoader because dataloader grabs individual indices of
    the dataset and calls cat (slow).
    Source: https://discuss.pytorch.org/t/dataloader-much-slower-than-manual-batching/27014/6
    """
    __slots__ = (
        "parquet_path", "variants_per_batch", "torch_device", "progress_logger", "shuffle",
        "random_generator", "_current_parquet_file", "process_executor",
        "_wanted_columns", "_parquet_files", "_temp_dir",
        "_num_rows", "_current_row", "_columns", "_variant_weights_column", "_variant_id_column",
        "_variant_columns", "_format_columns", "_sample_ids", "_property_names",
        "_tensor_baseline", "_tensor_scale",
        "_training_rows", "_validation_rows"
    )

    def __init__(
            self,
            parquet_path: Path,
            properties_scaling_json: Path,
            process_executor: concurrent.futures.ProcessPoolExecutor,
            variants_per_batch: int,
            torch_device_kind: training_utils.TorchDeviceKind,
            progress_logger: training_utils.ProgressLogger,
            temp_dir: Path = Default.temp_dir,
            shuffle: bool = Default.shuffle,
            scale_bool_values: bool = Default.scale_bool_properties,
            random_generator: common.GeneratorInit = Default.random_seed,
            excluded_properties: set[str] = Default.excluded_properties,
            _current_parquet_file: Optional[str] = None,
            _parquet_files: Optional[list[str]] = None,
            _current_row: int = 0
    ):
        self.parquet_path = tarred_properties_to_parquet.extract_tar_to_folder(parquet_path) \
            if parquet_path.suffix == ".tar" else parquet_path
        self.process_executor = process_executor
        self.variants_per_batch = variants_per_batch
        self.progress_logger = progress_logger
        self.torch_device = torch_device_kind.get_device(progress_logger=progress_logger)
        self.shuffle = shuffle
        self._temp_dir = self.__class__.get_temp_subdir(temp_dir)
        self.random_generator = \
            common.init_generator(generator_init=random_generator) if shuffle else None
        self._set_column_info(
            excluded_properties=excluded_properties,
            properties_scaling_json=properties_scaling_json,
            scale_bool_values=scale_bool_values
        )
        self._parquet_files = [] if _parquet_files is None else _parquet_files
        self._num_rows = tarred_properties_to_parquet.get_parquet_file_num_rows(self.parquet_path)
        self._current_row = _current_row
        self._current_parquet_file = _current_parquet_file

    @classmethod
    def get_temp_subdir(cls, temp_dir: Path) -> Path:
        parent_subdir = temp_dir / cls.__name__
        parent_subdir.mkdir(parents=True, exist_ok=True)
        temp_subdir = tempfile.mkdtemp(dir=parent_subdir)
        atexit.register(shutil.rmtree, temp_subdir)
        return Path(temp_subdir)

    def _get_pickle_folder(self) -> Path:
        return Path(tempfile.mkdtemp(dir=self._temp_dir))

    def _set_column_info(
            self,
            excluded_properties: set[str],
            properties_scaling_json: Path,
            scale_bool_values: bool
    ):
        f"""
        Set a bunch of column info so that it can be quickly retrieved later
          -self._wanted_columns is the list of raw (str) columns saved in the parquet file that we
           want to load, i.e. everything but the excluded properties
          -self._variant_columns is the unflattened (i.e. MultiIndex) columns that correspond to
           per-variant properties
          -self._format_columns is the unflattened columns that correspond to per-genotype / format
           properties
        Have a few annoying issues to work around:
          -Absolutely don't want {Keys.row} in any columns (it's the name of the row index, not a
           column).
          -Can only determine if a column is for an excluded property after unflattening
        """
        with open(properties_scaling_json, 'r') as f_in:
            properties_scaling = json.load(f_in)

        # get raw_columns and _columns, explicitly removing Keys.row, the name of the row index
        # -raw_columns are ALL the columns in the dataframe, flattened so columns are not
        #  multi-index (this is what's stored in the parquet file)
        # -self.columns is the unflattened multi-index (tuple of sample_id, property_name) columns
        #  that are easier to process and organize
        raw_columns, self._columns = zip(*(
            (flat_column, tarred_properties_to_parquet.unflatten_column_name(flat_column))
            for flat_column in tarred_properties_to_parquet.get_parquet_file_columns(
                self.parquet_path
            )
            if flat_column != Keys.row
        ))
        # get rid of all the columns we don't want (corresponding to excluded properties)
        self._wanted_columns, self._columns = zip(*(
            (raw_column, unflat_column)
            for raw_column, unflat_column in zip(raw_columns, self._columns)
            if unflat_column[1] not in excluded_properties
        ))
        self._variant_weights_column = None if Keys.variant_weights in excluded_properties \
            else next(
                index for index, column in enumerate(self._columns)
                if column[1] == Keys.variant_weights
            )
        self._variant_id_column = next(
            index for index, column in enumerate(self._columns) if column[1] == Keys.id
        )
        # get the DataFrame columns that correspond to per-variant / INFO columns
        self._variant_columns = numpy.array(
            [index for index, column in enumerate(self._columns)
             if not isinstance(column[0], str) and
             column[1] not in Default.non_ml_needed_for_training]
        )
        # get the DataFrame columns that correspond to per-sample / FORMAT columns
        self._format_columns = numpy.array(
            [index for index, column in enumerate(self._columns)
             if isinstance(column[0], str) and column[1] not in Default.non_ml_needed_for_training]
        )
        # get property names (corresponding to dimension 2 of output tensor)
        # organize with per-variant / INFO properties first
        self._property_names = [
            column[1] for column in self._columns
            if not isinstance(column[0], str) and
            column[1] not in Default.non_ml_needed_for_training
        ]
        # get sample IDs in order, and add per-sample / FORMAT properties to property names
        found_properties = set()
        found_sample_ids = set()
        self._sample_ids = []
        for sample_id, prop_name in self._columns:
            if prop_name in Default.non_ml_needed_for_training or not isinstance(sample_id, str):
                continue  # skip properties not in final tensor, or INFO properties
            if sample_id not in found_sample_ids:
                found_sample_ids.add(sample_id)
                self._sample_ids.append(sample_id)
            if prop_name not in found_properties:
                found_properties.add(prop_name)
                self._property_names.append(prop_name)

        # get baseline and scale to (approximately) z-score output tensor
        def _get_baseline(_prop_name) -> float:
            if scale_bool_values or Keys.positive_value not in properties_scaling[_prop_name]:
                return properties_scaling[_prop_name][Keys.baseline]
            else:
                return 0.0

        def _get_scale(_prop_name) -> float:
            if scale_bool_values or Keys.positive_value not in properties_scaling[_prop_name]:
                return properties_scaling[_prop_name][Keys.scale]
            else:
                return 1.0

        self._tensor_baseline = torch.tensor(
            [_get_baseline(prop_name) for prop_name in self._property_names],
            dtype=torch.float32, device=self.torch_device
        )
        self._tensor_scale = torch.tensor(
            [1.0 / _get_scale(prop_name) for prop_name in self._property_names],
            dtype=torch.float32, device=self.torch_device
        )

    def _set_parquet_files(self):
        self._parquet_files = list(self.parquet_path.glob("*.parquet"))
        if self.shuffle:
            # randomize order of visiting partitions
            self._parquet_files = self.random_generator.permutation(self._parquet_files).tolist()
        else:
            # reverse so we can visit in order by popping the files from the end
            self._parquet_files = self._parquet_files[::-1]

    @property
    def num_samples(self) -> int:
        """ return the number of samples """
        return len(self._sample_ids)

    @property
    def sample_ids(self):
        """ return array of sample IDs """
        return self._sample_ids

    @property
    def property_names(self) -> list[str]:
        """
        return the names of each property ((e.g. 3rd dimension in num-variants x num-samples
        x num-properties tensor)
        """
        return self._property_names

    @property
    def num_properties(self) -> int:
        """
        return the number of distinct properties (e.g. 3rd dimension in num-variants x num-samples
        x num-properties tensor)
        """
        return len(self.property_names)

    def __len__(self) -> int:
        """ return the number of remaining batches in the iterator """
        n_batches, extra_rows = divmod(self._num_rows - self._current_row, self.variants_per_batch)
        return n_batches + 1 if extra_rows else n_batches


class BatchPicklerBase:
    @classmethod
    def _load_buffer(
            cls,
            parquet_file: Path,
            wanted_columns: list[str],
            wanted_ids: Optional[Collection[str]] = None,
            variant_id_column: Optional[int] = None
    ) -> pandas.array:
        """ get the buffer that corresponds only to rows that are part of the remainders """
        if wanted_ids is None:
            filters = None
        else:
            if variant_id_column is None:
                raise ValueError("Must specify variant_id_column when specifying wanted_ids")
            _id_filter = (
                wanted_columns[variant_id_column], "in", wanted_ids
            )
            filters = [_id_filter]
        return pandas.read_parquet(parquet_file, columns=wanted_columns, filters=filters).values
