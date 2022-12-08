import glob
import os
import atexit
import shutil
import tempfile
import contextlib
import time
import json
import concurrent.futures
import pickle
import numpy
import pandas
import pandas.core.arrays
import torch
from collections import defaultdict
from typing import Optional, Union, Any, IO
from collections.abc import Iterator, Sequence
from gq_recalibrator import tarred_properties_to_parquet, training_utils
from sv_utils import get_truth_overlap, benchmark_variant_filter


ArrowStringArray = pandas.core.arrays.string_arrow.ArrowStringArray
ConfidentVariants = get_truth_overlap.ConfidentVariants
_debug_file_locks = False


class Keys:
    row = tarred_properties_to_parquet.Keys.row
    cuda = "cuda"
    cpu = "cpu"
    mps = "mps"  # Macbook Performance Shaders
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
    random_state = 0
    non_ml_properties = tarred_properties_to_parquet.non_ml_properties
    # don't exclude variant_weights: they're not fed into the neural net, but we use them in training:
    non_ml_needed_for_training = frozenset({Keys.id, Keys.variant_weights})
    excluded_properties = non_ml_properties.difference(non_ml_needed_for_training)
    scale_bool_properties = False,
    validation_proportion = 0.2
    temp_dir = tempfile.gettempdir()
    sleep_time_s = 0.1
    delete_after_read = True
    max_look_ahead_batches = 3


def get_torch_device(device_name: str, progress_logger: Optional[training_utils.ProgressLogger]) -> torch.device:
    if device_name == Keys.cuda:
        if torch.cuda.is_available():
            if progress_logger is not None:
                progress_logger.log(f"Using {Keys.cuda} for torch")
            torch.cuda.empty_cache()
            return torch.device("cuda:0")
        else:
            if progress_logger is not None:
                progress_logger.log(f"{Keys.cuda} not available, falling back to {Keys.cpu} for torch")
            return torch.device("cpu")
    elif device_name == Keys.mps:
        if hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
            if progress_logger is not None:
                progress_logger.log(f"Using {Keys.mps} for torch")
            torch.cuda.empty_cache()
            return torch.device("mps")
        else:
            if progress_logger is not None:
                progress_logger.log(f"{Keys.mps} not available, falling back to {Keys.cpu} for torch")
            return torch.device("cpu")

    else:
        if progress_logger is not None:
            progress_logger.log(f"Using {Keys.cpu} for torch")
        return torch.device("cpu")


if _debug_file_locks:
    @contextlib.contextmanager
    def locked_write(file_name: str, mode: str) -> IO:
        """
        Get context manager for writing to a file, generating lock file during write to indicate that reading should wait
        until the write operation is finished
        """
        lock_file = file_name + Keys.lock_suffix
        with open(lock_file, 'w') as f_lock:
            pass  # just create lock file
        with open(file_name, mode) as f_out:
            yield f_out
        print(f"wrote file {file_name}")
        os.remove(lock_file)
        print(f"unlocked file {file_name}")


    @contextlib.contextmanager
    def locked_read(
            file_name: str,
            mode: str,
            delete_after_read: bool = Default.delete_after_read,
            sleep_time_s: float = Default.sleep_time_s
    ) -> IO:
        # wait until file exists
        n_wait = 0
        while not os.path.isfile(file_name):
            n_wait += 1
            if n_wait % 100 == 0:
                print(f"waiting for {file_name} to exist")
            time.sleep(sleep_time_s)
        # wait until lock file is deleted
        lock_file = file_name + Keys.lock_suffix
        n_wait = 0
        while os.path.isfile(lock_file):
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
            os.remove(file_name)

else:
    @contextlib.contextmanager
    def locked_write(file_name: str, mode: str) -> IO:
        """
        Get context manager for writing to a file, generating lock file during write to indicate that reading should wait
        until the write operation is finished
        """
        lock_file = file_name + Keys.lock_suffix
        with open(lock_file, 'w') as f_lock:
            pass  # just create lock file
        with open(file_name, mode) as f_out:
            yield f_out
        os.remove(lock_file)


    @contextlib.contextmanager
    def locked_read(
            file_name: str,
            mode: str,
            delete_after_read: bool = Default.delete_after_read,
            sleep_time_s: float = Default.sleep_time_s
    ) -> IO:
        # wait until file exists
        while not os.path.isfile(file_name):
            time.sleep(sleep_time_s)
        # wait until lock file is deleted
        lock_file = file_name + Keys.lock_suffix
        while os.path.isfile(lock_file):
            time.sleep(sleep_time_s)
        # read
        with open(file_name, mode) as f_in:
            yield f_in
        if delete_after_read:
            os.remove(file_name)


class VcfTensorDataLoaderBase:
    """
    A DataLoader-like object for a set of tensors that can be much faster than
    TensorDataset + DataLoader because dataloader grabs individual indices of
    the dataset and calls cat (slow).
    Source: https://discuss.pytorch.org/t/dataloader-much-slower-than-manual-batching/27014/6
    """
    __slots__ = (
        "parquet_path", "variants_per_batch", "torch_device", "progress_logger", "shuffle", "random_state",
        "_current_parquet_file", "process_executor", "thread_executor", "_wanted_columns", "_parquet_files",
        "_buffer_index", "_buffer_matrix", "_temp_dir",
        "_num_rows", "_current_row", "_columns", "_variant_weights_column", "_variant_id_column", "_variant_columns",
        "_format_columns", "_sample_ids", "_property_names", "_futures", "_numpy_tensor", "_weights_tensor",
        "_ids_array", "_tensor_baseline", "_tensor_scale", "_training_rows", "_validation_rows"
    )

    def __init__(
            self,
            parquet_path: str,
            properties_scaling_json: str,
            process_executor: concurrent.futures.ProcessPoolExecutor,
            thread_executor: concurrent.futures.ThreadPoolExecutor,
            variants_per_batch: int,
            torch_device: str,
            progress_logger: training_utils.ProgressLogger,
            temp_dir: str = Default.temp_dir,
            shuffle: bool = Default.shuffle,
            scale_bool_values: bool = Default.scale_bool_properties,
            random_state: Union[int, numpy.random.RandomState, None] = Default.random_state,
            excluded_properties: set[str] = Default.excluded_properties,
            _current_parquet_file: Optional[str] = None
    ):
        self.parquet_path = tarred_properties_to_parquet.extract_tar_to_folder(parquet_path) \
            if parquet_path.endswith(".tar") else parquet_path
        self.process_executor = process_executor
        self.thread_executor = thread_executor
        self.variants_per_batch = variants_per_batch
        self.progress_logger = progress_logger
        self.torch_device = get_torch_device(torch_device, self.progress_logger)
        self.shuffle = shuffle
        self._temp_dir = self.__class__.get_temp_subdir(temp_dir)
        self.random_state = random_state if isinstance(random_state, numpy.random.RandomState) \
            else numpy.random.RandomState(random_state)
        self._set_column_info(
            excluded_properties=excluded_properties, properties_scaling_json=properties_scaling_json,
            scale_bool_values=scale_bool_values
        )
        self._numpy_tensor = numpy.empty(
            (self.variants_per_batch, self.num_samples, self.num_properties), dtype=numpy.float32
        )
        if Keys.variant_weights not in excluded_properties:
            self._weights_tensor = numpy.empty(self.variants_per_batch, dtype=numpy.float32)
        self._ids_array = pandas.array([""] * self.variants_per_batch, dtype="string[pyarrow]")
        self._buffer_index = 0
        self._buffer_matrix = pandas.DataFrame()
        self._parquet_files = []
        self._num_rows = tarred_properties_to_parquet.get_parquet_file_num_rows(self.parquet_path)
        self._current_row = 0
        self._current_parquet_file = _current_parquet_file
        self._futures = []

    @classmethod
    def get_temp_subdir(cls, temp_dir: str) -> str:
        parent_subdir = os.path.join(temp_dir, cls.__name__)
        os.makedirs(parent_subdir, exist_ok=True)
        temp_subdir = tempfile.mkdtemp(dir=parent_subdir)
        atexit.register(shutil.rmtree, temp_subdir)
        return temp_subdir

    def _set_column_info(self, excluded_properties: set[str], properties_scaling_json: str, scale_bool_values: bool):
        f"""
        Set a bunch of column info so that it can be quickly retrieved later
          -self._wanted_columns is the list of raw (str) columns saved in the parquet file that we want to load, i.e.
           everything but the excluded properties
          -self._variant_columns is the unflattened (i.e. MultiIndex) columns that correspond to per-variant properties
          -self._format_columns is the unflattened columns that correspond to per-genotype / format properties
        Have a few annoying issues to work around:
          -Absolutely don't want {Keys.row} in any columns (wanted_columns, format_columns). Also it can't be unflattened
          -Can only determine if a column is for an excluded property after unflattening
        """
        with open(properties_scaling_json, 'r') as f_in:
            properties_scaling = json.load(f_in)

        # raw_columns are ALL the columns in the dataframe, flattened so columns are not a multi-index
        raw_columns = [
            flat_column for flat_column in tarred_properties_to_parquet.get_parquet_file_columns(self.parquet_path)
            if flat_column != Keys.row
        ]
        # raw_columns are ALL the columns in the dataframe (except row), flattened so columns are not a multi-index
        # unflatten the columns to make mult-index (list of tuples) sample_id, property_name
        raw_columns, self._columns = zip(*(
            (flat_column, tarred_properties_to_parquet.unflatten_column_name(flat_column))
            for flat_column in tarred_properties_to_parquet.get_parquet_file_columns(self.parquet_path)
            if flat_column != Keys.row
        ))
        # get rid of all the columns we don't want
        self._wanted_columns, self._columns = zip(*(
            (raw_column, unflat_column) for raw_column, unflat_column in zip(raw_columns, self._columns)
            if unflat_column[1] not in excluded_properties
        ))
        self._variant_weights_column = None if Keys.variant_weights in excluded_properties \
            else next(index for index, column in enumerate(self._columns) if column[1] == Keys.variant_weights)
        self._variant_id_column = next(
            index for index, column in enumerate(self._columns) if column[1] == Keys.id
        )
        # get the DataFrame columns that correspond to per-variant / INFO columns
        self._variant_columns = numpy.array(
            [index for index, column in enumerate(self._columns)
             if not isinstance(column[0], str) and column[1] not in Default.non_ml_needed_for_training]
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
            if not isinstance(column[0], str) and column[1] not in Default.non_ml_needed_for_training
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
        self._parquet_files = glob.glob(f"{self.parquet_path}/*.parquet")
        if self.shuffle:
            # randomize order of visiting partitions
            self._parquet_files = self.random_state.permutation(self._parquet_files).tolist()
        else:
            # reverse so we can visit in order by popping the files from the end
            self._parquet_files = self._parquet_files[::-1]

    @staticmethod
    def _get_next_tensor_blocks(
            buffer_matrix: numpy.ndarray,
            variant_ids_column: int,
            variant_weights_column: Optional[int],
            variant_columns: numpy.ndarray,
            format_columns: numpy.ndarray,
            num_samples: int
    ) -> (ArrowStringArray, Optional[numpy.ndarray], numpy.ndarray, numpy.ndarray):
        """
        Convert DataFrame with multi-index columns (1st index corresponding to sample, with None for per-variant
        properties, 2nd index corresponding to property name) to 3-tensor (num_variants x num_samples x num_properties).
        NOTE:
        1) The per-variant properties will be duplicated during this process
        2) All properties will be converted to 32-bit floats
        3) It is assumed that Categorical properties are already one-hot encoded (by tarred_properties_to_parquet.py)

        Args:
            buffer_matrix: numpy.ndarray
                Numpy matrix (from DataFrame.values) of variant/sample properties with multi-index columns
            next_tensor: numpy.ndarray
                3-tensor of properties used for training or filtering
        Returns:
        """
        # fill num_rows x 1 x num_variant_properties matrix with per-variant properties. Note, the "1" dimension will
        # broadcast to fill the 2nd (num_samples) dim
        num_variant_properties = len(variant_columns)
        return (
            buffer_matrix.take(variant_ids_column, axis=1),
            None if variant_weights_column is None else buffer_matrix.take(variant_weights_column, axis=1),
            buffer_matrix.take(variant_columns, axis=1).reshape(buffer_matrix.shape[0], 1, num_variant_properties),
            buffer_matrix.take(format_columns, axis=1)
                .reshape(buffer_matrix.shape[0], num_samples, len(format_columns) // num_samples)
        )  # noqa E131   # ignore pep8 131 error, I think this is more legible

    @property
    def num_samples(self) -> int:
        return len(self._sample_ids)

    @property
    def sample_ids(self):
        return self._sample_ids

    @property
    def property_names(self) -> list[str]:
        return self._property_names

    @property
    def num_properties(self) -> int:
        return len(self.property_names)

    def __len__(self) -> int:
        # report the number of remaining batches in the iterator
        n_batches, extra_rows = divmod(self._num_rows - self._current_row, self.variants_per_batch)
        return n_batches + 1 if extra_rows else n_batches


class VcfFilterTensorDataLoader(VcfTensorDataLoaderBase):
    def __init__(
            self,
            parquet_path: str,
            properties_scaling_json: str,
            process_executor: concurrent.futures.ProcessPoolExecutor,
            thread_executor: concurrent.futures.ThreadPoolExecutor,
            variants_per_batch: int,
            torch_device: str,
            progress_logger: training_utils.ProgressLogger,
            temp_dir: str = Default.temp_dir,
            shuffle: bool = Default.shuffle,
            scale_bool_values: bool = Default.scale_bool_properties,
            random_state: Union[int, numpy.random.RandomState, None] = Default.random_state,
            excluded_properties: set[str] = Default.excluded_properties
    ):
        # basically it's the VcfTensorDataLoaderBase with variant weights property excluded and explicitly no
        # validation
        super().__init__(
            parquet_path=parquet_path, properties_scaling_json=properties_scaling_json,
            process_executor=process_executor, thread_executor=thread_executor,
            variants_per_batch=variants_per_batch, torch_device=torch_device, progress_logger=progress_logger,
            shuffle=shuffle, temp_dir=temp_dir, scale_bool_values=scale_bool_values, random_state=random_state,
            excluded_properties=excluded_properties.union((Keys.variant_weights,))
        )

    def __iter__(self):
        if self._buffer_future is None:
            # only executes on first iteration, so there's no remainder rows
            if not self._parquet_files:
                self._set_parquet_files()
            if self._current_parquet_file is None:
                self._current_parquet_file = self._parquet_files.pop()
            self._buffer_future = self.process_executor.submit(
                VcfFilterTensorDataLoader._parquet_file_to_buffer, parquet_file=self._current_parquet_file,
                wanted_columns=self._wanted_columns, random_state=self.random_state if self.shuffle else None
            )
        return self

    @staticmethod
    def _parquet_file_to_buffer(
            parquet_file: str,
            wanted_columns: list[str],
            random_state: Optional[numpy.random.RandomState],
            remainder_buffer: Optional[numpy.ndarray] = None
    ) -> (numpy.ndarray, Optional[numpy.random.RandomState]):
        buffer = pandas.read_parquet(parquet_file, columns=wanted_columns).values
        if random_state is not None:
            buffer = buffer.take(random_state.permutation(buffer.shape[0]), axis=0)
        return (buffer, random_state) if remainder_buffer is None \
            else (numpy.concatenate((remainder_buffer, buffer), axis=0), random_state)

    def __next__(self) -> (torch.Tensor, ArrowStringArray):
        """
        emit tensors for next batch
        Returns:
            properties_tensor: torch.Tensor
                num_variants x num_samples x num_properties 32-bit float tensor of training / filtering properties
            variant_ids: ArrowStringArray
                num_variants array of variant IDs
        """
        if self._current_row >= self._num_rows:
            # since batch size is unlikely to be an even divisor of the data set size, we can either exactly traverse
            # the data set (and have an uneven final batch) or have uniform batch size (but only approximately traverse
            # the data set each epoch).  Choose to use uniform batch size and approximate traversal:
            self._current_row -= self._num_rows
            raise StopIteration
        next_buffer_index = self._buffer_index + self.variants_per_batch
        if next_buffer_index > self._buffer_matrix.shape[0]:
            self._load_buffer()
            next_buffer_index = self._buffer_index + self.variants_per_batch

        num_variant_properties = len(self._variant_columns)
        (self._ids_array[:], __, self._numpy_tensor[:, :, :num_variant_properties],
         self._numpy_tensor[:, :, num_variant_properties:]) = VcfTensorDataLoaderBase._get_next_tensor_blocks(
            buffer_matrix=self._buffer_matrix[self._buffer_index:next_buffer_index, :],
            variant_ids_column=self._variant_id_column, variant_weights_column=None,
            variant_columns=self._variant_columns, format_columns=self._format_columns, num_samples=self.num_samples
        )
        self._buffer_index = next_buffer_index
        self._current_row += self.variants_per_batch

        return (
            self._tensor_scale * (
                torch.tensor(self._numpy_tensor, dtype=torch.float32, device=self.torch_device)
                - self._tensor_baseline
            ),
            self._ids_array
        )

    def _load_buffer(self):
        self._buffer_matrix, self.random_state = self._buffer_future.result()
        self._buffer_index = 0
        if not self._parquet_files:
            self._set_parquet_files()
        self._current_parquet_file = self._parquet_files.pop()
        remainder_rows = divmod(self._buffer_matrix.shape[0], self.variants_per_batch)[1]
        self._buffer_future = self.process_executor.submit(
            VcfFilterTensorDataLoader._parquet_file_to_buffer, parquet_file=self._current_parquet_file,
            wanted_columns=self._wanted_columns, random_state=self.random_state if self.shuffle else None,
            remainder_buffer=self._buffer_matrix[-remainder_rows:, :]
        )


class TrainingBatchIterator(Iterator):
    """
    Iterator class that approximately iterates at the appropriate validation proportion. Works over batches and doesn't
    allow rounding error to gradually accumulate excess training or validation variants in the remainder
    """
    __slots__ = (
        "num_training_variants", "num_validation_variants", "variants_per_batch", "pickle_folder", "parquet_file",
        "training_batch", "validation_batch", "previous_training_remainder_ids", "previous_validation_remainder_ids",
        "num_training_batches", "num_validation_batches", "_approx_validation_proportion", "num_training_remainder",
        "num_validation_remainder")
    __state_keys__ = (
        "num_training_variants", "num_validation_variants", "variants_per_batch", "pickle_folder", "parquet_file",
        "training_batch", "validation_batch", "previous_training_remainder_ids", "previous_validation_remainder_ids"
    )

    def __init__(
            self,
            num_training_variants: int,
            num_validation_variants: int,
            variants_per_batch: int,
            pickle_folder: str,
            parquet_file: str,
            training_batch: int = 0,
            validation_batch: int = 0,
            previous_training_remainder_ids: Sequence = (),
            previous_validation_remainder_ids: Sequence = ()
    ):
        self.num_training_variants = num_training_variants
        self.num_validation_variants = num_validation_variants
        self.variants_per_batch = variants_per_batch
        self.pickle_folder = pickle_folder
        self.parquet_file = parquet_file
        self.training_batch = training_batch
        self.validation_batch = validation_batch
        self.num_training_batches, self.num_training_remainder = divmod(
            num_training_variants, variants_per_batch
        )
        self.num_validation_batches, self.num_validation_remainder = divmod(
            num_validation_variants, variants_per_batch
        )
        self._approx_validation_proportion = \
            self.num_validation_batches / (self.num_validation_batches + self.num_training_batches)
        self.previous_training_remainder_ids = pandas.array(previous_training_remainder_ids,
                                                            dtype="string[pyarrow]")
        self.previous_validation_remainder_ids = pandas.array(previous_validation_remainder_ids,
                                                              dtype="string[pyarrow]")

    @property
    def is_training_batch(self) -> Optional[bool]:
        """ return True if next batch is training, False if it's validation, and None if it's time to stop iteration """
        if self.training_batch >= self.num_training_batches:
            return None if self.validation_batch >= self.num_validation_batches else False
        elif self.validation_batch >= self.num_validation_batches:
            return True
        return self.training_batch * self._approx_validation_proportion < \
               (self.validation_batch + 1) * (1.0 - self._approx_validation_proportion)

    @property
    def has_next(self) -> bool:
        return self.training_batch < self.num_training_batches or self.validation_batch < self.num_validation_batches

    @property
    def state_dict(self) -> dict[str, Any]:
        raw_state_dict = {key: getattr(self, key) for key in self.__state_keys__}
        return {
            key: value.tolist() if hasattr(value, "tolist") else value
            for key, value in raw_state_dict.items()
        }

    def to_json(self):
        json_file = os.path.join(self.pickle_folder, Keys.manifest_file)
        with locked_write(json_file, 'w') as f_out:
            json.dump(self.state_dict, f_out)

    @staticmethod
    def from_json(
            pickle_folder: str,
            delete_after_read: bool = Default.delete_after_read,
            sleep_time_s: float = Default.sleep_time_s
    ) -> "TrainingBatchIterator":
        json_file = os.path.join(pickle_folder, Keys.manifest_file)
        with locked_read(json_file, mode='r', delete_after_read=delete_after_read, sleep_time_s=sleep_time_s) as f_in:
            state_dict = json.load(f_in)
        return TrainingBatchIterator(**state_dict)

    def __iter__(self) -> "TrainingBatchIterator":
        return self

    def __next__(self) -> (bool, str, int):
        is_training_batch = self.is_training_batch
        if is_training_batch is None:
            raise StopIteration
        if is_training_batch:
            batch_number = self.training_batch
            self.training_batch += 1
        else:
            batch_number = self.validation_batch
            self.validation_batch += 1
        pickle_file = VcfTrainingTensorDataLoader.get_nth_pickle_file(
            pickle_folder=self.pickle_folder, n=batch_number, is_training=is_training_batch
        )
        return is_training_batch, pickle_file, batch_number


class VcfTrainingTensorDataLoader(VcfTensorDataLoaderBase):
    __slots__ = (
        "validation_proportion", "_training_batch", "_validation_batch", "batch_iterator", "previous_batch_iterator",
        "_genotype_is_good_dict", "_genotype_is_bad_dict", "_genotype_is_good_tensor", "_genotype_is_bad_tensor",
        "_trainable_variant_ids", "_gq_column", "_gq_tensor", "max_look_ahead_batches"
    )
    __state_keys__ = (
        "random_state", "_current_parquet_file", "_parquet_files", "_buffer_index", "_current_row",
        "_training_batch", "_validation_batch", "batch_iterator", "previous_batch_iterator"
    )

    def __init__(
            self,
            parquet_path: str,
            properties_scaling_json: str,
            truth_json: str,
            process_executor: concurrent.futures.ProcessPoolExecutor,
            thread_executor: concurrent.futures.ThreadPoolExecutor,
            variants_per_batch: int,
            torch_device: str,
            progress_logger: training_utils.ProgressLogger,
            temp_dir: str = Default.temp_dir,
            shuffle: bool = Default.shuffle,
            scale_bool_values: bool = Default.scale_bool_properties,
            random_state: Union[int, numpy.random.RandomState, None] = Default.random_state,
            excluded_properties: set[str] = Default.excluded_properties,
            validation_proportion: float = Default.validation_proportion,
            max_look_ahead_batches: int = Default.max_look_ahead_batches,
            batch_iterator: Optional[TrainingBatchIterator] = None,
            previous_batch_iterator: Optional[TrainingBatchIterator] = None
    ):
        # basically it's the VcfTensorDataLoaderBase with a train/test split
        super().__init__(
            parquet_path=parquet_path, properties_scaling_json=properties_scaling_json,
            process_executor=process_executor, thread_executor=thread_executor,
            variants_per_batch=variants_per_batch, torch_device=torch_device, progress_logger=progress_logger,
            temp_dir=temp_dir, shuffle=shuffle, scale_bool_values=scale_bool_values, random_state=random_state,
            excluded_properties=excluded_properties,
            _current_parquet_file=None if batch_iterator is None else batch_iterator.parquet_file
        )
        self.max_look_ahead_batches = max_look_ahead_batches
        self.batch_iterator = batch_iterator
        self.previous_batch_iterator = previous_batch_iterator
        if validation_proportion <= 0 or validation_proportion >= 1:
            raise ValueError("Validation proportion must be in open interval (0, 1)")
        self.validation_proportion = validation_proportion
        self._training_batch = 0
        self._validation_batch = 0
        self._genotype_is_good_dict, self._genotype_is_bad_dict = self._reorganize_confident_variants(
            benchmark_variant_filter.TruthData.load_confident_variants(truth_json)
        )
        self._trainable_variant_ids = frozenset(
            self._genotype_is_good_dict.keys()
        ).union(self._genotype_is_bad_dict.keys())
        self._genotype_is_good_tensor = numpy.empty((self.variants_per_batch, self.num_samples), dtype=numpy.float32)
        self._genotype_is_bad_tensor = numpy.empty((self.variants_per_batch, self.num_samples), dtype=numpy.float32)
        # get the DataFrame columns that correspond to per-sample / FORMAT columns
        self._gq_column = next(
            index for index, property_name in enumerate(self._property_names) if property_name == Keys.gq
        )
        self._gq_tensor = numpy.empty((self.variants_per_batch, self.num_samples), dtype=numpy.float32)

    def _reorganize_confident_variants(
            self,
            confident_variants: ConfidentVariants
    ) -> tuple[dict[str, numpy.ndarray], dict[str, numpy.ndarray]]:
        genotype_is_good = defaultdict(list)
        genotype_is_bad = defaultdict(list)
        sample_index_map = {sample_id: index for index, sample_id in enumerate(self.sample_ids)}
        for sample_id, sample_confident_variants in confident_variants.items():
            sample_index = sample_index_map[sample_id]
            for variant_id in sample_confident_variants.good_variant_ids:
                genotype_is_good[variant_id].append(sample_index)
            for variant_id in sample_confident_variants.bad_variant_ids:
                genotype_is_bad[variant_id].append(sample_index)

        int_dtype = numpy.min_scalar_type(self.num_samples)
        genotype_is_good = {
            variant_id: numpy.array(sample_indices, dtype=int_dtype)
            for variant_id, sample_indices in genotype_is_good.items()
        }
        genotype_is_bad = {
            variant_id: numpy.array(sample_indices, dtype=int_dtype)
            for variant_id, sample_indices in genotype_is_bad.items()
        }

        return genotype_is_good, genotype_is_bad

    def __iter__(self):
        while len(self._futures) < self.max_look_ahead_batches:
            # in order to truly be able to restore state, need to be able to restore future results at program start
            self._submit_next_batch()
        if self.batch_iterator is None:
            self._get_next_batch_iterator()
        return self

    def _get_next_batch_iterator(self):
        future, pickle_folder = self._futures.pop(0)
        future_exception = future.exception()
        if future_exception is not None:
            # need to shutdown other processes before we can propagate exception
            self.process_executor.shutdown(wait=False, cancel_futures=True)
            training_utils.kill_child_processes(parent_pid=os.getpid())
            raise RuntimeError(f"Pickling tensors in folder {pickle_folder}", future_exception)
        self.previous_batch_iterator = self.batch_iterator
        self.batch_iterator = TrainingBatchIterator.from_json(pickle_folder, delete_after_read=True)
        self._submit_next_batch()

    def _submit_next_batch(self):
        # get the next parquet file
        if not self._parquet_files:
            self._set_parquet_files()
        self._current_parquet_file = self._parquet_files.pop()
        # randomly get a (temporary) pickle folder
        pickle_folder = self._get_pickle_folder()
        # in order to truly be able to restore state, need to be able to restore future results at start time...
        __, previous_pickle_folder = self._futures[-1] if len(self._futures) > 0 else (None, None)
        future = self.process_executor.submit(
            VcfTrainingTensorDataLoader._pickle_batch_tensors,
            parquet_file=self._current_parquet_file, pickle_folder=pickle_folder,
            previous_pickle_folder=previous_pickle_folder,
            wanted_columns=self._wanted_columns, random_seed=self._next_random_seed,
            variants_per_batch=self.variants_per_batch, num_samples=self.num_samples,
            num_properties=self.num_properties, variant_weights_column=self._variant_weights_column,
            variant_columns=self._variant_columns, format_columns=self._format_columns,
            validation_proportion=self.validation_proportion, trainable_variant_ids=self._trainable_variant_ids,
            variant_id_column=self._variant_id_column
        )
        self._futures.append((future, pickle_folder))

    @property
    def _next_random_seed(self) -> Optional[int]:
        if self.shuffle:
            seed_dtype = numpy.uint32
            seed_dtype_info = numpy.iinfo(seed_dtype)
            return self.random_state.randint(seed_dtype_info.min, seed_dtype_info.max, dtype=seed_dtype)
        else:
            return None

    def __next__(self) -> (bool, torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor):
        """
        emit tensors for next batch
        Returns:
            properties_tensor: torch.Tensor
                num_variants x num_samples x num_properties 32-bit float tensor of training / filtering properties
            weights_tensor: torch.Tensor
                num_variants array of training weights that help balance variants
            variant_ids: ArrowStringArray
                num_variants array of variant IDs
        """
        if self._current_row >= self._num_rows:
            # since batch size is unlikely to be an even divisor of the data set size, we can either exactly traverse
            # the data set (and have an uneven final batch) or have uniform batch size (but only approximately traverse
            # the data set each epoch).  Choose to use uniform batch size and approximate traversal:
            self._current_row -= self._num_rows
            raise StopIteration

        if not self.batch_iterator.has_next:
            self._get_next_batch_iterator()

        is_training_batch, pickle_file, __ = next(self.batch_iterator)
        if is_training_batch:
            self._training_batch += 1
        else:
            self._validation_batch += 1
        self._current_row += self.variants_per_batch

        with locked_read(pickle_file, mode="rb", delete_after_read=True, sleep_time_s=Default.sleep_time_s) as f_in:
            ids_array, weights_tensor, numpy_tensor = pickle.load(f_in)
        self._gq_tensor[:] = numpy_tensor.take(self._gq_column, axis=2)
        VcfTrainingTensorDataLoader._fill_truth_tensor(ids_array, self._genotype_is_good_dict,
                                                       self._genotype_is_good_tensor)
        VcfTrainingTensorDataLoader._fill_truth_tensor(ids_array, self._genotype_is_bad_dict,
                                                       self._genotype_is_bad_tensor)
        return (
            is_training_batch,
            self._tensor_scale * (
                torch.tensor(numpy_tensor, dtype=torch.float32, device=self.torch_device)
                - self._tensor_baseline
            ).nan_to_num(),
            torch.tensor(weights_tensor, dtype=torch.float32, device=self.torch_device),
            torch.tensor(self._gq_tensor, dtype=torch.float32, device=self.torch_device),
            torch.tensor(self._genotype_is_good_tensor, dtype=torch.float32, device=self.torch_device),
            torch.tensor(self._genotype_is_bad_tensor, dtype=torch.float32, device=self.torch_device),
        )

    @staticmethod
    def _fill_truth_tensor(
            variant_ids: numpy.ndarray,
            variant_id_trainable_map: dict[str, numpy.ndarray],
            truth_tensor: numpy.ndarray
    ):
        truth_tensor[:] = 0.0
        for row, variant_id in enumerate(variant_ids):
            truth_indices = variant_id_trainable_map.get(variant_id, None)
            if truth_indices is not None:
                truth_tensor[row, :].put(truth_indices, 1.0)

    @staticmethod
    def _pickle_batch_tensors(
            parquet_file: str,
            pickle_folder: str,
            previous_pickle_folder: Optional[str],
            wanted_columns: list[str],
            random_seed: Optional[int],
            variants_per_batch: int,
            num_samples: int,
            num_properties: int,
            variant_weights_column: int,
            variant_columns: numpy.ndarray,
            format_columns: numpy.ndarray,
            validation_proportion: float,
            trainable_variant_ids: set[str],
            variant_id_column: int
    ) -> str:
        t0 = time.time()
        # # read the basic buffer into a numpy array, only loading trainable variant IDs
        id_filter = (wanted_columns[variant_id_column], 'in', trainable_variant_ids)
        buffer = pandas.read_parquet(parquet_file, columns=wanted_columns, filters=[id_filter]).values
        t1 = time.time()
        # split new buffer up so that the first block of rows is for validation, remaining for training,
        # and randomly permute the indices within training / permutation blocks if desired
        validation_index = round(buffer.shape[0] * validation_proportion)
        if random_seed is not None:
            random_state = numpy.random.RandomState(random_seed)
            perm_validation = random_state.permutation(validation_index)
            perm_training = validation_index + random_state.permutation(buffer.shape[0] - validation_index)
        else:
            perm_validation = numpy.arange(validation_index)
            perm_training = numpy.arange(validation_index, buffer.shape[0])
        t2 = time.time()
        # load remainder sizes
        if previous_pickle_folder is None:
            num_previous_training_remainder, num_previous_validation_remainder = 0, 0
        else:
            previous_remainder_sizes = os.path.join(previous_pickle_folder, Keys.remainder_sizes)
            with locked_read(previous_remainder_sizes, "rb") as f_in:
                num_previous_training_remainder, num_previous_validation_remainder = pickle.load(f_in)

        # construct an iterator to handle divvying up the variants into interleaved
        # training/validation batches
        training_batch_iterator = TrainingBatchIterator(
            num_training_variants=(
                buffer.shape[0] - validation_index + num_previous_training_remainder
            ),
            num_validation_variants=validation_index + num_previous_validation_remainder,
            variants_per_batch=variants_per_batch,
            pickle_folder=pickle_folder,
            parquet_file=parquet_file,
            previous_training_remainder_ids=(),
            previous_validation_remainder_ids=()
        )
        # save iterator so that the main tensor data loader can also iterate over the batches
        training_batch_iterator.to_json()
        # pickle remainder variants (variants that are not enough for a complete batch) if any
        # do this FIRST, because now there will be enough information for another process to proceed
        # onto the next batch of variants
        # first quickly output just the sizes, because that can be very quick and is enough to
        # enable planning
        remainder_sizes = os.path.join(pickle_folder, Keys.remainder_sizes)
        with locked_write(remainder_sizes, "wb") as f_out:
            pickle.dump(
                (
                    training_batch_iterator.num_training_remainder,
                    training_batch_iterator.num_validation_remainder
                ),
                f_out
            )
        t3 = time.time()
        # load previous remainder buffers and pickle new remainder buffers
        training_remainder_buffer = VcfTrainingTensorDataLoader._handle_new_and_old_remainder_buffers(
            is_training=True, previous_pickle_folder=previous_pickle_folder,
            pickle_folder=pickle_folder, buffer=buffer,
            training_batch_iterator=training_batch_iterator, buffer_indices=perm_training
        )
        validation_remainder_buffer = VcfTrainingTensorDataLoader._handle_new_and_old_remainder_buffers(
            is_training=False, previous_pickle_folder=previous_pickle_folder,
            pickle_folder=pickle_folder, buffer=buffer,
            training_batch_iterator=training_batch_iterator, buffer_indices=perm_validation
        )
        t4 = time.time()
        # construct numpy arrays that will hold tensors for each batch
        numpy_tensor = numpy.empty(
            (variants_per_batch, num_samples, num_properties), dtype=numpy.float32
        )
        weights_tensor = numpy.empty(variants_per_batch, dtype=numpy.float32)
        ids_array = pandas.array([""] * variants_per_batch, dtype="string[pyarrow]")

        # iterate over batches in order, and pickle the batch tensors
        num_variant_properties = len(variant_columns)
        training_buffer_end = 0
        validation_buffer_end = 0
        for is_training_batch, previous_remainder_sizes, batch_number in training_batch_iterator:
            has_remainder = batch_number == 0 and (
                training_remainder_buffer.shape[0] > 0 if is_training_batch
                else validation_remainder_buffer.shape[0] > 0
            )
            if has_remainder:
                # first batch has to handle remainders
                if is_training_batch:
                    remainder_buffer = training_remainder_buffer
                    num_remainder = remainder_buffer.shape[0]
                    training_buffer_begin = training_buffer_end
                    training_buffer_end += (variants_per_batch - num_remainder)
                    indices = perm_training[training_buffer_begin:training_buffer_end]
                else:
                    remainder_buffer = validation_remainder_buffer
                    num_remainder = remainder_buffer.shape[0]
                    validation_buffer_begin = validation_buffer_end
                    validation_buffer_end += (variants_per_batch - num_remainder)
                    indices = perm_validation[validation_buffer_begin:validation_buffer_end]
                (
                    ids_array[:num_remainder],
                    weights_tensor[:num_remainder],
                    numpy_tensor[:num_remainder, :, :num_variant_properties],
                    numpy_tensor[:num_remainder, :, num_variant_properties:]
                ) = VcfTensorDataLoaderBase._get_next_tensor_blocks(
                    buffer_matrix=remainder_buffer, variant_ids_column=variant_id_column,
                    variant_weights_column=variant_weights_column, variant_columns=variant_columns,
                    format_columns=format_columns, num_samples=num_samples
                )
                (
                    ids_array[num_remainder:],
                    weights_tensor[num_remainder:],
                    numpy_tensor[num_remainder:, :, :num_variant_properties],
                    numpy_tensor[num_remainder:, :, num_variant_properties:]
                ) = VcfTensorDataLoaderBase._get_next_tensor_blocks(
                    buffer_matrix=buffer.take(indices, axis=0), variant_ids_column=variant_id_column,
                    variant_weights_column=variant_weights_column, variant_columns=variant_columns,
                    format_columns=format_columns, num_samples=num_samples
                )
            else:
                # past the remainders, just take the needed variants from the buffer
                if is_training_batch:
                    training_buffer_begin = training_buffer_end
                    training_buffer_end += variants_per_batch
                    indices = perm_training[training_buffer_begin:training_buffer_end]
                else:
                    validation_buffer_begin = validation_buffer_end
                    validation_buffer_end += variants_per_batch
                    indices = perm_validation[validation_buffer_begin:validation_buffer_end]
                (
                    ids_array[:],
                    weights_tensor[:],
                    numpy_tensor[:, :, :num_variant_properties],
                    numpy_tensor[:, :, num_variant_properties:]
                ) = VcfTensorDataLoaderBase._get_next_tensor_blocks(
                    buffer_matrix=buffer.take(indices, axis=0),
                    variant_ids_column=variant_id_column, variant_weights_column=variant_weights_column,
                    variant_columns=variant_columns, format_columns=format_columns, num_samples=num_samples
                )
            with locked_write(previous_remainder_sizes, mode="wb") as f_out:
                pickle.dump((ids_array, weights_tensor, numpy_tensor), f_out, protocol=pickle.HIGHEST_PROTOCOL)
        t5 = time.time()
        print(
            f"time to load buffer: {t1 - t0}, make permutations: {t2 - t1}, plan batches: {t3 - t2}, "
            f"load and pickle remainders: {t4 - t3}, pickle batches: {t5 - t4}"
        )

        return parquet_file

    @staticmethod
    def load_remainder_buffer(
            pickle_folder: Optional[str],
            is_training: bool,
            num_columns: int,
            delete_after_read: bool = Default.delete_after_read,
            sleep_time_s: float = Default.sleep_time_s
    ) -> numpy.ndarray:
        if pickle_folder is None:
            return numpy.empty((0, num_columns), dtype=numpy.float32)
        pickle_file = VcfTrainingTensorDataLoader.get_nth_pickle_file(pickle_folder, n=Keys.remainder,
                                                                      is_training=is_training)
        with locked_read(pickle_file, "rb", delete_after_read=delete_after_read, sleep_time_s=sleep_time_s) as f_in:
            buffer = pickle.load(f_in)
        return buffer

    def _get_pickle_folder(self) -> str:
        return tempfile.mkdtemp(dir=self._temp_dir)

    @staticmethod
    def get_nth_pickle_file(pickle_folder: str, n: Union[int, str], is_training: bool) -> str:
        return os.path.join(pickle_folder,
                            f"{Keys.training if is_training else Keys.validation}_{n}{Keys.pickle_suffix}")

    @staticmethod
    def get_next_batch(training_batch: int, validation_batch: int, validation_proportion: float) -> (bool, int, int):
        is_training_batch = training_batch * validation_proportion < \
                            (validation_batch + 1) * (1.0 - validation_proportion)
        if is_training_batch:
            training_batch += 1
        else:
            validation_batch += 1
        return is_training_batch, training_batch, validation_batch

    @staticmethod
    def _handle_new_and_old_remainder_buffers(
            is_training: bool,
            previous_pickle_folder: str,
            pickle_folder: str,
            buffer: numpy.ndarray,
            training_batch_iterator: TrainingBatchIterator,
            buffer_indices: numpy.ndarray,
            delete_after_read: bool = Default.delete_after_read,
            sleep_time_s: float = Default.sleep_time_s
    ) -> numpy.ndarray:
        num_remainder = training_batch_iterator.num_training_remainder if is_training \
            else training_batch_iterator.num_validation_remainder
        # get file name for new remainder
        pickle_file = VcfTrainingTensorDataLoader.get_nth_pickle_file(
            pickle_folder, n=Keys.remainder, is_training=is_training
        )
        if num_remainder > buffer_indices.size:
            # need to copy over previous remainder to next remainder, so need to load previous
            # remainder first (this is slower for parallelization, so avoid if possible)
            previous_remainder_buffer = VcfTrainingTensorDataLoader.load_remainder_buffer(
                pickle_folder=previous_pickle_folder, is_training=is_training,
                num_columns=buffer.shape[1], delete_after_read=delete_after_read,
                sleep_time_s=sleep_time_s
            )
            new_remainder = numpy.concatenate(
                (previous_remainder_buffer, buffer.take(buffer_indices, axis=0)),
                axis=0
            )
            with locked_write(pickle_file, mode="wb") as f_out:
                pickle.dump(new_remainder, f_out, protocol=pickle.HIGHEST_PROTOCOL)
        else:
            # new remainder is a subset of this buffer. we can immediately output the new remainder
            # then load previous remainder (which we may have to wait for)
            # don't use negative range index because if num_remainder == 0 you would get
            # *ALL* indices instead of none of them
            end = len(buffer_indices)
            begin = end - num_remainder  # don't use negative range index becaues
            with locked_write(pickle_file, mode="wb") as f_out:
                pickle.dump(
                    buffer.take(buffer_indices[begin:end], axis=0),
                    f_out,
                    protocol=pickle.HIGHEST_PROTOCOL
                )
            previous_remainder_buffer = VcfTrainingTensorDataLoader.load_remainder_buffer(
                pickle_folder=previous_pickle_folder, is_training=is_training,
                num_columns=buffer.shape[1], delete_after_read=delete_after_read,
                sleep_time_s=sleep_time_s
            )
        return previous_remainder_buffer

    def __len__(self) -> int:
        # report the number of training batches in the iterator
        n_batches, extra_rows = divmod(round((1.0 - self.validation_proportion) * (self._num_rows - self._current_row)),
                                       self.variants_per_batch)
        return n_batches + 1 if extra_rows else n_batches

    @property
    def state_dict(self) -> dict[str, Any]:
        return {key: getattr(self, key) for key in VcfTrainingTensorDataLoader.__state_keys__}

    def load_state_dict(self, state_dict: dict[str, Any]):
        for key in VcfTrainingTensorDataLoader.__state_keys__:
            setattr(self, key, state_dict.get(key))
