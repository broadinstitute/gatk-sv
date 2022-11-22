import glob
import json
import concurrent.futures
import numpy
import pandas
import pandas.core.arrays
import torch
from collections import defaultdict
from typing import Optional, Union, Any
from gq_recalibrator import tarred_properties_to_parquet, training_utils
from sv_utils import get_truth_overlap, benchmark_variant_filter


ArrowStringArray = pandas.core.arrays.string_arrow.ArrowStringArray
ConfidentVariants = get_truth_overlap.ConfidentVariants


class Keys:
    cuda = "cuda"
    cpu = "cpu"
    id = tarred_properties_to_parquet.Keys.id
    variant_weights = tarred_properties_to_parquet.Keys.variant_weights
    positive_value = tarred_properties_to_parquet.Keys.positive_value
    negative_value = tarred_properties_to_parquet.Keys.negative_value
    baseline = tarred_properties_to_parquet.Keys.baseline
    scale = tarred_properties_to_parquet.Keys.scale
    gq = benchmark_variant_filter.Keys.gq


class Default:
    shuffle = True
    random_state = 0
    non_ml_properties = tarred_properties_to_parquet.non_ml_properties
    # don't exclude variant_weights: they're not fed into the neural net, but we use them in training:
    non_ml_needed_for_training = frozenset({Keys.id, Keys.variant_weights})
    excluded_properties = non_ml_properties.difference(non_ml_needed_for_training)
    scale_bool_properties = False,
    validation_proportion = 0.2


def get_torch_device(device_name: str, progress_logger: training_utils.ProgressLogger) -> torch.device:
    if device_name == Keys.cuda:
        if torch.cuda.is_available():
            progress_logger.log("Using cuda for torch")
            torch.cuda.empty_cache()
            return torch.device("cuda:0")
        else:
            progress_logger.log("Cuda not available, falling back to cpu for torch")
            return torch.device("cpu")
    else:
        progress_logger.log("Using cpu for torch")
        return torch.device("cpu")


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
        "_buffer_index", "_buffer_matrix",
        "_num_rows", "_current_row", "_columns", "_variant_weights_column", "_variant_id_column", "_variant_columns",
        "_format_columns", "_sample_ids", "_property_names", "_buffer_future", "_numpy_tensor", "_weights_tensor",
        "_ids_array", "_tensor_baseline", "_tensor_scale", "_training_rows", "_validation_rows"
    )

    def __init__(
            self,
            parquet_path: str,
            properties_scaling_json: str,
            process_executor: concurrent.futures.Executor,
            thread_executor: concurrent.futures.Executor,
            variants_per_batch: int,
            torch_device: str,
            progress_logger: training_utils.ProgressLogger,
            shuffle: bool = Default.shuffle,
            scale_bool_values: bool = Default.scale_bool_properties,
            random_state: Union[int, numpy.random.RandomState, None] = Default.random_state,
            excluded_properties: set[str] = Default.excluded_properties
    ):
        self.parquet_path = tarred_properties_to_parquet.extract_tar_to_folder(parquet_path) \
            if parquet_path.endswith(".tar") else parquet_path
        self.process_executor = process_executor
        self.thread_executor = thread_executor
        self.variants_per_batch = variants_per_batch
        self.progress_logger = progress_logger
        self.torch_device = get_torch_device(torch_device, self.progress_logger)
        self.shuffle = shuffle
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
        self._current_parquet_file = None
        self._buffer_future = None

    def _set_column_info(self, excluded_properties: set[str], properties_scaling_json: str, scale_bool_values: bool):
        """
        Set a bunch of column info so that it can be quickly retrieved later
          -self._wanted_columns is the list of raw (str) columns saved in the parquet file that we want to load, i.e.
           everything but the excluded properties
          -self._variant_columns is the unflattened (i.e. MultiIndex) columns that correspond to per-variant properties
          -self._format_columns is the unflattened columns that correspond to per-genotype / format properties
        Have a few annoying issues to work around:
          -Want "row" in wanted_columns, because it will be the index, but don't want it in variant_columns or
           format_columns because it's not in the columns of the final dataframe. Also it can't be unflattened
          -Can only determine if a column is for an excluded property after unflattening
        """
        with open(properties_scaling_json, 'r') as f_in:
            properties_scaling = json.load(f_in)

        raw_columns = tarred_properties_to_parquet.get_parquet_file_columns(self.parquet_path)
        self._columns = [
            raw_column if raw_column == "row" else tarred_properties_to_parquet.unflatten_column_name(raw_column)
            for raw_column in raw_columns
        ]
        self._wanted_columns, self._columns = zip(*(
            (raw_column, unflat_column) for raw_column, unflat_column in zip(raw_columns, self._columns)
            if unflat_column[1] not in excluded_properties
        ))
        # remove "row", as it will be an index
        self._columns = [column for column in self._columns if not isinstance(column, str)]
        if Keys.variant_weights not in excluded_properties:
            self._variant_weights_column = next(
                index for index, column in enumerate(self._columns) if column[1] == Keys.variant_weights
            )
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
        # get baseline and scale for final tensor
        self._property_names = [
            column[1] for column in self._columns
            if not isinstance(column[0], str) and column[1] not in Default.non_ml_needed_for_training
        ]
        found_properties = set()
        found_sample_ids = set()
        self._sample_ids = []
        for sample_id, prop_name in self._columns:
            if prop_name in Default.non_ml_needed_for_training or not isinstance(sample_id, str):
                continue
            if sample_id not in found_sample_ids:
                found_sample_ids.add(sample_id)
                self._sample_ids.append(sample_id)
            if prop_name not in found_properties:
                found_properties.add(prop_name)
                self._property_names.append(prop_name)

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
                .reshape(buffer_matrix.shape[0], len(format_columns) // num_samples, num_samples)
                .swapaxes(1, 2)
        )  # noqa E131

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
            process_executor: concurrent.futures.Executor,
            thread_executor: concurrent.futures.Executor,
            variants_per_batch: int,
            torch_device: str,
            progress_logger: training_utils.ProgressLogger,
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
            shuffle=shuffle, scale_bool_values=scale_bool_values, random_state=random_state,
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


class VcfTrainingTensorDataLoader(VcfTensorDataLoaderBase):
    __slots__ = (
        "validation_proportion", "_validation_end", "_validation_buffer_index", "_training_batch", "_validation_batch",
        "_genotype_is_good_dict", "_genotype_is_bad_dict", "_genotype_is_good_tensor", "_genotype_is_bad_tensor",
        "_trainable_variant_ids", "_gq_columns", "_gq_tensor"
    )
    __state_keys__ = (
        "random_state", "_current_parquet_file", "_parquet_files", "_buffer_index", "_current_row", "_validation_end",
        "_validation_buffer_index", "_training_batch", "_validation_batch"
    )

    def __init__(
            self,
            parquet_path: str,
            properties_scaling_json: str,
            truth_json: str,
            process_executor: concurrent.futures.Executor,
            thread_executor: concurrent.futures.Executor,
            variants_per_batch: int,
            torch_device: str,
            progress_logger: training_utils.ProgressLogger,
            shuffle: bool = Default.shuffle,
            scale_bool_values: bool = Default.scale_bool_properties,
            random_state: Union[int, numpy.random.RandomState, None] = Default.random_state,
            excluded_properties: set[str] = Default.excluded_properties,
            validation_proportion: float = Default.validation_proportion
    ):
        # basically it's the VcfTensorDataLoaderBase with a train/test split
        super().__init__(
            parquet_path=parquet_path, properties_scaling_json=properties_scaling_json,
            process_executor=process_executor, thread_executor=thread_executor,
            variants_per_batch=variants_per_batch, torch_device=torch_device, progress_logger=progress_logger,
            shuffle=shuffle, scale_bool_values=scale_bool_values, random_state=random_state,
            excluded_properties=excluded_properties
        )
        if validation_proportion <= 0 or validation_proportion >= 1:
            raise ValueError("Validation proportion must be in open interval (0, 1)")
        self.validation_proportion = validation_proportion
        self._validation_end = -1
        self._validation_buffer_index = 0
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
        self._gq_columns = numpy.array(
            [index for index, column in enumerate(self._columns)
             if isinstance(column[0], str) and column[1] == Keys.gq]
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
        if self._buffer_future is None:
            # only executes on first iteration, so there's no remainder rows
            if not self._parquet_files:
                self._set_parquet_files()
            if self._current_parquet_file is None:
                self._current_parquet_file = self._parquet_files.pop()
            self._buffer_future = self.process_executor.submit(
                VcfTrainingTensorDataLoader._parquet_file_to_buffer, parquet_file=self._current_parquet_file,
                wanted_columns=self._wanted_columns, random_state=self.random_state if self.shuffle else None,
                trainable_variant_ids=self._trainable_variant_ids, variant_id_column=self._variant_id_column,
                validation_proportion=self.validation_proportion
            )
            # self._buffer_future = self.process_executor.submit(
            #     VcfTrainingTensorDataLoader._load_parquet, parquet_file=self._current_parquet_file,
            #     wanted_columns=self._wanted_columns
            # )
        return self

    def _load_buffer(self):
        self._buffer_matrix, self.random_state, self._validation_end = self._buffer_future.result()

        # if self._buffer_matrix.shape[0] > 0:
        #     remainder_rows = divmod(self._buffer_matrix.shape[0] - self._validation_end, self.variants_per_batch)[1]
        #     remainder_validation_rows = divmod(self._validation_end, self.variants_per_batch)[1]
        #     remainder_buffer = self._buffer_matrix[-remainder_rows:, :]
        #     validation_remainder_buffer = self._buffer_matrix[
        #                                      self._validation_end - remainder_validation_rows:self._validation_end, :
        #                                  ]
        # else:
        #     remainder_buffer, validation_remainder_buffer = None, None
        # self._buffer_matrix, self.random_state, self._validation_end = VcfTrainingTensorDataLoader._form_next_buffer(
        #     buffer=self._buffer_future.result(), random_state=self.random_state if self.shuffle else None,
        #     trainable_variant_ids = self._trainable_variant_ids, variant_id_column=self._variant_id_column,
        #     validation_proportion=self.validation_proportion, remainder_buffer=None, validation_remainder_buffer=None
        # )
        self._buffer_index = self._validation_end
        self._validation_buffer_index = 0
        if not self._parquet_files:
            self._set_parquet_files()
        self._current_parquet_file = self._parquet_files.pop()
        remainder_rows = divmod(self._buffer_matrix.shape[0] - self._validation_end, self.variants_per_batch)[1]
        remainder_validation_rows = divmod(self._validation_end, self.variants_per_batch)[1]
        self._buffer_future = self.process_executor.submit(
            VcfTrainingTensorDataLoader._parquet_file_to_buffer, parquet_file=self._current_parquet_file,
            wanted_columns=self._wanted_columns, random_state=self.random_state if self.shuffle else None,
            trainable_variant_ids=self._trainable_variant_ids, variant_id_column=self._variant_id_column,
            validation_proportion=self.validation_proportion, remainder_buffer=self._buffer_matrix[-remainder_rows:, :],
            validation_remainder_buffer=self._buffer_matrix[
                                            self._validation_end - remainder_validation_rows:self._validation_end, :
                                        ]
        )
        # self._buffer_future = self.process_executor.submit(
        #     VcfTrainingTensorDataLoader._load_parquet, parquet_file=self._current_parquet_file,
        #     wanted_columns=self._wanted_columns
        # )

    @staticmethod
    def _load_parquet(
            parquet_file: str,
            wanted_columns: list[str],
    ) -> numpy.ndarray:
        return pandas.read_parquet(parquet_file, columns=wanted_columns).values

    @staticmethod
    def _form_next_buffer(
            buffer: numpy.ndarray,
            random_state: Optional[numpy.random.RandomState],
            validation_proportion: float,
            trainable_variant_ids: frozenset[str],
            variant_id_column: int,
            remainder_buffer: Optional[numpy.ndarray] = None,
            validation_remainder_buffer: Optional[numpy.ndarray] = None,
    ):
        # only take the trainable variants
        buffer = buffer.compress(
            [variant_id in trainable_variant_ids for variant_id in buffer.take(variant_id_column, axis=1)],
            axis=0
        )
        # split new buffer up so that the first section is for validation
        validation_index = round(buffer.shape[0] * validation_proportion)
        if remainder_buffer is None:
            # no remainder buffers
            if random_state is not None:
                # permute indices without mixing validation and training
                perm_validation = random_state.permutation(validation_index)
                perm_training = validation_index + random_state.permutation(buffer.shape[0] - validation_index)
                buffer = buffer.take(numpy.concatenate((perm_validation, perm_training)), axis=0)
            return buffer, random_state, validation_index
        else:
            # copy in remainder buffers after scrambling
            if random_state is None:
                buffer = numpy.concatenate((
                    validation_remainder_buffer,
                    buffer[:validation_index, :],
                    remainder_buffer,
                    buffer[validation_index:, :]
                ))
            else:
                perm_validation = random_state.permutation(validation_index)
                perm_training = validation_index + random_state.permutation(buffer.shape[0] - validation_index)
                buffer = numpy.concatenate((
                    validation_remainder_buffer,
                    buffer.take(perm_validation, axis=0),
                    remainder_buffer,
                    buffer.take(perm_training, axis=0)
                ))
            return buffer, random_state, validation_remainder_buffer.shape[0] + validation_index

    @staticmethod
    def _parquet_file_to_buffer(
            parquet_file: str,
            wanted_columns: list[str],
            random_state: Optional[numpy.random.RandomState],
            validation_proportion: float,
            trainable_variant_ids: set[str],
            variant_id_column: int,
            remainder_buffer: Optional[numpy.ndarray] = None,
            validation_remainder_buffer: Optional[numpy.ndarray] = None,
    ) -> (numpy.ndarray, Optional[numpy.random.RandomState], int):
        buffer = pandas.read_parquet(parquet_file, columns=wanted_columns).values
        # only take the trainable variants
        buffer = buffer.compress(
            [variant_id in trainable_variant_ids for variant_id in buffer.take(variant_id_column, axis=1)],
            axis=0
        )
        # split new buffer up so that the first section is for validation
        validation_index = round(buffer.shape[0] * validation_proportion)
        if remainder_buffer is None:
            # no remainder buffers
            if random_state is not None:
                # permute indices without mixing validation and training
                perm_validation = random_state.permutation(validation_index)
                perm_training = validation_index + random_state.permutation(buffer.shape[0] - validation_index)
                buffer = buffer.take(numpy.concatenate((perm_validation, perm_training)), axis=0)
            return buffer, random_state, validation_index
        else:
            # copy in remainder buffers after scrambling
            if random_state is None:
                buffer = numpy.concatenate((
                    validation_remainder_buffer,
                    buffer[:validation_index, :],
                    remainder_buffer,
                    buffer[validation_index:, :]
                ))
            else:
                perm_validation = random_state.permutation(validation_index)
                perm_training = validation_index + random_state.permutation(buffer.shape[0] - validation_index)
                buffer = numpy.concatenate((
                    validation_remainder_buffer,
                    buffer.take(perm_validation, axis=0),
                    remainder_buffer,
                    buffer.take(perm_training, axis=0)
                ))
            return buffer, random_state, validation_remainder_buffer.shape[0] + validation_index

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

        is_training_batch = self._training_batch * self.validation_proportion <= \
            (self._validation_batch + 1) * (1.0 - self.validation_proportion)

        num_variant_properties = len(self._variant_columns)
        if is_training_batch:
            # next batch to emit is for training
            buffer_begin = self._buffer_index
            buffer_end = buffer_begin + self.variants_per_batch
            if buffer_end > self._buffer_matrix.shape[0]:
                self._load_buffer()
                buffer_begin = self._buffer_index
                buffer_end = buffer_begin + self.variants_per_batch
            self._buffer_index = buffer_end
            self._training_batch += 1
        else:
            # next batch to emit is for validation
            buffer_begin = self._validation_buffer_index
            buffer_end = buffer_begin + self.variants_per_batch
            if buffer_end > self._buffer_matrix.shape[0]:
                self._load_buffer()
                buffer_begin = self._validation_buffer_index
                buffer_end = buffer_begin + self.variants_per_batch
            self._validation_buffer_index = buffer_end
            self._validation_batch += 1
        self._current_row += self.variants_per_batch

        (self._ids_array[:], self._weights_tensor[:], self._numpy_tensor[:, :, :num_variant_properties],
         self._numpy_tensor[:, :, num_variant_properties:]) = VcfTensorDataLoaderBase._get_next_tensor_blocks(
            buffer_matrix=self._buffer_matrix[buffer_begin:buffer_end, :],
            variant_ids_column=self._variant_id_column, variant_weights_column=self._variant_weights_column,
            variant_columns=self._variant_columns, format_columns=self._format_columns, num_samples=self.num_samples
        )
        self._gq_tensor[:] = self._buffer_matrix[buffer_begin:buffer_end, :].take(self._gq_columns, axis=1)
        VcfTrainingTensorDataLoader._fill_truth_tensor(self._ids_array, self._genotype_is_good_dict,
                                                       self._genotype_is_good_tensor)
        VcfTrainingTensorDataLoader._fill_truth_tensor(self._ids_array, self._genotype_is_bad_dict,
                                                       self._genotype_is_bad_tensor)

        return (
            is_training_batch,
            self._tensor_scale * (
                torch.tensor(self._numpy_tensor, dtype=torch.float32, device=self.torch_device)
                - self._tensor_baseline
            ).nan_to_num(),
            torch.tensor(self._weights_tensor, dtype=torch.float32, device=self.torch_device),
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
