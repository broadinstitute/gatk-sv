import glob
import json
import concurrent.futures
import numpy
import pandas
import pandas.core.arrays
import torch
from typing import Optional, Union
from gq_recalibrator import tarred_properties_to_parquet


ArrowStringArray = pandas.core.arrays.string_arrow.ArrowStringArray


class Keys:
    cuda = "cuda"
    cpu = "cpu"
    id = tarred_properties_to_parquet.Keys.id
    variant_weights = tarred_properties_to_parquet.Keys.variant_weights
    positive_value = tarred_properties_to_parquet.Keys.positive_value
    negative_value = tarred_properties_to_parquet.Keys.negative_value
    baseline = tarred_properties_to_parquet.Keys.baseline
    scale = tarred_properties_to_parquet.Keys.scale


class Default:
    shuffle = True
    random_state = 0
    non_ml_properties = tarred_properties_to_parquet.non_ml_properties
    # don't exclude variant_weights: they're not fed into the neural net, but we use them in training:
    non_ml_needed_for_training = frozenset({Keys.id, Keys.variant_weights})
    excluded_properties = non_ml_properties.difference(non_ml_needed_for_training)
    scale_bool_properties = False,
    validation_proportion = 0.2


def get_torch_device(device_name: str) -> torch.device:
    return torch.device("cuda:0" if device_name == Keys.cuda and torch.cuda.is_available() else "cpu")


class VcfTensorDataLoaderBase:
    """
    A DataLoader-like object for a set of tensors that can be much faster than
    TensorDataset + DataLoader because dataloader grabs individual indices of
    the dataset and calls cat (slow).
    Source: https://discuss.pytorch.org/t/dataloader-much-slower-than-manual-batching/27014/6
    """
    __slots__ = (
        "parquet_path", "variants_per_batch", "torch_device", "shuffle", "random_state",
        "_current_parquet_file", "executor", "_wanted_columns", "_parquet_files", "_buffer_index", "_buffer_matrix",
        "_num_rows", "_current_row", "_variant_weights_column", "_variant_id_column", "_variant_columns",
        "_format_columns", "_sample_ids", "_num_properties", "_buffer_future", "_numpy_tensor", "_weights_tensor",
        "_ids_array", "_tensor_baseline", "_tensor_scale", "_training_rows", "_validation_rows"
    )

    def __init__(
            self,
            parquet_path: str,
            properties_scaling_json: str,
            executor: concurrent.futures.Executor,
            variants_per_batch: int,
            torch_device: str,
            shuffle: bool = Default.shuffle,
            scale_bool_values: bool = Default.scale_bool_properties,
            random_state: Union[int, numpy.random.RandomState, None] = Default.random_state,
            excluded_properties: set[str] = Default.excluded_properties
    ):
        self.parquet_path = tarred_properties_to_parquet.extract_tar_to_folder(parquet_path) \
            if parquet_path.endswith(".tar") else parquet_path
        self.executor = executor
        self.variants_per_batch = variants_per_batch
        self.torch_device = get_torch_device(torch_device)
        self.shuffle = shuffle
        self.random_state = random_state if isinstance(random_state, numpy.random.RandomState) \
            else numpy.random.RandomState(random_state)
        self._set_column_info(
            excluded_properties=excluded_properties, properties_scaling_json=properties_scaling_json,
            scale_bool_values=scale_bool_values
        )
        num_format_properties = len(self._format_columns) // self.num_samples
        self._numpy_tensor = numpy.empty(
            (self.variants_per_batch, self.num_samples, len(self._variant_columns) + num_format_properties),
            dtype=numpy.float32
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
        unflat_columns = [
            raw_column if raw_column == "row" else tarred_properties_to_parquet.unflatten_column_name(raw_column)
            for raw_column in raw_columns
        ]
        self._wanted_columns, unflat_columns = zip(*(
            (raw_column, unflat_column) for raw_column, unflat_column in zip(raw_columns, unflat_columns)
            if unflat_column[1] not in excluded_properties
        ))
        # remove "row", as it will be an index
        unflat_columns = [column for column in unflat_columns if not isinstance(column, str)]
        if Keys.variant_weights not in excluded_properties:
            self._variant_weights_column = next(
                index for index, column in enumerate(unflat_columns) if column[1] == Keys.variant_weights
            )
        self._variant_id_column = next(
            index for index, column in enumerate(unflat_columns) if column[1] == Keys.id
        )
        print(f"variant_id_column: {self._variant_id_column}")
        # get the DataFrame columns that correspond to per-variant / INFO columns
        self._variant_columns = numpy.array(
            [index for index, column in enumerate(unflat_columns)
             if not isinstance(column[0], str) and column[1] not in Default.non_ml_needed_for_training]
        )
        print(f"variant_columns: {self._variant_columns}")
        # get the DataFrame columns that correspond to per-sample / FORMAT columns
        self._format_columns = numpy.array(
            [index for index, column in enumerate(unflat_columns)
             if isinstance(column[0], str) and column[1] not in Default.non_ml_needed_for_training]
        )
        # get baseline and scale for final tensor
        properties_tensor_order = [
            column[1] for column in unflat_columns
            if not isinstance(column[0], str) and column[1] not in Default.non_ml_needed_for_training
        ]
        found_properties = set()
        found_sample_ids = set()
        self._sample_ids = []
        for sample_id, prop_name in unflat_columns:
            if prop_name in Default.non_ml_needed_for_training or not isinstance(sample_id, str):
                continue
            if sample_id not in found_sample_ids:
                self._sample_ids.append(sample_id)
            if prop_name not in found_properties:
                found_properties.add(prop_name)
                properties_tensor_order.append(prop_name)
        print(f"properties_tensor_order: {properties_tensor_order}")

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

        self._num_properties = len(properties_tensor_order)

        self._tensor_baseline = torch.tensor(
            [_get_baseline(prop_name) for prop_name in properties_tensor_order],
            dtype=torch.float32, device=self.torch_device
        )
        self._tensor_scale = torch.tensor(
            [1.0 / _get_scale(prop_name) for prop_name in properties_tensor_order],
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
    def num_properties(self) -> int:
        return self._num_properties

    def __len__(self) -> int:
        # report the number of remaining batches in the iterator
        n_batches, extra_rows = divmod(self._num_rows - self._current_row, self.variants_per_batch)
        return n_batches + 1 if extra_rows else n_batches


class VcfFilterTensorDataLoader(VcfTensorDataLoaderBase):
    def __init__(
            self,
            parquet_path: str,
            properties_scaling_json: str,
            executor: concurrent.futures.Executor,
            variants_per_batch: int,
            torch_device: str,
            shuffle: bool = Default.shuffle,
            scale_bool_values: bool = Default.scale_bool_properties,
            random_state: Union[int, numpy.random.RandomState, None] = Default.random_state,
            excluded_properties: set[str] = Default.excluded_properties
    ):
        # basically it's the VcfTensorDataLoaderBase with variant weights property excluded and explicitly no
        # validation
        super().__init__(
            parquet_path=parquet_path, properties_scaling_json=properties_scaling_json, executor=executor,
            variants_per_batch=variants_per_batch, torch_device=torch_device, shuffle=shuffle,
            scale_bool_values=scale_bool_values, random_state=random_state,
            excluded_properties=excluded_properties.union((Keys.variant_weights,))
        )

    def __iter__(self):
        if self._buffer_future is None:
            # only executes on first iteration, so there's no remainder rows
            if not self._parquet_files:
                self._set_parquet_files()
            self._current_parquet_file = self._parquet_files.pop()
            self._buffer_future = self.executor.submit(
                self.__class__._parquet_file_to_df_buffer, parquet_file=self._current_parquet_file,
                wanted_columns=self._wanted_columns, random_state=self.random_state if self.shuffle else None
            )
        return self

    @staticmethod
    def _parquet_file_to_df_buffer(
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
        self._buffer_future = self.executor.submit(
            self.__class__._parquet_file_to_df_buffer, parquet_file=self._current_parquet_file,
            wanted_columns=self._wanted_columns, random_state=self.random_state if self.shuffle else None,
            remainder_buffer=self._buffer_matrix[-remainder_rows:, :]
        )


class VcfTrainingTensorDataLoader(VcfTensorDataLoaderBase):
    __slots__ = ("validation_proportion", "_validation_end", "_validation_buffer_index",
                 "_training_batch", "_validation_batch", "_validation_tensor", "_validation_weights_tensor",
                 "_validation_ids_array")

    def __init__(
            self,
            parquet_path: str,
            properties_scaling_json: str,
            executor: concurrent.futures.Executor,
            variants_per_batch: int,
            torch_device: str,
            shuffle: bool = Default.shuffle,
            scale_bool_values: bool = Default.scale_bool_properties,
            random_state: Union[int, numpy.random.RandomState, None] = Default.random_state,
            excluded_properties: set[str] = Default.excluded_properties,
            validation_proportion: float = Default.validation_proportion
    ):
        # basically it's the VcfTensorDataLoaderBase with a train/test split
        super().__init__(
            parquet_path=parquet_path, properties_scaling_json=properties_scaling_json, executor=executor,
            variants_per_batch=variants_per_batch, torch_device=torch_device, shuffle=shuffle,
            scale_bool_values=scale_bool_values, random_state=random_state, excluded_properties=excluded_properties
        )
        if validation_proportion <= 0 or validation_proportion >= 0.5:
            raise ValueError("Validation proportion must be in open interval (0, 0.5)")
        self.validation_proportion = validation_proportion
        self._validation_end = -1
        self._validation_buffer_index = 0
        self._training_batch = 0
        self._validation_batch = 0
        self._validation_tensor = numpy.empty(self._numpy_tensor.shape, self._numpy_tensor.dtype)
        self._validation_weights_tensor = numpy.empty(self._weights_tensor.shape, self._weights_tensor.dtype)
        self._validation_ids_array = self._ids_array.copy()

    def __iter__(self):
        if self._buffer_future is None:
            # only executes on first iteration, so there's no remainder rows
            if not self._parquet_files:
                self._set_parquet_files()
            self._current_parquet_file = self._parquet_files.pop()
            self._buffer_future = self.executor.submit(
                self.__class__._parquet_file_to_df_buffer, parquet_file=self._current_parquet_file,
                wanted_columns=self._wanted_columns, random_state=self.random_state if self.shuffle else None,
                validation_proportion=self.validation_proportion
            )
        return self

    def _load_buffer(self):
        self._buffer_matrix, self.random_state, self._validation_end = self._buffer_future.result()
        self._buffer_index = self._validation_end
        self._validation_buffer_index = 0
        if not self._parquet_files:
            self._set_parquet_files()
        self._current_parquet_file = self._parquet_files.pop()
        remainder_rows = divmod(self._buffer_matrix.shape[0] - self._validation_end, self.variants_per_batch)[1]
        remainder_validation_rows = divmod(self._validation_end, self.variants_per_batch)[1]
        self._buffer_future = self.executor.submit(
            self.__class__._parquet_file_to_df_buffer, parquet_file=self._current_parquet_file,
            wanted_columns=self._wanted_columns, random_state=self.random_state if self.shuffle else None,
            validation_proportion=self.validation_proportion,
            remainder_buffer=self._buffer_matrix[-remainder_rows:, :],
            validation_remainder_buffer=self._buffer_matrix[
                                            self._validation_end - remainder_validation_rows:self._validation_end, :
                                        ]
        )

    @staticmethod
    def phred_to_p(phred_matrix: numpy.ndarray) -> numpy.ndarray:
        phred_coef = numpy.log(10.0) / 10.0
        return 1.0 - 0.5 * numpy.exp(phred_coef * (1.0 - phred_matrix))


    @staticmethod
    def _parquet_file_to_df_buffer(
            parquet_file: str,
            wanted_columns: list[str],
            random_state: Optional[numpy.random.RandomState],
            validation_proportion: float,
            remainder_buffer: Optional[numpy.ndarray] = None,
            validation_remainder_buffer: Optional[numpy.ndarray] = None
    ) -> (numpy.ndarray, Optional[numpy.random.RandomState], int):
        buffer = pandas.read_parquet(parquet_file, columns=wanted_columns).values
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

    def __next__(self) -> (torch.Tensor, torch.Tensor, ArrowStringArray,
                           Optional[torch.Tensor], Optional[torch.Tensor], Optional[ArrowStringArray]):
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
        next_buffer_index = self._buffer_index + self.variants_per_batch
        if next_buffer_index > self._buffer_matrix.shape[0]:
            self._load_buffer()
            next_buffer_index = self._buffer_index + self.variants_per_batch

        num_variant_properties = len(self._variant_columns)
        (self._ids_array[:], self._weights_tensor[:], self._numpy_tensor[:, :, :num_variant_properties],
         self._numpy_tensor[:, :, num_variant_properties:]) = VcfTensorDataLoaderBase._get_next_tensor_blocks(
            buffer_matrix=self._buffer_matrix[self._buffer_index:next_buffer_index, :],
            variant_ids_column=self._variant_id_column, variant_weights_column=self._variant_weights_column,
            variant_columns=self._variant_columns, format_columns=self._format_columns, num_samples=self.num_samples
        )
        self._training_batch += 1
        self._buffer_index = next_buffer_index
        self._current_row += self.variants_per_batch
        if self._training_batch * self.validation_proportion >= \
                (self._validation_batch + 1) * (1.0 - self.validation_proportion):
            # time for the next validation batch
            next_validation_buffer_index = self._validation_buffer_index + self.variants_per_batch
            if next_validation_buffer_index > self._validation_end:
                self._load_buffer()
                next_validation_buffer_index = self._validation_buffer_index + self.variants_per_batch
            (self._validation_ids_array[:], self._validation_weights_tensor[:],
             self._validation_tensor[:, :, :num_variant_properties],
             self._validation_tensor[:, :, num_variant_properties:]) = \
                VcfTensorDataLoaderBase._get_next_tensor_blocks(
                    buffer_matrix=self._buffer_matrix[self._validation_buffer_index:next_validation_buffer_index, :],
                    variant_ids_column=self._variant_id_column, variant_weights_column=self._variant_weights_column,
                    variant_columns=self._variant_columns, format_columns=self._format_columns,
                    num_samples=self.num_samples
                )
            self._validation_batch += 1
            self._current_row += self.variants_per_batch
            self._validation_buffer_index = next_validation_buffer_index
            return (
                self._tensor_scale * (
                    torch.tensor(self._numpy_tensor, dtype=torch.float32, device=self.torch_device)
                    - self._tensor_baseline
                ),
                torch.tensor(self._weights_tensor, dtype=torch.float32, device=self.torch_device),
                self._ids_array,
                self._tensor_scale * (
                        torch.tensor(self._validation_tensor, dtype=torch.float32, device=self.torch_device)
                        - self._tensor_baseline
                ),
                torch.tensor(self._validation_weights_tensor, dtype=torch.float32, device=self.torch_device),
                self._validation_ids_array
            )
        else:
            # not time for a validation batch
            return (
                self._tensor_scale * (
                    torch.tensor(self._numpy_tensor, dtype=torch.float32, device=self.torch_device)
                    - self._tensor_baseline
                ),
                torch.tensor(self._weights_tensor, dtype=torch.float32, device=self.torch_device),
                self._ids_array,
                None,
                None,
                None
            )

    def __len__(self) -> int:
        # report the number of training batches in the iterator
        n_batches, extra_rows = divmod(round((1.0 - self.validation_proportion) * (self._num_rows - self._current_row)),
                                       self.variants_per_batch)
        return n_batches + 1 if extra_rows else n_batches
