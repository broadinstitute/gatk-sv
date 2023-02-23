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
import pandas.core.dtypes.api
import torch
from typing import Optional, Any, IO, TypeVar
from collections.abc import Iterator, Collection, Iterable
from gq_recalibrator import tarred_properties_to_parquet, training_utils
from sv_utils import common, get_truth_overlap, genomics_io

ArrowStringArray = pandas.core.arrays.string_arrow.ArrowStringArray
ConfidentVariants = get_truth_overlap.ConfidentVariants
DType = numpy.dtype | pandas.api.extensions.ExtensionDtype
PandasArray = pandas.core.arrays.PandasArray
ArrayType = TypeVar("ArrayType", numpy.ndarray, PandasArray, torch.Tensor)

_debug_file_locks = False


class Keys:
    row = tarred_properties_to_parquet.Keys.row
    positive_value = tarred_properties_to_parquet.Keys.positive_value
    negative_value = tarred_properties_to_parquet.Keys.negative_value
    baseline = tarred_properties_to_parquet.Keys.baseline
    scale = tarred_properties_to_parquet.Keys.scale
    id = genomics_io.Keys.id
    gq = genomics_io.Keys.gq
    allele_count = genomics_io.Keys.allele_count
    pickle_suffix = ".pickle"
    lock_suffix = ".lock"
    parquet_suffix = tarred_properties_to_parquet.Keys.parquet_suffix


class Default:
    scale_bool_properties = False
    temp_dir = Path(tempfile.gettempdir())
    sleep_time_s = 0.1
    delete_after_read = True
    max_look_ahead_batches = 3
    non_filtration_properties = tarred_properties_to_parquet.non_ml_properties


def _get_lock_path(file_path: Path) -> Path:
    """ Unified function for getting lock-file Path that corresponds to a target path """
    return file_path.with_suffix(file_path.suffix + Keys.lock_suffix)


if _debug_file_locks:
    # Use more complex functions that print out information about the locking/unlocking
    # process, only to be used if there is something going wrong with parallel pickling
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


def _get_baseline(
        property_name: str,
        properties_scaling: dict[str, dict[str, float]],
        scale_bool_values: bool
) -> float:
    """get baseline to (approximately) z-score output tensor"""
    if scale_bool_values or Keys.positive_value not in properties_scaling[property_name]:
        return properties_scaling[property_name][Keys.baseline]
    else:
        return 0.0


def _get_scale(
        property_name: str,
        properties_scaling: dict[str, dict[str, float]],
        scale_bool_values: bool
) -> float:
    """get multiplicative scale to (approximately) z-score output tensor"""
    if scale_bool_values or Keys.positive_value not in properties_scaling[property_name]:
        return 1.0 / properties_scaling[property_name][Keys.scale]
    else:
        return 1.0


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


class SupplementalPropertyGetter:
    """Holds info and methods for getting supplemental properties: properties that are needed for
    filtering and/or training but may not be fed into the neural net. e.g. variant id,
    variant weights, etc.
    """
    __slots__ = (
        "tensor_ind", "buffer_ind", "scale", "baseline", "dtype", "_is_variant", "_is_format"
    )

    def __init__(
            self,
            tensor_ind: Optional[torch.Tensor],
            buffer_ind: int | list[int],
            scale: Optional[torch.Tensor],
            baseline: Optional[torch.Tensor],
            dtype: DType
    ):
        self.tensor_ind = tensor_ind
        self.buffer_ind = buffer_ind
        self.scale = scale
        self.baseline = baseline
        self.dtype = dtype
        assert isinstance(buffer_ind, (int, list)), f"Bad buffer_ind type: {type(buffer_ind)}"
        self._is_variant = isinstance(buffer_ind, int)
        self._is_format = isinstance(buffer_ind, list)

    @staticmethod
    def build(
        property_name: str,
        unflatted_columns: list[tuple[Optional[str], str]],
        tensor_property_names: list[str],
        properties_scaling: dict[str, dict[str, float]],
        scale_bool_values: bool,
        dtypes: dict[str, DType]
    ) -> "SupplementalPropertyGetter":
        buffer_inds, sample_ids = zip(
            *(
                (ind, sample_id)
                for ind, (sample_id, _property_name) in enumerate(unflatted_columns)
                if _property_name == property_name
            )
        )
        if len(sample_ids) > 1 or sample_ids[0] is not None:
            # it's a format property
            buffer_ind = list(buffer_inds)
        else:
            # it's a variant property
            buffer_ind = buffer_inds[0]
        try:
            # if we can find it in the list of tensor properties, we can get it that way
            tensor_ind = torch.tensor(
                tensor_property_names.index(property_name), dtype=torch.int64
            )
        except ValueError:
            # we'll have to get it from the buffer and pass it back
            tensor_ind = None
        parquet_dtype = next(
            _dtype for (_sample_id, _property_name), _dtype in dtypes.items()
            if property_name == _property_name
        )
        if pandas.api.types.is_bool_dtype(parquet_dtype):
            dtype = numpy.dtype(bool)
        elif pandas.api.types.is_numeric_dtype(parquet_dtype):
            # it's a numeric property
            dtype = numpy.dtype(numpy.int64) if pandas.api.types.is_integer_dtype(parquet_dtype) \
                else numpy.dtype(numpy.float32)
        else:
            # it's a string
            dtype = pandas.core.dtypes.api.pandas_dtype("string[pyarrow]")
        if property_name in properties_scaling:
            scale = torch.tensor(
                _get_scale(
                    property_name=property_name,
                    properties_scaling=properties_scaling,
                    scale_bool_values=scale_bool_values
                ),
                dtype=torch.float32
            )
            baseline = torch.tensor(
                _get_baseline(
                    property_name=property_name,
                    properties_scaling=properties_scaling,
                    scale_bool_values=scale_bool_values
                ),
                dtype=torch.float32
            )
        else:
            scale = None
            baseline = None

        return SupplementalPropertyGetter(
            tensor_ind=tensor_ind,
            buffer_ind=buffer_ind,
            scale=scale,
            baseline=baseline,
            dtype=dtype
        )

    def to(self, device: torch.device) -> "SupplementalPropertyGetter":
        return SupplementalPropertyGetter(
            tensor_ind=None if self.tensor_ind is None else self.tensor_ind.detach().to(device),
            buffer_ind=self.buffer_ind,
            scale=None if self.scale is None else self.scale.detach().to(device),
            baseline=None if self.baseline is None else self.baseline.detach().to(device),
            dtype=self.dtype
        )

    @property
    def is_buffer(self) -> bool:
        return self.tensor_ind is None

    @property
    def is_tensor(self) -> bool:
        return self.tensor_ind is not None

    @property
    def is_variant(self) -> bool:
        return self._is_variant

    @property
    def is_format(self) -> bool:
        return self._is_format

    def get_from_tensor(self, tensor: ArrayType) -> ArrayType:
        return tensor[:, :, self.tensor_ind] if self._is_format else tensor[:, 0, self.tensor_ind]

    def get_from_buffer(self, buffer: ArrayType) -> ArrayType:
        return buffer.take(self.buffer_ind, axis=1)

    def scale_array(self, array: numpy.ndarray, device: torch.device) -> torch.tensor:
        return self.scale * (
            torch.tensor(array, dtype=torch.float32, device=device) - self.baseline
        )

    def get_empty_array(self, variants_per_batch: int) -> Optional[PandasArray | numpy.ndarray]:
        if self.is_tensor:
            return None
        shape: tuple[int, ...] = (variants_per_batch,) if self._is_variant \
            else (variants_per_batch, len(self.buffer_ind))
        if isinstance(self.dtype, numpy.dtype):
            return numpy.empty(shape, dtype=self.dtype)
        else:
            return pandas.array(
                numpy.full((int(numpy.prod(shape)),), "", dtype="object"),
                dtype=self.dtype
            )


class VcfTensorDataLoaderBase:
    """
    A DataLoader-like object for a set of tensors that can be much faster than
    TensorDataset + DataLoader because dataloader grabs individual indices of
    the dataset and calls cat (slow).
    Source: https://discuss.pytorch.org/t/dataloader-much-slower-than-manual-batching/27014/6
    """
    __slots__ = (
        "parquet_path", "variants_per_batch", "torch_device", "progress_logger", "shuffle",
        "keep_multiallelic", "keep_homref", "keep_homvar",
        "random_generator", "_current_parquet_file", "process_executor",
        "_wanted_columns", "_parquet_files", "_temp_dir",
        "_num_rows", "_current_row", "_columns", "_variant_columns", "_format_columns",
        "_supplemental_property_getters", "_sample_ids", "_property_names",
        "_tensor_baseline", "_tensor_scale", "_training_rows", "_validation_rows",
        "_variant_id_flat_column_name"
    )

    def __init__(
            self,
            parquet_path: Path,
            properties_scaling_json: Path,
            process_executor: concurrent.futures.ProcessPoolExecutor,
            variants_per_batch: int,
            keep_multiallelic: bool,
            keep_homref: bool,
            keep_homvar: bool,
            torch_device_kind: training_utils.TorchDeviceKind,
            progress_logger: training_utils.ProgressLogger,
            shuffle: bool,
            random_generator: common.GeneratorInit,
            supplemental_properties: set[str],
            non_filtration_properties: set[str] = Default.non_filtration_properties,
            temp_dir: Path = Default.temp_dir,
            scale_bool_values: bool = Default.scale_bool_properties,
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
        self.keep_multiallelic = keep_multiallelic
        self.keep_homref = keep_homref
        self.keep_homvar = keep_homvar
        self._temp_dir = self.__class__.get_temp_subdir(temp_dir)
        self.random_generator = \
            common.init_generator(generator_init=random_generator) if shuffle else None
        self._set_column_info(
            non_filtration_properties=non_filtration_properties,
            supplemental_properties=supplemental_properties,
            properties_scaling_json=properties_scaling_json,
            scale_bool_values=scale_bool_values
        )
        self._parquet_files = [] if _parquet_files is None else _parquet_files
        self._num_rows = tarred_properties_to_parquet.get_parquet_folder_num_rows(
            self.parquet_path
        )
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
            non_filtration_properties: set[str],
            supplemental_properties: set[str],
            properties_scaling_json: Path,
            scale_bool_values: bool
    ):
        f"""
        Set a bunch of column info so that it can be quickly retrieved later
          - self._wanted_columns is the list of raw (str) columns saved in the parquet file that we
            want to load, i.e. everything but the excluded properties
          - self._variant_columns is the unflattened (i.e. MultiIndex) columns that correspond to
            per-variant properties
          - self._format_columns is the unflattened columns that correspond to per-genotype /
            format properties
          - self._supplemental_property_getters is a dict from property name to
            SupplementalPropertyGetter to manage getting properties needed for training or
            filtering that are not fed into neural network
        Have a few annoying issues to work around:
          - Absolutely don't want {Keys.row} in any columns (it's the name of the row index, not a
            column).
          - Can only determine if a column is for an excluded property after unflattening
        Args:
            non_filtration_properties: properties that will not be used for training or filtering
            supplemental_properties: properties that may or may not be injested by the neural
                                     network, but are otherwise required for training or filtering,
                                     but are needed. e.g. variant_weights, variant_id, etc 
            properties_scaling_json: Path to JSON file with information about each property's
                                     baseline and scale, used to nominally center them on zero
                                     and in the rough range (-1, 1)
            scale_bool_values: If true, scale boolean values too. If False, render True -> 1 and
                               False -> 0. Generally a bad idea, because some boolean values are
                               extremely unbalanced.
        """
        with open(properties_scaling_json, 'r') as f_in:
            properties_scaling = json.load(f_in)

        # get flat_columns and _columns, explicitly removing Keys.row, the name of the row index
        # -flat_columns are ALL the columns in the dataframe, flattened so columns are not
        #  multi-index (this is what's stored in the parquet file)
        # -self.columns is the unflattened multi-index (tuple of sample_id, property_name) columns
        #  that are easier to process and organize
        flat_columns, self._columns = zip(*(
            (flat_column, tarred_properties_to_parquet.unflatten_column_name(flat_column))
            for flat_column in tarred_properties_to_parquet.get_parquet_folder_columns(
                self.parquet_path
            )
            if flat_column != Keys.row
        ))
        # need to get the flattened column name for variant ID, because it's used in filtering
        # which variants to open
        self._variant_id_flat_column_name = next(
            flat_column for flat_column, unflat_column in zip(flat_columns, self._columns)
            if unflat_column[1] == Keys.id
        )
        excluded_properties = non_filtration_properties.difference(supplemental_properties)
        # get rid of all the columns we don't want at all
        self._wanted_columns, self._columns = zip(*(
            (flat_column, unflat_column)
            for flat_column, unflat_column in zip(flat_columns, self._columns)
            if unflat_column[1] not in excluded_properties
        ))

        # get the DataFrame columns that correspond to per-variant / INFO columns
        self._variant_columns = numpy.array(
            [index for index, column in enumerate(self._columns)
             if not isinstance(column[0], str) and
             column[1] not in non_filtration_properties]
        )
        # get the DataFrame columns that correspond to per-sample / FORMAT columns
        self._format_columns = numpy.array(
            [index for index, column in enumerate(self._columns)
             if isinstance(column[0], str) and column[1] not in non_filtration_properties]
        )
        # get property names (corresponding to dimension 2 of output tensor)
        # organize with per-variant / INFO properties first
        self._property_names = [
            column[1] for column in self._columns
            if not isinstance(column[0], str) and
            column[1] not in non_filtration_properties
        ]
        # get sample IDs in order, and add per-sample / FORMAT properties to property names
        # NOTE we rely on tarred_properties_to_parquet to put the samples into a consistent order
        found_properties = set(self._property_names)
        found_sample_ids = set()
        self._sample_ids = []
        for sample_id, prop_name in self._columns:
            if prop_name in non_filtration_properties or not isinstance(sample_id, str):
                continue  # skip properties not in final tensor, or INFO properties
            if sample_id not in found_sample_ids:
                found_sample_ids.add(sample_id)
                self._sample_ids.append(sample_id)
            if prop_name not in found_properties:
                found_properties.add(prop_name)
                self._property_names.append(prop_name)

        # get the dtypes for each property, needed by supplemental property getters

        # get the supplemental property getters
        dtypes = tarred_properties_to_parquet.get_parquet_folder_dtypes(
            self.parquet_path, unflatten_column_names=True
        )
        self._supplemental_property_getters = {
            supplemental_property: SupplementalPropertyGetter.build(
                property_name=supplemental_property, unflatted_columns=self._columns,
                tensor_property_names=self._property_names, properties_scaling=properties_scaling,
                scale_bool_values=scale_bool_values, dtypes=dtypes
            ).to(device=self.torch_device)
            for supplemental_property in supplemental_properties
        }

        # get baseline and scale to (approximately) z-score output tensor
        self._tensor_baseline = torch.tensor(
            [
                _get_baseline(property_name=prop_name, properties_scaling=properties_scaling,
                              scale_bool_values=scale_bool_values)
                for prop_name in self._property_names
            ],
            dtype=torch.float32, device=self.torch_device
        )
        self._tensor_scale = torch.tensor(
            [
                _get_scale(property_name=prop_name, properties_scaling=properties_scaling,
                           scale_bool_values=scale_bool_values)
                for prop_name in self._property_names
            ],
            dtype=torch.float32, device=self.torch_device
        )

    def get_all_parquet_files(self) -> list[Path]:
        return sorted(
            self.parquet_path.glob(f"*{Keys.parquet_suffix}"),
            key=lambda p: int(p.name.split('.')[1])
        )

    def _set_parquet_files(self):
        self._parquet_files = list(self.get_all_parquet_files())
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
        """return the names of each property ((e.g. 3rd dimension in num-variants x num-samples
        x num-properties tensor)
        """
        return self._property_names

    @property
    def num_properties(self) -> int:
        """return the number of distinct properties (e.g. 3rd dimension in num-variants x
        num-samples x num-properties tensor)
        """
        return len(self.property_names)

    @property
    def supplemental_property_getters(self) -> dict[str, SupplementalPropertyGetter]:
        """Disconnect and move to cpu before pickling for work in remote process"""
        cpu = training_utils.TorchDeviceKind.cpu.get_device(progress_logger=None)
        return {
            property_name: getter.to(cpu)
            for property_name, getter in self._supplemental_property_getters.items()
        }

    def get_unscaled_property_tensor(
        self,
        property_name: str,
        filtration_properties_tensor: numpy.array,
        supplemental_properties_buffers: dict[str, Optional[numpy.ndarray | PandasArray]]
    ) -> torch.Tensor:
        getter = self._supplemental_property_getters[property_name]
        return torch.tensor(
            (
                getter.get_from_tensor(filtration_properties_tensor).astype(getter.dtype)
                if getter.is_tensor
                else supplemental_properties_buffers[property_name]
            ),
            device=self.torch_device
        )

    def get_scaled_property_tensor(
        self,
        property_name: str,
        scaled_filtration_properties_tensor: torch.Tensor,
        supplemental_properties_buffers: dict[str, Optional[numpy.ndarray | PandasArray]]
    ) -> torch.Tensor:
        getter = self._supplemental_property_getters[property_name]
        return (
            getter.get_from_tensor(scaled_filtration_properties_tensor) if getter.is_tensor
            else getter.scale_array(supplemental_properties_buffers[property_name],
                                    device=self.torch_device)
        )


class BatchPicklerBase:
    """Static class that is only used to organize helper methods for (process) parallel pickling
    batch tensors
    """
    @classmethod
    def _load_buffer(
            cls,
            parquet_file: Path,
            wanted_columns: list[str],
            wanted_ids: Optional[Collection[str]] = None,
            variant_id_flat_column_name: Optional[str] = "(nan, 'id')"
    ) -> numpy.ndarray:
        """get the buffer (2D numpy array) that contains desired rows and columns

        Args:
            parquet_file: Path to parquet file with DataFrame of properties
            wanted_columns: List of flattened column names to load (ignoring others). Order of
                            columns in buffer will correspond to order in this list.
            wanted_ids: Collection (e.g. set) of variant IDs that correspond to rows that should
                        be loaded. Order of rows will be the same as the order on disk, i.e.
                        wanted_ids does NOT set the order
            variant_id_flat_column_name: flat column name that corresponds to the column name for
                                         variant IDs. Typically "(nan, 'id')". Must be specified
                                         if wanted_ids is specified.
        """
        if wanted_ids is None:
            filters = None
        else:
            if variant_id_flat_column_name is None:
                raise ValueError(
                    "Must specify variant_id_flat_column_name when specifying wanted_ids"
                )
            _id_filter = (
                variant_id_flat_column_name, "in", wanted_ids
            )
            filters = [_id_filter]
        return pandas.read_parquet(parquet_file, columns=wanted_columns, filters=filters).values

    @classmethod
    def _fill_next_tensor_blocks(
        cls,
        buffer_matrix: numpy.ndarray,
        variant_columns: numpy.ndarray,
        format_columns: numpy.ndarray,
        num_samples: int,
        supplemental_property_getters: dict[str, SupplementalPropertyGetter],
        filtration_properties_tensor: numpy.ndarray,
        supplemental_property_buffers: dict[str, Optional[numpy.ndarray | PandasArray]],
        fill_start_index: Optional[int] = None,
        fill_end_index: Optional[int] = None,
    ) -> None:
        """
        Extract values from buffer_matrix, a 2d numpy array, and use them to fill a 3d numpy array,
        filtration_properties_tensor. Also fill buffers for any needed supplemental properties
        (properties needed for training or filtering that may not be fed to the neural network).

        NOTE:
        1) The per-variant properties will be duplicated during this process
        2) All filtration properties will be converted to 32-bit floats
        3) It is assumed that Categorical properties are already one-hot encoded (by
           tarred_properties_to_parquet.py)

        Args:
            buffer_matrix: 2D numpy array (from DataFrame.values) of variant/sample properties with
                           multi-index columns (1st index corresponding to sample, with None for
                           per-variant properties, 2nd index corresponding to property name)
        (num_variants x num_samples x num_properties)
            variant_columns: Array of indices to columns with per-variant (e.g. INFO) properties
            format_columns: Array of indices to columns with per-genotype (e.g. FORMAT) properties
            num_samples: Number of samples, i.e. genotypes per row
            supplemental_property_getters: dict from property name to getter object that can
                                           get supplemental (non-filtration) properties from the
                                           buffer if needed
            filtration_properties_tensor: num_variants x num_samples x num_properties tensor of
                                          properties that are fed to the neural network for
                                          training or filtering
            supplemental_property_buffers: dict from property name to array that will store
                                           supplemental (non-filtration) properties. If a property
                                           can be gotten from the filtration tensor, the buffer
                                           will be None and will not be filled
            fill_start_index: used for remainders / partial filling. The index in the buffer to
                              begin filling.
            fill_end_index: used for remainders / partial filling. The index in the buffer to
                            stop filling. Note in python that means this index will not be filled.
        """
        num_variant_properties = len(variant_columns)
        fill_range = slice(fill_start_index, fill_end_index)
        # fill the info / per-variant block of numpy_tensor, using broadcasting to repeat for every
        # sample
        filtration_properties_tensor[fill_range, :, :num_variant_properties] = buffer_matrix.take(
            variant_columns, axis=1
        ).reshape(buffer_matrix.shape[0], 1, num_variant_properties)
        # fill the format / per-genotype block of numpy_tensor
        filtration_properties_tensor[fill_range, :, num_variant_properties:] = buffer_matrix.take(
            format_columns, axis=1
        ).reshape(buffer_matrix.shape[0], num_samples, len(format_columns) // num_samples)
        # copy every supplemental property that's taken from the buffer
        for property_name, supplemental_property_getter in supplemental_property_getters.items():
            if supplemental_property_getter.is_buffer:
                supplemental_property_buffers[property_name][fill_range] = \
                    supplemental_property_getter.get_from_buffer(buffer_matrix)
