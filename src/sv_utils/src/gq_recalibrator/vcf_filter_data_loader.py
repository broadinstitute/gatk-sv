import attr
from pathlib import Path
import concurrent.futures
import pickle
import numpy
import pandas
import pandas.core.arrays
import torch
from typing import Optional, Any
from sv_utils import get_truth_overlap
from gq_recalibrator import training_utils
from gq_recalibrator.vcf_tensor_data_loader_base import (
    locked_read, locked_write,
    BatchIteratorBase, VcfTensorDataLoaderBase, BatchPicklerBase, Default, Keys
)

ArrowStringArray = pandas.core.arrays.string_arrow.ArrowStringArray
ConfidentVariants = get_truth_overlap.ConfidentVariants


class FilterBatchIterator(BatchIteratorBase):
    """
    Iterator class that approximately iterates at the appropriate validation proportion. Works over
    batches and doesn't allow rounding error to gradually accumulate excess training or validation
    variants in the remainder
    """
    __slots__ = (
        "num_variants", "batch", "num_batches", "num_variants_in_remainder_batch"
    )
    __state_keys__ = BatchIteratorBase.__state_keys__ + (
        "num_variants", "batch"
    )

    def __init__(
            self,
            pickle_folder: Path | str,
            variants_per_batch: int,
            parquet_file: Path | str,
            num_variants: int,
            batch: int = 0
    ):
        super().__init__(pickle_folder=pickle_folder, variants_per_batch=variants_per_batch,
                         parquet_file=parquet_file)
        self.num_variants = num_variants
        self.batch = batch
        self.num_batches, self.num_variants_in_remainder_batch = divmod(
            num_variants, variants_per_batch
        )
        if self.num_variants_in_remainder_batch > 0:
            self.num_batches += 1

    @property
    def has_next(self) -> bool:
        return self.batch < self.num_batches

    @property
    def has_remainder_batch(self) -> bool:
        return self.num_variants_in_remainder_batch > 0

    def __next__(self) -> (bool, str, int):
        if self.batch >= self.num_batches:
            raise StopIteration
        is_remainder_batch = self.batch == self.num_batches - 1 and self.has_remainder_batch
        num_variants_in_batch = self.num_variants_in_remainder_batch if is_remainder_batch \
            else self.variants_per_batch
        pickle_file = self.get_nth_pickle_file(n=self.batch)
        return is_remainder_batch, pickle_file, num_variants_in_batch

    def get_nth_pickle_file(self, n: int) -> Path:
        return self.pickle_folder / f"{n}{Keys.pickle_suffix}"


@attr.define(frozen=True, slots=True)
class FilterFuture:
    future: concurrent.futures.Future = attr.field()
    pickle_folder: Path = attr.field()


class VcfFilterTensorDataLoader(VcfTensorDataLoaderBase):
    __slots__ = ("max_look_ahead_batches", "_filter_futures",  "_batch_iterator")

    def __init__(
            self,
            parquet_path: Path,
            properties_scaling_json: Path,
            process_executor: concurrent.futures.ProcessPoolExecutor,
            variants_per_batch: int,
            torch_device_kind: training_utils.TorchDeviceKind,
            progress_logger: training_utils.ProgressLogger,
            temp_dir: Path = Default.temp_dir,
            scale_bool_values: bool = Default.scale_bool_properties,
            excluded_properties: set[str] = Default.excluded_properties,
            max_look_ahead_batches: int = Default.max_look_ahead_batches,
    ):
        # basically it's the VcfTensorDataLoaderBase with variant weights property excluded and
        # explicitly no validation
        super().__init__(
            parquet_path=parquet_path, properties_scaling_json=properties_scaling_json,
            process_executor=process_executor, variants_per_batch=variants_per_batch,
            torch_device_kind=torch_device_kind, progress_logger=progress_logger, shuffle=False,
            temp_dir=temp_dir, scale_bool_values=scale_bool_values,
            random_generator=None,
            excluded_properties=excluded_properties.union((Keys.variant_weights,))
        )
        self.max_look_ahead_batches = max_look_ahead_batches
        self._filter_futures = []
        self._batch_iterator = None

    def __iter__(self) -> "VcfFilterTensorDataLoader":
        self._current_row = 0
        self._set_parquet_files()
        while (
               len(self._parquet_files) > 0 and
               len(self._filter_futures) < self.max_look_ahead_batches
        ):
            self._submit_next_batch()
        self._get_next_batch_iterator()
        return self

    def _get_next_batch_iterator(self):
        """Called when we've finished iterating over the current batch of training and validation
        tensors and need to start on the next. Maintains current _batch_iterator and enough prior
        state info to allow for resume
        """
        if not self._filter_futures:
            raise StopIteration
        # Get next training future
        filter_future = self._filter_futures.pop(0)
        # Get and set current _batch_iterator to the one that describes the next batch, i.e. folder
        # of pickled training and validation tensors
        batch_iterator_state_dict = training_utils.get_result_or_raise(
            future=filter_future.future,
            executor=self.process_executor,
            context=f"Pickling tensors in folder {filter_future.pickle_folder}"
        )
        self._batch_iterator = FilterBatchIterator(**batch_iterator_state_dict)
        # Submit a new batch to be computed and pickled, so that when this batch is done iterating
        # hopefully the next batch will be ready to load
        self._submit_next_batch()

    def _submit_next_batch(self):
        # get the next parquet file
        if not self._parquet_files():
            return
        self._current_parquet_file = self._parquet_files.pop()
        pickle_folder = self._get_pickle_folder()
        future = self.process_executor.submit(
            FilterBatchPickler.pickle_batch_tensors,
            parquet_file=self._current_parquet_file, pickle_folder=pickle_folder,
            wanted_columns=self._wanted_columns, variants_per_batch=self.variants_per_batch,
            num_samples=self.num_samples, num_properties=self.num_properties,
            variant_columns=self._variant_columns, format_columns=self._format_columns,
        )
        self._filter_futures.append(
            FilterFuture(future=future, pickle_folder=pickle_folder)
        )

    @staticmethod
    @training_utils.reraise_with_stack
    def _parquet_file_to_buffer(
            parquet_file: Path,
            wanted_columns: list[str]
    ) -> (numpy.ndarray, Optional[numpy.random.Generator]):
        return pandas.read_parquet(parquet_file, columns=wanted_columns).values

    def __next__(self) -> torch.Tensor:
        """
        emit tensors for next batch
        Returns:
            properties_tensor: torch.Tensor
                num_variants x num_samples x num_properties 32-bit float tensor of training /
                filtering properties
            variant_ids: ArrowStringArray
                num_variants array of variant IDs
        """
        if not self._batch_iterator.has_next:
            self._get_next_batch_iterator()
        is_remainder_batch, pickle_file, num_variants_in_batch = next(self._batch_iterator)
        self._current_row += num_variants_in_batch

        with locked_read(
                pickle_file, mode="rb", delete_after_read=True, sleep_time_s=Default.sleep_time_s
        ) as f_in:
            numpy_tensor = pickle.load(f_in)

        return self._tensor_scale * (
                torch.tensor(numpy_tensor, dtype=torch.float32, device=self.torch_device)
                - self._tensor_baseline
            )


class FilterBatchPickler(BatchPicklerBase):
    @classmethod
    @training_utils.reraise_with_stack
    def pickle_batch_tensors(
            cls,
            parquet_file: Path,
            pickle_folder: Path,
            wanted_columns: list[str],
            variants_per_batch: int,
            num_samples: int,
            num_properties: int,
            variant_columns: numpy.ndarray,
            format_columns: numpy.ndarray,
    ) -> dict[str, Any]:
        # read the basic buffer into a pandas array
        buffer = cls._load_buffer(parquet_file=parquet_file, wanted_columns=wanted_columns)
        filter_batch_iterator = FilterBatchIterator(
            pickle_folder=pickle_folder, variants_per_batch=variants_per_batch,
            parquet_file=parquet_file, num_variants=buffer.shape[0]
        )

        # construct numpy array that will hold tensors for each batch
        numpy_tensor = numpy.empty(
            (variants_per_batch, num_samples, num_properties), dtype=numpy.float32
        )

        # iterate over batches in order, and pickle the batch tensors
        num_variant_properties = len(variant_columns)
        buffer_end = 0
        for is_remainder_batch, pickle_file, num_variants_in_batch in filter_batch_iterator:
            buffer_begin = buffer_end
            buffer_end += num_variants_in_batch
            if is_remainder_batch:
                (
                    numpy_tensor[:num_variants_in_batch, :, :num_variant_properties],
                    numpy_tensor[:num_variants_in_batch, :, num_variant_properties:]
                ) = cls._get_next_tensor_blocks(
                    buffer_matrix=buffer[buffer_begin:buffer_end, :],
                    variant_columns=variant_columns,
                    format_columns=format_columns,
                    num_samples=num_samples
                )
                with locked_write(pickle_file, mode="wb") as f_out:
                    pickle.dump(numpy_tensor[:num_variants_in_batch, :, :], f_out,
                                protocol=pickle.HIGHEST_PROTOCOL)
            else:
                (
                    numpy_tensor[:, :, :num_variant_properties],
                    numpy_tensor[:, :, num_variant_properties:]
                ) = cls._get_next_tensor_blocks(
                    buffer_matrix=buffer[buffer_begin:buffer_end, :],
                    variant_columns=variant_columns,
                    format_columns=format_columns,
                    num_samples=num_samples
                )
                with locked_write(pickle_file, mode="wb") as f_out:
                    pickle.dump(numpy_tensor, f_out, protocol=pickle.HIGHEST_PROTOCOL)

        return filter_batch_iterator.state_dict

    @classmethod
    def _get_next_tensor_blocks(
            cls,
            buffer_matrix: numpy.ndarray,
            variant_columns: numpy.ndarray,
            format_columns: numpy.ndarray,
            num_samples: int
    ) -> (numpy.ndarray, numpy.ndarray):
        """
        Convert DataFrame with multi-index columns (1st index corresponding to sample, with None
        for per-variant properties, 2nd index corresponding to property name) to 3-tensor
        (num_variants x num_samples x num_properties).
        NOTE:
        1) The per-variant properties will be duplicated during this process
        2) All properties will be converted to 32-bit floats
        3) It is assumed that Categorical properties are already one-hot encoded (by
           tarred_properties_to_parquet.py)

        Args:
            buffer_matrix: numpy.ndarray
                Numpy matrix (from DataFrame.values) of variant/sample properties with multi-index
                columns
            variant_columns: numpy.ndarray
                Array of indices to columns with per-variant (e.g. INFO) properties
            format_columns: numpy.ndarray
                Array of indices to columns with per-genotype (e.g. FORMAT) properties
            num_samples: int
                Number of samples / genotypes per row
        Returns:
            variant_properties: numpy.ndarray
                num_rows x 1 x num_variant_properties matrix of per-variant / INFO properties
            format_properties: numpy.ndarray
                num_rows x num_samples x num_format_properties matrix of per-genotype / FORMAT
                properties
        """
        # fill num_rows x 1 x num_variant_properties matrix with per-variant properties. Note, the
        # "1" dimension will broadcast to fill the 2nd (num_samples) dim
        num_variant_properties = len(variant_columns)
        return (
            buffer_matrix.take(variant_columns, axis=1).reshape(buffer_matrix.shape[0], 1,
                                                                num_variant_properties),
            buffer_matrix.take(format_columns, axis=1)
                .reshape(buffer_matrix.shape[0], num_samples, len(format_columns) // num_samples)
        )  # noqa E131   # ignore pep8 131 error, I think this is more legible
