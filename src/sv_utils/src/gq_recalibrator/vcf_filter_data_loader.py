import attr
from pathlib import Path
import concurrent.futures
import pickle
import numpy
import pandas.core.arrays
import torch
from typing import Any
from sv_utils import genomics_io, get_truth_overlap
from gq_recalibrator import (
    training_utils, tarred_properties_to_parquet, vcf_tensor_data_loader_base
)
from gq_recalibrator.vcf_tensor_data_loader_base import (
    locked_read, locked_write, ArrayType,
    BatchIteratorBase, VcfTensorDataLoaderBase, BatchPicklerBase, SupplementalPropertyGetter
)

ArrowStringArray = pandas.core.arrays.string_arrow.ArrowStringArray
ConfidentVariants = get_truth_overlap.ConfidentVariants


class Keys:
    pickle_suffix = vcf_tensor_data_loader_base.Keys.pickle_suffix
    id = tarred_properties_to_parquet.Keys.id
    is_multiallelic = tarred_properties_to_parquet.Keys.is_multiallelic
    gq = genomics_io.Keys.gq
    allele_count = genomics_io.Keys.allele_count


class Default:
    non_filtration_properties = vcf_tensor_data_loader_base.Default.non_filtration_properties
    supplemental_columns_for_filtering = frozenset((
        Keys.id, Keys.is_multiallelic, Keys.gq, Keys.allele_count
    ))
    scale_bool_properties = vcf_tensor_data_loader_base.Default.scale_bool_properties,
    temp_dir = vcf_tensor_data_loader_base.Default.temp_dir
    sleep_time_s = vcf_tensor_data_loader_base.Default.sleep_time_s
    delete_after_read = vcf_tensor_data_loader_base.Default.delete_after_read
    max_look_ahead_batches = vcf_tensor_data_loader_base.Default.max_look_ahead_batches


class FilterBatchIterator(BatchIteratorBase):
    """
    Iterator class that approximately iterates at the appropriate validation proportion. Works over
    batches and doesn't allow rounding error to gradually accumulate excess training or validation
    variants in the remainder
    """
    __slots__ = (
        "num_variants", "batch", "num_batches", "num_variants_in_remainder_batch"
    )
    __state_keys__ = BatchIteratorBase.__state_keys__ + ("num_variants",)

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
        self.batch += 1
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
            keep_multiallelic: bool,
            keep_homref: bool,
            keep_homvar: bool,
            torch_device_kind: training_utils.TorchDeviceKind,
            progress_logger: training_utils.ProgressLogger,
            temp_dir: Path = Default.temp_dir,
            scale_bool_values:  bool = Default.scale_bool_properties,
            non_filtration_properties: set[str] = Default.non_filtration_properties,
            supplemental_properties: set[str] = Default.supplemental_columns_for_filtering,
            max_look_ahead_batches: int = Default.max_look_ahead_batches,
    ):
        # basically it's the VcfTensorDataLoaderBase with variant weights property excluded and
        # explicitly no validation
        super().__init__(
            parquet_path=parquet_path, properties_scaling_json=properties_scaling_json,
            process_executor=process_executor, variants_per_batch=variants_per_batch,
            keep_multiallelic=keep_multiallelic, keep_homref=keep_homref, keep_homvar=keep_homvar,
            torch_device_kind=torch_device_kind, progress_logger=progress_logger, shuffle=False,
            temp_dir=temp_dir, scale_bool_values=scale_bool_values,
            random_generator=None, non_filtration_properties=non_filtration_properties,
            supplemental_properties=supplemental_properties
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

    @staticmethod
    def _batches_per_parquet_file(parquet_file: Path, variants_per_batch: int) -> int:
        num_rows = tarred_properties_to_parquet.get_parquet_file_num_rows(parquet_file)
        n_batches, extra_rows = divmod(num_rows, variants_per_batch)
        return n_batches + 1 if extra_rows > 0 else n_batches

    def __len__(self) -> int:
        """ return the number of remaining batches in the iterator """
        parquet_files = self.get_all_parquet_files() if not self._parquet_files \
            else self._parquet_files
        n_batches = sum(
            VcfFilterTensorDataLoader._batches_per_parquet_file(
                parquet_file=parquet_file, variants_per_batch=self.variants_per_batch
            ) for parquet_file in parquet_files
        )
        if self._batch_iterator is not None:
            n_batches += len(self._batch_iterator)
        return n_batches

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
        if not self._parquet_files:
            return
        self._current_parquet_file = self._parquet_files.pop()
        pickle_folder = self._get_pickle_folder()
        future = self.process_executor.submit(
            FilterBatchPickler.pickle_batch_tensors,
            parquet_file=self._current_parquet_file, pickle_folder=pickle_folder,
            wanted_columns=self._wanted_columns, variants_per_batch=self.variants_per_batch,
            num_samples=self.num_samples, num_properties=self.num_properties,
            variant_columns=self._variant_columns, format_columns=self._format_columns,
            supplemental_property_getters=self.supplemental_property_getters,
            keep_multiallelic=self.keep_multiallelic, keep_homref=self.keep_homref,
            keep_homvar=self.keep_homvar
        )
        self._filter_futures.append(
            FilterFuture(future=future, pickle_folder=pickle_folder)
        )

    def __next__(self) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        """
        emit tensors for next batch
        Returns:
            filtration_properties_tensor: num_variants x num_samples x num_properties 32-bit float
                                          tensor of properties to be injested by neural network
            original_gq: num_variants x num_samples torch tensor with GQ values from original
                         VCF
            is_filterable: num_variants x num_samples boolean array, set to True if the
                           corresponding genotype is filterable, False otherwise
        """
        if not self._batch_iterator.has_next:
            self._get_next_batch_iterator()
        is_remainder_batch, pickle_file, num_variants_in_batch = next(self._batch_iterator)
        self._current_row += num_variants_in_batch

        with locked_read(
                pickle_file, mode="rb", delete_after_read=True, sleep_time_s=Default.sleep_time_s
        ) as f_in:
            (
                filtration_properties_tensor, supplemental_properties_buffers, packed_is_filterable
            ) = pickle.load(f_in)

        # scale the filtration properties to be appropriate for the neural net
        scaled_filtration_properties_tensor = self._tensor_scale * (
                torch.tensor(
                    filtration_properties_tensor, dtype=torch.float32, device=self.torch_device
                ) - self._tensor_baseline
            )
        # unpack the is_filterable numpy array
        is_filterable = numpy.unpackbits(
            packed_is_filterable, count=num_variants_in_batch * self.num_samples
        ).astype(bool).reshape((num_variants_in_batch, self.num_samples))

        return (
            scaled_filtration_properties_tensor,
            self.get_unscaled_property_tensor(
                property_name=Keys.gq,
                filtration_properties_tensor=filtration_properties_tensor,
                supplemental_properties_buffers=supplemental_properties_buffers
            ),
            torch.tensor(is_filterable, dtype=torch.bool, device=self.torch_device)
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
            supplemental_property_getters: dict[str, SupplementalPropertyGetter],
            keep_multiallelic: bool,
            keep_homref: bool,
            keep_homvar: bool
    ) -> dict[str, Any]:
        # read the basic buffer into a pandas array
        buffer = cls._load_buffer(parquet_file=parquet_file, wanted_columns=wanted_columns)
        filter_batch_iterator = FilterBatchIterator(
            pickle_folder=pickle_folder, variants_per_batch=variants_per_batch,
            parquet_file=parquet_file, num_variants=buffer.shape[0]
        )

        # construct numpy array that will hold tensors for each batch
        filtration_properties_tensor = numpy.empty(
            (variants_per_batch, num_samples, num_properties), dtype=numpy.float32
        )
        supplemental_property_buffers = {
            property_name: property_getter.get_empty_array(variants_per_batch)
            for property_name, property_getter in supplemental_property_getters.items()
        }
        is_filterable = numpy.empty((variants_per_batch, num_samples), dtype=bool)

        # iterate over batches in order, and pickle the batch tensors
        buffer_end = 0
        for is_remainder_batch, pickle_file, num_variants_in_batch in filter_batch_iterator:
            buffer_begin = buffer_end
            buffer_end += num_variants_in_batch
            cls._fill_next_tensor_blocks(
                buffer_matrix=buffer[buffer_begin:buffer_end, :],
                variant_columns=variant_columns,
                format_columns=format_columns,
                num_samples=num_samples,
                supplemental_property_getters=supplemental_property_getters,
                filtration_properties_tensor=filtration_properties_tensor[:num_variants_in_batch],
                supplemental_property_buffers=supplemental_property_buffers,
                fill_end_index=num_variants_in_batch
            )
            packed_is_filterable = cls._pack_is_filterable(
                keep_multiallelic=keep_multiallelic,
                keep_homref=keep_homref,
                keep_homvar=keep_homvar,
                filtration_properties_tensor=filtration_properties_tensor[:num_variants_in_batch],
                is_filterable_matrix=is_filterable[:num_variants_in_batch],
                supplemental_property_getters=supplemental_property_getters,
                supplemental_property_buffers=supplemental_property_buffers
            )
            with locked_write(pickle_file, mode="wb") as f_out:
                pickle.dump(
                    (
                        filtration_properties_tensor[:num_variants_in_batch],
                        supplemental_property_buffers,
                        packed_is_filterable
                    ),
                    f_out,
                    protocol=pickle.HIGHEST_PROTOCOL
                )

        return filter_batch_iterator.state_dict

    @classmethod
    def _pack_is_filterable(
            cls,
            keep_multiallelic: bool,
            keep_homvar: bool,
            keep_homref: bool,
            filtration_properties_tensor: numpy.ndarray,
            is_filterable_matrix: numpy.ndarray,
            supplemental_property_getters: dict[str, SupplementalPropertyGetter],
            supplemental_property_buffers: dict[str, ArrayType]
    ) -> numpy.ndarray:
        is_filterable_matrix[:] = True
        num_rows = is_filterable_matrix.shape[0]

        if keep_multiallelic:
            is_multiallelic = supplemental_property_buffers[Keys.is_multiallelic][:num_rows]
            is_filterable_matrix[is_multiallelic, :] = False

        if keep_homref or keep_homvar:
            allele_count_getter = supplemental_property_getters[Keys.allele_count]
            allele_counts = (
                allele_count_getter.get_from_tensor(filtration_properties_tensor)
                if allele_count_getter.is_tensor
                else supplemental_property_buffers[Keys.allele_count][:num_rows]
            )
            if keep_homref:
                is_filterable_matrix[allele_counts == 0] = False
            if keep_homvar:
                is_filterable_matrix[allele_counts == 2] = False

        return numpy.packbits(is_filterable_matrix)
