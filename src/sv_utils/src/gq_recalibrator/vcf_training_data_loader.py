import attr
from pathlib import Path, PosixPath
import json
import concurrent.futures
import pickle
import numpy
import pandas
import pandas.core.arrays
import torch
from collections import defaultdict
from typing import Optional, Union, Any
from collections.abc import Sequence
from gq_recalibrator import (
    training_utils, tarred_properties_to_parquet, vcf_tensor_data_loader_base
)
from gq_recalibrator.vcf_tensor_data_loader_base import (
    locked_read, locked_write, BatchIteratorBase, VcfTensorDataLoaderBase, BatchPicklerBase,
    SupplementalPropertyGetter
)
from sv_utils import common, genomics_io, get_truth_overlap, benchmark_variant_filter

ArrowStringArray = pandas.core.arrays.string_arrow.ArrowStringArray
ConfidentVariants = get_truth_overlap.ConfidentVariants


class Keys:
    row = vcf_tensor_data_loader_base.Keys.row
    id = tarred_properties_to_parquet.Keys.id
    variant_weights = tarred_properties_to_parquet.Keys.variant_weights
    training = "training"
    validation = "validation"
    pickle_suffix = vcf_tensor_data_loader_base.Keys.pickle_suffix
    remainder = "remainder"
    gq = genomics_io.Keys.gq
    remainder_sizes = f"{remainder}_sizes{pickle_suffix}"
    training_remainder = f"{training}_{remainder}"
    validation_remainder = f"{validation}_{remainder}"


class Default:
    shuffle = True
    random_seed = 0
    non_filtration_properties = vcf_tensor_data_loader_base.Default.non_filtration_properties
    supplemental_columns_for_training = frozenset({Keys.id, Keys.variant_weights, Keys.gq})
    scale_bool_properties = vcf_tensor_data_loader_base.Default.scale_bool_properties,
    validation_proportion = 0.2
    temp_dir = vcf_tensor_data_loader_base.Default.temp_dir
    sleep_time_s = vcf_tensor_data_loader_base.Default.sleep_time_s
    delete_after_read = vcf_tensor_data_loader_base.Default.delete_after_read
    max_look_ahead_batches = vcf_tensor_data_loader_base.Default.max_look_ahead_batches


class TrainingBatchIterator(BatchIteratorBase):
    """
    Iterator class that approximately iterates at the appropriate validation proportion. Works over
    batches and doesn't allow rounding error to gradually accumulate excess training or validation
    variants in the remainder
    """
    __slots__ = (
        "num_training_variants", "num_validation_variants", "training_batch", "validation_batch",
        "training_remainder_ids", "validation_remainder_ids", "num_training_batches",
        "num_validation_batches", "_approx_validation_proportion", "num_training_remainder",
        "num_validation_remainder"
    )
    __state_keys__ = BatchIteratorBase.__state_keys__ + (
        "num_training_variants", "num_validation_variants", "training_remainder_ids",
        "validation_remainder_ids"
    )

    empty_ids = pandas.array([], dtype="string[pyarrow]")

    def __init__(
            self,
            pickle_folder: Path | str,
            variants_per_batch: int,
            parquet_file: Path | str,
            num_training_variants: int,
            num_validation_variants: int,
            training_batch: int = 0,
            validation_batch: int = 0,
            training_remainder_ids: Sequence[str] = empty_ids,
            validation_remainder_ids: Sequence[str] = empty_ids
    ):
        super().__init__(pickle_folder=pickle_folder, variants_per_batch=variants_per_batch,
                         parquet_file=parquet_file)
        self.num_training_variants = num_training_variants
        self.num_validation_variants = num_validation_variants
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
        self.training_remainder_ids = pandas.array(training_remainder_ids, dtype="string[pyarrow]")
        self.validation_remainder_ids = pandas.array(validation_remainder_ids,
                                                     dtype="string[pyarrow]")

    @property
    def is_training_batch(self) -> Optional[bool]:
        """
        return True if next batch is training,
               False if it's validation,
               None if it's time to stop iteration
        """
        if self.training_batch >= self.num_training_batches:
            return None if self.validation_batch >= self.num_validation_batches else False
        elif self.validation_batch >= self.num_validation_batches:
            return True
        return (
            self.training_batch * self._approx_validation_proportion <
            (self.validation_batch + 1) * (1.0 - self._approx_validation_proportion)
        )

    @property
    def has_next(self) -> bool:
        return (
            self.training_batch < self.num_training_batches or
            self.validation_batch < self.num_validation_batches
        )

    @property
    def has_remainders(self) -> bool:
        return self.num_training_remainder > 0 or self.num_validation_remainder > 0

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
        pickle_file = self.get_nth_pickle_file(n=batch_number, is_training=is_training_batch)
        return is_training_batch, pickle_file, batch_number

    @property
    def remainder_ids(self) -> set[str]:
        return set(self.training_remainder_ids).union(self.validation_remainder_ids)

    def get_nth_pickle_file(self, n: Union[int, str], is_training: bool) -> Path:
        return self.pickle_folder / (
            f"{Keys.training if is_training else Keys.validation}_{n}{Keys.pickle_suffix}"
        )


@attr.define(frozen=True, slots=True)
class TrainingFuture:
    future: concurrent.futures.Future = attr.field()
    pickle_folder: Path = attr.field()
    generator_state: Optional[common.BitGeneratorState] = attr.field()


class VcfTrainingTensorDataLoader(VcfTensorDataLoaderBase):
    __slots__ = (
        "validation_proportion", "_training_batch", "_validation_batch", "_batch_iterator",
        "_previous_batch_iterators", "_genotype_is_good_dict", "_genotype_is_bad_dict",
        "_genotype_is_good_tensor", "_genotype_is_bad_tensor", "_trainable_variant_ids",
        "max_look_ahead_batches", "_resume_generator_state", "_training_futures",
    )
    __state_keys__ = (
        "_resume_generator_state", "_current_parquet_file", "_parquet_files", "_current_row",
        "_training_batch", "_validation_batch", "_previous_batch_iterators"
    )

    def __init__(
            self,
            parquet_path: Path,
            properties_scaling_json: Path,
            truth_json: Path,
            process_executor: concurrent.futures.ProcessPoolExecutor,
            variants_per_batch: int,
            keep_multiallelic: bool,
            keep_homref: bool,
            keep_homvar: bool,
            torch_device_kind: training_utils.TorchDeviceKind,
            progress_logger: training_utils.ProgressLogger,
            temp_dir: Path = Default.temp_dir,
            shuffle: bool = Default.shuffle,
            scale_bool_values: bool = Default.scale_bool_properties,
            random_generator: common.GeneratorInit = Default.random_seed,
            non_filtration_properties: set[str] = Default.non_filtration_properties,
            supplemental_properties: set[str] = Default.supplemental_columns_for_training,
            validation_proportion: float = Default.validation_proportion,
            max_look_ahead_batches: int = Default.max_look_ahead_batches,
            _resume_generator_state: Optional[common.BitGeneratorState] = None,
            _current_parquet_file: Optional[Path] = None,
            _parquet_files: Optional[list[Path]] = None,
            _current_row: int = 0,
            _training_batch: int = 0,
            _validation_batch: int = 0,
            _batch_iterator: Optional[TrainingBatchIterator] = None,
            _previous_batch_iterators: Optional[list[TrainingBatchIterator]] = None,
    ):
        # basically it's the VcfTensorDataLoaderBase with a train/test split
        super().__init__(
            parquet_path=parquet_path, properties_scaling_json=properties_scaling_json,
            process_executor=process_executor, variants_per_batch=variants_per_batch,
            torch_device_kind=torch_device_kind, progress_logger=progress_logger,
            temp_dir=temp_dir, shuffle=shuffle, keep_multiallelic=keep_multiallelic,
            keep_homref=keep_homref, keep_homvar=keep_homvar, scale_bool_values=scale_bool_values,
            random_generator=(
                random_generator if _resume_generator_state is None else _resume_generator_state
            ),
            non_filtration_properties=non_filtration_properties,
            supplemental_properties=supplemental_properties,
            _current_parquet_file=_current_parquet_file,
            _parquet_files=_parquet_files,
            _current_row=_current_row
        )
        self._training_futures = []
        self._resume_generator_state = self.random_generator
        self.max_look_ahead_batches = max_look_ahead_batches
        self._batch_iterator = None
        if _previous_batch_iterators is not None and len(_previous_batch_iterators) > 0:
            # need to get a new pickle_folder, the old one is irrelevant now
            _previous_pickle_folder = self._get_pickle_folder()
            self._previous_batch_iterators = [
                TrainingBatchIterator(
                    **_previous_batch_iterator.state_dict
                ).set_pickle_folder(_previous_pickle_folder)
                for _previous_batch_iterator in _previous_batch_iterators
            ]
            # create a process to write the remainder files for the batch before the one we'll
            # resume with
            self.progress_logger("readied restart_pickle_previous_remainders")
            future = self.process_executor.submit(
                TrainingBatchPickler.restart_pickle_previous_remainders,
                previous_batch_iterators=self._previous_batch_iterators,
                wanted_columns=self._wanted_columns,
                variant_id_column=self._variant_id_flat_column_name,
                supplemental_property_getters=self.supplemental_property_getters
            )
            self._training_futures.append(
                TrainingFuture(
                    future=future,
                    pickle_folder=_previous_pickle_folder,
                    generator_state=self.random_generator.bit_generator.state
                )
            )
        else:
            self.progress_logger("skipped restart_pickle_previous_remainders")
            self._previous_batch_iterators = []

        if validation_proportion <= 0 or validation_proportion >= 1:
            raise ValueError("Validation proportion must be in open interval (0, 1)")
        self.validation_proportion = validation_proportion
        self._training_batch = _training_batch
        self._validation_batch = _validation_batch
        self._genotype_is_good_dict, self._genotype_is_bad_dict = \
            self._reorganize_confident_variants(
                benchmark_variant_filter.TruthData.load_confident_variants(f"{truth_json}")
            )
        self._trainable_variant_ids = frozenset(
            self._genotype_is_good_dict.keys()
        ).union(self._genotype_is_bad_dict.keys())
        self._genotype_is_good_tensor = numpy.empty((self.variants_per_batch, self.num_samples),
                                                    dtype=numpy.float32)
        self._genotype_is_bad_tensor = numpy.empty((self.variants_per_batch, self.num_samples),
                                                   dtype=numpy.float32)

    def __len__(self) -> int:
        """ return the number of remaining batches in the iterator """
        n_batches, extra_rows = divmod(self._num_rows - self._current_row, self.variants_per_batch)
        return n_batches + 1 if extra_rows else n_batches

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

    def __iter__(self) -> "VcfTrainingTensorDataLoader":
        """ initialize for iteration over next full epoch """
        while len(self._training_futures) < self.max_look_ahead_batches:
            # while we haven't started prepping far enough into the future, submit more batches
            self._submit_next_batch()
        while self._batch_iterator is None:
            # if we don't have a current batch iterator, get it from futures. This should only
            # happen for the first iteration epoch
            self._get_next_batch_iterator()
        return self

    def _get_next_batch_iterator(self):
        """Called when we've finished iterating over the current batch of training and validation
        tensors and need to start on the next. Maintains current _batch_iterator and enough prior
        state info to allow for resume
        """
        # Maintain enough previous batch iterators to resume iteration if the program is
        # interrupted
        self._append_previous_batch_iterator(self._batch_iterator)
        # Get next training future
        training_future = self._training_futures.pop(0)
        # Maintain previous random generator state, so that iteration can be resumed if the
        #            program is interrupted
        self._resume_generator_state = training_future.generator_state
        # Get and set current _batch_iterator to the one that describes the next batch, i.e. folder
        # of pickled training and validation tensors
        batch_iterator_state_dict = training_utils.get_result_or_raise(
            future=training_future.future,
            executor=self.process_executor,
            context=f"Pickling tensors in folder {training_future.pickle_folder}"
        )
        self._batch_iterator = None if batch_iterator_state_dict is None \
            else TrainingBatchIterator(**batch_iterator_state_dict)
        # Submit a new batch to be computed and pickled, so that when this batch is done iterating
        # hopefully the next batch will be ready to load
        self._submit_next_batch()

    def _append_previous_batch_iterator(self, batch_iterator: Optional[TrainingBatchIterator]):
        """Maintain enough previous batch iterators to resume iteration if the program is
        interrupted
        """
        if batch_iterator is None:
            # just starting, nothing to do
            return
        num_needed_training_ids = len(batch_iterator.training_remainder_ids)
        num_needed_validation_ids = len(batch_iterator.validation_remainder_ids)
        if num_needed_training_ids == 0 and num_needed_validation_ids == 0:
            # no remainders at all
            self._previous_batch_iterators = []
            return
        self._previous_batch_iterators.append(batch_iterator)
        # we'll need to carry over some or all existing previous iterators
        for ind in range(len(self._previous_batch_iterators) - 1, 1, -1):
            # see if we have enough remainders to stop early
            previous_batch_iterator = self._previous_batch_iterators[ind]
            num_needed_training_ids -= len(previous_batch_iterator.training_remainder_ids)
            num_needed_validation_ids -= len(previous_batch_iterator.validation_remainder_ids)
            if num_needed_training_ids <= 0 and num_needed_validation_ids <= 0:
                # keep this and all later batch iterators
                self._previous_batch_iterators = self._previous_batch_iterators[ind:]
                return

    def _submit_next_batch(self):
        """Create/launch concurrent future object to pickle next batch of tensors"""
        # get the next parquet file
        if not self._parquet_files:
            self._set_parquet_files()
        self._current_parquet_file = self._parquet_files.pop()
        # randomly get a (temporary) pickle folder
        pickle_folder = self._get_pickle_folder()
        # in order to truly be able to restore state, need to be able to restore future results at
        # start time...
        previous_pickle_folder = (
            self._training_futures[-1].pickle_folder
            if self._training_futures
            else None
        )
        future = self.process_executor.submit(
            TrainingBatchPickler.pickle_batch_tensors,
            parquet_file=self._current_parquet_file, pickle_folder=pickle_folder,
            previous_pickle_folder=previous_pickle_folder,
            wanted_columns=self._wanted_columns, random_seed=self._next_random_seed,
            variants_per_batch=self.variants_per_batch, num_samples=self.num_samples,
            num_properties=self.num_properties,
            variant_columns=self._variant_columns, format_columns=self._format_columns,
            validation_proportion=self.validation_proportion,
            trainable_variant_ids=self._trainable_variant_ids,
            variant_id_flat_column_name=self._variant_id_flat_column_name,
            supplemental_property_getters=self.supplemental_property_getters
        )
        self._training_futures.append(
            TrainingFuture(future=future, pickle_folder=pickle_folder,
                           generator_state=self.random_generator.bit_generator.state)
        )

    @property
    def _next_random_seed(self) -> Optional[int]:
        """
        To avoid transmission costs and make resume easier, send random seeds to child processes,
        rather than entire random states
        This function generates the next random seed to send (or returns None if shuffling is not
        requested)
        """
        if self.shuffle:
            seed_dtype = numpy.uint32
            seed_dtype_info = numpy.iinfo(seed_dtype)
            return self.random_generator.integers(seed_dtype_info.min, seed_dtype_info.max,
                                                  dtype=seed_dtype)
        else:
            return None

    def __next__(self) -> (
            bool, torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor
    ):
        """
        emit tensors for next training batch
        Returns:
            properties_tensor: torch.Tensor
                num_variants x num_samples x num_properties 32-bit float tensor of training /
                filtering properties
            weights_tensor: torch.Tensor
                num_variants array of training weights that help balance variants
            variant_ids: ArrowStringArray
                num_variants array of variant IDs
        """
        if self._current_row >= self._num_rows:
            # since batch size is unlikely to be an even divisor of the data set size, we can
            # either exactly traverse the data set (and have an uneven final batch) or have uniform
            # batch size (but only approximately traverse the data set each epoch).  Choose to use
            # uniform batch size and approximate traversal:
            self._current_row -= self._num_rows
            raise StopIteration

        if not self._batch_iterator.has_next:
            self._get_next_batch_iterator()

        is_training_batch, pickle_file, __ = next(self._batch_iterator)
        if is_training_batch:
            self._training_batch += 1
        else:
            self._validation_batch += 1
        self._current_row += self.variants_per_batch

        with locked_read(
                pickle_file, mode="rb", delete_after_read=True, sleep_time_s=Default.sleep_time_s
        ) as f_in:
            filtration_properties_tensor, supplemental_properties_buffers = pickle.load(f_in)
        # fill good and bad truth tensors. Note variant ID is always a non-filtration property,
        # so we know the buffer was passed back
        VcfTrainingTensorDataLoader._fill_truth_tensor(
            supplemental_properties_buffers[Keys.id],
            self._genotype_is_good_dict,
            self._genotype_is_good_tensor
        )
        VcfTrainingTensorDataLoader._fill_truth_tensor(
            supplemental_properties_buffers[Keys.id],
            self._genotype_is_bad_dict,
            self._genotype_is_bad_tensor
        )
        # scale filtration properties for learning
        scaled_filtration_properties_tensor = self._tensor_scale * (
            torch.tensor(
                filtration_properties_tensor, dtype=torch.float32, device=self.torch_device
            ) - self._tensor_baseline
        ).nan_to_num()
        return (
            is_training_batch,
            scaled_filtration_properties_tensor,
            self.get_unscaled_property_tensor(
                property_name=Keys.variant_weights,
                filtration_properties_tensor=filtration_properties_tensor,
                supplemental_properties_buffers=supplemental_properties_buffers
            ),
            self.get_scaled_property_tensor(
                property_name=Keys.gq,
                scaled_filtration_properties_tensor=scaled_filtration_properties_tensor,
                supplemental_properties_buffers=supplemental_properties_buffers
            ),
            torch.tensor(self._genotype_is_good_tensor, dtype=torch.float32,
                         device=self.torch_device),
            torch.tensor(self._genotype_is_bad_tensor, dtype=torch.float32,
                         device=self.torch_device)
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
    def get_nth_pickle_file(pickle_folder: Path, n: Union[int, str], is_training: bool) -> Path:
        return pickle_folder / (
            f"{Keys.training if is_training else Keys.validation}_{n}{Keys.pickle_suffix}"
        )

    @staticmethod
    def get_next_batch(
            training_batch: int, validation_batch: int, validation_proportion: float
    ) -> (bool, int, int):
        is_training_batch = training_batch * validation_proportion < \
                            (validation_batch + 1) * (1.0 - validation_proportion)
        if is_training_batch:
            training_batch += 1
        else:
            validation_batch += 1
        return is_training_batch, training_batch, validation_batch

    @property
    def state_dict(self) -> dict[str, Any]:
        return {key: getattr(self, key) for key in VcfTrainingTensorDataLoader.__state_keys__}

    @property
    def state_dict_json_str(self) -> str:
        """ useful for debugging, get the current state_dict in pretty json form """
        return json.dumps(
            {
                key: (
                    value.state_dict if isinstance(value, TrainingBatchIterator)
                    else [f"{path}" for path in value] if key == "_parquet_files"
                    else f"{value}" if isinstance(value, Path | PosixPath)
                    else value
                )
                for key, value in self.state_dict.items()
            },
            indent="  "
        )


class TrainingBatchPickler(BatchPicklerBase):
    @classmethod
    @training_utils.reraise_with_stack
    def pickle_batch_tensors(
            cls,
            parquet_file: Path,
            pickle_folder: Path,
            previous_pickle_folder: Optional[Path],
            wanted_columns: list[str],
            random_seed: Optional[int],
            variants_per_batch: int,
            num_samples: int,
            num_properties: int,
            variant_columns: numpy.ndarray,
            format_columns: numpy.ndarray,
            validation_proportion: float,
            trainable_variant_ids: set[str],
            variant_id_flat_column_name: str,
            supplemental_property_getters: dict[str, SupplementalPropertyGetter],
            keep_multiallelic: bool,
            keep_homref: bool,
            keep_homvar: bool
    ) -> dict[str, Any]:
        # read the basic buffer into a numpy array, only loading trainable variant IDs
        buffer = cls._load_buffer(
            parquet_file=parquet_file,
            wanted_columns=wanted_columns,
            wanted_ids=trainable_variant_ids,
            variant_id_flat_column_name=variant_id_flat_column_name
        )
        # split new buffer up so that the first block of rows is for validation, remaining for
        # training, and randomly permute the indices within training / permutation blocks if
        # desired
        validation_index = round(buffer.shape[0] * validation_proportion)
        if random_seed is not None:
            random_generator = common.init_generator(generator_init=random_seed)
            perm_validation = random_generator.permutation(validation_index)
            perm_training = validation_index + random_generator.permutation(
                buffer.shape[0] - validation_index
            )
        else:
            perm_validation = numpy.arange(validation_index)
            perm_training = numpy.arange(validation_index, buffer.shape[0])
        # load remainder sizes
        if previous_pickle_folder is None:
            num_previous_training_remainder, num_previous_validation_remainder = 0, 0
        else:
            previous_remainder_sizes = previous_pickle_folder / Keys.remainder_sizes
            with locked_read(previous_remainder_sizes, "rb") as f_in:
                num_previous_training_remainder, num_previous_validation_remainder = \
                    pickle.load(f_in)

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
            training_remainder_ids=TrainingBatchIterator.empty_ids,
            validation_remainder_ids=TrainingBatchIterator.empty_ids
        )
        # pickle remainder variants (variants that are not enough for a complete batch) if any
        # do this FIRST, because now there will be enough information for another process to
        # proceed onto the next batch of variants
        # first quickly output just the sizes, because that can be very quick and is enough to
        # enable planning
        cls._write_remainder_sizes(
            pickle_folder=pickle_folder,
            num_training_remainder=training_batch_iterator.num_training_remainder,
            num_validation_remainder=training_batch_iterator.num_validation_remainder
        )
        # load previous remainder buffers and pickle new remainder buffers
        training_remainder_buffer, training_batch_iterator.training_remainder_ids = \
            cls._handle_new_and_old_remainder_buffers(
                is_training=True, previous_pickle_folder=previous_pickle_folder,
                num_previous_remainder=num_previous_training_remainder,
                pickle_folder=pickle_folder,
                num_remainder=training_batch_iterator.num_training_remainder,
                buffer=buffer, variant_id_getter=supplemental_property_getters[Keys.id],
                buffer_indices=perm_training
            )
        validation_remainder_buffer, training_batch_iterator.validation_remainder_ids = \
            cls._handle_new_and_old_remainder_buffers(
                is_training=False, previous_pickle_folder=previous_pickle_folder,
                num_previous_remainder=num_previous_validation_remainder,
                pickle_folder=pickle_folder,
                num_remainder=training_batch_iterator.num_validation_remainder,
                buffer=buffer, variant_id_getter=supplemental_property_getters[Keys.id],
                buffer_indices=perm_validation
            )
        # construct numpy arrays that will hold tensors for each batch
        filtration_properties_tensor = numpy.empty(
            (variants_per_batch, num_samples, num_properties), dtype=numpy.float32
        )

        supplemental_property_buffers = {
            property_name: property_getter.get_empty_array(variants_per_batch)
            for property_name, property_getter in supplemental_property_getters.items()
        }

        # iterate over batches in order, and pickle the batch tensors
        training_buffer_end = 0
        validation_buffer_end = 0
        for is_training_batch, pickle_file, batch_number in training_batch_iterator:
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
                cls._fill_next_tensor_blocks(
                    buffer_matrix=remainder_buffer,
                    variant_columns=variant_columns,
                    format_columns=format_columns,
                    num_samples=num_samples,
                    supplemental_property_getters=supplemental_property_getters,
                    filtration_properties_tensor=filtration_properties_tensor,
                    supplemental_property_buffers=supplemental_property_buffers,
                    fill_end_index=num_remainder
                )
                cls._fill_next_tensor_blocks(
                    buffer_matrix=buffer.take(indices, axis=0),
                    variant_columns=variant_columns,
                    format_columns=format_columns,
                    num_samples=num_samples,
                    supplemental_property_getters=supplemental_property_getters,
                    filtration_properties_tensor=filtration_properties_tensor,
                    supplemental_property_buffers=supplemental_property_buffers,
                    fill_start_index=num_remainder
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
                cls._fill_next_tensor_blocks(
                    buffer_matrix=buffer.take(indices, axis=0),
                    variant_columns=variant_columns,
                    format_columns=format_columns,
                    num_samples=num_samples,
                    supplemental_property_getters=supplemental_property_getters,
                    filtration_properties_tensor=filtration_properties_tensor,
                    supplemental_property_buffers=supplemental_property_buffers
                )
            with locked_write(pickle_file, mode="wb") as f_out:
                pickle.dump((filtration_properties_tensor, supplemental_property_buffers), f_out,
                            protocol=pickle.HIGHEST_PROTOCOL)
        return training_batch_iterator.state_dict

    @classmethod
    @training_utils.reraise_with_stack
    def restart_pickle_previous_remainders(
            cls,
            previous_batch_iterators: list[TrainingBatchIterator],
            wanted_columns: list[str],
            variant_id_flat_column_name: str,
            supplemental_property_getters: dict[str, SupplementalPropertyGetter],
            keep_multiallelic: bool,
            keep_homref: bool,
            keep_homvar: bool
    ) -> None:
        """
        For restarting training after interruption, recreate the remainder tensor for the batch
        prior to the resumed batch
        """
        for ind, batch_iterator in enumerate(previous_batch_iterators):
            print(f"{ind}: {json.dumps(batch_iterator.json_dict, indent='  ')}")
        last_batch_iterator = previous_batch_iterators[-1]
        if not last_batch_iterator.has_remainders:
            raise ValueError("Should not _restart_pickle_previous_remainders with no remainders")
        cls._write_remainder_sizes(
            pickle_folder=last_batch_iterator.pickle_folder,
            num_training_remainder=last_batch_iterator.num_training_remainder,
            num_validation_remainder=last_batch_iterator.num_validation_remainder
        )

        # concatenate the remainder buffers together
        buffer = numpy.concatenate(
            tuple(
                cls._load_buffer(
                    parquet_file=_batch_iterator.parquet_file, wanted_columns=wanted_columns,
                    wanted_ids=_batch_iterator.remainder_ids,
                    variant_id_flat_column_name=variant_id_flat_column_name
                )
                for _batch_iterator in previous_batch_iterators
            )
        )

        # get the buffer's variant ids in order
        variant_ids = supplemental_property_getters[Keys.id].get_from_buffer(buffer)
        # sort variant IDs, keeping track of the original indices
        buffer_indices = numpy.argsort(variant_ids)
        variant_ids = variant_ids.take(buffer_indices)

        def _pickle_last_remainder(is_training: bool):
            """ function to pickle training or validation remainder"""
            use_case = "training" if is_training else "validation"
            num_remainder = getattr(last_batch_iterator, f"num_{use_case}_remainder")
            if num_remainder > 0:
                # write the previous training/validation remainder buffer
                # a. get file name for new remainder
                pickle_file = last_batch_iterator.get_nth_pickle_file(
                    n=Keys.remainder, is_training=is_training
                )
                # b. get the complete set of new training remainder IDs
                remainder_ids = numpy.concatenate(
                    tuple(
                        getattr(_batch_iterator, f"{use_case}_remainder_ids")
                        for _batch_iterator in previous_batch_iterators
                    )
                )[-num_remainder:]
                # c. get the indices into the buffer of the permuted remainder ids
                permuted_indices = buffer_indices.take(
                    numpy.searchsorted(variant_ids, remainder_ids)
                )
                # c. write the appropriate buffer rows in the correct permuted order
                with locked_write(pickle_file, mode="wb") as f_out:
                    pickle.dump(
                        buffer.take(permuted_indices, axis=0),
                        f_out,
                        protocol=pickle.HIGHEST_PROTOCOL
                    )
                print(f"pickled {pickle_file}")

        # pickle remainders
        _pickle_last_remainder(is_training=True)
        _pickle_last_remainder(is_training=False)
        # don't return this, because it's a phony batch iterator now. It was just used to recreate
        # a pickle file for the first real batch iterator
        # return last_batch_iterator.state_dict
        return None

    @staticmethod
    def _write_remainder_sizes(
            pickle_folder: Path,
            num_training_remainder: int,
            num_validation_remainder: int
    ):
        """"Write remainder sizes so that parallel pickling can happen quickly"""
        remainder_sizes_file = pickle_folder / Keys.remainder_sizes
        with locked_write(remainder_sizes_file, "wb") as f_out:
            pickle.dump(
                (num_training_remainder, num_validation_remainder),
                f_out
            )

    @classmethod
    def _handle_new_and_old_remainder_buffers(
            cls,
            is_training: bool,
            previous_pickle_folder: Path,
            num_previous_remainder: int,
            pickle_folder: Path,
            num_remainder: int,
            buffer: numpy.ndarray,
            variant_id_getter: SupplementalPropertyGetter,
            buffer_indices: numpy.ndarray,
            delete_after_read: bool = Default.delete_after_read,
            sleep_time_s: float = Default.sleep_time_s
    ) -> (numpy.ndarray, ArrowStringArray):
        has_previous_remainder = num_previous_remainder > 0
        # get file name for new remainder
        pickle_file = VcfTrainingTensorDataLoader.get_nth_pickle_file(
            pickle_folder, n=Keys.remainder, is_training=is_training
        )
        if num_remainder > buffer_indices.size:
            # need to copy over previous remainder to next remainder, so need to load previous
            # remainder first (this is slower for parallelization, so avoid if possible)
            previous_remainder_buffer = cls._load_remainder_buffer(
                pickle_folder=previous_pickle_folder, is_training=is_training,
                num_columns=buffer.shape[1], delete_after_read=delete_after_read,
                sleep_time_s=sleep_time_s
            ) if has_previous_remainder else cls._get_empty_buffer(
                num_columns=buffer.shape[1]
            )
            buffer = buffer.take(buffer_indices, axis=0)
            remainder_variant_ids = variant_id_getter.get_from_buffer(buffer)
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
            if num_remainder > 0:
                end = len(buffer_indices)
                begin = end - num_remainder
                buffer = buffer.take(buffer_indices[begin:end], axis=0)
                remainder_variant_ids = variant_id_getter.get_from_buffer(buffer)
                with locked_write(pickle_file, mode="wb") as f_out:
                    pickle.dump(buffer, f_out, protocol=pickle.HIGHEST_PROTOCOL)
            else:
                remainder_variant_ids = TrainingBatchIterator.empty_ids
            previous_remainder_buffer = cls._load_remainder_buffer(
                pickle_folder=previous_pickle_folder, is_training=is_training,
                num_columns=buffer.shape[1], delete_after_read=delete_after_read,
                sleep_time_s=sleep_time_s
            ) if has_previous_remainder else cls._get_empty_buffer(
                num_columns=buffer.shape[1]
            )
        return previous_remainder_buffer, remainder_variant_ids

    @classmethod
    def _load_remainder_buffer(
            cls,
            pickle_folder: Optional[Path],
            is_training: bool,
            num_columns: int,
            delete_after_read: bool = Default.delete_after_read,
            sleep_time_s: float = Default.sleep_time_s
    ) -> numpy.ndarray:
        if pickle_folder is None:
            return cls._get_empty_buffer(num_columns=num_columns)
        pickle_file = VcfTrainingTensorDataLoader.get_nth_pickle_file(
            pickle_folder, n=Keys.remainder, is_training=is_training
        )
        with locked_read(
                pickle_file, "rb", delete_after_read=delete_after_read, sleep_time_s=sleep_time_s
        ) as f_in:
            buffer = pickle.load(f_in)
        return buffer

    @staticmethod
    def _get_empty_buffer(num_columns: int) -> numpy.ndarray:
        return numpy.empty((0, num_columns), dtype=numpy.float32)
