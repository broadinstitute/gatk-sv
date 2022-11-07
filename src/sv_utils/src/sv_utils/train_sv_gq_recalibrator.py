#!/usr/bin/env python

import sys
import tempfile
import glob
import argparse
import warnings
import time
import json
from tqdm.auto import tqdm as tqdm
import numpy.random
import pandas
import pandas.core.arrays
import torch
import concurrent.futures

from sv_utils import common, tarred_properties_to_parquet
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import dask.dataframe
from typing import Optional, Union

tqdm.monitor_interval = 0

DaskDataFrame = dask.dataframe.DataFrame
PandasDataFrame = pandas.DataFrame
DataFrame = Union[DaskDataFrame, PandasDataFrame]


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
    temp_dir = tempfile.gettempdir()
    variants_per_batch = 100
    shuffle = True
    random_state = 0
    non_ml_properties = tarred_properties_to_parquet.non_ml_properties
    # don't exclude variant_weights: they're not fed into the neural net, but we use them in training:
    non_ml_needed_for_training = frozenset({Keys.id, Keys.variant_weights})
    excluded_properties = non_ml_properties.difference(non_ml_needed_for_training)
    torch_device = Keys.cpu
    scale_bool_properties = False


def __parse_arguments(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Train filter to recalibrate genotype quality from extracted VCF properties",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        prog=argv[0]
    )
    parser.add_argument("--properties", "-p", type=str, required=True,
                        help="full path to tarred parquet file with variant properties")
    parser.add_argument("--properties-scaling-json", "-j", type=str, required=True,
                        help="full path to JSON with baseline and scales needed for training / filtering properties")
    parser.add_argument("--output-model", "-o", type=str, required=True,
                        help="path to output file with trained model")
    parser.add_argument("--temp-dir", "-t", type=str, default=Default.temp_dir,
                        help="full path to preferred temporary directory")
    parser.add_argument("--variants-per-batch", type=int, default=Default.variants_per_batch,
                        help="number of variants used in each training batch")
    parser.add_argument("--shuffle", type=common.argparse_bool, default=Default.shuffle,
                        help="if True, shuffle order of variants while training")
    parser.add_argument("--random-seed", type=int, default=Default.random_state,
                        help="initial random seed")
    parser.add_argument("--excluded-properties", type=str, default=','.join(Default.excluded_properties),
                        help="comma-separated list of properties to not use in tensor")
    parser.add_argument("--torch-device", type=str, default=Default.torch_device, choices=[Keys.cpu, Keys.cuda],
                        help="Device on which to perform training")
    parser.add_argument("--scale-bool-properties", type=common.argparse_bool, default=Default.scale_bool_properties,
                        help="If true, z-score bool variables too. If false, leave as 0.0 or 1.0. Scaling can result in"
                             "large values if the bool variable is mostly true or mostly false.")
    return parser.parse_args(argv[1:] if len(argv) > 1 else ["--help"])


def main(argv: Optional[list[str]] = None):
    args = __parse_arguments(sys.argv if argv is None else argv)
    train_sv_gq_recalibrator(
        properties_path=args.properties,
        properties_scaling_json=args.properties_scaling_json,
        output_model_path=args.output_model,
        temp_dir=args.temp_dir,
        variants_per_batch=args.variants_per_batch,
        shuffle=args.shuffle,
        random_state=args.random_seed,
        excluded_properties=set(args.excluded_properties.split(',')),
        torch_device=args.torch_device
    )


def train_sv_gq_recalibrator(
        properties_path: str,
        properties_scaling_json: str,
        output_model_path: str,
        temp_dir: str = Default.temp_dir,
        variants_per_batch: int = Default.variants_per_batch,
        shuffle: bool = Default.shuffle,
        random_state: Union[int, numpy.random.RandomState, None] = Default.random_state,
        excluded_properties: set[str] = Default.excluded_properties,
        torch_device: str = Default.torch_device,
        scale_bool_values: bool = Default.scale_bool_properties
):
    t0 = time.time()
    num_batches = 0
    with concurrent.futures.ProcessPoolExecutor(max_workers=1) as executor:
        vcf_tensor_data_loader = VcfTensorDataLoader(
            parquet_path=properties_path, properties_scaling_json=properties_scaling_json,
            scale_bool_values=scale_bool_values, variants_per_batch=variants_per_batch, shuffle=shuffle,
            random_state=random_state, excluded_properties=excluded_properties, torch_device=torch_device,
            executor=executor
        )
        for batch_tensor, batch_weights, batch_ids in tqdm(
                vcf_tensor_data_loader, desc="batch", mininterval=0.5, maxinterval=float('inf'), smoothing=0
        ):
            if not isinstance(batch_tensor, torch.Tensor):
                raise ValueError(f"batch {num_batches} is a {type(batch_tensor)}")
            num_batches += 1
    t1 = time.time()
    print(f"Got {num_batches} batches in {t1 - t0} s")
    print(f"batch_ids:\n{batch_ids}")
    print(f"batch_weights:\n{batch_weights}")
    print(f"batch_tensor_shape:\n{batch_tensor.shape}")


def get_torch_device(device_name: str) -> torch.device:
    return torch.device("cuda:0" if device_name == Keys.cuda else "cpu")


class VcfTensorDataLoader:
    """
    A DataLoader-like object for a set of tensors that can be much faster than
    TensorDataset + DataLoader because dataloader grabs individual indices of
    the dataset and calls cat (slow).
    Source: https://discuss.pytorch.org/t/dataloader-much-slower-than-manual-batching/27014/6
    """
    __slots__ = (
        "parquet_path", "variants_per_batch", "shuffle", "random_state", "torch_device",
        "executor", "_wanted_columns", "_parquet_files", "_buffer_index", "_buffer_matrix", "_num_rows", "_current_row",
        "_variant_weights_column", "_variant_id_column", "_variant_columns", "_format_columns", "_num_samples",
        "_num_properties", "_buffer_future", "_numpy_tensor", "_weights_tensor", "_ids_array", "_tensor_baseline",
        "_tensor_scale"
    )

    def __init__(
            self,
            parquet_path: str,
            properties_scaling_json: str,
            executor: concurrent.futures.Executor,
            variants_per_batch: int = Default.variants_per_batch,
            shuffle: bool = Default.shuffle,
            scale_bool_values: bool = Default.scale_bool_properties,
            random_state: Union[int, numpy.random.RandomState, None] = Default.random_state,
            excluded_properties: set[str] = Default.excluded_properties,
            torch_device: str = Default.torch_device
    ):
        """
        Initialize a FastTensorDataLoader.
        :param *tensors: tensors to store. Must have the same length @ dim 0.
        :param batch_size: batch size to load.
        :param shuffle: if True, shuffle the data *in-place* whenever an
            iterator is created out of this object.
        :returns: A FastTensorDataLoader.
        """
        self.parquet_path = tarred_properties_to_parquet.extract_tar_to_folder(parquet_path) \
            if parquet_path.endswith(".tar") else parquet_path
        self.executor = executor
        self.variants_per_batch = variants_per_batch
        self.shuffle = shuffle
        self.random_state = random_state if isinstance(random_state, numpy.random.RandomState) \
            else numpy.random.RandomState(random_state)
        self.torch_device = get_torch_device(torch_device)
        self._set_column_info(
            excluded_properties=excluded_properties, properties_scaling_json=properties_scaling_json,
            scale_bool_values=scale_bool_values
        )
        num_format_properties = len(self._format_columns) // self._num_samples
        self._numpy_tensor = numpy.empty(
            (self.variants_per_batch, self._num_samples, len(self._variant_columns) + num_format_properties),
            dtype=numpy.float32
        )
        self._weights_tensor = numpy.empty(self.variants_per_batch, dtype=numpy.float32)
        self._ids_array = pandas.array([""] * self.variants_per_batch, dtype="string[pyarrow]")
        self._buffer_index = 0
        self._buffer_matrix = pandas.DataFrame()
        self._parquet_files = []
        self._num_rows = tarred_properties_to_parquet.get_parquet_file_num_rows(self.parquet_path)
        self._current_row = 0
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
        self._variant_weights_column = next(
            index for index, column in enumerate(unflat_columns) if column[1] == Keys.variant_weights
        )
        print(f"variant_weights_column: {self._variant_weights_column}")
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
        self._num_samples = len({column[0] for column in unflat_columns if isinstance(column[0], str)})
        # get baseline and scale for final tensor
        properties_tensor_order = [
            column[1] for column in unflat_columns
            if not isinstance(column[0], str) and column[1] not in Default.non_ml_needed_for_training
        ]
        found_properties = set()
        for sample_id, prop_name in unflat_columns:
            if prop_name in Default.non_ml_needed_for_training or not isinstance(sample_id, str):
                continue
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

    def _load_buffer(self):
        self._buffer_matrix, self.random_state = self._buffer_future.result()
        self._buffer_index = 0
        if not self._parquet_files:
            self._set_parquet_files()
        parquet_file = self._parquet_files.pop()
        remainder_rows = divmod(self._buffer_matrix.shape[0], self.variants_per_batch)[1]
        self._buffer_future = self.executor.submit(
            VcfTensorDataLoader._parquet_file_to_df_buffer_s, parquet_file=parquet_file,
            wanted_columns=self._wanted_columns, random_state=self.random_state if self.shuffle else None,
            remainder_buffer=self._buffer_matrix[-remainder_rows:, :]
        )

    def _parquet_file_to_df_buffer(self, parquet_file: str) -> pandas.DataFrame:
        df_buffer = pandas.read_parquet(parquet_file, columns=self._wanted_columns)
        tarred_properties_to_parquet.unflatten_columns(df_buffer)
        return df_buffer.iloc[self.random_state.permutation(df_buffer.shape[0])] if self.shuffle else df_buffer

    @staticmethod
    def _parquet_file_to_df_buffer_s(
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

    @staticmethod
    def _one_hot_property(df_buffer: pandas.DataFrame, property_name: str) -> pandas.DataFrame:
        """
        Should probably DELETE THIS once this script is functional: one-hot encoding best done in
        tarred_paroperties_to_parquet.py
        """
        df_property = df_buffer.loc[:, (slice(None), property_name)]
        is_variant_property = df_property.shape[1] == 1 and \
                              (df_property.columns[0][0] is None or numpy.isnan(df_property.columns[0][0]))
        categories = {category for c in df_property.columns for category in df_property[c].cat.categories}
        is_string_set = any(',' in category for category in categories)
        if is_string_set:
            categories = {category for string_set in categories for category in string_set.split(',')}
        categories = sorted(categories)

        new_columns = pandas.MultiIndex.from_tuples(
            [(None, f"{property_name}={category}") for category in categories] if is_variant_property
            else [(column[0], f"{property_name}={category}")
                  for column in df_property.columns for category in categories],
            names=["sample", "property"]
        )
        new_df_property = df_property.astype("object").apply(
            lambda s: pandas.Series(
                [e == category for e in s for category in categories] if is_variant_property
                else [category in set(e.split(',')) for e in s for category in categories],
                index=new_columns, dtype=bool
            ),
            axis=1
        )
        return pandas.concat((df_buffer.drop(columns=df_property.columns), new_df_property), axis=1)

    def _set_parquet_files(self):
        self._parquet_files = glob.glob(f"{self.parquet_path}/*.parquet")
        if self.shuffle:
            # randomize order of visiting partitions
            self._parquet_files = self.random_state.permutation(self._parquet_files).tolist()
        else:
            # reverse so we can visit in order by popping the files from the end
            self._parquet_files = self._parquet_files[::-1]

    def __iter__(self):
        if self._buffer_future is None:
            if not self._parquet_files:
                self._set_parquet_files()
            parquet_file = self._parquet_files.pop()
            self._buffer_future = self.executor.submit(
                VcfTensorDataLoader._parquet_file_to_df_buffer_s, parquet_file=parquet_file,
                wanted_columns=self._wanted_columns, random_state=self.random_state if self.shuffle else None
            )
        return self

    def __next__(self) -> (torch.Tensor, torch.Tensor, pandas.core.arrays.string_arrow.ArrowStringArray):
        """
        emit tensors for next batch
        Returns:
            properties_tensor: torch.Tensor
                num_variants x num_samples x num_properties 32-bit float tensor of training / filtering properties
            weights_tensor: torch.Tensor
                num_variants array of training weights that help balance variants
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
         self._numpy_tensor[:, :, num_variant_properties:]) = VcfTensorDataLoader._get_next_tensor_blocks(
            buffer_matrix=self._buffer_matrix[self._buffer_index:next_buffer_index, :],
            variant_ids_column=self._variant_id_column, variant_weights_column=self._variant_weights_column,
            variant_columns=self._variant_columns, format_columns=self._format_columns, num_samples=self._num_samples
        )
        self._buffer_index = next_buffer_index
        self._current_row += self.variants_per_batch

        return (
            self._tensor_scale * (
                torch.tensor(self._numpy_tensor, dtype=torch.float32, device=self.torch_device)
                - self._tensor_baseline
            ),
            torch.tensor(self._weights_tensor, dtype=torch.float32, device=self.torch_device),
            self._ids_array
        )

    @staticmethod
    def _get_next_tensor_blocks(
            buffer_matrix: numpy.ndarray,
            variant_ids_column: int,
            variant_weights_column: int,
            variant_columns: numpy.ndarray,
            format_columns: numpy.ndarray,
            num_samples: int
    ) -> (pandas.core.arrays.string_arrow.ArrowStringArray, numpy.ndarray, numpy.ndarray, numpy.ndarray):
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
            buffer_matrix.take(variant_weights_column, axis=1),
            buffer_matrix.take(variant_columns, axis=1).reshape(buffer_matrix.shape[0], 1, num_variant_properties),
            buffer_matrix.take(format_columns, axis=1)
                .reshape(buffer_matrix.shape[0], len(format_columns) // num_samples, num_samples)
                .swapaxes(1, 2)
        )  # noqa E131

    def __len__(self) -> int:
        # report the number of remaining batches in the iterator
        n_batches, extra_rows = divmod(self._num_rows - self._current_row, self.variants_per_batch)
        return n_batches + 1 if extra_rows else n_batches


if __name__ == "__main__":
    main()
