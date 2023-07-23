import ast
from pathlib import Path
import tempfile
import shutil
import attr
import numpy
import pandas
import dask.dataframe
import dask.array
from dask.delayed import Delayed
from typing import Callable, ClassVar, Optional, TypeVar
from sv_utils import common

DaskDataFrame = dask.dataframe.DataFrame
DaskSeries = dask.dataframe.Series
DaskArray = dask.array.Array
DataFrame = DaskDataFrame | pandas.DataFrame


class Keys:
    parquet_suffix = ".parquet"
    numpy_suffix = ".npy"
    sample_id = "sample_id"
    property = "property"


class Default:
    compression_algorithm = "snappy"
    # I'd like to just include this, but it causes circular imports. Need to rearchitect the
    # organization of all modules, or just encode this information twice.
    #column_levels = genomics_io.Default.column_levels
    column_levels = (Keys.sample_id, Keys.property)


@attr.frozen(slots=False, weakref_slot=False)
class TempPersister:
    """Class to help speed up dask computations by "persisting" a set of lazy calculations to disk.
    This can help if there are multiple calculation steps and values from some stage may be used
    multiple times: then you can persist and access those results quickly, rather than lazily
    computing them multiple times. Generated temporary files and folders will be removed when the
    relevant dataframes or arrays go out of scope.

    Attributes:
        temp_dir: path to store temporary folders
        engine: engine for writing parquet files for dask DataFrames
        compression_algorithm: algorithm for compressing parquet files for dask DataFrames
    """
    temp_dir: Path
    engine: str = "pyarrow"
    compression_algorithm: str = Default.compression_algorithm
    # DataType: persist_to_disk passes original object type through after persisting it.
    ArrayDataFrameOrNone: ClassVar[TypeVar] = TypeVar(
        "ArrayDataFrameOrNone", DaskDataFrame, DaskSeries, DaskArray, None
    )
    persist_folder_field: ClassVar[str] = "temp_persist_folder"

    @staticmethod
    def _add_hook(
        dask_object: ArrayDataFrameOrNone, temp_persist_folder: Path
    ) -> ArrayDataFrameOrNone:
        """Actually set the persist folder"""
        if not hasattr(dask.dataframe.DataFrame, "__del__"):
            dask.dataframe.DataFrame.__del__ = TempPersister.remove_temp_persist_folder
            dask.dataframe.Series.__del__ = TempPersister.remove_temp_persist_folder
            dask.array.Array.__del__ = TempPersister.remove_temp_persist_folder

        setattr(dask_object, TempPersister.persist_folder_field, temp_persist_folder)
        return dask_object

    def persist_to_disk(self, dask_object: ArrayDataFrameOrNone) -> ArrayDataFrameOrNone:
        """Persist dask object to disk, monkey-patch to remove temporary files when it goes out of
        scope, and the return the persisted object
        """
        if dask_object is None or hasattr(dask_object, TempPersister.persist_folder_field):
            return dask_object  # either None, or already persisted
        elif isinstance(dask_object, DaskDataFrame):
            return self._persist_dataframe_to_disk(df=dask_object)
        elif isinstance(dask_object, DaskSeries):
            return self._persist_series_to_disk(s=dask_object)
        else:
            return self._persist_array_to_disk(arr=dask_object)

    def _persist_series_to_disk(self, s: DaskSeries) -> DaskSeries:
        """Persist a dask Series. Use _persist_dataframe_to_disk and do a bit of cleanup to avoid
        duplicating the somewhat finicky dask code
        """
        df = self._persist_dataframe_to_disk(s.to_frame())
        new_s = df[s.name]
        self._add_hook(new_s, temp_persist_folder=getattr(df, TempPersister.persist_folder_field))
        # remove the df hooks so the temp folder isn't removed
        delattr(df, TempPersister.persist_folder_field)
        return new_s

    def _persist_dataframe_to_disk(self, df: DaskDataFrame) -> DaskDataFrame:
        """Persist a dask DataFrame"""
        temp_persist_folder = Path(
            tempfile.mkdtemp(
                suffix=f"{Keys.parquet_suffix}", dir=self.temp_dir
            )
        )
        need_flatten = df.columns.nlevels > 1
        if need_flatten:
            df = flatten_columns(df)
        index_name = df.index.name
        divisions = df.divisions
        df.to_parquet(
            f"{temp_persist_folder}", engine=self.engine, compression=self.compression_algorithm,
            write_metadata_file=None, write_index=(index_name is not None)
        )
        new_df = dask.dataframe.read_parquet(
            f"{temp_persist_folder}", engine=self.engine, index=index_name
        )
        if not all(division is None for division in divisions):
            new_df.divisions = divisions
        if need_flatten:
            new_df = unflatten_columns(df)
        return self._add_hook(new_df, temp_persist_folder=temp_persist_folder)

    def _persist_array_to_disk(self, arr: DaskArray) -> DaskArray:
        """Persist a dask Array"""
        temp_persist_folder = Path(
            tempfile.mkdtemp(
                suffix=f"{Keys.numpy_suffix}", dir=self.temp_dir
            )
        )
        dask.array.to_npy_stack(temp_persist_folder, arr)
        new_arr = dask.array.from_npy_stack(f"{temp_persist_folder}")
        return self._add_hook(new_arr, temp_persist_folder=temp_persist_folder)

    @staticmethod
    def remove_temp_persist_folder(
        dask_object: DaskDataFrame | DaskSeries | DaskArray,
        persist_folder_field: str = persist_folder_field
    ) -> None:
        """Remote temporary folders backing this dask object. Don't invoke directly."""
        persist_folder = getattr(dask_object, persist_folder_field, None)
        if persist_folder is not None:
            print(f"removing {persist_folder}")
            shutil.rmtree(persist_folder, ignore_errors=True)


def compute_if_unknown(dask_object) -> int:
    return dask_object.compute() if hasattr(dask_object, "compute") else dask_object


def dask_shape(dask_object) -> tuple[int, ...]:
    """Get dask shape as integers, computing them if needed"""
    return tuple(
        compute_if_unknown(s) for s in dask_object.shape
    )


def flatten_columns(df: DataFrame) -> DataFrame:
    old_columns = df.columns
    if df.columns.nlevels > 1:
        df.columns = pandas.Index(
            [f"{column}" for column in old_columns.to_flat_index()]
        )
    return df


def unflatten_column_name(flat_column_name: str) -> tuple[Optional[str], str]:
    try:
        return ast.literal_eval(flat_column_name.replace("(nan,", "(None,"))
    except ValueError as value_error:
        common.add_exception_context(value_error, f"unflattening '{flat_column_name}'")
        raise


def unflatten_columns(df: DataFrame) -> DataFrame:
    old_columns = df.columns
    try:
        multicolumns = pandas.MultiIndex.from_tuples(
            [unflatten_column_name(flat_column_name) for flat_column_name in old_columns]
        )
    except (SyntaxError, ValueError):
        # couldn't unflatten, so probably not a flattened DataFrame. just return df unchanged
        multicolumns = None
    if multicolumns is not None and multicolumns.nlevels > 1:
        if multicolumns.nlevels == 2:
            multicolumns.names = Default.column_levels
        df.columns = multicolumns
    return df


def pandas_values_to_value_counts(pandas_obj: pandas.DataFrame | pandas.Series) -> pandas.Series:
    try:
        return (
            pandas.Series([], dtype=numpy.int64) if pandas_obj.empty
            else numpy_values_to_value_counts(
                pandas_obj.values.ravel() if isinstance(pandas_obj, pandas.DataFrame)
                else pandas_obj.values
            )
        )
    except Exception as exception:
        if isinstance(pandas_obj, pandas.DataFrame):
            dtypes = {f"{pandas_obj[col].dtype}" for col in pandas_obj.columns}
            context = f"getting value counts from dataframe with dtypes: {', '.join(dtypes)}"
        else:
            context = f"getting value counts from series with dtype: {pandas_obj.dtype}"
        common.add_exception_context(exception, context)
        raise


def get_dtype_info(dt: type) -> numpy.iinfo | numpy.finfo:
    """Get int or float dtype info, whichever is appropriate"""
    numpy_dt = dt.numpy_dtype if hasattr(dt, "numpy_dtype") else dt
    return numpy.iinfo(numpy_dt) if numpy.issubdtype(numpy_dt, numpy.integer) \
        else numpy.finfo(numpy_dt)


def numpy_values_to_value_counts(values: numpy.ndarray) -> pandas.Series:
    """Get unique set of all values, and count of number of times each occurs"""
    if values.size == 0:
        return pandas.Series([], dtype=numpy.int64)
    elif numpy.issubdtype(values.dtype, numpy.integer) and get_dtype_info(values.dtype).min >= 0:
        # if original scores data was strictly non-negative ints, we can use this faster routine
        _value_counts = pandas.Series(numpy.bincount(values))
        return _value_counts.loc[_value_counts > 0]
    else:
        # use general purpose routine
        # first get rid of null scores (set to minimum allowed value by the data type)
        values = values.compress(values != numpy.iinfo(values.dtype).min)
        return pandas.Series(values).value_counts(dropna=False).sort_index()


def wanted_numpy_values_to_value_counts(
    values_block: numpy.ndarray, value_is_wanted_block: numpy.ndarray
) -> pandas.Series:
    """Return value of counts of all the values that are wanted (corresponding entry in
    mask array value_is_wanted == True). Efficiently manages dask parallelism"""
    return numpy_values_to_value_counts(
        values_block.ravel().compress(value_is_wanted_block.ravel())
    )


def wanted_dask_values_to_value_counts(
    values: DaskArray, value_is_wanted: DaskArray
) -> pandas.Series:
    """Return value of counts of all the values that are wanted (corresponding entry in
    mask array value_is_wanted == True). Efficiently manages dask parallelism"""
    block_value_counts = [
        dask.delayed(wanted_numpy_values_to_value_counts)(
            values_block, wanted_block
        )
        for values_block, wanted_block in zip(
            values.blocks, value_is_wanted.blocks
        )
    ]
    return reduce_add_value_counts(block_value_counts).compute()


TypeToReduce = TypeVar("TypeToReduce")


def reduce_delayed(
        block_values: list[TypeToReduce],
        reduction_function: Callable[[TypeToReduce, TypeToReduce], TypeToReduce],
        length: Optional[int] = None,
        empty_block_values_result: Optional[TypeToReduce] = None
) -> TypeToReduce:
    """Take list of values from different blocks/ partitions / groups / etc, and combine
    them into one values Series by using the provided reduction function.
    If the length is known in advance (and passed), this can work if the list is itself delayed,
    otherwise the list must be a regular list of delayed objects.
    """
    # ensure the reduction_function is delayed too
    try:
        if not isinstance(reduction_function, Delayed):
            reduction_function = dask.delayed(reduction_function)
    except Exception as exception:
        common.add_exception_context(
            exception, f"type(reduction_function): {type(reduction_function)}"
        )
        raise
    # this function is expected to be acting on dask Delayed objects, so rather than simply
    # summing down the list, break into tree of operations that can be conducted in parallel
    if length is None:  # didn't pass length hint, will only work if the list is not delayed
        length = len(block_values)
    if length == 0:
        if empty_block_values_result is None:
            raise RuntimeError("Got empty list but no empty_block_values_result was provided")
        else:
            return empty_block_values_result
    while length > 1:
        _start_idx = length % 2
        block_values = block_values[:_start_idx] + [
            reduction_function(block_values[_idx], block_values[_idx + 1])
            for _idx in range(_start_idx, length, 2)
        ]
        length = (length + _start_idx) // 2
    # return sole remaining value counts Series
    return block_values[0]


def add_block_value_counts(
        value_counts_1: pandas.Series, value_counts_2: pandas.Series
) -> pandas.Series:
    """Pairwise combine two sets of value counts from different blocks. Shared values have
    their counts added while unique values effectively are merged unchanged (the series
    missing that value has its counts set to 0)
    """
    return value_counts_1.add(value_counts_2, fill_value=0).sort_index().astype(numpy.int64)


def reduce_add_value_counts(
    block_value_counts: list[pandas.Series], length: Optional[int] = None
) -> pandas.Series:
    """Take list of value counts from different blocks/ partitions / groups / etc, and combine
    them into one value_counts Series by adding the counts for common values (and transferring
    new values unchanged)
    """
    return reduce_delayed(
        block_values=block_value_counts,
        reduction_function=add_block_value_counts,
        length=length,
        empty_block_values_result=dask.delayed(pandas.Series([], dtype=numpy.int64))
    )


def series_a_in_series_b(series_a: DaskSeries, series_b: DaskSeries) -> DaskSeries:
    """Return boolean DaskSeries with same index and divisions as series_a:
        True where that element of series_a is anywhere in series_b, and False otherwise
    """
    # noinspection PyTypeChecker
    return dask.dataframe.from_delayed(
        [
            reduce_delayed(
                block_values=[
                    dask.delayed(pandas.Series.isin)(_a_partition, _b_partition)
                    for _b_partition in series_b.partitions
                ],
                reduction_function=numpy.logical_or,
                length=series_b.npartitions,
                empty_block_values_result=dask.delayed(pandas.Series([], dtype=bool))
            )
            for _a_partition in series_a.partitions
        ],
        meta=pandas.Series([], dtype=bool),
        divisions=series_a.divisions
    )


Column = tuple[Optional[str], str]
MinMax = tuple[int, int]
StatsDict = dict[Column, MinMax]


def reconcile_dtypes(delayed_partitions: list[pandas.DataFrame]) -> list[pandas.DataFrame]:
    def _stats(_p: pandas.DataFrame) -> StatsDict:
        _s = {
            col: (_p[col].min(), _p[col].max())
            for col in _p.columns
            if pandas.api.types.is_integer_dtype(_p[col].dtype)
        }
        return _s

    def _reduce_stats(_s1: StatsDict, _s2: StatsDict) -> StatsDict:
        _s: StatsDict = {}
        for _k, _tup1 in _s1.items():
            _tup2 = _s2.get(_k, None)
            if _tup2 is None:
                _s[_k] = _tup1
            else:
                _s[_k] = (min(_tup1[0], _tup2[0]), max(_tup1[1], _tup2[1]))
        return _s

    overall_stats = reduce_delayed(
        [dask.delayed(_stats)(partition) for partition in delayed_partitions],
        _reduce_stats,
        length=len(delayed_partitions),
        empty_block_values_result={}
    ).compute()
    dtypes = {
        col: numpy.min_scalar_type(_max if _min >= 0 else -max(-_min, _max))
        for col, (_min, _max) in overall_stats.items()
    }

    def _convert_dtypes(
        _partition: pandas.DataFrame, _dtypes: dict[Column, numpy.dtype]
    ) -> pandas.DataFrame:
        for col, dtype in _dtypes.items():
            _partition[col] = _partition[col].astype(dtype)
        return _partition
    return [
        dask.delayed(_convert_dtypes)(partition, dtypes) for partition in delayed_partitions
    ]

