import pandas
import numpy
import warnings
from typing \
    import Optional, Callable, Any, Union, Mapping, List, Tuple, Iterator, Collection, Text, Dict, TypeVar, Iterable
from types import MappingProxyType
from tqdm.auto import tqdm
with warnings.catch_warnings():
    # pympler has some deprecated matrix stuff, I don't want warnings cluttering things up
    warnings.simplefilter("ignore")
    from pympler import asizeof

from sv_utils import common, parallel_tools, genomics_io

# type definitions used for type hints
PandasObject = Union[pandas.DataFrame, pandas.Series]
Record = Union[numpy.record, Tuple]
OverlapOutput = TypeVar("OverlapOutput")  # anything that can be jammed into pandas.DataFrame, but must be consistent
KWArgs = Mapping[Text, Any]
OverlapFuncWithoutKWArgs = Callable[[Record, pandas.DataFrame], OverlapOutput]
OverlapFuncWithKWargs = Callable[[Record, pandas.DataFrame, KWArgs], OverlapOutput]
OverlapFunc = Union[OverlapFuncWithoutKWArgs, OverlapFuncWithKWargs]
# SimpleIntervalReducer takes overlap_results, simple_intervals and reduces them to results on primary interval
SimpleIntervalReducer = Callable[[pandas.DataFrame, pandas.DataFrame, KWArgs], OverlapOutput]
EvalSlice = slice
OtherSlice = Optional[slice]
IndexChunk = Tuple[numpy.array, Optional[numpy.array]]
SliceChunk = Tuple[OtherSlice, EvalSlice]
Chunk = Union[IndexChunk, SliceChunk]
IndexIntType = numpy.int64


class Keys:
    id = genomics_io.Keys.id
    contig = genomics_io.Keys.contig
    begin = genomics_io.Keys.begin
    end = genomics_io.Keys.end
    svtype = genomics_io.Keys.svtype
    svlen = genomics_io.Keys.svlen
    bnd_contig_2 = genomics_io.Keys.bnd_contig_2
    bnd_end_2 = genomics_io.Keys.bnd_end_2
    source = genomics_io.Keys.source
    cpx_intervals = "cpx_intervals"
    primary_id = "primary_id"
    is_single_interval = "is_single_interval"


# properties that you need to compute interval overlaps
_overlap_required_properties = (Keys.contig, Keys.begin, Keys.end)


class Default:
    breakend_types = frozenset({"BND", "CTX"})
    expand_non_point_types = frozenset({"DUP"})
    ins_types = frozenset({"INS", "INS:ME", "DUP"})
    expand_point_svs_bp = 300  # when checking overlap, expand point SVs by this many BP in each direction
    point_sv_scale_factor = 1.0  # when checking overlap, expand point SVs this * SVLEN in each direction
    expand_non_point_svs_bp = 75   # when checking overlap, expand non-point SVs by this many BP in each direction
    non_point_sv_scale_factor = 0.25  # when checking overlap, expand non-point SVs this * SVLEN in each direction
    max_overlap_task_size = 1.0e7      # maximum task size for interval overlap tasks
    n_jobs = common.num_physical_cpus  # by default, use all available cores for interval overlap tasks
    missing_value = genomics_io.Default.missing_value
    stat_reducers = MappingProxyType({
        "min": min,
        "max": max,
        "mean": numpy.mean,
        "median": numpy.median
    })
    # properties that you need for computing *variant* overlaps (i.e. includes properties used for breaking down
    # variants into simple intervals
    wanted_properties = \
        _overlap_required_properties + (
            Keys.id, Keys.svlen, Keys.svtype, Keys.bnd_contig_2, Keys.bnd_end_2, Keys.cpx_intervals, Keys.source
        )


def fix_variants(
        variants: pandas.DataFrame,
        expand_point_svs_bp: int = Default.expand_point_svs_bp,
        expand_non_point_svs_bp: int = Default.expand_non_point_svs_bp,
        point_sv_scale_factor: float = Default.point_sv_scale_factor,
        non_point_sv_scale_factor: float = Default.non_point_sv_scale_factor,
        expand_non_point_types: Collection[str] = Default.expand_non_point_types,
        breakend_types: Collection[str] = Default.breakend_types,
        needed_columns: Iterable[str] = ()
) -> pandas.DataFrame:
    f"""
    "Fix" variants in preparation for computing overlaps:
        1) preprocess BNDs/complex SVs with multiple intervals, breaking them up into many simple intervals (but labeled
           so that the results can be re-combined afterwards)
        2) expand intervals (especially intervals that are points on the reference, e.g. INS)
        3) return only those columns needed for overlap calculations or merging the simple interval results back onto
           the original complex SVs
    Args:
        variants: pandas.DataFrame
            Table of SV intervals to "fix" for overlap calculations
        expand_point_svs_bp: int (default={Default.expand_point_svs_bp})
            Number of base pairs to expand length=0 intervals (in each direction)
        expand_non_point_svs_bp: int (default={Default.point_sv_scale_factor})
            Number of base pairs to expand length > 0 intervals (in each direction)
        point_sv_scale_factor: float (default={Default.point_sv_scale_factor})
            Proportion of SVLEN to expand length=0 intervals (in each direction)
        non_point_sv_scale_factor: float (default={Default.non_point_sv_scale_factor})
            Proportion of SVLEN to expand length > 0 intervals (in each direction)
        breakend_types: Collection[str] (default={Default.breakend_types})
            SVTYPEs that should be interpreted as breakends
        needed_columns: Iterable[str] (default=())
            Any extra columns that should be left in the fixed DataFrame, in case of bespoke overlap calculations
    Returns:
        fixed_df: pandas.DataFrame
            Table of simple intervals ready for basic overlap calculations
    """
    # preprocess BNDs and complex variants into simple (one-interval) variants
    variants = _preprocess_bnd_complex(variants, breakend_types=breakend_types)

    # for variants that are points on the reference, expand them (in each direction) by
    #    a) expand_point_svs_bp
    #    a) svlen * insertion scale factor if svtype == "INS" and svlen is available
    svlen = variants[Keys.svlen].copy().fillna(0) if Keys.svlen in variants else None
    is_point_sv = variants[Keys.end] - variants[Keys.begin] <= 1
    if (svlen is not None and point_sv_scale_factor > 0) or expand_point_svs_bp > 0:
        if svlen is None:
            point_sv_half_w = expand_point_svs_bp
        else:
            point_sv_half_w = (
                expand_point_svs_bp + point_sv_scale_factor * svlen[is_point_sv].abs()
            ).round().astype(numpy.int32)
        variants.loc[is_point_sv, Keys.begin] = \
            numpy.maximum(variants.loc[is_point_sv, Keys.begin], point_sv_half_w + 1) - point_sv_half_w
        variants.loc[is_point_sv, Keys.end] = variants.loc[is_point_sv, Keys.end] + point_sv_half_w
    if (svlen is not None and non_point_sv_scale_factor > 0) or expand_non_point_svs_bp > 0:
        # expand non-point SVs too
        is_expand_non_point = ~is_point_sv
        if Keys.svtype in variants and expand_non_point_types is not None:
            is_expand_non_point = is_expand_non_point & variants[Keys.svtype].isin(expand_non_point_types)
        if svlen is None:
            non_point_sv_half_w = expand_non_point_svs_bp
        else:
            non_point_sv_half_w = (
                expand_non_point_svs_bp + non_point_sv_scale_factor * svlen[is_expand_non_point].abs()
            ).round().astype(numpy.int32)
        variants.loc[is_expand_non_point, Keys.begin] = \
            numpy.maximum(variants.loc[is_expand_non_point, Keys.begin], non_point_sv_half_w + 1) \
            - non_point_sv_half_w
        variants.loc[is_expand_non_point, Keys.end] = variants.loc[is_expand_non_point, Keys.end] + non_point_sv_half_w

    #   keep only needed columns for overlap detection and reduction to original intervals
    needed_columns = {Keys.contig, Keys.begin, Keys.end, Keys.svtype, Keys.primary_id, Keys.is_single_interval}.union(
        needed_columns
    )
    keep_columns = [c for c in variants.columns if c in needed_columns]
    return variants.loc[:, keep_columns]


def _simple_interval_summary_stats_results_reducer(
        overlap_results: pandas.DataFrame,
        simple_intervals: pandas.DataFrame,
        breakend_types: Collection[str] = Default.breakend_types,
        stat_reducers: Mapping[str, Callable[[pandas.Series], float]] = Default.stat_reducers
) -> Dict[str, float]:
    f"""
    Ingest overlap results on simple intervals, and combine onto a single complex/BND SV by "reducing" them via
    supplied reducers (e.g. max, mean, etc)
    Args:
        overlap_results: pandas.DataFrame
            Table of raw overlap results (e.g. reciprocal overlap, overlap support) for simple intervals corresponding
            to one original SV
        simple_intervals: pandas.DataFrame
            Table of simple intervals used to perform overlap calculations, with extra columns indicating the SVs that
            the simple intervals originate in.
        breakend_types: Collection[str] (default={Default.breakend_types})
            Collection of SVTYPEs that will be interpreted as breakends
        stat_reducers: Mapping[str, Callable[[pandas.Series], float]] (default={Default.stat_reducers})
            Mapping of reducer name to reducer func. The reducers ingest a pandas.Series (column of overlap values) and
            emit a summary statistic (e.g. mean, max, etc)
    Returns:
        summary_overlap_stats: Dict[str, float]
            Dict with keys being summary statistic descriptions, and values being floats with overlap statistics
            appropriate for the entire SV
    """
    primary_id = simple_intervals.iloc[0][Keys.primary_id]
    svtype = simple_intervals.loc[primary_id, Keys.svtype]
    overlap_types = overlap_results.columns
    if svtype in breakend_types:
        # BND field_type, no primary interval, but there are two breakends to reduce
        primary_stats = {
            f"primary_{overlap_type}": numpy.nan
            for overlap_type in overlap_types
        }
        reduced_stats = {
            f"{key}_{overlap_type}": func(overlap_results[overlap_type])
            for key, func in stat_reducers.items()
            for overlap_type in overlap_types
        }
    else:
        # there is a primary interval, but also extra intervals to reduce
        primary_stats = {
            f"primary_{overlap_type}": overlap_results.loc[primary_id, overlap_type]
            for overlap_type in overlap_types
        }
        reduced_stats = {
            f"{key}_{overlap_type}": func(overlap_results[overlap_type])
            for key, func in stat_reducers.items()
            for overlap_type in overlap_types
        }

    return {**primary_stats, **reduced_stats}


def postprocess_bnd_complex(
        simple_interval_overlaps: pandas.DataFrame,
        simple_variants: pandas.DataFrame,
        original_index: pandas.Index,
        reduce_single_intervals: bool = True,
        simple_interval_reducer: SimpleIntervalReducer = _simple_interval_summary_stats_results_reducer,
        reducer_kwargs: KWArgs = MappingProxyType({})
) -> pandas.DataFrame:
    f"""
    Ingest overlap results on table of simple intervals, and combine onto original complex/BND SVs by "reducing" them
    via statistical reducers (e.g. max, mean, etc)
    Args:
        simple_interval_overlaps: pandas.DataFrame
            Overlap results from simple intervals
        simple_variants: pandas.DataFrame
            Table of simple variant intervals
        original_index: pandas.index
            Index of the original SVs, before being broken down into simple intervals.
        reduce_single_intervals: bool (default=True)
            If True apply statistical reduction onto all SVs.
            If False, only apply reduction to SVs with more than one simple interval. Mainly useful for saving compute
            time in situations where the reducing function is simple and produces unchanged results if there is only
            one interval.
        simple_interval_reducer: SimpleIntervalReducer (default={_simple_interval_summary_stats_results_reducer})
            Function used to perform the reduction.
        reducer_kwargs: KWArgs (default={{}})
            kwargs to pass to simple_interval_reducer
    Returns:
        svs_overlap_results: pandas.DataFrame
            Overlap results for original SVs
    """
    if reduce_single_intervals:
        is_single_interval = None
    else:
        is_single_interval = simple_variants[Keys.is_single_interval]
        if is_single_interval.all():
            # no need to reduce anything, just put in correct order
            return simple_interval_overlaps.loc[original_index]

    # reconstruct
    complex_variants = pandas.DataFrame.from_dict(
        {
            primary_id: simple_interval_reducer(
                simple_interval_overlaps.loc[variant_intervals.index], variant_intervals, **reducer_kwargs
            )
            for primary_id, variant_intervals in tqdm(
                (simple_variants if reduce_single_intervals else simple_variants.loc[~is_single_interval]).groupby(
                    Keys.primary_id, sort=False, as_index=False, group_keys=False
                ), desc="reducing simple intervals", mininterval=5.0
            )
        },
        orient="index"
    )
    return complex_variants.loc[original_index] if reduce_single_intervals else \
        genomics_io.vcat_with_categoricals(
            (complex_variants, simple_interval_overlaps.loc[simple_variants.loc[is_single_interval].index]),
        ).loc[original_index]  # noqa E123


def apply_interval_overlap_func(
        func: OverlapFunc,
        intervals_df: pandas.DataFrame,
        other_intervals_df: Optional[pandas.DataFrame] = None,
        evidence_is_unsorted: bool = True,
        exclude_self_overlap: bool = True,
        n_jobs: Optional[int] = Default.n_jobs,
        description: Optional[str] = None,
        property_names: Union[Collection[Text], Text, None] = None,
        max_task_size: float = Default.max_overlap_task_size,
        required_worker_memory: Optional[float] = None,
        required_master_memory: Optional[float] = None,
        require_physical_cpus: bool = True,
        kwargs: KWArgs = MappingProxyType({})
) -> pandas.DataFrame:
    f"""
    For each interval in intervals_df, find subset of intervals that overlap it,
    and evaluate
        func(eval_interval, overlapping_intervals, **kwargs).

    This function is finds other-overlaps (e.g. Evidence that
    overlaps with Variants) as opposed to self-overlaps (e.g. Evidence that
    overlaps with Evidence).
    Args:
        func: Callable
            It should implement
                wanted_properties = func(interval, overlap_intervals)
            Note: if a single property is desired but that property is a
            Container, func should return (property,) to disambiguate from the
            case of returning multiple samples_data.
        intervals_df: pandas.DataFrame
            Table of GenomeIntervals to evaluate (i.e. columns must contain "contig", "begin", and "end").
        other_intervals_df: Optional[pandas.DataFrame]
            Table of GenomeIntervals to search for overlaps with intervals_df. If not provided, use intervals_df.
        evidence_is_unsorted: bool (Default=True)
            Evidence is assumed to be unsorted unless this flag is changed.
        exclude_self_overlap: bool (Default=True)
            If other_intervals_df is not provided (i.e. searching for overlaps between elements of intervals_df),
            do not consider an individual interval to overlap with itself.
        n_jobs: int or None (Default=None)
            Number of parallel workers. If n_jobs in [-1, None] then use maximum
            number supported by system (subject to constraints on cluster usage
            and memory)
        description: str or None (Default=None)
            Description of job for progress bar. If None, defaults to
            'overlap of [function prop_name]'
        property_names: list-like or None
            List of names of samples_data returned by func. If None, will be
            inferred if func returns
                mapping (e.g. dict)
                pandas.Series
                object with ._fields  (e.g. namedtuple)
                numpy.recarray
        max_task_size: float (Default={Default.max_overlap_task_size})
            Break up tasks so that no one task is larger than this
        required_worker_memory: float (Default=None)
            Used to constrain number of parallel workers based on available
            memory. Amount of memory needed by master process.
            If empty,
        required_master_memory: float (Default=None)
            Used to constrain number of parallel workers based on available
            memory. Amount of memory needed by each worker process.
        require_physical_cpus: bool (Default=True)
            If True, limit number of jobs to number of physical cores in system, not number of hyperthreads.
        kwargs: KWArgs: (Default=empty dict())
            list of keyword arguments passed to func
    Returns:
        results: pandas.DataFrame
            table of results, with rows corresponding to intervals_df, and
            columns corresponding to samples_data
    """
    # signal to take action to avoid self-overlaps only if it's a possibility
    if len(intervals_df) == 0:
        raise ValueError("Input intervals are empty")
    missing_columns = set(_overlap_required_properties).difference(intervals_df.columns)
    if missing_columns:
        raise ValueError(f"Input intervals are missing required overlap column(s): {','.join(missing_columns)}")
    exclude_self_overlap = exclude_self_overlap and other_intervals_df is None
    interval_sorter = _IntervalSorter(intervals_df, evidence_needs_sort=evidence_is_unsorted)
    intervals_df = interval_sorter.sorted
    other_intervals_df = intervals_df if other_intervals_df is None \
        else _IntervalSorter(other_intervals_df, evidence_needs_sort=evidence_is_unsorted).sorted
    missing_columns = set(_overlap_required_properties).difference(other_intervals_df.columns)
    if missing_columns:
        raise ValueError(f"Input other intervals are missing required overlap column(s): {','.join(missing_columns)}")

    tasks, task_sizes, max_task_intervals, max_task_other_intervals = _make_interval_overlap_tasks(
        intervals_df, other_intervals_df,
        max_task_size=max_task_size, require_physical_cpus=require_physical_cpus, n_jobs=n_jobs
    )

    func_kwargs = {
        "func": func, "kwargs": kwargs, "property_names": property_names, "exclude_self_overlap": exclude_self_overlap
    }
    # place reasonable guesses on memory usage to limit n_jobs
    if required_master_memory is None:
        num_props = \
            1 if property_names is None or isinstance(property_names, Text) \
            else len(property_names)
        # assume num_props floats + one 64-bit index
        prop_size = (num_props + 1) * 8.0 / 2.0 ** 30
        # required final memory is size of samples_data (in GB) * _num_rows
        required_master_memory = prop_size * len(intervals_df)
    if required_worker_memory is None:
        # each worker holds a fraction of the data set, final results,
        # and func_kwargs
        kwargs_size = asizeof.asizeof(func_kwargs)
        # put in 5x safety factor
        mem_use_fraction = max_task_intervals / max(1, len(intervals_df))
        mem_use_other = max_task_other_intervals / max(1, len(other_intervals_df))
        required_worker_memory = (
            kwargs_size +
            mem_use_fraction * (
                required_master_memory +
                intervals_df.memory_usage(index=True, deep=True).sum()
            ) +
            mem_use_other * (
                + other_intervals_df.memory_usage(index=True, deep=True).sum()
            )
        ) / 2.0 ** 30

    if description is None:
        description = 'overlap of %s' % func.__name__  # pragma: no cover
    results_itr = parallel_tools.parmap(
        _overlap_eval_func, tasks, task_sizes=task_sizes, ordered=False,
        kwargs=func_kwargs, starmap=True, description=description, n_jobs=n_jobs,
        required_master_memory=required_master_memory,
        required_worker_memory=required_worker_memory,
        require_physical_cpus=require_physical_cpus
    )
    results = pandas.concat(results_itr, axis=0).loc[intervals_df.index]

    return interval_sorter.unsort(results)


def get_connected_components_from_contig_intervals(
        contig_intervals_df: pandas.DataFrame,
        evidence_needs_sort: bool = True,
        allowed_gap: int = 0,
) -> Iterator[pandas.DataFrame]:
    """
    From sorted list of intervals, find connected components (maximal lists of
    intervals that cannot be divided without separating overlapping intervals)
    Args:
        contig_intervals_df: pandas.DataFrame
            Table of intervals from one contig (with index, begin, end)
        evidence_needs_sort: bool (Default=True)
            evidence_df must be sorted key=(begin, end) for this algorithm. If it
            is already sorted, this flag can be set to False for speed-up.
        allowed_gap: int (Default=0)
            minimum gap between end of one interval and beginning of another
            that does not connect component
    Yields:
        connected_components: pandas.DataFrame
            intervals in each connected component sorted by (begin, end)
    """
    if len(contig_intervals_df) == 0:
        return

    if evidence_needs_sort:
        contig_intervals_df = contig_intervals_df.sort_values([Keys.begin, Keys.end])

    max_end = numpy.maximum.accumulate(contig_intervals_df[Keys.end].values)[:-1]
    start_inds = 1 + numpy.nonzero(contig_intervals_df[Keys.begin].values[1:] >= max_end + allowed_gap)[0]
    if start_inds.size:
        yield contig_intervals_df.iloc[:start_inds[0]]
        for n in range(len(start_inds) - 1):
            start = start_inds[n]
            end = start_inds[n + 1]
            yield contig_intervals_df.iloc[start:end]
        yield contig_intervals_df.iloc[start_inds[-1]:]
    else:
        yield contig_intervals_df


"""
####################################### Only private classes and methods below #########################################
"""


class _IntervalSorter:
    """
    The overlap functions here work by iterating through appropriately sorted DataFrames, however, users will expect
    results in the original order of the input data, and with the same index. This helper class handles sorting input
    data and then restoring the original order and index for the output.
    """
    __slots__ = ("df", "index", "unsort_ind")
    df: pandas.DataFrame
    index: pandas.Index
    unsort_ind: Union[numpy.ndarray, slice]

    def __init__(self, df: pandas.DataFrame, evidence_needs_sort: bool = True):
        self.index = df.index
        df = df.reset_index(drop=True)
        if evidence_needs_sort:
            df = df.sort_values([Keys.contig, Keys.begin, Keys.end])
            sort_ind = numpy.array(df.index)
            self.df = df.reset_index(drop=True)
            self.unsort_ind = numpy.empty_like(sort_ind)
            self.unsort_ind.put(sort_ind, numpy.arange(sort_ind.size))
        else:
            self.unsort_ind = slice(None)
            self.df = df

    @property
    def sorted(self) -> pandas.DataFrame:
        return self.df

    def unsort(self, other_df: pandas.DataFrame) -> pandas.DataFrame:
        other_df = other_df.iloc[self.unsort_ind]
        other_df.index = self.index
        return other_df


def _check_categorical(new_value: str, series: pandas.Series) -> str:
    """
    Check if a new value is a previously-seen category value for given Series. If not, add it to the Series' categories
    Args:
        new_value: str
        series: pandas.Series
    Returns:
        new_value
    """
    cat = series.cat
    if new_value not in cat.categories:
        cat.add_categories(new_value)
    return new_value


def _get_variants_distal_source(variants: pandas.DataFrame, source: pandas.Series, id_suffix: str) -> pandas.DataFrame:
    """
    Extract new genomic intervals from a source field (generally "source" or "cpx_intervals"). Get the interval location
    from the source, and other info from the original variants. For IDs for the new intervals with id_suffix.
    Args:
        variants: pandas.DataFrame
            table of original variant info
        source: pandas.Series
            series with text that can be parsed to give genomic loci
        id_suffix: str
            suffix to add to variant IDs to produce new interval IDs
    Returns:
        distal_src: pandas.DataFrame
            table of new intervals extracted from source
    """
    _one = numpy.int32(1)

    distal_src = variants.copy()

    def _split_source(_s):
        _s0 = _s
        try:
            _svtype, _s = _s.split("_", 1)
            _contig, _s = _s.split(":", 1)
            _begin, _end = _s.split("-", 1)
        except ValueError as value_error:
            common.add_exception_context(value_error, f'error processing source "{_s0}"')
            raise
        _check_categorical(_contig, variants[Keys.contig])
        _check_categorical(_svtype, variants[Keys.svtype])
        return _svtype, _contig, max(_one, numpy.int32(_begin)), numpy.int32(_end)

    distal_src[Keys.svtype], distal_src[Keys.contig], distal_src[Keys.begin], distal_src[Keys.end] = zip(
        *(_split_source(_s) for _s in source)
    )
    distal_src.index = distal_src.index + id_suffix
    distal_src[Keys.contig] = distal_src[Keys.contig].astype("category")
    distal_src[Keys.svtype] = distal_src[Keys.svtype].astype("category")
    _check_valid_locations(distal_src)
    return distal_src


def _is_not_missing(series: pandas.Series, missing_value: str) -> numpy.ndarray:
    """
    return numpy array of bools that corresponds to whether each element of series is present (True) or missing (False)
    """
    return numpy.logical_not(genomics_io.is_missing(series, missing_value=missing_value))


def _check_valid_locations(
        variants: pandas.DataFrame,
        non_null_columns: Iterable[str] = (Keys.begin, Keys.end, Keys.contig, Keys.svtype, Keys.svlen)
):
    """ Validate that variants table has all required fields defined with no missing values """
    for column in non_null_columns:
        if column not in variants:
            raise ValueError(f"{column} not in variants")
        if variants[column].hasnans:
            raise ValueError(f"{column} has missing values")


def _preprocess_bnd_complex(
        variants: pandas.DataFrame,
        breakend_types: Collection[str] = Default.breakend_types,
        missing_value: str = Default.missing_value
) -> pandas.DataFrame:
    f"""
    Divide SVs with multiple intervals into individual simple intervals, and add fields {Keys.primary_id}, and
    {Keys.is_single_interval} to enable combining overlap results from the different simple intervals to get an overall
    result for the entire SV.
    Args:
        variants: pandas.DataFrame
            Raw table of variant info for SVs
        breakend_types: Collection[str] (default={Default.breakend_types})
            SVTYPEs that correspond to breakends
        missing_value: str (Default={Default.missing_value}
            Code for missing categorical/string values
    Returns:
        simple_variant_intervals: pandas.DataFrame
            Table of variant simple intervals
    """
    if variants.columns.nlevels > 1:
        # variants has multi-level columns. Overlap calculations operate at the variant level, not per-sample genotypes.
        # drop sample info and proceed
        variants = variants.xs(None, axis=1, level=genomics_io.Keys.sample_id)
    _check_valid_locations(variants)
    # add primary id to variants and move it to be the first column
    variants = variants.assign(**{Keys.primary_id: variants.index})
    variants = variants[[variants.columns[-1]] + variants.columns[:-1].to_list()]
    # ensure that certain important columns are categoricals
    for key in (Keys.contig, Keys.svtype):
        if not pandas.api.types.is_categorical_dtype(variants[key]):
            variants[key] = variants[key].astype("category")

    is_bnd = variants[Keys.svtype].isin(breakend_types).values if Keys.svtype in variants \
        else common.false(len(variants))
    is_cpx = _is_not_missing(variants[Keys.cpx_intervals], missing_value) if Keys.cpx_intervals in variants \
        else common.false(len(variants))
    is_src = _is_not_missing(variants[Keys.source], missing_value) & ~is_cpx if Keys.source in variants \
        else common.false(len(variants))
    try:
        is_single_interval = ~(is_bnd | is_src | is_cpx)
    except ValueError as value_error:
        context = \
            f"is_bnd.shape={is_bnd.shape}, " \
            f"svtype.shape={variants[Keys.svtype].shape if Keys.svtype in variants else None}, "\
            f"is_cpx.shape={is_cpx.shape}, "\
            f"cpx_intervals.shape={variants[Keys.cpx_intervals].shape if Keys.cpx_intervals in variants else None}, "\
            f"is_src.shape={is_src.shape}, "\
            f"source.shape={variants[Keys.source].shape if Keys.source in variants else None}"
        common.add_exception_context(value_error, context)
        raise
    variants = variants.assign(**{Keys.is_single_interval: is_single_interval})
    if is_single_interval.all():
        return variants

    distal_bnd = variants.loc[is_bnd].copy()
    if is_bnd.any():
        distal_bnd[Keys.contig] = [_check_categorical(contig, variants[Keys.contig])
                                   for contig in distal_bnd[Keys.bnd_contig_2]]
        distal_bnd[Keys.contig] = distal_bnd[Keys.contig].astype("category")
        if Keys.bnd_end_2 in distal_bnd:
            distal_bnd[Keys.end] = distal_bnd[Keys.bnd_end_2]
        else:
            variants.loc[is_bnd, Keys.end] = variants.loc[is_bnd, Keys.begin] + 1
        distal_bnd[Keys.begin] = distal_bnd[Keys.end] - 1
        distal_bnd.index = distal_bnd.index.astype("object") + "_distal"
        variants.loc[is_bnd, Keys.begin] = variants.loc[is_bnd, Keys.begin]
        variants.loc[is_bnd, Keys.end] = variants.loc[is_bnd, Keys.end]
        _check_valid_locations(distal_bnd)
        concat_intervals = (variants, distal_bnd)
    else:
        concat_intervals = (variants,)

    if is_src.any():
        concat_intervals = concat_intervals + (
            _get_variants_distal_source(
                variants.loc[is_src], variants.loc[is_src, Keys.source], "_source"
            ),
        )

    if is_cpx.any():
        cpx_index = 0
        cpx_intervals = variants[Keys.cpx_intervals].copy().astype("object")
        is_cpx_this = is_cpx.copy()
        while is_cpx_this.any():
            cpx_intervals_this = cpx_intervals[is_cpx_this].apply(lambda i: i[0])
            concat_intervals = concat_intervals + (
                _get_variants_distal_source(
                    variants.loc[is_cpx_this], cpx_intervals_this, f"_cpx_{cpx_index}"
                ),
            )
            cpx_index += 1
            cpx_intervals = cpx_intervals.apply(
                lambda i: i[1:] if isinstance(i, tuple) and len(i) > 1 else missing_value
            )
            is_cpx_this = _is_not_missing(cpx_intervals, missing_value)

    # join the various simple intervals together
    # sort by contig, begin, end
    # then drop unnecessary columns
    return genomics_io.sort_intervals_table(
        genomics_io.vcat_with_categoricals(concat_intervals), drop_index=False, inplace=False
    ).drop({Keys.bnd_contig_2, Keys.bnd_end_2, Keys.source, Keys.cpx_intervals}.intersection(variants.columns),
           inplace=False, axis=1)


def _make_interval_overlap_tasks(
        intervals_df: pandas.DataFrame,
        other_intervals_df: pandas.DataFrame,
        max_task_size: float,
        n_jobs: Optional[int] = Default.n_jobs,
        require_physical_cpus: bool = True
) -> (Iterator[Tuple[pandas.DataFrame, Optional[pandas.DataFrame]]], Tuple[int], int, int):
    f"""
    Return generator of tasks for parallel execution, and list of estimated task
    sizes. Tasks will typically be connected components (maximal lists of
    intervals that cannot be divided without separating overlapping intervals)
    from intervals_df; however components may be lumped together if they are very
    small, or split into multiple tasks if they are very large.

    Args:
        intervals_df: pandas.DataFrame
            Table of GenomeIntervals to evaluate.
        other_intervals_df: pandas.DataFrame
            Table of GenomeIntervals to search for overlaps with intervals_df.
        max_task_size: float
            Break up tasks so that no one task is larger than this
        n_jobs: int or None (Default={Default.n_jobs})
            Number of parallel workers. If n_jobs in [-1, None] then use maximum
            number supported by system (subject to constraints on cluster usage
            and memory)
        require_physical_cpus: bool (Default=True)
            If True, limit number of jobs to number of physical cores in system, not number of hyperthreads.
    Returns:
        tasks: generator[tuples]
            Generator of tasks (check_intervals, eval_intervals) to evaluate
        task_sizes: list[int]
            List of estimated task evaluation time (arbitrary units)
        max_num_intervals: int
            Maximum number of intervals passed for any task
        max_num_other_intervals: int
            Maximum number of other_intervals passed for any task
    """
    num_contigs = len(intervals_df[Keys.contig].unique())
    if num_contigs == 0:
        return (([], []) for __ in []), [], 0, 0

    needed_columns = [Keys.contig, Keys.begin, Keys.end]

    contig_iter: Iterator[Tuple[pandas.DataFrame, pandas.DataFrame]] = (
        (
            contig_intervals.drop(Keys.contig, axis=1),
            other_intervals_df.loc[
                other_intervals_df[Keys.contig] == contig, [Keys.begin, Keys.end]
            ]
        )
        for contig, contig_intervals in intervals_df[needed_columns].groupby(
            Keys.contig, sort=False, as_index=False, group_keys=False
        )
    )

    check_slices, eval_slices, task_sizes = zip(
        *parallel_tools.parmap(
            _merge_overlap_task_chunks,
            tasks=contig_iter, num_tasks=num_contigs, starmap=True, flatmap=True,
            ordered=False, description="making tasks",
            kwargs={"max_task_size": max_task_size},
            n_jobs=n_jobs, require_physical_cpus=require_physical_cpus
        )
    )

    max_num_other_intervals = max(s.stop - s.start + 1 for s in check_slices)
    max_num_intervals = max(s.stop - s.start + 1 for s in eval_slices)

    tasks = ((other_intervals_df.loc[check_slice].drop(Keys.contig, axis=1),
              intervals_df.loc[eval_slice].drop(Keys.contig, axis=1))
             for check_slice, eval_slice in zip(check_slices, eval_slices))
    return tasks, task_sizes, max_num_intervals, max_num_other_intervals


class _SliceFactory:
    """
    Helper class that builds the minimal index slice that contains all the required indices (which are added via the
    .append() function)
    """
    __slots__ = ("start", "stop")
    start: int
    stop: int

    def __init__(self, *args, require_contiguous: bool = True):
        if len(args) < 1:
            # null slice
            self.start = 1
            self.stop = 0
        elif isinstance(args[0], numpy.ndarray):
            if len(args) > 1:
                raise TypeError("SliceFactory was passed an extra argument after indices")
            indices = args[0]
            if indices.size == 0:
                # empty indices = null slice
                self.start = 1
                self.stop = 0
            else:
                if require_contiguous:
                    _SliceFactory.assert_contiguous(indices)
                self.start = indices[0]
                self.stop = indices[-1]
        else:
            if len(args) > 2:
                raise TypeError("SliceFactory was passed more than 2 arguments")
            self.start = args[0]
            self.stop = args[1]

    def __len__(self):
        return max(self.stop - self.start + 1, 0)

    def __str__(self) -> str:
        return f"SliceFactory({self.start}, {self.stop})"

    def __repr__(self) -> Dict[str, int]:
        return {"start": self.start, "stop": self.stop}

    @staticmethod
    def assert_contiguous(indices: numpy.ndarray):
        """
        Throw a ValueError if indices are not contiguous
        Args:
            indices: numpy.ndarray
        """
        if len(indices) == 0:
            return
        try:
            if not numpy.array_equal(indices, numpy.arange(indices[0], indices[0] + len(indices))):
                raise ValueError("indices not contiguous")
        except Exception as err:
            common.add_exception_context(err, f"indices is {indices}")
            raise

    @property
    def is_empty(self) -> bool:
        return self.stop < self.start

    def is_contiguous(self, indices: numpy.ndarray) -> bool:
        """
        Return true if the new indices would be contiguous with the existing factory
        """
        return indices.size == 0 or self.is_empty or (indices[0] <= self.stop + 1 and indices[-1] >= self.start - 1)

    def append(self, indices: numpy.array, require_contiguous: bool = True):
        if indices.size == 0:
            # no index was added
            return
        elif self.is_empty:
            self.start = indices[0]
            self.stop = indices[-1]
        else:
            if require_contiguous and (indices[0] > self.stop + 1 or indices[-1] < self.start - 1):
                raise ValueError(
                    f"index is not sequential with previous chunk\n{indices}\nslice({self.start, self.stop})"
                )
            self.start = min(self.start, indices[0])
            self.stop = max(self.stop, indices[-1])

    @property
    def slice(self) -> slice:
        return slice(self.start, self.stop)


def _merge_overlap_task_chunks(
        contig_intervals: pandas.DataFrame,
        other_intervals: pandas.DataFrame,
        max_task_size: float
) -> Tuple[Tuple[slice, slice, int], ...]:
    f"""
    Join very small chunks into a larger chunk. Split very large chunks into smaller chunks.
    This should aid in load balancing for parallel execution.
    Args:
        contig_intervals: pandas.DataFrame
            Table of intervals on this contig with columns {Keys.begin} and {Keys.end}. These are the intervals of
            interest, and the final output results should have the same index.
        other_intervals: pandas.DataFrame
            Table of intervals to search for overlaps with connected_components.
        max_task_size: float (Default=1.0e6)
            Break up tasks so that no one task is larger than this
    Returns:
        tasks: Tuple[Tuple[slice, slice, int], ...]:
            Tuple of tasks, each task being a slice into other intervals , a slice into contig_intervals, and an
            associated task size roughly proportional to execution time. slices are passed to DataFrame.loc[] to select
            data that should be passed to workers for parallel evaluation.
    """
    connected_components = get_connected_components_from_contig_intervals(
        contig_intervals, evidence_needs_sort=False
    )

    # store check/other chunk, eval chunk, and task size in list of tuples
    slice_chunks: List[Tuple[_SliceFactory, _SliceFactory, int]] = []
    prev_eval_chunk = _SliceFactory()
    prev_check_chunk = _SliceFactory()
    other_index = other_intervals.index.values
    other_begin = other_intervals[Keys.begin].values
    other_end = other_intervals[Keys.end].values
    prev_size = 0
    for component in connected_components:
        overlapper_coverage = _get_overlapper_coverage(other_intervals, component)
        task_size = _get_overlap_task_size(other_intervals, component, overlapper_coverage=overlapper_coverage)
        if not prev_eval_chunk.is_empty:
            # previous chunk + this one is still small, group them together
            if task_size + prev_size < max_task_size:
                check_component_index = other_index.take(_get_other_check_chunk(component, other_begin, other_end))
                # Note: no way to guarantee check chunks are contiguous.
                prev_check_chunk.append(check_component_index, require_contiguous=False)
                prev_eval_chunk.append(component.index.values)
                prev_size += task_size
                continue
            else:
                # previous small chunk can't be grouped with this one, add it now
                slice_chunks.append((prev_check_chunk, prev_eval_chunk, prev_size))
        n_chunk_tasks = int(numpy.ceil(task_size / max_task_size))
        if n_chunk_tasks <= 1:
            # this chunk is small, save it for later to see if it can be grouped
            check_component_index = other_index.take(_get_other_check_chunk(component, other_begin, other_end))
            prev_check_chunk = _SliceFactory(check_component_index, require_contiguous=False)
            prev_eval_chunk = _SliceFactory(component.index.values)
            prev_size = task_size
        else:
            # be conservative
            n_chunk_tasks = min(
                int(numpy.ceil(task_size / (0.8 * max_task_size))),
                component.size
            )
            prev_eval_chunk = _SliceFactory()
            prev_check_chunk = _SliceFactory()
            prev_size = 0

            for eval_chunk in _split_chunk(component, n_chunk_tasks):
                other_chunk_index = other_index.take(
                    _get_other_check_chunk(eval_chunk, other_begin, other_end)
                )
                slice_chunks.append(
                    (_SliceFactory(other_chunk_index, require_contiguous=False), _SliceFactory(eval_chunk.index.values),
                     _get_overlap_task_size(other_intervals, eval_chunk, overlapper_coverage=overlapper_coverage))
                )

    if not prev_eval_chunk.is_empty:
        slice_chunks.append((prev_check_chunk, prev_eval_chunk, prev_size))

    return tuple(
        (other_slice.slice, eval_slice.slice, task_size) for other_slice, eval_slice, task_size in slice_chunks
    )


def _get_overlapper_coverage(other_intervals: pandas.DataFrame, component: pandas.DataFrame) -> float:
    """ Estimate the average depth of coverage that other_intervals have over the connected component """
    if len(other_intervals) == 0:
        return 0.0
    component_begin = component[Keys.begin].values[0]
    component_end = component[Keys.end].values[-1]
    other_overlap = numpy.nonzero(
        numpy.logical_and(
            other_intervals[Keys.end].values > component_begin,
            other_intervals[Keys.begin].values < component_end
        )
    )[0]
    if other_overlap.size == 0:
        return 0.0
    # note: other intervals are sorted by (begin, end). no guarantee on where max end is
    other_begin = numpy.maximum(other_intervals[Keys.begin].values.take(other_overlap), component_begin)
    other_end = numpy.minimum(other_intervals[Keys.end].values.take(other_overlap), component_end)
    return (other_end - other_begin).astype(numpy.uint64).sum() / max(float(component_end - component_begin), 1.0)


def _get_overlap_task_size(
        other_intervals: pandas.DataFrame,
        component: pandas.DataFrame,
        overlapper_coverage: Optional[float] = None
) -> int:
    """
    Estimated task duration (arbitrary units) to evaluate self-overlap function
    on this component/chunk.
    Args:
        other_intervals: numpy.recarray
            recarray of intervals that may align with intervals in component
        component: numpy.recarray
            recarray of intervals to evaluate
    Return:
        task_size: int
            Estimated time (arbitrary units) to calculate self-overlap on this
            component
    """
    if overlapper_coverage is None:
        overlapper_coverage = _get_overlapper_coverage(other_intervals, component)
    avg_num_overlap = overlapper_coverage * numpy.mean(component[Keys.end].values - component[Keys.begin].values)
    return max(int(round(avg_num_overlap * len(component))), 1)


def _split_chunk(chunk: pandas.DataFrame, n_chunk_tasks: int) -> Iterator[pandas.DataFrame]:
    """
    Split chunk into requested number of smaller chunks for evaluation. Yield
    info for each sub-chunk sequentially.
    Args:
        chunk: list(tuple(int))
            List of intervals (begin, end, index)
        n_chunk_tasks: int
            number of sub-chunks to split chunk into (for evaluation)
    Yields:
        eval_chunk: list(tuple(int))
            List of intervals to evaluate in this sub-chunk
    """
    _SliceFactory.assert_contiguous(chunk.index.values)

    chunk_len = len(chunk)
    chunk_divs = numpy.arange(chunk_len + 1) if n_chunk_tasks >= chunk_len \
        else numpy.round(numpy.linspace(0, chunk_len, n_chunk_tasks + 1)).astype(int)
    for k in range(chunk_divs.size - 1):
        div_start = chunk_divs[k]
        div_stop = chunk_divs[k + 1]
        yield chunk.iloc[div_start:div_stop]


def _get_other_check_chunk(
        eval_chunk: pandas.DataFrame,
        other_begin: numpy.ndarray,
        other_end: numpy.ndarray
) -> numpy.ndarray:
    """
    Find appropriate check_chunk from other_begin: subset that overlaps the interval [min(eval_chunk), max(eval_chunk))
    NOTE: there's no way to sort other intervals to guarantee that the overlapping indices are contiguous
    Params:
        eval_chunk: pandas.DataFrame
            Table of intervals to evaluate
        other_begin: numpy.ndarray
            array of .begin from intervals that may overlap with eval_chunk.
        other_end: numpy.ndarray
            array of .end from intervals that may overlap with eval_chunk.
    Returns:
         check_chunk: numpy.array(int64)
            Index of intervals in other_begin that may overlap with intervals in eval_chunk
    """
    eval_begin = eval_chunk[Keys.begin].values[0]
    eval_end = eval_chunk[Keys.end].values.max()
    check_end = other_begin.searchsorted(eval_end, side="left")
    wanted = numpy.nonzero(other_end[:check_end] > eval_begin)[0]
    return wanted


def _overlap_eval_func(
        check_intervals: pandas.DataFrame,
        eval_intervals: pandas.DataFrame,
        func: OverlapFunc,
        exclude_self_overlap: bool,
        property_names: Union[Collection[Text], Text, None] = None,
        kwargs: KWArgs = MappingProxyType({})
) -> pandas.DataFrame:
    """
    Helper function run by workers. For each interval in eval_intervals, it
    finds the subset of intervals in check_intervals that overlap it, evaluating
        func(eval_interval, overlapping_intervals, **kwargs).
    Args:
        check_intervals: pandas.DataFrame
            Table of Intervals that may overlap Intervals of interest.
            Should be sorted by (begin, end)
        eval_intervals: pandas.DataFrame
            Table of Intervals to evaluate.
            Should be sorted by (begin, end)
        func: Callable
            Function that implements
            wanted_props = func(interval, overlapping_intervals, **kwargs)
        exclude_self_overlap: bool
            If other_intervals_df is not provided (i.e. searching for overlaps between elements of intervals_df),
            do not consider an individual interval to overlap with itself.
        property_names: list-like, str, or None
            List of names of samples_data returned by func. If None, will be
            inferred if func returns pandas.Series object, or object with
            ._fields
        kwargs: KWArgs (default = empty dict())
            Keyword args to be passed to func.
    Returns:
        results: pandas.DataFrame
            table of results, with rows corresponding to eval_intervals, and
            columns corresponding to samples_data
    """
    results_itr = _gen_overlap_results(
        check_intervals, eval_intervals, func, exclude_self_overlap=exclude_self_overlap, kwargs=kwargs
    )
    return _pack_overlap_results(results_itr, eval_intervals.index, property_names)


def _gen_overlap_results(
        check_intervals: pandas.DataFrame,
        eval_intervals: pandas.DataFrame,
        func: OverlapFunc,
        exclude_self_overlap: bool,
        kwargs: Mapping = MappingProxyType({})
) -> OverlapOutput:
    """
    Iterate through eval_intervals and find intervals from check_intervals that
    overlap with them. Yield the result of func(interval, overlappers) for
    each interval in eval_intervals.
    Args:
        check_intervals: pandas.DataFrame
            Table of Intervals that may overlap Intervals of interest.
            Should be sorted by (begin, end)
        eval_intervals: pandas.DataFrame
            Table of Intervals to evaluate.
            Should be sorted by (begin, end)
        func: Callable
            Function that implements
            wanted_props = func(interval, overlapping_intervals)
        exclude_self_overlap: bool
            If other_intervals_df is not provided (i.e. searching for overlaps between elements of intervals_df),
            do not consider an individual interval to overlap with itself.

    Yields:
        func_result: OverlapOutput
            Properties obtained from
                func(interval, overlappers, **kwargs)
    """
    if len(check_intervals) == 0:
        # no overlapping check intervals, just evaluate func over empty check_intervals for each eval_interval:
        for eval_interval in eval_intervals.itertuples():
            yield func(eval_interval, check_intervals, **kwargs)
        return

    check_begin = check_intervals[Keys.begin].values
    check_end = check_intervals[Keys.end].values
    check_index = check_intervals.index.values if exclude_self_overlap else None
    last_begin = -1
    # keep track of indices into check_intervals:
    # -left is index of 1st check_interval with .end > eval_interval.begin
    # -mid is index of 1st check_interval with .begin >= eval_interval.begin
    #   (hence all subsequent intervals have .end > eval_interval.begin)
    # -right is index of 1st check_interval with .begin >= eval_interval.end
    # -front is index array of check_intervals that overlap eval_inteval
    #  but have .begin < eval_interval.begin
    left = numpy.int64(0)
    mid = numpy.int64(0)
    right = numpy.int64(0)
    front = None
    for eval_interval in eval_intervals.itertuples():
        # noinspection PyUnresolvedReferences
        i_begin = eval_interval.begin
        if i_begin > last_begin:
            # this interval has a larger .begin than the previous
            last_begin = i_begin  # update last_begin
            # .end may have moved forward OR backward, so check for updates to
            # right starting from previous location of .begin, which is safe:
            # noinspection PyUnresolvedReferences
            right = mid + check_begin[mid:].searchsorted(eval_interval.end, side='left')
            if right > mid:
                # mid is between previous value and new right
                mid += check_begin[mid:right + 1].searchsorted(i_begin, side='left')
            if mid > left:
                # check for updates to left if mid is greater than left: find first check end bigger than i_begin
                delta_left = numpy.nonzero(check_end[left:mid] > i_begin)[0]
                if delta_left.size:
                    left += delta_left[0]
                else:
                    left = mid
            # front is array of indices to check_intervals with .begin
            front = numpy.arange(left, mid).compress(
                check_end[left:mid] > i_begin
            )
        else:
            # this interval has same .begin as previous, thus only .end can have increased.
            # Check for updates to right:
            # noinspection PyUnresolvedReferences
            right += check_begin[right:].searchsorted(
                eval_interval.end, side='left'
            )
        if exclude_self_overlap:
            # noinspection PyUnresolvedReferences
            self_ind = check_index.searchsorted(eval_interval.Index, side="left")
            overlap_inds = numpy.concatenate((
                front, numpy.arange(mid, self_ind), numpy.arange(self_ind + 1, right)
            ))
        else:
            overlap_inds = numpy.concatenate((
                front, numpy.arange(mid, right)
            ))
        yield func(eval_interval, check_intervals.iloc[overlap_inds], **kwargs)


def _pack_overlap_results(
        results: Iterator[OverlapOutput],
        index: pandas.Index,
        property_names: Union[Collection[Text], Text, None]
) -> pandas.DataFrame:
    """
    Consume generated results from one chunk and pack into pandas.DataFrame
    Args:
        results: Iterator
            Iterator that yields results of overlap func on sequential intervals
        index: pandas.core.indexes.base.Index
            Index into original intervals DataFrame
        property_names: str or list-like[str]
            Names of samples_data calculated by func.
    Returns:
        packed_results: pandas.DataFrame
            samples_data from overlap func in tabular form
    """
    results = common.PeekableIter(results)
    first_result = results.peek_next()
    results_are_mapping = isinstance(first_result, Mapping)
    if property_names is None or not property_names:
        multi_column_output = hasattr(first_result, '__len__') and len(first_result) > 1
        if results_are_mapping:
            property_names = list(first_result.keys())
        elif isinstance(first_result, pandas.Series):
            property_names = first_result.index
        elif hasattr(first_result, '_fields'):
            # noinspection PyProtectedMember
            property_names = first_result._fields
        elif hasattr(first_result, 'dtype') and hasattr(first_result.dtype, 'names'):
            property_names = first_result.dtype.names
        else:
            if multi_column_output:
                property_names = ['overlap_property_%d' % n for n in range(len(first_result))]
            else:
                property_names = ['overlap_property']
    else:
        if isinstance(property_names, str):
            property_names = [property_names]
        multi_column_output = len(property_names) > 1

    packed_results = pandas.DataFrame(
        results if results_are_mapping or not multi_column_output else (tuple(r) for r in results),
        index=index
    )
    if len(property_names) == len(packed_results.columns):
        # should basically always be the case, but could be the user didn't supply property names and the result is
        # returning inconsistently-named properties
        packed_results.columns = property_names

    return packed_results
