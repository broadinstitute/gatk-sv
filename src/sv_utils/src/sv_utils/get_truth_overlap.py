#!/usr/bin/env python

import sys
import os
import argparse
import json
import warnings
import dill
import numpy
import pandas
import scipy
import scipy.stats
from sv_utils import common, genomics_io, pedigree_tools, interval_overlaps
from typing import List, Text, Optional, Dict, Tuple, TypeVar, Union, Collection, Callable, Iterable, Set, Mapping
from types import MappingProxyType

EncodedVcfField = Union[int, str]
PandasObj = TypeVar("PandasObj", pandas.Series, pandas.DataFrame)
TestTruthOverlapStats = Dict[str, pandas.DataFrame]  # map from sample ID to table of per-variant overlap statistics
SvTypeCutoffSelector = Callable[[pandas.DataFrame], pandas.Series]
Array = Union[numpy.ndarray, pandas.array]


class Keys:
    id = genomics_io.Keys.id
    contig = genomics_io.Keys.contig
    begin = genomics_io.Keys.begin
    end = genomics_io.Keys.end
    gt = genomics_io.Keys.gt
    cn = genomics_io.Keys.cn
    rd_cn = genomics_io.Keys.rd_cn
    svtype = genomics_io.Keys.svtype
    svlen = genomics_io.Keys.svlen
    allele_frequency = genomics_io.Keys.allele_frequency
    bnd_contig_2 = genomics_io.Keys.bnd_contig_2
    bnd_end_2 = genomics_io.Keys.bnd_end_2
    gq = genomics_io.Keys.gq
    vapor_read_scores = genomics_io.Keys.vapor_read_scores
    vapor_p_non_ref = "p_non_ref"
    father = "father"
    mother = "mother"
    child = "child"
    reciprocal_overlap = "reciprocal_overlap"
    overlap_support = "overlap_support"
    inverse_distance = "inverse_distance"
    sample_ids = genomics_io.Keys.sample_ids
    is_called = "is_called"
    all_overlappers = "all"
    sample_id = genomics_io.Keys.sample_id
    property = genomics_io.Keys.property
    bnd_types = "BND types"


class Default:
    wanted_properties = interval_overlaps.Default.wanted_properties + (
        Keys.allele_frequency, Keys.gt, Keys.cn, Keys.rd_cn
    )
    f_beta = 1.0
    min_overlap_cutoff_precision = 0.99
    min_vapor_precision = 0.999
    inheritance_af_rareness = 0.05
    # when checking overlap, pseudo-size of variants with 0 / >0 reference length = scale_factor * SVLEN + expand_bp:
    expand_point_svs_bp = interval_overlaps.Default.expand_point_svs_bp
    point_sv_scale_factor = interval_overlaps.Default.point_sv_scale_factor
    expand_non_point_svs_bp = interval_overlaps.Default.expand_non_point_svs_bp
    non_point_sv_scale_factor = interval_overlaps.Default.non_point_sv_scale_factor
    dev_module_03_dir = os.path.dirname(os.path.realpath(__file__))
    project_dir = os.path.realpath(os.path.join(dev_module_03_dir, "..", ".."))
    check_some_overlap_key = f"max_{Keys.all_overlappers}_{Keys.overlap_support}"
    sv_selector_size_ranges = MappingProxyType(
        {"off-size": (-numpy.inf, 50), "small": (50, 500), "large": (500, numpy.inf)}
    )
    breakend_types = interval_overlaps.Default.breakend_types
    num_threads = common.num_physical_cpus  # by default, use all available cores for interval overlap tasks
    use_copy_number = genomics_io.Default.use_copy_number
    p_misaligned_long_read = 0.1
    vapor_good_alt_reads_threshold = 2
    vapor_bad_alt_reads_threshold = 0
    vapor_bad_cov_threshold = 5


class SvTypeCutoffInfo:
    """Class for holding data for optimal overlap cutoffs"""
    __slots__ = ("selector",
                 "stat_good", "num_good", "cutoff_good", "f_score_good", "precision_good",
                 "stat_bad", "num_bad", "cutoff_bad", "f_score_bad", "precision_bad")
    selector: SvTypeCutoffSelector  # function for selecting variants appropriate for these cutoffs
    stat_good: str                  # column prop_name of for statistic used for selecting "good" variants
    num_good: int                   # number of good variants used to derive cutoff
    cutoff_good: float              # optimal threshold for this variant for selecting "good" variants
    f_score_good: float             # f-score at "good" threshold
    precision_good: float           # precision at "good" threshold
    stat_bad: str                   # column prop_name of for statistic used for selecting "bad" variants
    num_bad: int                    # number of bad variants used to derive cutoff
    cutoff_bad: float               # optimal threshold for this variant for selecting "bad" variants
    f_score_bad: float              # f-score at "bad" threshold
    precision_bad: float            # precision at "bad" threshold

    def __init__(self, selector: SvTypeCutoffSelector,
                 stat_good: str, num_good: int, cutoff_good: float, f_score_good: float, precision_good: float,
                 stat_bad: str, num_bad: int, cutoff_bad: float, f_score_bad: float, precision_bad: float):
        # the explicit casting below *shouldn't* be necessary, but I've had weird errors with somehow getting numpy
        # types instead of python built-ins
        self.selector = selector
        self.stat_good = stat_good
        self.num_good = int(num_good)
        self.cutoff_good = float(cutoff_good)
        self.f_score_good = float(f_score_good)
        self.precision_good = float(precision_good)
        self.stat_bad = stat_bad
        self.num_bad = int(num_bad)
        self.cutoff_bad = float(cutoff_bad)
        self.f_score_bad = float(f_score_bad)
        self.precision_bad = float(precision_bad)

    def __str__(self):
        return json.dumps(
            {slot: getattr(self, slot) for slot in self.__slots__ if slot != "selector"},
            indent="  "
        )


CutoffInfo = Dict[str, SvTypeCutoffInfo]


class SampleConfidentVariants:
    __slots__ = ("good_variant_ids", "bad_variant_ids")
    good_variant_ids: Tuple[str, ...]
    bad_variant_ids: Tuple[str, ...]

    def __init__(self, good_variant_ids: Iterable[str], bad_variant_ids: Iterable[str]):
        self.good_variant_ids = tuple(good_variant_ids)
        self.bad_variant_ids = tuple(bad_variant_ids)

    @property
    def __dict__(self) -> Dict[str, Tuple[str, ...]]:
        return {"good_variant_ids": self.good_variant_ids, "bad_variant_ids": self.bad_variant_ids}

    @property
    def all_confident_variant_ids(self) -> Set[str]:
        return set(self.good_variant_ids).union(self.bad_variant_ids)

    def select(self, selector: SvTypeCutoffSelector, metrics: pandas.DataFrame) -> "SampleConfidentVariants":
        selected_ids = set(metrics.index[selector(metrics)])
        return SampleConfidentVariants(selected_ids.intersection(self.good_variant_ids),
                                       selected_ids.intersection(self.bad_variant_ids))


ConfidentVariants = Mapping[str, SampleConfidentVariants]  # mapping from sample ID to SampleConfidentVariants


class PrecisionRecallCurve:
    """ Class to form, manipulate, and query data for precision vs recall curves """
    __slots__ = ("dataframe", "num_good", "num_bad", "f_beta", "is_high_cutoff")
    dataframe: pandas.DataFrame
    num_good: int
    num_bad: int
    f_beta: float
    is_high_cutoff: bool

    f_score_key = "f"
    recall_key = "recall"
    precision_key = "precision"
    threshold_key = "threshold"

    def __init__(
            self,
            dataframe: pandas.DataFrame,
            num_good: int,
            num_bad: int,
            is_sorted: bool = False,
            f_beta: float = Default.f_beta,
            is_high_cutoff: bool = True
    ):
        self.dataframe = dataframe if is_sorted else dataframe.sort_index()
        self.num_good = num_good
        self.num_bad = num_bad
        self.f_beta = f_beta
        self.is_high_cutoff = is_high_cutoff
        if PrecisionRecallCurve.f_score_key not in self.dataframe:
            self._set_f_score()

    @staticmethod
    def from_arrays(
            thresholds: numpy.ndarray,
            precision: numpy.ndarray,
            recall: numpy.ndarray,
            num_good: int,
            num_bad: int,
            f_score: Optional[numpy.ndarray] = None,
            is_sorted: bool = False,
            f_beta: float = Default.f_beta,
            is_high_cutoff: bool = True
    ) -> "PrecisionRecallCurve":
        dataframe = pandas.DataFrame(
            {PrecisionRecallCurve.precision_key: precision, PrecisionRecallCurve.recall_key: recall},
            index=pandas.Index(thresholds, name=PrecisionRecallCurve.threshold_key)
        )
        if f_score is not None:
            dataframe[PrecisionRecallCurve.f_score_key] = f_score
        return PrecisionRecallCurve(dataframe=dataframe, num_good=num_good, num_bad=num_bad,
                                    is_sorted=is_sorted, f_beta=f_beta, is_high_cutoff=is_high_cutoff)

    @staticmethod
    def from_good_bad_thresholds(
            good_thresholds: numpy.ndarray,
            bad_thresholds: numpy.ndarray,
            f_beta: float = Default.f_beta,
            is_high_cutoff: bool = True
    ) -> "PrecisionRecallCurve":
        # map NaN thresholds to -inf/+inf (depending on high/low cutoff), because they will never be selected
        nan_replacement = -numpy.inf if is_high_cutoff else numpy.inf
        thresholds = numpy.nan_to_num(numpy.concatenate((good_thresholds, bad_thresholds)), nan=nan_replacement)
        # construct precision-recall curve (vs threshold)
        n_good = len(good_thresholds)
        n_bad = len(bad_thresholds)
        if n_good == 0 or n_bad == 0:
            # noinspection PyTypeChecker
            return PrecisionRecallCurve.empty_curve

        sort_ind = numpy.argsort(thresholds)[::-1] if is_high_cutoff else numpy.argsort(thresholds)
        thresholds = thresholds.take(sort_ind)
        recall = numpy.concatenate(
            (numpy.full(n_good, 1.0 / n_good), numpy.zeros(n_bad))
        ).take(sort_ind).cumsum()
        precision = numpy.concatenate(
            (numpy.ones(n_good), numpy.zeros(n_bad))
        ).take(sort_ind).cumsum() \
            / numpy.arange(1, len(sort_ind) + 1)

        if is_high_cutoff:
            # de-duplicate thresholds, taking the *LAST* match to a threshold, since it's all-or-nothing if multiple
            # variants match
            __, de_dup_index = numpy.unique(thresholds[::-1], return_index=True)

            def _de_dup(_arr):
                return _arr[::-1].take(de_dup_index)[::-1]
        else:
            # de-duplicate thresholds, taking the *FIRST* match to a threshold, since it's all-or-nothing if multiple
            # variants match
            __, de_dup_index = numpy.unique(thresholds, return_index=True)

            def _de_dup(_arr):
                return _arr.take(de_dup_index)

        thresholds, recall, precision = _de_dup(thresholds), _de_dup(recall), _de_dup(precision)
        # return as a DataFrame for organization purposes
        return PrecisionRecallCurve.from_arrays(thresholds=thresholds, precision=precision, recall=recall,
                                                num_good=n_good, num_bad=n_bad, is_sorted=True, f_beta=f_beta,
                                                is_high_cutoff=is_high_cutoff)

    def loc(self, indexer: pandas.api.indexers.BaseIndexer) -> "PrecisionRecallCurve":
        # note: don't know that the indexer leaves dataframe sorted, but presumably the point of doing the indexer is
        # to manipulate the underlying data
        return PrecisionRecallCurve(dataframe=self.dataframe.loc[indexer], num_good=self.num_good, num_bad=self.num_bad,
                                    is_sorted=True, f_beta=self.f_beta)

    @property
    def threshold(self) -> numpy.ndarray:
        return self.dataframe.index.to_numpy()

    @property
    def precision(self) -> numpy.ndarray:
        return self.dataframe[PrecisionRecallCurve.precision_key]

    @property
    def recall(self) -> numpy.ndarray:
        return self.dataframe[PrecisionRecallCurve.recall_key]

    @property
    def f_score(self) -> numpy.ndarray:
        return self.dataframe[PrecisionRecallCurve.f_score_key]

    @staticmethod
    def get_f_score(
            precision: numpy.ndarray,
            recall: numpy.ndarray,
            f_beta: float = Default.f_beta
    ) -> numpy.ndarray:
        f_beta_2 = f_beta * f_beta
        return numpy.nan_to_num((1.0 + f_beta_2) * precision * recall / (f_beta_2 * precision + recall), nan=0.0)

    @f_score.setter
    def f_score(self, value):
        self.dataframe[PrecisionRecallCurve.f_score_key] = value

    def _set_f_score(self):
        self.f_score = self.get_f_score(precision=self.precision, recall=self.recall, f_beta=self.f_beta)

    @property
    def is_empty(self) -> bool:
        """ return True if there is no data that could be selected with a finite threshold, False otherwise """
        return self.dataframe.empty or (
            self.dataframe.shape[0] == 1 and not numpy.isfinite(self.dataframe.index[0])
        )

    def get_point_at_max_f_score(self) -> Optional[pandas.Series]:
        if self.is_empty:
            raise ValueError("No max f-score: PrecisionRecallCurve is empty")
        return self.dataframe.loc[self.dataframe[PrecisionRecallCurve.f_score_key].idxmax()]

    def get_point_at_max_precision(self) -> Optional[pandas.Series]:
        if self.is_empty:
            raise ValueError("No max precision: PrecisionRecallCurve is empty")
        return self.dataframe.loc[self.dataframe[PrecisionRecallCurve.precision_key].idxmax()]

    def get_point_at_threshold(self, threshold: float) -> pandas.Series:
        """
        Get precision, recall, f-score at point on curve corresponding to given threshold. This will be the first row
        with curve_threshold >= requested_threshold (for high cutoffs) or curve_threshold <= requested_threshold (for
        low cutoffs)
        Args:
            threshold: float
                Requested threshold on curve
        Returns:
            prc_point: pandas.Series
                Series with threshold, precision, and recall. Due to how pandas.Series work, threshold will be accessed
                as prc_point.name
        """
        if self.is_empty:
            warnings.warn("No point at threshold: PrecisionRecallCurve is empty")
            return PrecisionRecallCurve.empty_point(threshold)
        row_index = len(self.dataframe) - 1 - self.dataframe.index[::-1].searchsorted(threshold, side="left") \
            if self.is_high_cutoff else self.dataframe.index.searchsorted(threshold, side="left")

        if row_index < 0 or row_index >= len(self.dataframe):
            _log(str(self.dataframe.index))
            warnings.warn(f"No data at requested threshold ({threshold})")
            return PrecisionRecallCurve.empty_point(threshold)
        return self.dataframe.iloc[row_index].rename(threshold)

    # noinspection PyMethodParameters
    @common.classproperty
    def empty_curve(cls: type) -> "PrecisionRecallCurve":
        return cls(
            pandas.DataFrame(
                None, columns=[PrecisionRecallCurve.precision_key, PrecisionRecallCurve.recall_key,
                               PrecisionRecallCurve.f_score_key],
                index=pandas.Index([], name=PrecisionRecallCurve.threshold_key, dtype=numpy.float64)
            ),
            num_good=0, num_bad=0, is_sorted=True
        )

    @classmethod
    def empty_point(cls: type, threshold: float = numpy.nan) -> pandas.Series:
        return pandas.Series(numpy.nan, index=[PrecisionRecallCurve.precision_key, PrecisionRecallCurve.recall_key,
                                               PrecisionRecallCurve.f_score_key],
                             name=threshold)


def _log(out: str, end='\n'):
    print(out, flush=True, end=end, file=sys.stderr)


def split_vcf_dataframe(
        variants: pandas.DataFrame,
        guarantee_allele_frequency: bool = False,
        use_copy_number: bool = Default.use_copy_number
) -> (pandas.DataFrame, pandas.DataFrame):
    f"""
    Split variant table into a location/properties table and a genotypes table
    Args:
        variants: pandas.DataFrame
            Table of variants loaded from VCF
        guarantee_allele_frequency: bool (default=False)
            If true and {Keys.allele_frequency} is not in variants, add it to variant_properties by counting the
            proportion of called alleles that are non-ref.
        use_copy_number: bool (default={Default.use_copy_number})
            Where genotype is insufficient, use copy number for estimating allele frequency and carrier status
    Returns:
        variant_properties: pandas.DataFrame
            Table of variant location and properties
        variant_carrier_status: pandas.DataFrame
            Num variants x num samples table of carrier status:
               0: is not a carrier
               1: is a carrier (may have some no-calls, but at least one called non-ref allele)
               pandas.NA: carrier status is unknown (due to no-calls)
    """
    # variant properties are columns where sample = None (as opposed to sample properties)
    variant_properties = variants.xs(None, level=Keys.sample_id, axis=1)

    for prop in [Keys.contig, Keys.begin, Keys.end, Keys.svtype, Keys.svlen]:
        if prop not in variant_properties:
            raise ValueError(f"{prop} not present in variants")
        if variant_properties[prop].hasnans:
            raise ValueError(f"{prop} has missing values")

    if guarantee_allele_frequency and Keys.allele_frequency not in variant_properties:
        variant_properties[Keys.allele_frequency] = \
            genomics_io.get_or_estimate_allele_frequency(variants, use_copy_number=use_copy_number)

    return variant_properties, genomics_io.get_carrier_status(variants, use_copy_number=use_copy_number)


def load_and_split_vcf(
        vcf: str,
        wanted_properties: Optional[Collection[str]],
        restrict_samples: Optional[Iterable[str]],
        guarantee_allele_frequency: bool = False,
        use_copy_number: bool = Default.use_copy_number
) -> (pandas.DataFrame, pandas.DataFrame):
    """
    load vcf data from file and split into non-genotype and genotype columns
    Args:
        vcf: str
            path to vcf
        wanted_properties: Optional[Collection[str]]
            Collection of desired VCF properties (note: property names from vcf standard are specified in
            genomics_io.Keys, and other properties are converted to lower case, and should be specified thusly)
        restrict_samples: Iterable[str] or None
            if not None, only take genotype data corresponding to the requested samples
        guarantee_allele_frequency: bool (default=False)
            If true and {Keys.allele_frequency} is not in variants, add it to variant_properties by counting the
            proportion of called alleles that are non-ref.
        use_copy_number: bool (default={Default.use_copy_number})
            Where genotype is insufficient, use copy number for estimating allele frequency and carrier status
    Returns:
        variant_properties: pandas.DataFrame
            Table of variant location and properties
        variant_carrier_status: pandas.DataFrame
            Num variants x num samples table of carrier status:
               0: is not a carrier
               1: is a carrier (may have some no-calls, but at least one called non-ref allele)
               pandas.NA: carrier status is unknown (due to no-calls)    """
    try:
        _log(f"load {vcf} ...", end='')
        # get all the properties we can, but don't error out HERE if some aren't present (most likely allele frequency)
        return split_vcf_dataframe(
            genomics_io.vcf_to_pandas(vcf, samples=restrict_samples, wanted_properties=wanted_properties,
                                      missing_properties_action=genomics_io.ErrorAction.Ignore,
                                      missing_samples_action=genomics_io.ErrorAction.Ignore),
            guarantee_allele_frequency=guarantee_allele_frequency, use_copy_number=use_copy_number
        )
    except Exception as err:
        common.add_exception_context(err, f"Error loading {vcf}")
        raise


def _indices_where_carrier(carrier_status: pandas.Series) -> pandas.Index:
    """ return indices (generally variant IDs) where the sample is a carrier """
    # noinspection PyTypeChecker
    return carrier_status[carrier_status].index


def _is_called_where_not_ref(carrier_status: pandas.Series) -> pandas.Series:
    """ return Series containing only non-REF GTs, with bool values True if the GT is called, and False if no-call """
    # noinspection PyTypeChecker
    return (~carrier_status[carrier_status | carrier_status.isna()].isna()).rename(Keys.is_called)


def quantify_overlap(
        interval: interval_overlaps.Record,
        overlappers: pandas.DataFrame,
) -> Dict[str, float]:
    f"""
    Compute overlap statistics between test interval and supplied overlapping intervals. NOTE that this function is not
    used directly (at least in the main get_truth_overlap algorithm), but is called by quantify_overlap_by_svtype() with
    different subsets of overlapping intervals.
    Parameters
    ----------
    interval: interval_overlaps.Record
        Test interval to check for overlaps
    overlappers: pandas.DataFrame
        Table of overlapping intervals with columns {Keys.begin} and {Keys.end}. Overlap is presumed, not checked.
    Returns
    -------
    overlap_statistics: Dict[str, float]
        Dict with three statistics:
            -reciprocal_overlap: maximum reciprocal overlap between interval and and overlapper
            -overlap_support: proportion of interval that is covered by any overlapper
            -inverse_distance: 1 / minimum distance from interval begin/end to overlapper begin/end
    """
    if len(overlappers) == 0:
        return {Keys.reciprocal_overlap: 0.0, Keys.overlap_support: 0.0, Keys.inverse_distance: 0.0}

    o_end = overlappers[Keys.end].values
    o_begin = overlappers[Keys.begin].values

    max_o_end = numpy.maximum.accumulate(o_end)
    break_indices = numpy.flatnonzero(o_begin[1:] > max_o_end[:-1])
    if len(break_indices) == 0:
        support = numpy.minimum(max_o_end[-1], interval.end) - numpy.maximum(o_begin[0], interval.begin)
    else:
        ends = max_o_end.take(numpy.concatenate((break_indices, [-1])))
        begins = o_begin.take(numpy.concatenate(([0], break_indices + 1)))
        support = (numpy.minimum(ends, interval.end) - numpy.maximum(begins, interval.begin)).sum()
    support_proportion = support / (interval.end - interval.begin)

    reciprocal_overlap = (
        (numpy.minimum(o_end, interval.end) - numpy.maximum(o_begin, interval.begin)) /
        (numpy.maximum(o_end - o_begin, interval.end - interval.begin))
    ).max()

    distance = numpy.min(
        numpy.abs(
            [o_begin - interval.begin, o_begin - interval.end, o_end - interval.begin, o_end - interval.end]
        )
    )

    return {
        Keys.reciprocal_overlap: reciprocal_overlap.max(),
        Keys.overlap_support: support_proportion,
        Keys.inverse_distance: 1.0 / (1.0 + distance)
    }


def quantify_overlap_by_svtype(
        interval: interval_overlaps.Record,
        overlappers: pandas.DataFrame,
) -> Dict[str, float]:
    """
    Compute overlap statistics between test interval and supplied overlapping intervals.
    Consider multiple subsets of overlappers (these restrictions are combinitorial):
        -restricted to only called genotypes, and unrestricted by genotype
        -restricted to on specific SVTYPE, and unrestricted
    For each subset, call quantify_overlap() to get reciprocal_overlap, overlap_support, and inverse_distance

    Parameters
    ----------
    interval: interval_overlaps.Record
        Test interval to check for overlaps
    overlappers: pandas.DataFrame
        Table of overlapping intervals with columns {Keys.begin} and {Keys.end}. Overlap is presumed, not checked.
    Returns
    -------
    overlap_statistics: Dict[str, float]
        Dict with multiple statistics (three for each subset of overlappers):
            -reciprocal_overlap: maximum reciprocal overlap between interval and and overlapper
            -overlap_support: proportion of interval that is covered by any overlapper
            -inverse_distance: 1 / minimum distance from interval begin/end to overlapper begin/end
        Keys are of form {subset_name}_{stat_name}
    """
    # measure reciprocal overlap and overlap support
    # keep "called" and "all" versions of all stats (based on whether we allow overlaps with no-call variants)
    # if there are overlappers, break down stats by svtype
    called_overlappers = overlappers.loc[overlappers[Keys.is_called]]

    return {
        **{
            f"{prefix}_{k}": v
            for check_overlappers, prefix in zip((called_overlappers, overlappers),
                                                 (Keys.is_called, Keys.all_overlappers))
            for k, v in quantify_overlap(interval, check_overlappers).items()
        },
        **{
            f"{prefix}_{svtype}_{k}": v
            for check_overlappers, prefix in zip((called_overlappers, overlappers),
                                                 (Keys.is_called, Keys.all_overlappers))
            for svtype in check_overlappers[Keys.svtype].unique()
            for k, v in quantify_overlap(interval,
                                         check_overlappers.loc[check_overlappers[Keys.svtype] == svtype]).items()
        }
    }


def get_test_truth_overlap_stats(
        test_variant_locations: pandas.DataFrame,
        truth_variant_locations: pandas.DataFrame,
        expand_point_svs_bp: int = Default.expand_point_svs_bp,
        point_sv_scale_factor: float = Default.point_sv_scale_factor,
        expand_non_point_svs_bp: int = Default.expand_non_point_svs_bp,
        non_point_sv_scale_factor: float = Default.non_point_sv_scale_factor,
        breakend_types: Collection[str] = Default.breakend_types,
        overlap_func: interval_overlaps.OverlapFunc = quantify_overlap_by_svtype,
        num_threads: int = Default.num_threads
) -> Optional[pandas.DataFrame]:
    f"""
    Get detailed statistics of overlaps between variants in test data set and truth data set. This routine uses
    all test and truth locations, so generally get_test_truth_non_ref_variant_overlap_stats() is invoked, which iterates
    over samples and filters variants by carrier status for each sample.
    Args:
        test_variant_locations: pandas.DataFrame
            Table with genomic location data (contig, begin, end, svlen, svtype) and optionally {Keys.allele_frequency}
        truth_variant_locations: pandas.DataFrame
            Table with genomic location data (contig, begin, end, svlen, svtype)
        expand_point_svs_bp: int (default = {Default.expand_point_svs_bp})
            How much to expand SVs with size 0 on the reference (in each direction), in BP
        point_sv_scale_factor: float (default = {Default.point_sv_scale_factor})
            How much to expand SVs with size 0 on the reference (in each direction), as a proportion of svlen,
        expand_non_point_svs_bp: int (default = {Default.expand_non_point_svs_bp})
            How much to expand SVs with size > 0 on the reference (in each direction), in BP
        non_point_sv_scale_factor: float (default = {Default.non_point_sv_scale_factor})
            How much to expand SVs with size > 0 on the reference (in each direction), as a proportion of svlen
        breakend_types: Collection[str] (default={Default.breakend_types})
            SVTYPEs that should be interpreted as breakends
        overlap_func: interval_overlaps.OverlapFunc (default={quantify_overlap_by_svtype})
            Function that computes overlap statistics, given a test variant interval and table of overlapping truth
            variant intervals
        num_threads: int (default = {Default.num_threads})
            Number of threads to use for interval overlap calculations
    Returns
        overlap_stats: pandas.DataFrame
            Table with variants as rows, and overlap properties as columns
    """
    if len(truth_variant_locations) == 0 or len(test_variant_locations) == 0:
        return None  # no overlaps
    truth_variant_locations = interval_overlaps.fix_variants(
        truth_variant_locations, expand_point_svs_bp=expand_point_svs_bp, point_sv_scale_factor=point_sv_scale_factor,
        expand_non_point_svs_bp=expand_non_point_svs_bp, non_point_sv_scale_factor=non_point_sv_scale_factor,
        breakend_types=breakend_types, needed_columns={Keys.is_called}
    )
    original_index = test_variant_locations.index.copy()
    test_variant_locations = interval_overlaps.fix_variants(
        test_variant_locations, expand_point_svs_bp=expand_point_svs_bp, point_sv_scale_factor=point_sv_scale_factor,
        expand_non_point_svs_bp=expand_non_point_svs_bp, non_point_sv_scale_factor=non_point_sv_scale_factor,
        breakend_types=breakend_types
    )

    # quantify to what degree each simple interval from a variant is a "hit" or a "miss"
    simple_test_overlaps = interval_overlaps.apply_interval_overlap_func(
        overlap_func, test_variant_locations, truth_variant_locations,
        description="getting test variant overlaps",
        evidence_is_unsorted=True, n_jobs=num_threads
    )

    if not isinstance(simple_test_overlaps, pandas.DataFrame):
        raise ValueError(f"simple_test_overlaps of field_type {type(simple_test_overlaps)}")

    # apply mathematical reduction functions to map overlaps on the simple intervals to the original variants
    return interval_overlaps.postprocess_bnd_complex(
        simple_test_overlaps, test_variant_locations, original_index
    )


def get_test_truth_non_ref_variant_overlap_stats(
        test_variant_locations: pandas.DataFrame,
        test_carrier_status: pandas.DataFrame,
        truth_variant_locations: pandas.DataFrame,
        truth_carrier_status: pandas.DataFrame,
        expand_point_svs_bp: int = Default.expand_point_svs_bp,
        point_sv_scale_factor: float = Default.point_sv_scale_factor,
        expand_non_point_svs_bp: int = Default.expand_non_point_svs_bp,
        non_point_sv_scale_factor: float = Default.non_point_sv_scale_factor,
        breakend_types: Collection[str] = Default.breakend_types,
        overlap_func: interval_overlaps.OverlapFunc = quantify_overlap_by_svtype,
        num_threads: int = Default.num_threads
) -> Dict[str, pandas.DataFrame]:
    f"""
    Get detailed statistics of overlaps between variants in test data set and truth data set for all common samples,
    (where each sample is not called as REF)
    Args:
        test_variant_locations: pandas.DataFrame
            Table with genomic location data (contig, begin, end, svlen, svtype) and optionally {Keys.allele_frequency}
        test_carrier_status: pandas.DataFrame
            Num variants x num samples table of carrier status for "test" data set:
               0: is not a carrier
               1: is a carrier (may have some no-calls, but at least one called non-ref allele)
               pandas.NA: carrier status is unknown (due to no-calls)
        truth_variant_locations: pandas.DataFrame
            Table with genomic location data (contig, begin, end, svlen, svtype)
        truth_carrier_status: pandas.DataFrame
            Num variants x num samples table of carrier status for "truth" data set:
               0: is not a carrier
               1: is a carrier (may have some no-calls, but at least one called non-ref allele)
               pandas.NA: carrier status is unknown (due to no-calls)
        expand_point_svs_bp: int (default = {Default.expand_point_svs_bp})
            How much to expand SVs with size 0 on the reference (in each direction), in BP
        point_sv_scale_factor: float (default = {Default.point_sv_scale_factor})
            How much to expand SVs with size 0 on the reference (in each direction), as a proportion of svlen,
        expand_non_point_svs_bp: int (default = {Default.expand_non_point_svs_bp})
            How much to expand SVs with size > 0 on the reference (in each direction), in BP
        non_point_sv_scale_factor: float (default = {Default.non_point_sv_scale_factor})
            How much to expand SVs with size > 0 on the reference (in each direction), as a proportion of svlen
        breakend_types: Collection[str] (default={Default.breakend_types})
            SVTYPEs that should be interpreted as breakends
        overlap_func: interval_overlaps.OverlapFunc (default={quantify_overlap_by_svtype})
            Function that computes overlap statistics, given a test variant interval and table of overlapping truth
            variant intervals
        num_threads: int (default = {Default.num_threads})
            Number of threads to use for interval overlap calculations
    Returns:
        overlap_stats: Dict[str, pandas.DataFrame]
            Mapping from sample_id to sample_overlap_stast, a table with called variants as rows, and overlap properties
            as columns
    """
    if Keys.allele_frequency not in test_variant_locations.columns:
        # If allele_frequency is not present, estimate it from cohort carrier status. It's better if this can be done
        # at least before converting carrier status to genotypes
        test_variant_locations = test_variant_locations.assign(
            **{Keys.allele_frequency: 1.0 - (1.0 - test_carrier_status.mean(axis=1, skipna=True)) ** 0.5}
        )
    overlap_stats = {}

    for sample_id in truth_carrier_status.columns:
        if sample_id not in test_carrier_status.columns:
            continue  # no overlap possible
        # get indices of variants that are present for this sample
        test_carrier_index = _indices_where_carrier(test_carrier_status[sample_id])

        # find variant calls in test set that have high degree of overlap with good calls in the truth set,
        # or the overlap is unknown (i.e. the truth set is no-call)
        _log(f"Finding overlap stats for {sample_id}")
        truth_is_called_where_not_ref = _is_called_where_not_ref(truth_carrier_status[sample_id])
        sample_overlap_stats = get_test_truth_overlap_stats(
            truth_variant_locations=pandas.concat(
                (truth_variant_locations.loc[truth_is_called_where_not_ref.index], truth_is_called_where_not_ref),
                axis=1
            ),
            test_variant_locations=test_variant_locations.loc[test_carrier_index],
            expand_point_svs_bp=expand_point_svs_bp, point_sv_scale_factor=point_sv_scale_factor,
            expand_non_point_svs_bp=expand_non_point_svs_bp, non_point_sv_scale_factor=non_point_sv_scale_factor,
            breakend_types=breakend_types, overlap_func=overlap_func, num_threads=num_threads
        )
        if sample_overlap_stats is None:
            _log(f"No overlaps for {sample_id}")
            continue
        overlap_stats[sample_id] = sample_overlap_stats
        overlap_stats[sample_id][Keys.svtype] = test_variant_locations.loc[test_carrier_index, Keys.svtype]
        overlap_stats[sample_id][Keys.svlen] = test_variant_locations.loc[test_carrier_index, Keys.svlen]
        overlap_stats[sample_id][Keys.allele_frequency] = \
            test_variant_locations.loc[test_carrier_index, Keys.allele_frequency]

        _log(f"Done with {sample_id}")

    return overlap_stats


def _get_common_sample_ids(test_vcfs: Union[str, Iterable[str]],
                           truth_vcfs: Union[str, Iterable[str]]) -> Optional[Set[str]]:
    """
    Find samples that are present in both the test and truth vcfs
    Args:
        test_vcfs: str or Iterable[str]
            Path to test vcf(s)
        truth_vcfs: str or Iterable[str]
            Paths to truth vcf(s)
    Returns:
        common_samples: Set[str]
            Sample IDs that are present in both test and truth
    """
    if isinstance(test_vcfs, str):
        test_vcfs = (test_vcfs,)
    test_samples = set()
    for test_vcf in test_vcfs:
        test_samples.update(genomics_io.get_vcf_sample_ids(test_vcf))
    truth_samples = set()
    for truth_vcf in truth_vcfs:
        truth_samples.update(genomics_io.get_vcf_sample_ids(truth_vcf))
    return test_samples.intersection(truth_samples)


def get_test_truth_overlap_stats_from_vcfs(
        test_vcfs: Union[str, Iterable[str]],
        truth_vcfs: Union[str, Iterable[str]],
        wanted_properties: Optional[Collection[str]] = Default.wanted_properties,
        expand_point_svs_bp: int = Default.expand_point_svs_bp,
        point_sv_scale_factor: float = Default.point_sv_scale_factor,
        expand_non_point_svs_bp: int = Default.expand_non_point_svs_bp,
        non_point_sv_scale_factor: float = Default.non_point_sv_scale_factor,
        breakend_types: Collection[str] = Default.breakend_types,
        overlap_func: interval_overlaps.OverlapFunc = quantify_overlap_by_svtype,
        use_copy_number: bool = Default.use_copy_number,
        num_threads: int = Default.num_threads
) -> Dict[str, pandas.DataFrame]:
    """
    Get detailed statistics of overlaps between variants in test data set and truth data set for all common samples
    Args:
        test_vcfs: Union[str, Iterable[str]]
            Full path to test VCF (or multiple test VCFs)
        truth_vcfs: Union[str, Iterable[str]]
            Full path to VCFs containing "truth" data (or multiple VCFs)
        wanted_properties: Optional[Collection[str] (default=None)
            If None, get all available properties in the VCF (determined by call to _get_all_properties()).
            Otherwise get the requested properties (which should be lower-case, and standard properties are elements of
            Keys). If extra property names are requested that are not present in the VCF, behavior is controlled by the
            value of missing_properties_action
        expand_point_svs_bp: int (default = {Default.expand_point_svs_bp})
            How much to expand SVs with size 0 on the reference (in each direction), in BP
        point_sv_scale_factor: float (default = {Default.point_sv_scale_factor})
            How much to expand SVs with size 0 on the reference (in each direction), as a proportion of svlen,
        expand_non_point_svs_bp: int (default = {Default.expand_non_point_svs_bp})
            How much to expand SVs with size > 0 on the reference (in each direction), in BP
        non_point_sv_scale_factor: float (default = {Default.non_point_sv_scale_factor})
            How much to expand SVs with size > 0 on the reference (in each direction), as a proportion of svlen
        breakend_types: Collection[str] (default={Default.breakend_types})
            SVTYPEs that should be interpreted as breakends
        overlap_func: interval_overlaps.OverlapFunc (default={quantify_overlap_by_svtype})
            Function that computes overlap statistics, given a test variant interval and table of overlapping truth
            variant intervals
        use_copy_number: bool (default={Default.use_copy_number})
            Where genotype is insufficient, use copy number for estimating allele frequency and carrier status
        num_threads: int (default = {Default.num_threads})
            Number of threads to use for interval overlap calculations
    Returns:
        overlap_stats: Dict[str, pandas.DataFrame]
            Mapping from sample_id to overlap_stats, a table with called variants as rows, and overlap properties as
            columns
    """
    overlap_stats = {}
    if isinstance(test_vcfs, str):
        test_vcfs = (test_vcfs,)
    for test_vcf in test_vcfs:
        restrict_samples = _get_common_sample_ids(test_vcf, truth_vcfs)
        if not restrict_samples:
            warnings.warn(f"Skipping test VCF {test_vcf} because it contains no samples in truth set")
            continue

        test_variant_locations, test_carrier_status = load_and_split_vcf(
            test_vcf, wanted_properties=wanted_properties, restrict_samples=restrict_samples,
            guarantee_allele_frequency=True, use_copy_number=use_copy_number
        )

        if isinstance(truth_vcfs, str):
            truth_vcfs = (truth_vcfs,)
        for truth_vcf in truth_vcfs:
            if not restrict_samples.intersection(genomics_io.get_vcf_sample_ids(truth_vcf)):
                warnings.warn(f"Skipping truth VCF {truth_vcf} because it contains no samples in test set")
                continue
            truth_variant_locations, truth_carrier_status = load_and_split_vcf(
                truth_vcf, wanted_properties=wanted_properties, restrict_samples=restrict_samples,
                guarantee_allele_frequency=False, use_copy_number=use_copy_number
            )
            try:
                overlap_stats.update(
                    get_test_truth_non_ref_variant_overlap_stats(
                        test_variant_locations, test_carrier_status, truth_variant_locations, truth_carrier_status,
                        expand_point_svs_bp=expand_point_svs_bp, point_sv_scale_factor=point_sv_scale_factor,
                        expand_non_point_svs_bp=expand_non_point_svs_bp,
                        non_point_sv_scale_factor=non_point_sv_scale_factor, breakend_types=breakend_types,
                        overlap_func=overlap_func, num_threads=num_threads
                    )
                )
            except Exception as err:
                common.add_exception_context(err, f"Processing {truth_vcf}")
                raise
    return overlap_stats


def _choose_cutoff(
        decision_curve: PrecisionRecallCurve,
        default_cutoff: float,
        min_overlap_cutoff_precision: float
) -> pandas.Series:
    """
    Return point on precision vs recall curve that is an optimal cutoff
    Args:
        decision_curve: PrecisionRecallCurve
            Object with data for precision vs recall curves formed by considering possible cutoffs vs the resulting
            heritability
        default_cutoff: float
            Cutoff to apply if there is insufficient data to estimate precision/recall
        min_overlap_cutoff_precision: float = {Default.min_overlap_cutoff_precision}
            Minimum allowed proportion of variants with overlap that will pass "good" or "bad" threshold that will be
            supported by inheritance or pacbio support. i.e. 1 - min_overlap_cutoff_precision are confident but wrong.
    Returns:
        optimal_decision_curve_point: pandas.Series
            Data at optimal cutoff (cutoff, precision, recall)
    """
    if decision_curve.is_empty:  # no data, return default with bad f-score & precision
        return pandas.Series({PrecisionRecallCurve.precision_key: 0, PrecisionRecallCurve.recall_key: 0,
                              PrecisionRecallCurve.f_score_key: 0}, name=default_cutoff)
    _min_precision_curve = decision_curve.loc(decision_curve.precision >= min_overlap_cutoff_precision)
    if _min_precision_curve.is_empty:  # can't meet required precision
        return decision_curve.get_point_at_max_precision()
    else:
        return _min_precision_curve.get_point_at_max_f_score()


def _get_stat_optimal_overlap_cutoffs(
        trios_overlap_info: Tuple[Mapping[str, pandas.DataFrame], ...],
        vapor_info: Mapping[str, pandas.DataFrame],
        selector: SvTypeCutoffSelector,
        stat_good: str,
        stat_bad: str,
        f_beta: float = Default.f_beta,
        min_overlap_cutoff_precision: float = Default.min_overlap_cutoff_precision,
        inheritance_af_rareness: float = Default.inheritance_af_rareness,
        check_some_overlap_key: str = Default.check_some_overlap_key
) -> SvTypeCutoffInfo:
    f"""
    Find the optimal "good" and "bad" cutoffs for a specified overlap statistic
    Args:
        trios_overlap_info: Tuple[Mapping[str, pandas.DataFrame], ...],
            Tuple with a mapping for each trio from family members (keys: {Keys.father}, {Keys.mother}, and
            {Keys.child}) to a pandas Dataframe with overlap info
        vapor_info: Mapping[str, pandas.DataFrame]
            Mapping from sample ID to table with overlap stats and the estimated probability that the variant is "good"
            (non-ref) from VaPoR
        selector: SvTypeCutoffSelector
            Function that consumes variant properties and returns boolean index that selects for desired variants (those
            that are intended to use this cutoff).
        stat_good: str
            Column prop_name of statistic to find cutoffs for good overlaps
        stat_bad: str
            Column prop_name of statistic to find cutoffs for bad overlaps
        f_beta: float = {Default.f_beta}
            Relative importance of recall, as compared to precision.
        min_overlap_cutoff_precision: float = {Default.min_overlap_cutoff_precision}
            Minimum allowed proportion of variants with overlap that will pass "good" or "bad" threshold that will be
            supported by inheritance or pacbio support. i.e. 1 - min_overlap_cutoff_precision are confident but wrong.
        inheritance_af_rareness: float = {Default.inheritance_af_rareness}
            Maximum allele frequency where inheritance will be assessed. i.e. if AF > inheritance_af_rareness, don't
            attempt to determine if a variant is mendelian, ignore it for the purposes of finding cutoffs.
        check_some_overlap_key: str = {Default.check_some_overlap_key}
            Field in overlap info to check if non-zero to ensure there is at least some overlap. Any overlap-related
            field will do, this is arbitrary and unlikely to affect results unless changed to something inappropriate
            (i.e. not related to overlap)
    Returns:
        stat_optimal_overlap_cutoffs: SvTypeCutoffInfo
            info necessary to find confidently "good"/"bad" overlaps with truth data
    """
    # restrict trios_overlap_info and vapor info to only the selected SVs
    trios_overlap_info = tuple(
        {
            family_member_key: info.loc[selector(info)]
            for family_member_key, info in trio_overlap_info.items()
        }
        for trio_overlap_info in trios_overlap_info
    )
    vapor_info = {
        sample_id: info.loc[selector(info)]
        for sample_id, info in vapor_info.items()
    }

    def _stat_or_nan(_stat: str, _overlaps: pandas.DataFrame, _v_ids: Union[pandas.Series, Set[str]]) -> Array:
        """ If the sample overlaps contain this stat, return the stat values for the requested variant IDs.
            otherwise, return appropriately-sized array of NaN """
        if _stat in _overlaps:
            if isinstance(_v_ids, pandas.Series):
                _v_ids = _v_ids.loc[_v_ids].index
            return _overlaps.loc[_v_ids, _stat].values
        else:
            if isinstance(_v_ids, set):
                _num_locs = sum(_i in _v_ids for _i in _overlaps.index)
            else:
                _num_locs = sum(_v_ids)
            return numpy.full(_num_locs, numpy.nan)

    good_cutoffs = numpy.zeros(0)
    bad_cutoffs = numpy.zeros(0)
    for trio in trios_overlap_info:
        father, mother, child = trio[Keys.father], trio[Keys.mother], trio[Keys.child]
        # good SVs:
        #  AF < "rareness" threshold AND child has it and at least one parent has it AND both have non-zero overlap
        #  -> then the SV is good and any sample with non-zero overlap has a good cutoff
        af = child[Keys.allele_frequency]
        maybe_good_child_svs = set(
            child.loc[(child[check_some_overlap_key] > 0) & (af < inheritance_af_rareness)].index
        )
        maybe_good_parent_svs = set(father.loc[father[check_some_overlap_key] > 0].index)\
            .union(mother.loc[mother[check_some_overlap_key] > 0].index)
        good_svs = maybe_good_child_svs.intersection(maybe_good_parent_svs)
        # note: still have to check parent overlap, because the OTHER parent could have non-zero overlap and make
        #       a good SV, but this one could have no overlap
        good_father_svs = father.loc[good_svs.intersection(father.index), check_some_overlap_key] > 0
        good_mother_svs = mother.loc[good_svs.intersection(mother.index), check_some_overlap_key] > 0
        good_cutoffs = numpy.concatenate((
            good_cutoffs,
            _stat_or_nan(stat_good, father, good_father_svs),
            _stat_or_nan(stat_good, mother, good_mother_svs),
            _stat_or_nan(stat_good, child, good_svs)
        ))

        # bad SVs
        # #  child has it but parent doesn't OR child doesn't have it and parent has zero overlap
        #  child has it but parent doesn't
        bad_svs = set(child.index).difference(set(father.index).union(mother.index))
        bad_cutoffs = numpy.concatenate((
            bad_cutoffs,
            _stat_or_nan(stat_bad, child, bad_svs)
        ))

    for vapor_variants in vapor_info.values():
        sv_is_good = vapor_variants[Keys.vapor_p_non_ref] > min_overlap_cutoff_precision
        sv_is_bad = vapor_variants[Keys.vapor_p_non_ref] < 1.0 - min_overlap_cutoff_precision
        good_cutoffs = numpy.concatenate((
            good_cutoffs, _stat_or_nan(stat_good, vapor_variants, sv_is_good)
        ))
        bad_cutoffs = numpy.concatenate((
            bad_cutoffs, _stat_or_nan(stat_bad, vapor_variants, sv_is_bad)
        ))

    # set nan (missing) values to 0, since "missing" only happens if a certain SV type isn't present in the overlap
    # function, corresponding to 0 overlap support and reciprocal overlap
    numpy.nan_to_num(good_cutoffs, nan=0.0, copy=False).sort()
    numpy.nan_to_num(bad_cutoffs, nan=0.0, copy=False).sort()

    # make typical PRC, balancing precision vs recall for finding good overlaps
    decision_curve_good = PrecisionRecallCurve.from_good_bad_thresholds(
        good_thresholds=good_cutoffs, bad_thresholds=bad_cutoffs, f_beta=f_beta
    )
    # make opposite PRC, trying to find *bad* overlaps
    decision_curve_bad = PrecisionRecallCurve.from_good_bad_thresholds(
        good_thresholds=bad_cutoffs, bad_thresholds=good_cutoffs, f_beta=f_beta, is_high_cutoff=False
    )

    good_cutoff_point = _choose_cutoff(decision_curve_good, default_cutoff=0.5,
                                       min_overlap_cutoff_precision=min_overlap_cutoff_precision)
    # don't always use for training cutoffs, but final cutoff will always toss things with no overlap
    bad_cutoff_point = _choose_cutoff(decision_curve_bad, default_cutoff=0.0,
                                      min_overlap_cutoff_precision=min_overlap_cutoff_precision)

    # noinspection PyTypeChecker
    return SvTypeCutoffInfo(
        selector=selector,
        stat_good=stat_good, num_good=len(good_cutoffs), cutoff_good=good_cutoff_point.name,
        f_score_good=good_cutoff_point[PrecisionRecallCurve.f_score_key],
        precision_good=good_cutoff_point[PrecisionRecallCurve.precision_key],
        stat_bad=stat_bad, num_bad=len(bad_cutoffs), cutoff_bad=bad_cutoff_point.name,
        f_score_bad=bad_cutoff_point[PrecisionRecallCurve.f_score_key],
        precision_bad=bad_cutoff_point[PrecisionRecallCurve.precision_key]
    )


def get_trios_overlap_info(
        overlap_stats: Dict[str, pandas.DataFrame],
        ped_files: Union[str, Iterable[str]]
) -> Tuple[Dict[str, pandas.DataFrame], ...]:
    """ Organize overlap info into trios, as specified by the supplied pedigree files """
    if isinstance(ped_files, str):
        ped_files = (ped_files,)

    return tuple(
        {
            Keys.father: overlap_stats[pedigree_line.father_id],
            Keys.mother: overlap_stats[pedigree_line.mother_id],
            Keys.child: overlap_stats[pedigree_line.proband_id]
        }
        for ped_file in ped_files
        for pedigree_line in pedigree_tools.PedigreeFileInfo(ped_file).pedigree_lines
        if not pedigree_line.any_unknown and all(
            sample_id in overlap_stats
            for sample_id in (pedigree_line.father_id, pedigree_line.mother_id, pedigree_line.proband_id)
        )
    )


def get_vapor_p_non_ref_old(vapor_variants: pandas.DataFrame) -> pandas.DataFrame:
    """ Given table of vapor data, return a one-column table of probabilities that each vapor variant is non-REF
     based on GQ """
    p_gt_bad = 10 ** -vapor_variants[Keys.gq]
    return pandas.DataFrame(
        numpy.where(
            vapor_variants[Keys.gt] == "0/0", p_gt_bad, 1.0 - p_gt_bad
        ),
        columns=[Keys.vapor_p_non_ref],
        index=vapor_variants.index
    )


def get_vapor_p_non_ref(
        vapor_variants: pandas.DataFrame,
        p_misaligned: float = Default.p_misaligned_long_read
) -> pandas.DataFrame:
    """ Given table of vapor data, return a one-column table of probabilities that each vapor variant is non-REF
     using a negative binomial model of supporting reads """
    vapor_read_scores = vapor_variants[Keys.vapor_read_scores].apply(
        lambda scores: [float(score) for score in scores.split(',')] if scores else []
    )
    num_reads = vapor_read_scores.apply(len)
    num_alt_reads = vapor_read_scores.apply(lambda _scores: sum(1 for _score in _scores if _score > 0))
    # model het likelihood as binomial distribution with 50/50 chance a read is ALT or REF
    likelihood_het = scipy.stats.binom.pmf(num_alt_reads, num_reads, 0.5)
    # model homref likelihood as binomial distribution with p_misaligned chance that read looks ALT
    likelihood_ref = scipy.stats.binom.pmf(num_alt_reads, num_reads, p_misaligned)
    # model homvar likelihood as binomial distribution with p_misaligned chance that read looks REF
    likelihood_homvar = scipy.stats.binom.pmf(num_reads - num_alt_reads, num_reads, p_misaligned)
    # set flat priors on true genotype, just allow likelihoods to dominate prediction
    return pandas.DataFrame(
        1.0 - likelihood_ref / (likelihood_ref + likelihood_het + likelihood_homvar),
        columns=[Keys.vapor_p_non_ref],
        index=vapor_variants.index
    )


def get_vapor_p_non_ref_threshold(
        vapor_variants: pandas.DataFrame,
        good_alt_reads_threshold: int = Default.vapor_good_alt_reads_threshold,
        bad_alt_reads_threshold: int = Default.vapor_bad_alt_reads_threshold,
        bad_cov_threshold: int = Default.vapor_bad_cov_threshold
) -> pandas.DataFrame:
    """ Given table of vapor data, return a one-column table of probabilities that each vapor variant is non-REF
     based on thresholds of the number of total and alt-supporting reads """
    vapor_read_scores = vapor_variants[Keys.vapor_read_scores].apply(
        lambda scores: [float(score) for score in scores.split(',')] if scores else []
    )
    num_alt_reads = vapor_read_scores.apply(lambda _scores: sum(1 for _score in _scores if _score > 0))
    num_reads = vapor_read_scores.apply(lambda _scores: len(_scores))
    return pandas.DataFrame(
        numpy.where(
            numpy.logical_and(num_alt_reads <= bad_alt_reads_threshold, num_reads >= bad_cov_threshold),
            0.0,
            numpy.where(num_alt_reads >= good_alt_reads_threshold, 1.0, 0.5)
        ),
        columns=[Keys.vapor_p_non_ref],
        index=vapor_variants.index
    )


def get_sv_selectors(
        all_sv_types: Set[str],
        bnd_types: Iterable[str] = interval_overlaps.Default.breakend_types,
        size_ranges: Mapping[str, Tuple[float, float]] = Default.sv_selector_size_ranges
) -> Dict[str, SvTypeCutoffSelector]:
    """
    Return Dict with key = SV category and value = SVTypeCutoffSelector (a function that selects only the variants from
    the desired category)
    """
    def _get_selector(_svtype: str, _low: float, _high: float) -> SvTypeCutoffSelector:
        return lambda info: (info[Keys.svtype] == _svtype) & (info[Keys.svlen] >= _low) & (info[Keys.svlen] < _high)

    # don't break down BNDs (and similar types) by size because it doesn't make sense to
    # don't break down INV by size because there are too few
    size_divided_types = all_sv_types.difference(bnd_types).difference({"INV"})
    return {
        **{
            Keys.bnd_types: lambda info: info[Keys.svtype].isin(bnd_types)
        },
        **{
            "INV": lambda info: info[Keys.svtype] == "INV"
        },
        **{
            f"{size_name} {svtype}": _get_selector(svtype, low, high)
            for svtype in size_divided_types for size_name, (low, high) in size_ranges.items()
        }
    }


def get_optimal_overlap_cutoffs(
        overlap_stats: Mapping[str, pandas.DataFrame],
        ped_files: Optional[Union[str, Iterable[str]]] = None,
        vapor_files: Optional[Dict[str, str]] = None,
        f_beta: float = Default.f_beta,
        min_overlap_cutoff_precision: float = Default.min_overlap_cutoff_precision,
        inheritance_af_rareness: float = Default.inheritance_af_rareness,
        bnd_types: Iterable[str] = interval_overlaps.Default.breakend_types,
        size_ranges: Mapping[str, Tuple[float, float]] = Default.sv_selector_size_ranges
) -> CutoffInfo:
    f"""
    Find optimal thresholds for declaring confident "good" or "bad" overlaps with truth data based on
        -Mendelian inheritance in trios    AND/OR
        -support in pacbio reads from VaPoR.
    Args:
        overlap_stats: Dict[str, pandas.DataFrame]
            Map from sample ID to table of overlap statistics (rows indexed by variant ID, columns by statistic)
        ped_files: Optional[Union[str, Iterable[str]]] = None
            Pedigree files with trios in overlap stats
        vapor_files: Optional[Dict[str, str]] = None
            Map from sample ID to vapor file (providing estimation of pacbio support for each short-read variant)
        f_beta: float = {Default.f_beta},
            Relative importance of recall, as compared to precision.
        min_overlap_cutoff_precision: float = {Default.min_overlap_cutoff_precision}
            Minimum allowed proportion of variants with overlap that will pass "good" or "bad" threshold that will be
            supported by inheritance or pacbio support. i.e. 1 - min_overlap_cutoff_precision are confident but wrong.
        inheritance_af_rareness: float = {Default.inheritance_af_rareness}
            Maximum allele frequency where inheritance will be assessed. i.e. if AF > inheritance_af_rareness, don't
            attempt to determine if a variant is mendelian, ignore it for the purposes of finding cutoffs.
        bnd_types: Iterable[str] {interval_overlaps.Default.breakend_types}
            SV types that are considered breakends.
        size_ranges: Mapping[str, Tuple[float, float]] {Default.sv_selector_size_ranges}
            Ranges of svlen to break SV types into
    Returns:
        optimal_overlap_cutoffs: CutoffInfo
            Mapping from SV category to SvTypeCutoffInfo (class holding info necessary to find confidently "good"/"bad"
            overlaps with truth data)
    """
    if ped_files is None:
        ped_files = ()
    elif isinstance(ped_files, str):
        ped_files = (ped_files,)
    if vapor_files is None:
        vapor_files = {}
    if not ped_files and not vapor_files:
        raise ValueError("Need one or more ped files and/or one or more vapor files to compute optimal overlaps")

    trios_overlap_info = tuple(
        {
            Keys.father: overlap_stats[pedigree_line.father_id],
            Keys.mother: overlap_stats[pedigree_line.mother_id],
            Keys.child: overlap_stats[pedigree_line.proband_id]
        }
        for pedigree_line in pedigree_tools.PedigreeFileInfo.load(pedigree_files=ped_files).pedigree_lines
        if not pedigree_line.any_unknown and all(
            sample_id in overlap_stats
            for sample_id in (pedigree_line.father_id, pedigree_line.mother_id, pedigree_line.proband_id)
        )
    )
    print(f"Found {len(trios_overlap_info)} trios")

    # augment overlap stats with probability the variant is non-ref, as estimated by VaPoR:
    vapor_info = {
        sample_id:
            get_vapor_p_non_ref(genomics_io.vapor_to_pandas(vapor_file))
            .join(overlap_stats[sample_id], how="inner")
        for sample_id, vapor_file in vapor_files.items()
        if sample_id in overlap_stats
    }
    print(f"Got {len(vapor_info)} vapor samples")
    if not trios_overlap_info and not vapor_info:
        raise ValueError("No trios from ped files or vapor files were present in the overlap info")

    all_sv_types = set(
        svtype
        for sample_overlap_stats in overlap_stats.values()
        for svtype in sample_overlap_stats[Keys.svtype].unique()
    )
    sv_selectors = get_sv_selectors(all_sv_types=all_sv_types, bnd_types=bnd_types, size_ranges=size_ranges)

    def _get_optimal_cutoff(_cutoffs: List[SvTypeCutoffInfo], _good: bool) -> SvTypeCutoffInfo:
        def _precision(_cutoff: SvTypeCutoffInfo) -> float:
            return _cutoff.precision_good if _good else _cutoff.precision_bad

        def _fscore(_cutoff: SvTypeCutoffInfo) -> float:
            return _cutoff.f_score_good if _good else _cutoff.f_score_bad

        def _cutoff_sort_key(_cutoff: SvTypeCutoffInfo) -> float:
            # choose in priority:
            #  a) if anything meets the precision goal, whatever has the highest f-score
            #  b) otherwise, whatever has the highest precision
            return min_overlap_cutoff_precision + (1.0 - min_overlap_cutoff_precision) * _fscore(_cutoff) \
                if _precision(_cutoff) >= min_overlap_cutoff_precision \
                else _precision(_cutoff)

        _cutoffs.sort(key=_cutoff_sort_key)
        return _cutoffs[-1]

    def _remove_called_status(_stat: str) -> str:
        return '_'.join(_stat.split(f"_{Keys.is_called}_", 1)) if f"_{Keys.is_called}_" in _stat else \
            '_'.join(_stat.split(f"_{Keys.all_overlappers}_", 1)) if f"_{Keys.all_overlappers}_" in _stat \
            else _stat

    # get all the stats that could be used for overlap thresholds
    all_stats = set.union(
        *(set(sample_overlap_stats.columns) for sample_overlap_stats in overlap_stats.values())
    ).difference({Keys.svtype, Keys.svlen, Keys.allele_frequency})
    # we want to use called variants for good overlaps (confident there is an overlapper) and any non-ref for bad
    # overlaps (confident there's no overlapper). So loop through and link them together
    all_stats_linked = {}
    for stat in all_stats:
        base_stat = _remove_called_status(stat)
        if Keys.is_called in stat:
            all_stats_linked[base_stat] = (stat, ) + all_stats_linked.get(base_stat, (None, None))[1:]
        else:
            all_stats_linked[base_stat] = all_stats_linked.get(base_stat, (None, None))[:1] + (stat,)
    overlap_cutoffs = {}
    for selector_name, selector in sv_selectors.items():
        stats_overlap_cutoff_info = [
            _get_stat_optimal_overlap_cutoffs(
                trios_overlap_info, vapor_info, selector, stat_good, stat_bad, f_beta=f_beta,
                min_overlap_cutoff_precision=min_overlap_cutoff_precision,
                inheritance_af_rareness=inheritance_af_rareness
            )
            for (stat_good, stat_bad) in all_stats_linked.values()
        ]
        n_good = next((c_i.num_good for c_i in stats_overlap_cutoff_info if c_i.num_good > 0), 0)
        n_bad = next((c_i.num_bad for c_i in stats_overlap_cutoff_info if c_i.num_bad > 0), 0)
        print(f"{selector_name}: {n_good} good, {n_bad} bad")
        for c_i in stats_overlap_cutoff_info:
            print(f"\t{c_i.stat_good}: good_c={c_i.cutoff_good:.3g}, f={c_i.f_score_good:.3g},"
                  f" prec={c_i.precision_good:.3g}; bad_c={c_i.cutoff_bad:.3g}, f={c_i.f_score_bad:.3g},"
                  f" prec={c_i.precision_bad:.3g}")

        good_cutoff_stats = _get_optimal_cutoff(stats_overlap_cutoff_info, _good=True)
        bad_cutoff_stats = _get_optimal_cutoff(stats_overlap_cutoff_info, _good=False)
        if good_cutoff_stats.f_score_good == 0 and bad_cutoff_stats.f_score_bad == 0:
            # no power, don't bother using this selector
            continue
        overlap_cutoffs[selector_name] = SvTypeCutoffInfo(
            selector=selector,
            stat_good=good_cutoff_stats.stat_good, num_good=good_cutoff_stats.num_good,
            cutoff_good=good_cutoff_stats.cutoff_good, f_score_good=good_cutoff_stats.f_score_good,
            precision_good=good_cutoff_stats.precision_good,
            stat_bad=bad_cutoff_stats.stat_bad, num_bad=bad_cutoff_stats.num_bad,
            cutoff_bad=bad_cutoff_stats.cutoff_bad, f_score_bad=bad_cutoff_stats.f_score_bad,
            precision_bad=bad_cutoff_stats.precision_bad
        )

    return overlap_cutoffs


def save_optimal_overlap_cutoffs(
        optimal_overlap_cutoffs: CutoffInfo,
        optimal_overlap_cutoffs_file: str
):
    """ Save optimal cutoffs object by pickling it (with dill module, so it's possible to pickle functions) """
    with open(optimal_overlap_cutoffs_file, "wb") as f_out:
        dill.dump(optimal_overlap_cutoffs, f_out)


def load_optimal_overlap_cutoffs(optimal_overlap_cutoffs_file: str) -> CutoffInfo:
    """ Load optimal cutoffs object by unpickling it (with dill module, so it's possible to pickle functions) """
    with open(optimal_overlap_cutoffs_file, "rb") as f_in:
        return dill.load(f_in)


def display_optimal_overlap_cutoffs(optimal_overlap_cutoffs: CutoffInfo):
    """ Display optimal cutoffs in relatively human-readable format """
    for sv_category, sv_cutoff in optimal_overlap_cutoffs.items():
        print(f"{sv_category}: {sv_cutoff}")


def apply_overlap_cutoffs(
        sample_overlap_stats: pandas.DataFrame, cutoff_info: CutoffInfo
) -> SampleConfidentVariants:
    """ Given overlap statistics and overlap cutoffs, apply the cutoffs and obtain confident good/bad variants """
    is_good = pandas.Series(numpy.zeros(len(sample_overlap_stats), dtype=bool), index=sample_overlap_stats.index)
    is_bad = pandas.Series(numpy.zeros(len(sample_overlap_stats), dtype=bool), index=sample_overlap_stats.index)
    for selector_name, selection_cutoff_info in cutoff_info.items():
        svtype_overlap_stats = sample_overlap_stats.loc[selection_cutoff_info.selector(sample_overlap_stats)]
        is_good[svtype_overlap_stats.index] = \
            (svtype_overlap_stats[selection_cutoff_info.stat_good]
             if selection_cutoff_info.stat_good in svtype_overlap_stats else 0.0) >= selection_cutoff_info.cutoff_good
        is_bad[svtype_overlap_stats.index] = \
            (svtype_overlap_stats[selection_cutoff_info.stat_bad]
             if selection_cutoff_info.stat_bad in svtype_overlap_stats else 0.0) <= selection_cutoff_info.cutoff_bad

    return SampleConfidentVariants(
        good_variant_ids=sorted(is_good[is_good].index),
        bad_variant_ids=sorted(is_bad[is_bad].index)
    )


def select_confident_vapor_variants(
        vapor_file: str,
        valid_variant_ids: Set[str],
        strategy: str,
        read_strategy_good_support_threshold: int,
        read_strategy_bad_support_threshold: int,
        read_strategy_bad_cov_threshold: int,
        precision: float = Default.min_vapor_precision
) -> SampleConfidentVariants:
    f"""
    Use provided the VaPoR file to select confident good/bad variants
    Args:
        vapor_file: str
            Full path to VaPoR file
        valid_variant_ids: Set[str]
            Set of variants for which we want to estimate the probability that they are non-REF
        precision: float (default={Default.min_vapor_precision})
            Threshold for deciding that that non-REF probability is high/low enough to be a good/bad variant.
    Returns:
        sample_confident_variants: SampleConfidentVariants
            Object holding the confident variants for this sample
    """
    vapor_data = genomics_io.vapor_to_pandas(vapor_file)
    if strategy == 'GQ':
        vapor_p_non_ref = get_vapor_p_non_ref_old(vapor_data)
    elif strategy == 'GT':
        vapor_p_non_ref = get_vapor_p_non_ref(vapor_data)
    elif strategy == 'READS':
        vapor_p_non_ref = get_vapor_p_non_ref_threshold(vapor_data,
                                                        read_strategy_good_support_threshold,
                                                        read_strategy_bad_support_threshold,
                                                        read_strategy_bad_cov_threshold)
    else:
        raise ValueError("Unsupported strategy: {}".format(strategy))
    return SampleConfidentVariants(
        good_variant_ids=sorted(
            valid_variant_ids.intersection(
                vapor_p_non_ref.loc[vapor_p_non_ref[Keys.vapor_p_non_ref] > precision].index
            )
        ),
        bad_variant_ids=sorted(
            valid_variant_ids.intersection(
                vapor_p_non_ref.loc[vapor_p_non_ref[Keys.vapor_p_non_ref] < 1.0 - precision].index
            )
        )
    )


def reconcile_confident_variants(
        overlap_variants: Optional[SampleConfidentVariants],
        vapor_variants: Optional[SampleConfidentVariants]
) -> SampleConfidentVariants:
    """ Given confident variants obtained by VaPoR and overlaps, reconcile any differences by favoring VaPoR """
    if overlap_variants is None:
        return vapor_variants
    elif vapor_variants is None:
        return overlap_variants
    else:
        # return union of both methods, but when there are disagreements, prefer VaPoR
        return SampleConfidentVariants(
            good_variant_ids=sorted(
                set(overlap_variants.good_variant_ids)
                .difference(vapor_variants.bad_variant_ids)
                .union(vapor_variants.good_variant_ids)
            ),
            bad_variant_ids=sorted(
                set(overlap_variants.bad_variant_ids)
                .difference(vapor_variants.good_variant_ids)
                .union(vapor_variants.bad_variant_ids)
            )
        )


def output_confident_variants(
        overlap_info: ConfidentVariants,
        output_file: str
):
    """ Output the confident variants, either to a json file or stdout """
    if output_file == '-':
        json.dump(overlap_info, sys.stdout, indent=2, default=lambda x: x.__dict__)
        sys.stdout.flush()
    else:
        _log(f"Saving overlaps to {output_file} ...", end='')
        with open(output_file, 'w') as f_out:
            json.dump(overlap_info, f_out, default=lambda x: x.__dict__)
        _log(" okay")


def _is_file(filename: Optional[str]) -> bool:
    """ return boolean if the file path is a real file on the system """
    return filename is not None and os.path.isfile(filename)


def get_truth_overlap(
        test_vcfs: Union[str, Iterable[str]],
        truth_vcfs: Union[str, Iterable[str]],
        ped_files: Optional[Union[str, Iterable[str]]],
        optimal_overlap_cutoffs_file: str,
        vapor_files: Optional[Dict[str, str]] = None,
        expand_point_svs_bp: int = Default.expand_point_svs_bp,
        point_sv_scale_factor: float = Default.point_sv_scale_factor,
        expand_non_point_svs_bp: int = Default.expand_non_point_svs_bp,
        non_point_sv_scale_factor: float = Default.non_point_sv_scale_factor,
        f_beta: float = Default.f_beta,
        min_overlap_cutoff_precision: float = Default.min_overlap_cutoff_precision,
        min_vapor_precision: float = Default.min_vapor_precision,
        inheritance_af_rareness: float = Default.inheritance_af_rareness,
        breakend_types: Collection[str] = Default.breakend_types,
        overlap_func: interval_overlaps.OverlapFunc = quantify_overlap_by_svtype,
        use_copy_number: bool = Default.use_copy_number,
        num_threads: int = Default.num_threads
) -> ConfidentVariants:
    f"""
    Given input test VCF(s), truth VCF(s), and pedigree file(s) find subset of variants that are good/bad with high
    confidence.
    1) Find multiple statistics describing overlap between test VCF variants and truth VCF variants, for each sample
       that is present in both the test and truth data
    2) If not optimal overlap cutoffs are supplied, data-mine to find new optimal overlap cutoffs
       a) good cutoffs that admit confident good variants that are nearly always inherited
       b) bad cutoffs that admit confident bad variants that are almost never inherited
       Then save those cutoffs
    3) Apply optimal overlap cutoffs to select confident variants from the overlap statistics
    4) If VaPoR data is available, use it to select confident variants, then reconcile any differences
    6) Output the confident variants.
    Args:
        test_vcfs: Union[str, Iterable[str]]
            Full path(s) to one or more VCFs with variants of unknown quality
        truth_vcfs: Union[str, Iterable[str]]
            Full path(s) to one or more VCFs with variants presumed to be true
        ped_files: Union[str, Iterable[str]]
            Full path(s) to one or more pedigree files describing the family structure of participants in these data
        optimal_overlap_cutoffs_file: str
            Full path to file with pickled optimal overlap cutoffs. If a file is present at call time, existing cutoffs
            will be used. If not, new cutoffs will be calculated and the result will be saved to this path.
            If new cutoffs must be calculated, at least one of ped_files or vapor_files must be provided, so that there
            is a source of truth data for selecting optimal cutoffs.
        vapor_files: Optional[Dict[str, str]] = None
            Map from sample ID to vapor file (providing estimation of pacbio support for each short-read variant)
        expand_point_svs_bp: int (default = {Default.expand_point_svs_bp})
            How much to expand SVs with size 0 on the reference (in each direction), in BP
        point_sv_scale_factor: float (default = {Default.point_sv_scale_factor})
            How much to expand SVs with size 0 on the reference (in each direction), as a proportion of svlen,
        expand_non_point_svs_bp: int (default = {Default.expand_non_point_svs_bp})
            How much to expand SVs with size > 0 on the reference (in each direction), in BP
        non_point_sv_scale_factor: float (default = {Default.non_point_sv_scale_factor})
            How much to expand SVs with size > 0 on the reference (in each direction), as a proportion of svlen
        f_beta: float = {Default.f_beta},
            Relative importance of recall, as compared to precision.
        min_overlap_cutoff_precision: float (default={Default.min_overlap_cutoff_precision})
            Minimum allowed proportion of variants with overlap that will pass "good" or "bad" threshold that will be
            supported by inheritance or pacbio support. i.e. 1 - min_overlap_cutoff_precision are confident but wrong.
        min_vapor_precision: float (default={Default.min_vapor_precision})
            Minimum allowed precision for selecting good or bad variants from VaPoR
        inheritance_af_rareness: float = {Default.inheritance_af_rareness}
            Maximum allele frequency where inheritance will be assessed. i.e. if AF > inheritance_af_rareness, don't
            attempt to determine if a variant is mendelian, ignore it for the purposes of finding cutoffs.
        breakend_types: Collection[str] (default={Default.breakend_types})
            SVTYPEs that should be interpreted as breakends
        overlap_func: interval_overlaps.OverlapFunc (default={quantify_overlap_by_svtype})
            Function that computes overlap statistics, given a test variant interval and table of overlapping truth
            variant intervals
        use_copy_number: bool (default={Default.use_copy_number})
            Where genotype is insufficient, use copy number for estimating allele frequency and carrier status
        num_threads: int (default = {Default.num_threads})
            Number of threads to use for interval overlap calculations
    Returns:
        confident_variants: ConfidentVariants
            Mapping from sample ID to SampleConfidentVariants, an object that holds variant IDs that are good/bad with
            high confidence
    """
    if ped_files is None and vapor_files is None:
        if not _is_file(optimal_overlap_cutoffs_file):
            raise ValueError("Must provide pedigree and/or VaPoR file(s) to discover optimal truth overlap, or provide "
                             "a save file for previously calculated optimal overlap cutoffs.")

    overlap_stats = get_test_truth_overlap_stats_from_vcfs(
        test_vcfs=test_vcfs, truth_vcfs=truth_vcfs, expand_point_svs_bp=expand_point_svs_bp,
        point_sv_scale_factor=point_sv_scale_factor, expand_non_point_svs_bp=expand_non_point_svs_bp,
        non_point_sv_scale_factor=non_point_sv_scale_factor, breakend_types=breakend_types, overlap_func=overlap_func,
        use_copy_number=use_copy_number, num_threads=num_threads
    )
    if not overlap_stats:
        raise ValueError("There were no samples in both the test and truth sets")

    if _is_file(optimal_overlap_cutoffs_file):
        optimal_overlap_cutoffs = load_optimal_overlap_cutoffs(optimal_overlap_cutoffs_file)
    else:
        optimal_overlap_cutoffs = get_optimal_overlap_cutoffs(
            overlap_stats, ped_files=ped_files, vapor_files=vapor_files,
            min_overlap_cutoff_precision=min_overlap_cutoff_precision,
            f_beta=f_beta, inheritance_af_rareness=inheritance_af_rareness
        )
        display_optimal_overlap_cutoffs(optimal_overlap_cutoffs)
        if optimal_overlap_cutoffs_file is not None:
            save_optimal_overlap_cutoffs(optimal_overlap_cutoffs, optimal_overlap_cutoffs_file)

    # get per-sample confident variants using overlaps and optimal cutoffs
    confident_variants_overlaps = {
        sample_id: apply_overlap_cutoffs(sample_overlap_stats, optimal_overlap_cutoffs)
        for sample_id, sample_overlap_stats in overlap_stats.items()
    }
    # get per-sample confident variants using VaPoR
    valid_variant_ids = {v_id for overlap_df in overlap_stats.values() for v_id in overlap_df.index}
    confident_variants_vapor = {} if vapor_files is None else {
        sample_id: select_confident_vapor_variants(vapor_file=vapor_file, valid_variant_ids=valid_variant_ids,
                                                   precision=min_vapor_precision)
        for sample_id, vapor_file in vapor_files.items()
    }
    # combine results, reconcile any differences by preferring VaPoR
    confident_variants = {
        sample_id: reconcile_confident_variants(
            overlap_variants=confident_variants_overlaps.get(sample_id, None),
            vapor_variants=confident_variants_vapor.get(sample_id, None)
        )
        for sample_id in set(confident_variants_overlaps.keys()).union(confident_variants_vapor.keys())
    }
    return confident_variants


def get_vapor_files(vapor_json: Optional[str]) -> Optional[Dict[str, str]]:
    """ Load map from sample ID to corresponding vapor file, from input JSON """
    if vapor_json is None:
        return None
    if not os.path.isfile(vapor_json):
        raise ValueError(f"vapor json ({vapor_json}) is not a path to a valid file")
    with open(vapor_json, 'r') as f_in:
        vapor_files = json.load(f_in)
    for sample_id, vapor_file in vapor_files.items():
        if not isinstance(sample_id, str) or not isinstance(vapor_file, str) or not os.path.isfile(vapor_file):
            raise ValueError(f"vapor_json ({vapor_json}) must consist of sample_id, VaPoR file pairs. Error with pair "
                             f"{sample_id}, {vapor_file}")
    return vapor_files


def __parse_arguments(argv: List[Text]) -> argparse.Namespace:
    # noinspection PyTypeChecker
    parser = argparse.ArgumentParser(
        description="Find overlapping variants between test VCF and truth VCF",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        prog=argv[0]
    )
    parser.add_argument("--test-vcf", '-v', type=str, action="extend", nargs="+", required=True,
                        help="VCF with variants of unknown quality to be compared to truth")
    parser.add_argument("--truth-vcf", '-t', type=str, action="extend", nargs='+', required=True,
                        help="VCF with presumably true variants")
    parser.add_argument("--ped-file", "-p", type=str, action="extend", nargs='*',
                        help="Pedigree file for samples. Specify to calculate optimal overlap cutoffs, or omit to use "
                             "previously saved optimal cutoffs.")
    parser.add_argument("--vapor-json", "-j", type=str,
                        help="Json file with mapping from sample ID to corresponding VaPoR file. Specify to calculate "
                             "optimal overlap cutoffs, or omit to use previously saved optimal cutoffs.")
    parser.add_argument("--optimal-overlap-cutoffs-file", "-c", type=str,
                        help="File with pickled optimal overlap cutoffs. If an existing file is provided, it will be "
                             "used, otherwise new cutoffs will be calculated. If a path is provided that is not an "
                             "existing file, the new cutoffs will be saved to that path.")
    parser.add_argument("--output", "-O", type=str, default="-",
                        help="File to output results to. If omitted or set to '-', print to stdout")
    parser.add_argument("--expand-point-svs-bp", type=int, default=Default.expand_point_svs_bp,
                        help="When checking overlap, expand SVs with size 0 on the reference by this many BP in each"
                             " direction")
    parser.add_argument("--point-sv-scale-factor", type=float, default=Default.point_sv_scale_factor,
                        help="When checking overlap, expand SVs with size 0 on the reference by scale-factor * SVLEN in"
                             " each direction")
    parser.add_argument("--expand-non-point-svs-bp", type=int, default=Default.expand_non_point_svs_bp,
                        help="When checking overlap, expand SVs with size > 0 on the reference by this many BP in each"
                             " direction")
    parser.add_argument("--non-point-sv-scale-factor", type=float, default=Default.non_point_sv_scale_factor,
                        help="When checking overlap, expand SVs with size > 0 on the reference by scale-factor * SVLEN"
                             " in each direction")
    parser.add_argument("--min-overlap-cutoff-precision", type=float, default=Default.min_overlap_cutoff_precision,
                        help="Minimum allowed precision for selecting cutoffs for good or bad variants from overlap")
    parser.add_argument("--min-vapor-precision", type=float, default=Default.min_vapor_precision,
                        help="Minimum allowed precision for selecting good or bad variants from VaPoR")
    parser.add_argument("--f-beta", type=float, default=Default.f_beta,
                        help="beta factor for f-score, weighting importance of recall relative to precision")
    parser.add_argument("--inheritance-af-rareness", type=float, default=Default.inheritance_af_rareness,
                        help="Maximum allele frequency for a variant to use trio inheritance as a truth signal.")
    parser.add_argument("--use-copy-number", type=bool, default=Default.use_copy_number,
                        help="Where genotype is insufficient, use copy number for estimating allele frequency and "
                             "carrier status")
    parser.add_argument("--num_threads", "-@", type=int, default=Default.num_threads,
                        help="number of threads for compressing output vcf")
    parsed_arguments = parser.parse_args(argv[1:] if len(argv) > 1 else ["--help"])
    if parsed_arguments.test_vcf is None:
        raise ValueError("Must supply one or more --test-vcf")
    if parsed_arguments.truth_vcf is None:
        raise ValueError("Must supply one or more --truth-vcf")
    if parsed_arguments.ped_file is None and parsed_arguments.vapor_json is None:
        if not _is_file(parsed_arguments.optimal_overlap_cutoffs_file):
            raise ValueError("Provide --optimal-overlap-cutoffs-file to use previously computed optimal cutoffs, *OR* "
                             "--ped-file to calculate new optimal cutoffs.")
    return parsed_arguments


def main(argv: Optional[List[Text]] = None) -> ConfidentVariants:
    arguments = __parse_arguments(sys.argv if argv is None else argv)
    confident_variants = get_truth_overlap(
        test_vcfs=arguments.test_vcf,
        truth_vcfs=arguments.truth_vcf,
        ped_files=None if arguments.ped_file is None else arguments.ped_file,
        vapor_files=get_vapor_files(arguments.vapor_json),
        optimal_overlap_cutoffs_file=arguments.optimal_overlap_cutoffs_file,
        expand_point_svs_bp=arguments.expand_point_svs_bp,
        point_sv_scale_factor=arguments.point_sv_scale_factor,
        expand_non_point_svs_bp=arguments.expand_non_point_svs_bp,
        non_point_sv_scale_factor=arguments.non_point_sv_scale_factor,
        f_beta=arguments.f_beta,
        min_overlap_cutoff_precision=arguments.min_overlap_cutoff_precision,
        min_vapor_precision=arguments.min_vapor_precision,
        inheritance_af_rareness=arguments.inheritance_af_rareness,
        use_copy_number=arguments.use_copy_number,
        num_threads=arguments.num_threads
    )
    output_confident_variants(confident_variants, output_file=arguments.output)
    return confident_variants


if __name__ == "__main__":
    main()
