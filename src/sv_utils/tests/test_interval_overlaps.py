import os
import glob
import numpy
import pandas
import pytest
from typing import Tuple, Union, Dict, Collection

from sv_utils import interval_overlaps, genomics_io, get_truth_overlap
import common_test_utils

Record = Union[numpy.record, Tuple]


class Default:
    num_check_itrs = 2
    num_contig_1a = 50
    num_contig_2a = 100
    contig_2a_width = 21
    num_contig_1b = 45
    num_contig_2b = 35
    contig_1b_width = 10
    # if True, randomize the order of intervals in each table (in a way consistent
    # with necessary correlations, e.g. intervals and expected values scrambled
    # identically). (Don't alter index sort order.)
    scramble_data_order = True
    # if True, alter the index of each table to create gaps, as though a subset of a
    # true "initial" table was passed
    simulate_subset = True
    # if True, reorder the index to be out of sorted order (don't alter data order)
    scramble_index_order = True
    resources_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources")
    small_vcfs_dir = os.path.join(resources_dir, "small_vcfs")
    small_vcfs = tuple(sorted(glob.glob(os.path.join(small_vcfs_dir, "*.vcf.gz"))))
    small_vcf = small_vcfs[0]
    wanted_properties = interval_overlaps.Default.wanted_properties


class Keys:
    contig = genomics_io.Keys.contig
    begin = genomics_io.Keys.begin
    end = genomics_io.Keys.end
    value = "value"
    num_overlap = "num_overlap"
    value_sum = "value_sum"
    primary_id = interval_overlaps.Keys.primary_id
    is_single_interval = interval_overlaps.Keys.is_single_interval


def get_interval_overlaps_test_data(
        num_contig_1a: int = Default.num_contig_1a,
        num_contig_2a: int = Default.num_contig_2a,
        contig_2a_width: int = Default.contig_2a_width,
        scramble_data_order: bool = Default.scramble_data_order,
        simulate_subset: bool = Default.simulate_subset,
        scramble_index_order: bool = Default.scramble_index_order,
) -> (pandas.DataFrame, pandas.DataFrame):
    """
    This seemingly crazy function is designed to produce intervals that have non-trivial overlaps that nevertheless have
    a closed-form mathematical solution for the overlap function. This function produces the intervals and the solution
    to the self-overlap problem.
    So this function tests the mathematical correctness of the overlap code, without testing its ability to load VCF or
    BED files.
    Args:
        num_contig_1a: int
            number of intervals to generate in contig 1
        num_contig_2a: int
            number of intervals to generate in contig 2
        contig_2a_width: int
            how wide to make the intervals in contig 2 (contig 1 intervals grow in width as they progress in position)
        scramble_data_order: bool
            if True, scramble the order of rows in the data set. (tests if overlap code can handle out-of-order data)
        simulate_subset: bool
            if True, alter the index so that it has gaps in it, rather than being numbered contiguously. (tests if the
            code somehow relies on sequential aspects of the index, which it shouldn't)
        scramble_index_order: bool
            if True, alter the index so that it is out of order, independent of data order. (tests if the code requires
            the index to be in order, which it shouldn't)
    Returns:
        intervals_a: pandas.DataFrame
            table of genomic intervals that have non-trivial overlaps that have a closed-form solution
        self_overlap_values: pandas.DataFrame
            table of values that describe the closed-form solution to overlaps between intervals in intervals_a
    """
    intervals_a = pandas.DataFrame(
        [
            {Keys.contig: '1', Keys.begin: x, Keys.end: 2 * x + 1, Keys.value: x} for x in range(num_contig_1a)
        ] +
        [
            {Keys.contig: '2', Keys.begin: x, Keys.end: x + contig_2a_width, Keys.value: x}
            for x in range(num_contig_2a)
        ]
    )
    intervals_a[Keys.contig] = intervals_a[Keys.contig].astype("category")
    num_a = len(intervals_a)

    # solve the math problem: how many self-overlaps should each interval in intervals_a have, and what should the sum
    # of overlapping values be
    def _self_overlap_1(x: int) -> Dict[str, int]:
        left_start = int(numpy.ceil(x / 2.0))
        right_stop = min(2 * x + 1, num_contig_1a)
        num_overlap = (x - left_start) + (right_stop - x - 1)
        value_sum = sum(v for v in range(left_start, x)) + sum(v for v in range(x + 1, right_stop))

        return {Keys.num_overlap: num_overlap, Keys.value_sum: value_sum}

    def _self_overlap_2(x: int) -> Dict[str, int]:
        left_start = max(0, x - contig_2a_width + 1)
        right_stop = min(x + contig_2a_width, num_contig_2a)
        num_overlap = (x - left_start) + (right_stop - x - 1)
        value_sum = sum(v for v in range(left_start, x)) + sum(v for v in range(x + 1, right_stop))

        return {Keys.num_overlap: num_overlap, Keys.value_sum: value_sum}

    self_overlap_values = pandas.DataFrame(
        [_self_overlap_1(x) for x in range(num_contig_1a)] +
        [_self_overlap_2(x) for x in range(num_contig_2a)]
    )

    # Scramble data and mess up index in various ways. The algorithm should not rely on any aspect of ordering of the
    # index or data, nor on it being a single block. These checks have caught bugs before, and should not be
    # streamlined away.
    if scramble_data_order:
        perm = numpy.random.permutation(num_a)
        intervals_a = intervals_a.iloc[perm].reset_index().drop("index", axis=1)
        self_overlap_values = self_overlap_values.iloc[perm].reset_index().drop("index", axis=1)

    if simulate_subset:
        new_index = numpy.sort(numpy.random.permutation(2 * num_a)[:num_a])
        intervals_a.index = new_index
        self_overlap_values.index = new_index

    if scramble_index_order:
        new_index = numpy.random.permutation(intervals_a.index)
        intervals_a.index = new_index
        self_overlap_values.index = new_index

    return intervals_a, self_overlap_values


def get_interval_overlaps_other_data(
        num_contig_1a: int = Default.num_contig_1a,
        num_contig_2a: int = Default.num_contig_2a,
        contig_2a_width: int = Default.contig_2a_width,
        num_contig_1b: int = Default.num_contig_1b,
        num_contig_2b: int = Default.num_contig_2b,
        contig_1b_width: int = Default.contig_1b_width,
        scramble_data_order: bool = Default.scramble_data_order,
        simulate_subset: bool = Default.simulate_subset,
        scramble_index_order: bool = Default.scramble_index_order
) -> (pandas.DataFrame, pandas.DataFrame, pandas.DataFrame):
    """
    This seemingly crazy function is designed to produce intervals that have non-trivial overlaps that nevertheless have
    a closed-form mathematical solution for the overlap function. This function produces two sets of intervals, and the
    solution to the problem of overlaps between them.
    So this function tests the mathematical correctness of the overlap code, without testing its ability to load VCF or
    BED files.
    Args:
        num_contig_1a: int
            number of 'a' intervals to generate in contig 1
        num_contig_2a: int
            number of 'a' intervals to generate in contig 2
        contig_2a_width: int
            how wide to make the 'a' intervals in contig 2 (contig 1 intervals grow in width as they progress in
            position)
        num_contig_1b: int
            number of 'b' intervals to generate in contig 1
        num_contig_2b: int
            number of 'b' intervals to generate in contig 2
        contig_1b_width: int
            how wide to make the 'b' intervals in contig 1 (contig 2 intervals grow in width as they progress in
            position)
        scramble_data_order: bool
            if True, scramble the order of rows in the data set. (tests if overlap code can handle out-of-order data)
        simulate_subset: bool
            if True, alter the index so that it has gaps in it, rather than being numbered contiguously. (tests if the
            code somehow relies on sequential aspects of the index, which it shouldn't)
        scramble_index_order: bool
            if True, alter the index so that it is out of order, independent of data order. (tests if the code requires
            the index to be in order, which it shouldn't)
    Returns:
        intervals_a: pandas.DataFrame
            table of genomic intervals that have non-trivial overlaps that have a closed-form solution
        intervals_b: pandas.DataFrame
            table of genomic intervals that have non-trivial overlaps that have a closed-form solution, distinct from
            intervals_a
        other_overlap_values: pandas.DataFrame
            table of values that describe the closed-form solution to overlaps between intervals in intervals_a and
            intervals_b

    """
    intervals_a = get_interval_overlaps_test_data(
        num_contig_1a=num_contig_1a, num_contig_2a=num_contig_2a,
        contig_2a_width=contig_2a_width, scramble_data_order=False,
        simulate_subset=False, scramble_index_order=False
    )[0]
    num_a = len(intervals_a)

    intervals_b = pandas.DataFrame(
        [
            {Keys.contig: '1', Keys.begin: x, Keys.end: x + contig_1b_width, Keys.value: x}
            for x in range(num_contig_1b)
        ] +
        [
            {Keys.contig: '2', Keys.begin: 2 * x, Keys.end: 3 * x + 1, Keys.value: x} for x in range(num_contig_2b)
        ]
    )
    intervals_b[Keys.contig] = intervals_b[Keys.contig].astype("category")
    num_b = len(intervals_b)

    # solve the math problem: how many intervals_b overlaps should each interval in intervals_a have, and what should
    # the sum of overlapping values be
    def _other_overlap_1(x: int) -> Dict[str, int]:
        left_start = max(0, x - contig_1b_width + 1)
        right_stop = min(2 * x + 1, num_contig_1b)
        num_overlap = right_stop - left_start
        value_sum = sum(v for v in range(left_start, right_stop))

        return {Keys.num_overlap: num_overlap, Keys.value_sum: value_sum}

    def _other_overlap_2(x: int) -> (int, int):
        left_start = int(numpy.ceil(x / 3.0))
        right_stop = min(int(numpy.ceil((x + contig_2a_width) / 2.0)),
                         num_contig_2b)
        num_overlap = right_stop - left_start
        value_sum = sum(v for v in range(left_start, right_stop))

        return {Keys.num_overlap: num_overlap, Keys.value_sum: value_sum}

    other_overlap_values = pandas.DataFrame(
        [_other_overlap_1(x) for x in range(num_contig_1a)] +
        [_other_overlap_2(x) for x in range(num_contig_2a)]
    )

    # Scramble data and mess up index in various ways. The algorithm should not rely on any aspect of ordering of the
    # index or data, nor on it being a single block. These checks have caught bugs before, and should not be
    # streamlined away.
    if scramble_data_order:
        perm = numpy.random.permutation(num_a)
        intervals_a = intervals_a.iloc[perm].reset_index(drop=True)
        other_overlap_values = other_overlap_values.iloc[perm].reset_index(drop=True)
        perm = numpy.random.permutation(num_b)
        intervals_b = intervals_b.iloc[perm].reset_index(drop=True)

    if simulate_subset:
        new_index = numpy.sort(numpy.random.permutation(2 * num_a)[:num_a])
        intervals_a.index = new_index
        other_overlap_values.index = new_index
        new_index = numpy.sort(numpy.random.permutation(2 * num_b)[:num_b])
        intervals_b.index = new_index

    if scramble_index_order:
        new_index = numpy.random.permutation(intervals_a.index)
        intervals_a.index = new_index
        other_overlap_values.index = new_index

        intervals_b.index = numpy.random.permutation(intervals_b.index)

    return intervals_a, intervals_b, other_overlap_values


def check_overlap_values_match_expected(
        capsys,
        overlap_values: pandas.DataFrame,
        expected_overlap_values: pandas.DataFrame,
        intervals_a: pandas.DataFrame,
        intervals_b: pandas.DataFrame
):
    """
    Helper function to automate the task of checking if the call to compute overlap values produced identical results to
    those derived by closed-form solution.
    """
    is_self_overlap = id(intervals_a) == id(intervals_b)
    assert len(overlap_values) == len(expected_overlap_values), "wrong length"
    for (r_ind, properties_row), (er_ind, expected_properties_row) in zip(
            overlap_values.iterrows(),
            expected_overlap_values.iterrows()
    ):
        # if things don't match, spit out warning detail before failing
        if not properties_row.equals(expected_properties_row):
            with capsys.disabled():
                print(f"row {r_ind}: found != expected")
                print(f"\t{str(properties_row)}")
                print(f"\t{str(expected_properties_row)}")
                interval_a = intervals_a.loc[r_ind]
                print(f"\t{interval_a}")
                # manually find all overlaps (slow for large datasets, otherwise we'd just do this!)
                # noinspection PyTypeChecker
                b_overlaps: pandas.Series = (intervals_b[Keys.contig] == interval_a[Keys.contig]) & \
                                            (intervals_b[Keys.begin] < interval_a[Keys.end]) & \
                                            (intervals_b[Keys.end] > interval_a[Keys.begin])
                # exclude self-overlaps
                if is_self_overlap:
                    b_overlaps.loc[r_ind] = False
                all_overlappers_str = intervals_b.loc[b_overlaps].sort_values([Keys.begin, Keys.end]).reset_index()\
                    .to_string().replace('\n', "\n\t")
                print(f"\t{all_overlappers_str}")
                print(f"\tfound {b_overlaps.sum()} overlaps")
        assert properties_row.equals(expected_properties_row)


# noinspection PyUnusedLocal
def _overlap_func(
        interval: Record,
        overlappers: pandas.DataFrame,
) -> (int, int):
    f"""
    Calculate overlap function for tst
    Args:
        interval: numpy.record
            Record with info corresponding to a piece of Evidence (not used in this function)
        overlappers: numpy.recarray
            Recarray with Evidence info for evidence that overlaps interval
    Returns:
        num_overlap: int
            number of intervals that overlap the test interval
        overlap_sum: int
            sum of {Keys.value} of intervals that overlap the test interval
    """
    num_overlap = len(overlappers)
    overlap_sum = overlappers[Keys.value].sum()
    return num_overlap, overlap_sum


def test_self_overlap(capsys, num_check_itrs: int = Default.num_check_itrs):
    for itr in range(num_check_itrs):
        overlap_props = [Keys.num_overlap, Keys.value_sum]
        intervals_a, expected_self_overlap_values = get_interval_overlaps_test_data()
        self_overlap_values = interval_overlaps.apply_interval_overlap_func(
            _overlap_func, intervals_a, description="finding self overlaps", property_names=overlap_props
        )

        check_overlap_values_match_expected(
            capsys, self_overlap_values, expected_self_overlap_values, intervals_a, intervals_a
        )


def test_other_overlap(capsys, num_check_itrs: int = Default.num_check_itrs):
    for itr in range(num_check_itrs):
        intervals_a, intervals_b, expected_other_overlap_values \
            = get_interval_overlaps_other_data()
        expected_other_overlap_values = expected_other_overlap_values.loc[
            intervals_a.index
        ]
        other_overlap_values = interval_overlaps.apply_interval_overlap_func(
            _overlap_func, intervals_a, intervals_b,
            description="finding other overlaps",
            property_names=('num_overlap', 'value_sum')
        )

        check_overlap_values_match_expected(
            capsys, other_overlap_values, expected_other_overlap_values, intervals_a, intervals_b
        )


@pytest.mark.integration_test
def test_real_variants_overlap(small_vcf: str = Default.small_vcf,
                               wanted_properties: Collection[str] = Default.wanted_properties):
    # Load actual variants from small multi-sample VCF that's been curated to hold a variety of SV types, including
    # BND and CPX. In this case, we don't care about genotypes, so only load variant properties necessary for computing
    # overlaps
    variants = genomics_io.vcf_to_pandas(small_vcf, wanted_properties=wanted_properties)
    # Break down multi-interval variants into table of simple intervals
    variant_simple_intervals = interval_overlaps.fix_variants(variants)
    # Compute overlaps amongst the simple intervals
    simple_interval_overlaps = interval_overlaps.apply_interval_overlap_func(
        get_truth_overlap.quantify_overlap, variant_simple_intervals
    )
    # Apply mathematical reduction operations (min, max, mean, etc) to reduce simple interval overlaps onto variants
    variant_overlaps = interval_overlaps.postprocess_bnd_complex(simple_interval_overlaps=simple_interval_overlaps,
                                                                 simple_variants=variant_simple_intervals,
                                                                 original_index=variants.index)
    # now test some stuff:
    # * every variant should be represented as a primary_id in variant_simple_intervals, with no extras
    assert not set(variant_simple_intervals[Keys.primary_id].unique()).symmetric_difference(variants.index)
    # * there should be no negative simple overlaps
    assert not (simple_interval_overlaps < 0).any(axis=None)
    # * the index for simple_interval_overlaps should be identical to the index for variant_simple_intervals
    common_test_utils.assert_indices_equal(simple_interval_overlaps.index, variant_simple_intervals.index,
                                           context="simple_interval_overlaps.index != variant_simple_intervals.index")
    # * there should be no negative variant overlaps
    assert not (variant_overlaps < 0).any(axis=None)
    # * there should be some non-trivial variant overlaps
    assert (variant_overlaps > 0).any(axis=None)
    # * the index for variant_overlaps should be identical to the index for variants
    common_test_utils.assert_indices_equal(variant_overlaps.index, variants.index,
                                           context="variant_overlaps.index != variants.index")
