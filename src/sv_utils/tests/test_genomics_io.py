import os
import glob
import numpy
import pandas
from typing import Sequence, Set, Tuple

from sv_utils import common, genomics_io
import pysam
import common_test_utils


class Keys:
    id = genomics_io.Keys.id
    contig = genomics_io.Keys.contig
    begin = genomics_io.Keys.begin
    end = genomics_io.Keys.end


class Default:
    num_scrambled_trials = 10
    resources_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources")
    small_vcfs_dir = os.path.join(resources_dir, "small_vcfs")
    small_vcfs = tuple(sorted(glob.glob(os.path.join(small_vcfs_dir, "*.vcf.gz"))))
    small_vcf = small_vcfs[0]
    small_beds_dir = os.path.join(resources_dir, "small_beds")
    small_beds = tuple(sorted(glob.glob(os.path.join(small_beds_dir, "chr*.bed.gz"))))
    small_bed = small_beds[0]
    test_bed_literal_eval_columns = frozenset({"ac", "algorithms", "alts", "cpx_intervals", "evidence", "filter"})
    num_subset_samples = 2
    num_subset_properties = 2
    num_subset_variants = 2


def get_nontrivial_permutation(permutation_length: int) -> numpy.ndarray:
    """
    Get a non-identity permutation. For longer permutations, it should be non-identity basically every time, but this
    makes it guaranteed
    Args:
        permutation_length: int
            Number of items to permute
    Returns:
        permutation: numpy.ndarray
            Random non-identity permutation
    """
    permutation = numpy.random.permutation(permutation_length)
    while numpy.array_equal(permutation, numpy.arange(permutation_length)):
        permutation = numpy.random.permutation(permutation_length)
    return permutation


def test_intervals_sort(num_scrambled_trials: int = Default.num_scrambled_trials):
    # to hit all branches, do this test three ways
    for test_method in range(3):
        sorted_intervals = pandas.DataFrame(
            [
                ['a', "chr1", 100, 200],
                ['b', "chr1", 100, 10000],
                ['c', "chr1", 101, 200],
                ['d', "chr1", 101, 201],
                ['e', "chr2", 50, 100],
                ['f', "chr2", 100, 200],
                ['g', "chr10", 100, 200],
                ['h', "chrX", 100, 200],
                ['i', "chrY", 100, 200],
                ['j', "chrM", 100, 200],
                ['k', "chrUn_KI270302v1", 100, 200]
            ],
            columns=[Keys.id, Keys.contig, Keys.begin, Keys.end],
        ).set_index(Keys.id)

        if test_method == 0:  # for test_method == 0 just leave the sorted_intervals as-is
            method_str = "object contig"
        else:   # make contig categorical
            categories = sorted_intervals[Keys.contig].unique()  # <- delivers in order of appearance, so sorted
            if test_method == 2:  # scramble the order so the categories are not in sorted order
                categories = categories[get_nontrivial_permutation(len(categories))]
                method_str = "unordered category contig"
            else:
                method_str = "ordered category contig"
            sorted_intervals[Keys.contig] = pandas.Categorical(
                sorted_intervals[Keys.contig].values, categories=categories, ordered=True
            )

        for trail_ind in range(num_scrambled_trials):
            # get a non-identity permutation (should be basically every time, but this makes it deterministic):
            permutation = get_nontrivial_permutation(len(sorted_intervals))

            scrambled_intervals = sorted_intervals.iloc[permutation]
            assert not scrambled_intervals.equals(sorted_intervals)  # the scramble should have done something
            common_test_utils.assert_dataframes_equal(
                genomics_io.sort_intervals_table(scrambled_intervals, inplace=False, drop_index=False),
                sorted_intervals,
                f"{method_str}, inplace=False, drop_index=False"
            )
            common_test_utils.assert_dataframes_equal(
                genomics_io.sort_intervals_table(scrambled_intervals, inplace=False, drop_index=True),
                sorted_intervals.reset_index(drop=True, inplace=False),
                f"{method_str}, inplace=False, drop_index=True"
            )

            genomics_io.sort_intervals_table(scrambled_intervals, inplace=True, drop_index=False)
            common_test_utils.assert_dataframes_equal(scrambled_intervals, sorted_intervals,
                                                      f"{method_str}, inplace=True, drop_index=False")

            # re-scramble and try again with drop_index=True
            scrambled_intervals = sorted_intervals.iloc[permutation]
            genomics_io.sort_intervals_table(scrambled_intervals, inplace=True, drop_index=True)
            common_test_utils.assert_dataframes_equal(
                scrambled_intervals, sorted_intervals.reset_index(drop=True, inplace=False),
                f"{method_str}, inplace=True, drop_index=True"
            )


def test_get_sample_ids(
        small_vcf: str = Default.small_vcf
):
    # get the sample IDs by grepping the header
    samples_grep = common.command_results(f"zgrep -m1 ^#[^#] {small_vcf} | cut -f 10- | less").split()
    # get the sample IDs by parsing the VCF header
    sample_ids_parsed_header = genomics_io.get_vcf_sample_ids(small_vcf)
    # get the sample IDs by loading the variants and getting the samples from the table
    sample_ids_table = genomics_io.get_sample_ids_from_variants(
        genomics_io.vcf_to_pandas(small_vcf, wanted_properties=(genomics_io.Keys.gt,))
    )
    assert samples_grep == sample_ids_parsed_header
    # because of the call to unique, the sample_ids_table could be out of order
    assert sorted(samples_grep) == sorted(sample_ids_table)


def test_subset_vcf(
        vcf: str = Default.small_vcf,
        num_subset_samples: int = Default.num_subset_samples,
        num_subset_properties: int = Default.num_subset_properties,
):
    # load all the variants
    whole_variants = genomics_io.vcf_to_pandas(vcf)

    # while we're at it, test get_vcf_variant_ids
    variant_ids = genomics_io.get_vcf_variant_ids(vcf)
    assert whole_variants.index.equals(variant_ids)

    # get subset of samples
    subset_samples = numpy.random.choice(genomics_io.get_vcf_sample_ids(vcf), num_subset_samples, replace=False)\
        .tolist()
    # get a subset of variant properties (no associated sample) and sample properties (contained in format)
    variant_properties = whole_variants.loc[:, (None, slice(None))].columns\
        .get_level_values(genomics_io.Keys.property).to_list()
    sample_properties = whole_variants.loc[:, (subset_samples[0], slice(None))].columns\
        .get_level_values(genomics_io.Keys.property).to_list()
    subset_properties = numpy.concatenate((
        numpy.random.choice(variant_properties, num_subset_properties),
        numpy.random.choice(sample_properties, num_subset_properties)
    )).tolist()

    for check_subset_samples in [True, False]:
        for check_subset_properties in [True, False]:
            if not (check_subset_samples or check_subset_properties):
                continue
            load_variants = genomics_io.vcf_to_pandas(
                vcf, wanted_properties=subset_properties if check_subset_properties else None,
                samples=subset_samples if check_subset_samples else None
            )
            samp_select = [None] + subset_samples if check_subset_samples else slice(None)
            prop_select = subset_properties if check_subset_properties else slice(None)
            manual_variants = whole_variants.loc[:, (samp_select, prop_select)]
            common_test_utils.assert_dataframes_equal(
                load_variants, manual_variants,
                context=f"subset_samples:{check_subset_samples}, subset_properties:{check_subset_properties}",
                check_column_order=False
            )


def _get_uncompressed_text(filename: str) -> Tuple[str, ...]:
    with pysam.BGZFile(filename, "rb") as f_in:
        return tuple(f_in)


def test_read_write_bed(
        tmpdir,
        small_vcf: str = Default.small_vcf,
        small_bed: str = Default.small_bed,
        test_bed_literal_eval_columns: Set[str] = Default.test_bed_literal_eval_columns
):
    temp_out_dir = tmpdir.mkdir("test_write_bed")
    # check that reading in the same variants as BED files or VCFs yields the same variants
    # for the VCF, restrict to variant-level data (no genotype / FORMAT data)
    # also: don't try to deal with breaking complex/BND into multiple intervals, just take the primary interval
    small_vcf_variants = genomics_io.vcf_to_pandas(small_vcf, samples={}, drop_missing_columns=True,
                                                   drop_trivial_multi_index=True)
    small_bed_variants = genomics_io.bed_to_pandas(small_bed, literal_eval_columns=test_bed_literal_eval_columns)
    # note that the order of columns is likely to be different
    assert not set(small_vcf_variants.columns).symmetric_difference(small_bed_variants.columns)
    common_test_utils.assert_dataframes_equal(small_bed_variants, small_vcf_variants,
                                              context=f"vcf={small_vcf}, bed={small_bed}", check_column_order=False)
    # check that writing the variants to a bed file is identical (when uncompressed) to the existing bed file. check the
    # re-ordering mechanism works to create a valid bed file
    temp_out_bed = os.path.join(temp_out_dir, os.path.basename(small_bed))
    genomics_io.pandas_to_bed(temp_out_bed, small_vcf_variants)
    assert _get_uncompressed_text(temp_out_bed) == _get_uncompressed_text(small_bed)
    # Ensure this isn't a trivial test (empty files, text getter somehow doesn't work). VCF and BED should be different
    assert _get_uncompressed_text(temp_out_bed) != _get_uncompressed_text(small_vcf)


def test_vcat_with_categoricals(small_vcfs: Sequence[str] = Default.small_vcfs,
                                num_subset_samples: int = Default.num_subset_samples):
    # We have a lot of samples. To save time just test with a few, if a few work, they'll all work
    subset_samples = numpy.random.choice(genomics_io.get_vcf_sample_ids(small_vcfs[0]), num_subset_samples,
                                         replace=False)
    small_dataframes = [
        genomics_io.vcf_to_pandas(small_vcf, samples=subset_samples) for small_vcf in small_vcfs
    ]
    # note the main purpose of this method (vs just calling pandas.concat) is managing concatenation of potentially
    # incompatible categoricals. The small vcfs were chosen to have non-overlapping contigs, so the main "test" is that
    # this runs with no exception
    joint_dataframe = genomics_io.vcat_with_categoricals(small_dataframes)
    assert sum(len(small_dataframe) for small_dataframe in small_dataframes) == len(joint_dataframe)
    for small_vcf, small_dataframe in zip(small_vcfs, small_dataframes):
        common_test_utils.assert_dataframes_equal(small_dataframe, joint_dataframe.loc[small_dataframe.index],
                                                  f"inclusion of {small_vcf} variants")
