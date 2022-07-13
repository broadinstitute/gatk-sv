import os
import glob
import pytest
import numpy
from typing import List, Set

from sv_utils import make_cross_validation_vcfs, genomics_io, pedigree_tools


class Default:
    resources_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources")
    small_vcfs_dir = os.path.join(resources_dir, "small_vcfs")
    small_vcfs = tuple(sorted(glob.glob(os.path.join(small_vcfs_dir, "*.vcf.gz"))))
    small_vcf = small_vcfs[0]
    ped_file = os.path.join(resources_dir, "1KGP_2504_and_698_with_GIAB.ped")
    num_splits = make_cross_validation_vcfs.Default.num_splits
    index_output_vcf = True
    num_threads = make_cross_validation_vcfs.Default.num_threads


def _select_unrelated_sample(candidate_samples: Set[str], pedigree_file_info: pedigree_tools.PedigreeFileInfo,
                             selected_samples: List[str]):
    # select a random sample from the set of candidates
    next_sample = numpy.random.choice(list(candidate_samples))
    selected_samples.append(next_sample)
    # remove all related participants from the set of candidates
    related = {
        participant_id
        for index in pedigree_file_info.participant_indices[next_sample]
        for participant_id in pedigree_file_info.pedigree_lines[index].trio_ids
    }
    candidate_samples.difference_update(related)


def _choose_random_unrelated(sample_ids: List[str], num_selected: int,
                             pedigree_file_info: pedigree_tools.PedigreeFileInfo) -> List[str]:
    candidate_sample_ids = set(sample_ids)
    truth_sample_ids = []
    for __ in range(num_selected):
        _select_unrelated_sample(candidate_sample_ids, pedigree_file_info, truth_sample_ids)
    return truth_sample_ids


@pytest.mark.integration_test
def test_make_cross_validation_vcfs(
        tmpdir,
        vcf: str = Default.small_vcf,
        ped_file: str = Default.ped_file,
        num_splits: int = Default.num_splits,
        index_output_vcf: bool = Default.index_output_vcf,
        num_threads: int = Default.num_threads
):
    # get the sample_ids from the vcf. cross-validation is over samples, not variants, so we'll need to check these
    sample_ids = genomics_io.get_vcf_sample_ids(vcf)
    # also get the variant IDs and properties, to check that they're not changed
    variant_ids = genomics_io.get_vcf_variant_ids(vcf)
    properties = sorted(genomics_io.get_vcf_properties(vcf))
    # make a temp dir, for writing out cross-validation vcfs, and also to make a "truth samples" file
    temp_out_dir = tmpdir.mkdir("test_make_cross_validation_vcfs")

    # load trio data
    pedigree_file_info = pedigree_tools.PedigreeFileInfo.load(ped_file)
    # noinspection PyTypeChecker
    num_trios = pedigree_file_info.subset_participants(sample_ids).num_trios
    expected_num_test_trios = num_trios / num_splits

    # randomly select num_splits samples to be "truth samples" that should be evenly distributed in the folds
    # Normally you'd pass any samples that had truth data, e.g. pacbio data. NOTE that make_cross_validation_vcfs
    # does not want to break up trios, so make an effort to randomly select unrelated individuals to be "truth" data
    truth_sample_ids = _choose_random_unrelated(sample_ids, num_splits, pedigree_file_info=pedigree_file_info)
    # write the truth samples to the file
    truth_samples_file = os.path.join(temp_out_dir, "truth_samples.txt")
    with open(truth_samples_file, 'w') as f_out:
        for sample_id in truth_sample_ids:
            f_out.write(sample_id + "\n")

    # call the main make_cross_validation_function
    cross_validation_vcfs = make_cross_validation_vcfs.make_cross_validation_vcfs(
        vcf=vcf,
        ped_file=ped_file,
        truth_samples_file=truth_samples_file,
        output_folder=temp_out_dir,
        num_splits=num_splits,
        index_output_vcf=index_output_vcf,
        num_threads=num_threads
    )
    # iterate over folds and ensure validity
    for train_vcf, test_vcf in cross_validation_vcfs:
        # first check that all the variants and properties are present
        if index_output_vcf:
            assert os.path.isfile(train_vcf + ".tbi")
            assert os.path.isfile(test_vcf + ".tbi")
        assert genomics_io.get_vcf_variant_ids(train_vcf).equals(variant_ids)
        assert genomics_io.get_vcf_variant_ids(test_vcf).equals(variant_ids)
        assert sorted(genomics_io.get_vcf_properties(train_vcf)) == properties
        assert sorted(genomics_io.get_vcf_properties(test_vcf)) == properties
        # next check that all the sample IDs are present, and in roughly the right ratios
        train_sample_ids = genomics_io.get_vcf_sample_ids(train_vcf)
        test_sample_ids = genomics_io.get_vcf_sample_ids(test_vcf)
        #     -all the sample IDs are present, and there are no extras
        assert sorted(train_sample_ids + test_sample_ids) == sorted(sample_ids)
        #     -they're divided more-or-less at the right ratios
        sample_split_ratio = len(test_sample_ids) / len(sample_ids)
        assert 1.0 / (num_splits + 0.5) <= sample_split_ratio <= 1.0 / (num_splits - 0.5)
        #     -there's one truth sample in the test VCF
        assert len(set(test_sample_ids).intersection(truth_sample_ids)) == 1
        #     -the trios are divided up appropriately
        # noinspection PyTypeChecker
        num_test_trios = pedigree_file_info.subset_participants(test_sample_ids).num_trios
        assert expected_num_test_trios / 2 <= num_test_trios <= expected_num_test_trios * 2
