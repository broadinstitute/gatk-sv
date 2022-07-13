import os
import glob
import pytest
import json
from typing import Iterable, Dict, Sequence
import pysam
import pandas

from sv_utils import genomics_io, get_truth_overlap


class Default:
    resources_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources")
    small_vcfs_dir = os.path.join(resources_dir, "small_vcfs")
    short_reads_vcfs = tuple(sorted(glob.glob(os.path.join(small_vcfs_dir, "*.vcf.gz"))))
    long_reads_dir = os.path.join(resources_dir, "pacbio")
    long_reads_vcfs = tuple(sorted(glob.glob(os.path.join(long_reads_dir, "*.vcf.gz"))))
    vapor_dir = os.path.join(resources_dir, "vapor")
    vapor_files = tuple(sorted(glob.glob(os.path.join(vapor_dir, "*.bed.gz"))))
    ped_file = os.path.join(resources_dir, "1KGP_2504_and_698_with_GIAB.ped")
    optimal_overlap_cutoffs_file = os.path.join(resources_dir, "optimal_overlap_cutoffs.pickle")
    pacbio_no_vapor_sample_id = "HG00733"
    vapor_no_pacbio_sample_id = "HG00512"
    pacbio_and_vapor_sample_id = "HG00514"


@pytest.fixture(scope="function")
def short_reads_vcf(tmpdir, short_reads_vcfs: Sequence[str] = Default.short_reads_vcfs):
    """
    yield temporary vcf that is all short-reads VCFs concatenated together
    It will be deleted when the test using it ends
    """
    concat_vcf = os.path.join(tmpdir, "short-reads.vcf.gz")
    newline = '\n'.encode("utf-8")
    comment = '#'.encode("utf-8")
    with pysam.BGZFile(concat_vcf, "wb") as f_out:
        with pysam.BGZFile(short_reads_vcfs[0], "rb") as f_in:
            for line in f_in:
                f_out.write(line + newline)
        for vcf in short_reads_vcfs[1:]:
            with pysam.BGZFile(vcf, "rb") as f_in:
                for line in f_in:
                    if line.startswith(comment):
                        continue
                    f_out.write(line + newline)

    yield concat_vcf
    if os.path.isfile(concat_vcf):
        os.remove(concat_vcf)


def get_vapor_dict(vapor_files: Iterable[str]) -> Dict[str, str]:
    """
    yield temporary JSON file for vapor data that will be deleted when the test function ends
    """
    def _vapor_file_to_sample_id(vapor_file: str) -> str:
        return os.path.basename(vapor_file).split('-', maxsplit=1)[0]

    return {_vapor_file_to_sample_id(vapor_file): vapor_file for vapor_file in vapor_files}


def test_read_vapor_dict(
    tmpdir: str,
    vapor_files: Iterable[str] = Default.vapor_files
):
    vapor_dict = get_vapor_dict(vapor_files)
    vapor_json = os.path.join(tmpdir, "vapor.json")
    with open(vapor_json, 'w') as f_out:
        json.dump(vapor_dict, f_out)

    assert get_truth_overlap.get_vapor_files(vapor_json) == vapor_dict


def _num_confident_variants(sample_confident_variants: get_truth_overlap.SampleConfidentVariants) -> int:
    return len(sample_confident_variants.good_variant_ids) + len(sample_confident_variants.bad_variant_ids)


def assert_has_confident_variants(confident_variants: get_truth_overlap.ConfidentVariants, sample_id: str):
    assert sample_id in confident_variants, f"{sample_id} is not a key in confident_variants"
    assert _num_confident_variants(confident_variants[sample_id]) > 0


def assert_has_no_confident_variants(confident_variants: get_truth_overlap.ConfidentVariants, sample_id: str):
    assert sample_id not in confident_variants or _num_confident_variants(confident_variants[sample_id]) == 0, \
        f"{sample_id} is not a key in confident_variants"


def assert_all_variant_ids_valid(confident_variants: get_truth_overlap.ConfidentVariants,
                                 valid_variant_ids: pandas.Index):
    for sample_id, sample_confident_variants in confident_variants.items():
        invalid_ids = sample_confident_variants.all_confident_variant_ids.difference(valid_variant_ids)
        assert not invalid_ids, f"In sample_confident_variants for {sample_id}, invalid variants: {invalid_ids}"


@pytest.mark.integration_test
def test_get_truth_overlap(
        capsys,
        short_reads_vcf,
        vapor_files: Iterable[str] = Default.vapor_files,
        long_reads_vcfs: Iterable[str] = Default.long_reads_vcfs,
        ped_file: str = Default.ped_file,
        optimal_overlap_cutoffs_file: str = Default.optimal_overlap_cutoffs_file,
        pacbio_no_vapor_sample_id: str = Default.pacbio_no_vapor_sample_id,
        vapor_no_pacbio_sample_id: str = Default.vapor_no_pacbio_sample_id,
        pacbio_and_vapor_sample_id: str = Default.pacbio_and_vapor_sample_id
):
    valid_variant_ids = genomics_io.get_vcf_variant_ids(short_reads_vcf)

    vapor_dict = get_vapor_dict(vapor_files)
    confident_variants = get_truth_overlap.get_truth_overlap(
        test_vcfs=(short_reads_vcf,), truth_vcfs=long_reads_vcfs, ped_files=ped_file, vapor_files=vapor_dict,
        optimal_overlap_cutoffs_file=optimal_overlap_cutoffs_file
    )

    # The mechanism of determining confident variants is different for each of these cases. Ensure it's working:
    assert_has_confident_variants(confident_variants, pacbio_no_vapor_sample_id)
    assert_has_confident_variants(confident_variants, vapor_no_pacbio_sample_id)
    assert_has_confident_variants(confident_variants, pacbio_and_vapor_sample_id)
    assert_all_variant_ids_valid(confident_variants, valid_variant_ids)

    # make sure that everything works without vapor
    no_vapor_confident_variants = get_truth_overlap.get_truth_overlap(
        test_vcfs=(short_reads_vcf,), truth_vcfs=long_reads_vcfs, ped_files=ped_file, vapor_files=None,
        optimal_overlap_cutoffs_file=optimal_overlap_cutoffs_file
    )
    assert_has_confident_variants(no_vapor_confident_variants, pacbio_no_vapor_sample_id)
    assert_has_no_confident_variants(no_vapor_confident_variants, vapor_no_pacbio_sample_id)
    assert_has_confident_variants(no_vapor_confident_variants, pacbio_and_vapor_sample_id)
    assert_all_variant_ids_valid(no_vapor_confident_variants, valid_variant_ids)
    with capsys.disabled():
        get_truth_overlap.output_confident_variants(no_vapor_confident_variants, "-")
