#!/usr/bin/env python

import sys
import os
import warnings
import argparse
import tarfile
import pysam
import numpy
from typing import List, Text, Dict, Mapping, Iterable, Iterator, Callable, Tuple

# define type for function that matches sample IDs between cohort VCF and benchmark tar file
SampleIdMatcher = Callable[[List[str], List[str]], Iterator[Tuple[str, str]]]


def __parse_arguments(argv: List[Text]) -> argparse.Namespace:

    # noinspection PyTypeChecker
    parser = argparse.ArgumentParser(
        description="Get TSV file with mapping from sample IDs in external benchmarking set to corresponding sample IDs"
                    "in cohort VCF",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        prog=argv[0]
    )
    parser.add_argument("benchmark_tar", type=str, help="Tar file with bed files used for per-sample benchmarking")
    parser.add_argument("cohort_vcf", type=str, help="Cohort VCF to be used in MainVcfQc")
    parser.add_argument("save_tsv_file", type=str, help="File to save samples map in TSV format")

    parsed_arguments = parser.parse_args(argv[1:] if len(argv) > 1 else ["--help"])
    return parsed_arguments


def match_benchmark_contains_vcf_case_insensitive(
        benchmark_sample_ids: Iterable[str], vcf_sample_ids: Iterable[str],
) -> Iterator[Tuple[str, str]]:
    """ SampleIdMatcher that assumes benchmark sample IDs contain VCF sample IDs, possibly in a different case """
    def _make_searchable(_sample_ids: Iterable[str]) -> (List[str], List[str]):
        # given iterable of sample IDs, make version that is searchable (lower-case, with underscores removed), and sort
        # both the searchable sample IDs and originals according to the order of the searchable list
        return zip(
            *sorted([(_sample_id, _sample_id.lower().strip('_')) for _sample_id in _sample_ids],
                    key=lambda tup: tup[1])
        )
    benchmark_sample_ids, searchable_benchmark_ids = _make_searchable(benchmark_sample_ids)
    vcf_sample_ids, searchable_vcf_sample_ids = _make_searchable(vcf_sample_ids)
    for vcf_ind, (benchmark_ind, searchable_benchmark_id) in zip(
            numpy.searchsorted(searchable_vcf_sample_ids, searchable_benchmark_ids, side="left"),
            enumerate(searchable_benchmark_ids)
    ):
        if searchable_vcf_sample_ids[vcf_ind] in searchable_benchmark_id:
            yield benchmark_sample_ids[benchmark_ind], vcf_sample_ids[vcf_ind]
        elif vcf_ind > 0 and searchable_vcf_sample_ids[vcf_ind - 1] in searchable_benchmark_id:
            yield benchmark_sample_ids[benchmark_ind], vcf_sample_ids[vcf_ind - 1]
        else:
            warnings.warn(f"{searchable_benchmark_id} has no match")


def main(arguments: List[str]):
    """ Get command-line arguments, create map from benchmark sample ID to cohort VCF sample ID, save as TSV """
    options = __parse_arguments(arguments)
    samples_map = get_rename_benchmark_samples_map(
        benchmark_tar=options.benchmark_tar,
        cohort_vcf=options.cohort_vcf,
        get_sample_id_matches=match_benchmark_contains_vcf_case_insensitive
    )
    save_samples_map_tsv(samples_map=samples_map, save_tsv_file=options.save_tsv_file)


def get_rename_benchmark_samples_map(
        benchmark_tar: str,
        cohort_vcf: str,
        get_sample_id_matches: SampleIdMatcher = match_benchmark_contains_vcf_case_insensitive
) -> Dict[str, str]:
    """
    Get mapping from benchmark sample IDs to cohort sample IDs. This function is defined with passable matching function
    in case future cohort VCFs / benchmark data sets require a different logic.
    Args:
        benchmark_tar: str
            path to tar file with benchmarking BED files
        cohort_vcf: str
            path to cohort VCF
        get_sample_id_matches: SampleIdMatcher (default=match_benchmark_startswith_vcf_lower)
            Callable that takes list of sample IDs from the benchmark set and a list of sample IDs from the cohort VCF,
            and yields tuples of matching (benchmark_sample_id, cohort sample_id) pairs
    Returns:
        benchmark_samples_map: Dict[str, str]
            Mapping from benchmark sample ID to corresponding cohort VCF sample ID
    """
    benchmark_sample_ids = get_benchmark_sample_ids(benchmark_tar)
    vcf_sample_ids = get_vcf_sample_ids(cohort_vcf)
    return {
        benchmark_id: vcf_id for benchmark_id, vcf_id in get_sample_id_matches(benchmark_sample_ids, vcf_sample_ids)
    }


def get_vcf_sample_ids(vcf: str) -> List[str]:
    """ Extract sample IDs from VCF """
    with pysam.VariantFile(vcf, 'r') as f_in:
        return list(f_in.header.samples)


def get_benchmark_sample_ids(benchmark_tar: str) -> List[str]:
    """ Extract sample IDs from benchmark tar file """
    return [os.path.basename(bed_file).split('.', 1)[0] for bed_file in get_tar_files(benchmark_tar)
            if bed_file.endswith(".bed.gz")]


def get_tar_files(path_to_tar_file: str) -> List[str]:
    """ Get contents of tar file """
    with tarfile.open(path_to_tar_file, 'r') as tar_in:
        return [tar_info.name for tar_info in tar_in.getmembers()]


def save_samples_map_tsv(
        samples_map: Mapping[str, str],
        save_tsv_file: str
):
    """ Save mapping from benchmark sample ID to cohort sample ID in TSV format """
    with open(save_tsv_file, 'w') as f_out:
        for benchmark_sample_id, vcf_sample_id in samples_map.items():
            f_out.write(f"{benchmark_sample_id}\t{vcf_sample_id}\n")


if __name__ == "__main__":
    main(sys.argv)
