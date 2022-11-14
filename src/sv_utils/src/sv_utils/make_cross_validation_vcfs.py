#!/usr/bin/env python

import sys
import os
from typing import Tuple, List, Text, Optional, Iterator, Collection, Set
import argparse
import numpy
import pysam
from sv_utils import common, genomics_io, pedigree_tools


class Default:
    num_splits = 5
    num_threads = common.num_physical_cpus
    output_folder = os.getcwd()
    index_output_vcf = False


def _divy_family_samples(family_sample_ids: List[Set[Text]], batch_ind: int, num_splits: int) -> Iterator[Text]:
    # try to evenly distribute
    num_divy = num_splits - batch_ind
    num_divied_samples = sum(len(family) for family in family_sample_ids) // num_divy
    while num_divied_samples > 0:
        divied_family = family_sample_ids.pop()
        num_divied_samples -= len(divied_family)
        yield from divied_family


def _divy_isolated_samples(isolated_sample_ids: Set[Text], batch_ind: int, num_splits: int) -> Iterator[Text]:
    # try to evenly distribute
    num_divy = num_splits - batch_ind
    num_divied_samples = len(isolated_sample_ids) // num_divy
    while num_divied_samples > 0:
        yield isolated_sample_ids.pop()
        num_divied_samples -= 1


def _iter_sample_splits(
        sample_ids: numpy.ndarray,
        pedigree_file_info: Optional[pedigree_tools.PedigreeFileInfo] = None,
        truth_samples: Optional[Collection[str]] = None,
        num_splits: int = Default.num_splits
) -> Iterator[Tuple[numpy.ndarray, numpy.ndarray]]:
    """
    Divide up samples into num_splits batches, trying to keep trios and truth samples evenly divided
    Args:
        sample_ids:
        pedigree_file_info:
        truth_samples:
        num_splits:

    Returns:

    """
    sample_order = {
        sample_id: index for index, sample_id in enumerate(sample_ids)
    }
    truth_samples = {} if truth_samples is None else set(truth_samples)
    # get tuple of sample_ids in complete trio structures, don't break up families if they exist (e.g. quads)

    family_sample_ids = [] if pedigree_file_info is None else [
        family_participant_ids
        for __, family_participant_ids in
        pedigree_file_info.subset_participants(sample_ids).iter_family_participant_ids()
    ]

    isolated_sample_ids = set(sample_ids).difference(sample_id for family in family_sample_ids for sample_id in family)
    isolated_truth_samples = truth_samples.intersection(isolated_sample_ids)
    isolated_non_truth_samples = isolated_sample_ids.difference(isolated_truth_samples)
    truth_family_sample_ids = [
        family for family in family_sample_ids if any(sample_id in truth_samples for sample_id in family)
    ]
    non_truth_family_sample_ids = [
        family for family in family_sample_ids if not any(sample_id in truth_samples for sample_id in family)
    ]
    # now evenly divy up into batches, in order of decreasing desirability
    # 1) truth_family_sample_ids (samples in trios with truth information)
    # 2) non_truth_family_sample_ids (samples in trios with no truth information)
    # 3) isolated_truth_samples (samples not in a trio that have truth information)
    # 4) isolated_non_truth_samples (samples not in a trio with no truth information)
    batches = [set() for _ in range(num_splits)]
    if truth_family_sample_ids:
        for batch_ind, batch in enumerate(batches):
            batch.update(_divy_family_samples(truth_family_sample_ids, batch_ind, num_splits))
    if non_truth_family_sample_ids:
        for batch_ind, batch in enumerate(batches):
            batch.update(_divy_family_samples(non_truth_family_sample_ids, batch_ind, num_splits))
    if isolated_truth_samples:
        for batch_ind, batch in enumerate(batches):
            batch.update(_divy_isolated_samples(isolated_truth_samples, batch_ind, num_splits))
    if isolated_non_truth_samples:
        for batch_ind, batch in enumerate(batches):
            batch.update(_divy_isolated_samples(isolated_non_truth_samples, batch_ind, num_splits))

    # Each batch serves as a test set. For each test set, glue other batches into a training set
    for batch_ind, test_batch in enumerate(batches):
        training_batch = {
            sample_id
            for batch in (batches[:batch_ind] + batches[batch_ind + 1:])
            for sample_id in batch
        }

        # yield the samples in the same order as the original VCF, in case that somehow matters
        yield sample_ids.take(sorted(sample_order[sample_id] for sample_id in training_batch)), \
            sample_ids.take(sorted(sample_order[sample_id] for sample_id in test_batch))


def subsample_vcf(
        vcf: str,
        sample_ids: numpy.ndarray,
        out_vcf: str,
        index_output_vcf: bool = Default.index_output_vcf,
        num_threads: int = Default.num_threads
) -> str:
    print(f"Creating {out_vcf} with {len(sample_ids)} samples ...", end='', flush=True, file=sys.stderr)
    with pysam.VariantFile(vcf, 'r', threads=num_threads) as f_in:
        # noinspection PyTypeChecker
        f_in.subset_samples(sample_ids.tolist())
        with pysam.VariantFile(out_vcf, 'w', threads=num_threads, header=f_in.header,) as f_out:
            for record in f_in.fetch():
                f_out.write(record)
    if index_output_vcf:
        pysam.tabix_index(out_vcf, force=True, preset="vcf")
    print(" done", flush=True, file=sys.stderr)
    return out_vcf


def make_cross_validation_vcfs(
        vcf: Text,
        ped_file: Optional[Text] = None,
        truth_samples_file: Optional[Text] = None,
        output_folder: str = Default.output_folder,
        num_splits: int = Default.num_splits,
        index_output_vcf: bool = Default.index_output_vcf,
        num_threads: int = Default.num_threads
) -> Tuple[Tuple[str, str], ...]:
    pedigree_file_info = None if ped_file is None else pedigree_tools.PedigreeFileInfo.load(ped_file)
    if truth_samples_file is None:
        truth_samples = {}
    else:
        with open(truth_samples_file, 'r') as f_in:
            truth_samples = {line.strip() for line in f_in}
    sample_ids = numpy.array(genomics_io.get_vcf_sample_ids(vcf))
    return tuple(
        (
            subsample_vcf(vcf, train_sample_ids, out_vcf=os.path.join(output_folder, f"train_{split_ind}.vcf.gz"),
                          index_output_vcf=index_output_vcf, num_threads=num_threads),
            subsample_vcf(vcf, test_sample_ids, out_vcf=os.path.join(output_folder, f"test_{split_ind}.vcf.gz"),
                          index_output_vcf=index_output_vcf, num_threads=num_threads)
        )
        for split_ind, (train_sample_ids, test_sample_ids) in enumerate(
            _iter_sample_splits(sample_ids, pedigree_file_info=pedigree_file_info, truth_samples=truth_samples,
                                num_splits=num_splits)
        )
    )


def __parse_arguments(argv: List[Text]) -> argparse.Namespace:
    # noinspection PyTypeChecker
    parser = argparse.ArgumentParser(
        description="Subsample a VCF (by samples) for k-fold cross-validation",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        prog=argv[0]
    )
    parser.add_argument("vcf", type=str, help="Multisample VCF to be sub-sampled")
    parser.add_argument("--ped-file", "-p", type=str, required=False)
    parser.add_argument("--truth-samples-file", "-t", type=str, required=False,
                        help="Text file with one sample ID per line. These sample IDs will be evenly distributed "
                             "between folds of the cross-validation.")
    parser.add_argument("--output-folder", '-o', type=str, default=Default.output_folder,
                        help="Path to save output VCFs")
    parser.add_argument("--num_splits", '-k', type=int, default=Default.num_splits,
                        help="Number of splits for cross-validation")
    parser.add_argument("--index-output-vcf", type=bool, default=Default.index_output_vcf,
                        help="if true, create tabix index for output vcf")
    parser.add_argument("--num-threads", type=int, default=Default.num_threads,
                        help="Number of threads for reading/writing compressed VCFs")

    parsed_arguments = parser.parse_args(argv[1:] if len(argv) > 1 else ["--help"])
    if not os.path.isfile(parsed_arguments.vcf):
        raise ValueError("vcf must be a path to a valid VCF file")
    return parsed_arguments


def main(argv: Optional[List[Text]] = None) -> Tuple[Tuple[str, str], ...]:
    arguments = __parse_arguments(sys.argv if argv is None else argv)
    return make_cross_validation_vcfs(
        vcf=arguments.vcf,
        ped_file=arguments.ped_file,
        truth_samples_file=arguments.truth_samples_file,
        output_folder=arguments.output_folder,
        num_splits=arguments.num_splits,
        index_output_vcf=arguments.index_output_vcf,
        num_threads=arguments.num_threads
    )


if __name__ == "__main__":
    main()
