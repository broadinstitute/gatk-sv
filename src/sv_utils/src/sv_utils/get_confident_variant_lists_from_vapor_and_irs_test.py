#!/usr/bin/env python

import sys
import argparse

import pandas as pd

from pysam import VariantFile, VariantRecord

from sv_utils import common, genomics_io, interval_overlaps, get_truth_overlap
from typing import List, Text, Optional, Dict, Tuple, TypeVar, Union, Collection, Callable, Iterable, Set, Mapping


def get_confident_variants_vapor(vapor_files: Optional[Dict[str, str]],
                                 precision: float,
                                 valid_variant_ids: set,
                                 strategy: str,
                                 read_strategy_good_support_threshold: int,
                                 read_strategy_bad_support_threshold: int,
                                 read_strategy_bad_cov_threshold: int
                                 ) -> get_truth_overlap.ConfidentVariants:
    vapor_info = {
        sample_id:
            get_truth_overlap.select_confident_vapor_variants(vapor_file=vapor_file,
                                                              valid_variant_ids=valid_variant_ids,
                                                              precision=precision,
                                                              strategy=strategy,
                                                              read_strategy_good_support_threshold=read_strategy_good_support_threshold,
                                                              read_strategy_bad_support_threshold=read_strategy_bad_support_threshold,
                                                              read_strategy_bad_cov_threshold=read_strategy_bad_cov_threshold)
        for sample_id, vapor_file in vapor_files.items()
    }
    return vapor_info


# non veriant genotypes (taken from svtk.utils)
NULL_GT = [(0, 0), (None, None), (0, ), (None, )]


def get_called_samples(record: VariantRecord) -> set:
    samples = set()
    for sample in record.samples.keys():
        if record.samples[sample]['GT'] not in NULL_GT:
            samples.add(sample)
    return samples


def get_irs_sample_confident_variants(vcf: str,
                                      valid_irs_variant_ids: set,
                                      samples_list_to_report_mapping: Mapping[Set[str], Tuple[Set[str], Tuple[Set[str]]]]) \
        -> get_truth_overlap.ConfidentVariants:
    irs_confident_variants = {}
    with VariantFile(vcf) as vcf:
        for record in vcf:
            if record.id in valid_irs_variant_ids:
                called_samples = get_called_samples(record)
                for sample_list in samples_list_to_report_mapping:
                    if record.id in samples_list_to_report_mapping[sample_list][0]:
                        for sample in called_samples:
                            if sample in sample_list:
                                if sample not in irs_confident_variants:
                                    irs_confident_variants[sample] = \
                                        get_truth_overlap.SampleConfidentVariants(good_variant_ids={record.id},
                                                                                  bad_variant_ids=set())
                                else:
                                    new_good_ids = set(irs_confident_variants[sample].__dict__['good_variant_ids'])
                                    new_good_ids.add(record.id)
                                    irs_confident_variants[sample] = \
                                        get_truth_overlap.SampleConfidentVariants(
                                            good_variant_ids=new_good_ids,
                                            bad_variant_ids=set(irs_confident_variants[sample].__dict__['bad_variant_ids']))
                    if record.id in samples_list_to_report_mapping[sample_list][1]:
                        for sample in called_samples:
                            if sample in sample_list:
                                if sample not in irs_confident_variants:
                                    irs_confident_variants[sample] = \
                                        get_truth_overlap.SampleConfidentVariants(good_variant_ids=set(),
                                                                                  bad_variant_ids={record.id})
                                else:
                                    new_bad_ids = set(irs_confident_variants[sample].__dict__['bad_variant_ids'])
                                    new_bad_ids.add(record.id)
                                    irs_confident_variants[sample] = \
                                        get_truth_overlap.SampleConfidentVariants(
                                            good_variant_ids=set(irs_confident_variants[sample].__dict__['good_variant_ids']),
                                            bad_variant_ids=new_bad_ids)
    return irs_confident_variants


def get_confident_variant_ids_from_irs_report(irs_test_report: str,
                                              irs_good_pvalue_threshold: float,
                                              irs_bad_pvalue_threshold: float,
                                              min_probes: int,
                                              valid_irs_variant_ids: set) -> Tuple[Set[str], Set[str]]:
    irs_results = genomics_io.tsv_to_pandas(data_file=irs_test_report, require_header=True, header_start='')
    irs_results.set_index("ID", inplace=True)
    irs_good_variant_ids = valid_irs_variant_ids.intersection(
        irs_results.loc[(irs_results["NPROBES"] >= min_probes) &
                        (pd.notna(irs_results["PVALUE"])) &
                        (irs_results["PVALUE"] <= irs_good_pvalue_threshold)].index
    )
    irs_bad_variant_ids = valid_irs_variant_ids.intersection(
        irs_results.loc[(irs_results["NPROBES"] >= min_probes) &
                        (pd.notna(irs_results["PVALUE"])) &
                        (irs_results["PVALUE"] >= irs_bad_pvalue_threshold)].index
    )

    return irs_good_variant_ids, irs_bad_variant_ids


def read_list_file(path: str) -> Iterable[str]:
    with open(path, 'r') as f:
        items = [line for line in f.read().splitlines() if line]
    if len(items) == 0:
        raise ValueError("list empty")
    return items


def __parse_arguments(argv: List[Text]) -> argparse.Namespace:
    # noinspection PyTypeChecker
    parser = argparse.ArgumentParser(
        description="Get lists of confident variant calls per sample given VaPoR and Array IRS Test inputs",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        prog=argv[0]
    )
    parser.add_argument("--vcf", "-v", type=str,
                        help="GATK-SV VCF file.", required=True)
    parser.add_argument("--vapor-json", "-j", type=str,
                        help="Json file with mapping from sample ID to corresponding VaPoR file.", required=True)
    parser.add_argument("--vapor-min-precision", type=float, default=get_truth_overlap.Default.min_vapor_precision,
                        help="Minimum allowed precision for selecting good or bad variants from VaPoR")
    parser.add_argument("--vapor-strategy", type=str, default="READS",
                        help="GQ (just based on VaPoR_GQ), GT (based on re-genotyped model), or READS (based on read support threshold)")
    parser.add_argument("--vapor-read-support-pos-thresh", type=int, default=2,
                        help="Min Number of supporting vapor reads required for positive example")
    parser.add_argument("--vapor-read-support-neg-thresh", type=int, default=0,
                        help="MaNumber of supporting vapor reads required for neg example")
    parser.add_argument("--vapor-read-support-neg-cov-thresh", type=int, default=5,
                        help="MaNumber of covering vapor reads required for neg example")
    parser.add_argument("--output", "-O", type=str, default="-",
                        help="File to output results to. If omitted or set to '-', print to stdout")
    parser.add_argument("--vapor-max-cnv-size", type=int, default="5000",
                        help="Maximum size CNV to trust vapor results for")
    parser.add_argument("--irs-sample-batch-lists", type=str,
                        help="list of lists of samples used in each IRS test batch")
    parser.add_argument("--irs-contigs-file", type=str,
                        help="list of contigs to restrict IRS variants to")
    parser.add_argument("--irs-test-report-list", type=str,
                        help="list of IRS results files")
    parser.add_argument("--irs-good-pvalue-threshold", type=float, default=0.001,
                        help="Maximum pvalue to choose a good record from the IRS report")
    parser.add_argument("--irs-bad-pvalue-threshold", type=float, default=0.5,
                        help="Minimum pvalue to choose a bad record from the IRS report")
    parser.add_argument("--irs-min-probes", type=int, default=4,
                        help="IRS results file")
    parser.add_argument("--irs-min-cnv-size", type=int, default="50000",
                        help="Minimum size CNV to trust IRS results for")

    parsed_arguments = parser.parse_args(argv[1:] if len(argv) > 1 else ["--help"])
    return parsed_arguments


def main(argv: Optional[List[Text]] = None) -> get_truth_overlap.ConfidentVariants:
    arguments = __parse_arguments(sys.argv if argv is None else argv)

    valid_vapor_variant_ids = set()
    valid_irs_variant_ids = set()

    if arguments.irs_contigs_file is not None:
        irs_contigs = frozenset([line.split("\t")[0] for line in read_list_file(arguments.irs_contigs_file)])
    else:
        irs_contigs = None

    # scan the vcf to get a list of valid variants
    with VariantFile(arguments.vcf) as vcf:
        for record in vcf:
            svtype = record.info['SVTYPE']
            if not (svtype == 'DEL' or svtype == 'DUP') or record.info['SVLEN'] <= arguments.vapor_max_cnv_size:
                valid_vapor_variant_ids.add(record.id)
            if (svtype == 'DEL' or svtype == 'DUP') and record.info['SVLEN'] >= arguments.irs_min_cnv_size \
                    and (irs_contigs is None or record.contig in irs_contigs):
                valid_irs_variant_ids.add(record.id)

    vapor_files = get_truth_overlap.get_vapor_files(arguments.vapor_json)

    vapor_confident_variants = get_confident_variants_vapor(
        vapor_files=vapor_files,
        precision=arguments.vapor_min_precision,
        valid_variant_ids=valid_vapor_variant_ids,
        strategy=arguments.vapor_strategy,
        read_strategy_good_support_threshold=arguments.vapor_read_support_pos_thresh,
        read_strategy_bad_support_threshold=arguments.vapor_read_support_neg_thresh,
        read_strategy_bad_cov_threshold=arguments.vapor_read_support_neg_cov_thresh
    )

    if arguments.irs_sample_batch_lists is not None:
        sample_list_file_to_report_file_mapping = zip(read_list_file(arguments.irs_sample_batch_lists),
                                                      read_list_file(arguments.irs_test_report_list))

        samples_list_to_confident_irs_variant_ids_mapping = {
            frozenset(read_list_file(sample_list)):
                get_confident_variant_ids_from_irs_report(report_list,
                                                          arguments.irs_good_pvalue_threshold,
                                                          arguments.irs_bad_pvalue_threshold,
                                                          arguments.irs_min_probes,
                                                          valid_irs_variant_ids)
            for sample_list, report_list in sample_list_file_to_report_file_mapping
        }

        # for each variant in the IRS table that passes filters as good,
        # find all non ref samples and add variant ID to good list
        irs_sample_confident_variants = get_irs_sample_confident_variants(arguments.vcf,
                                                                          valid_irs_variant_ids,
                                                                          samples_list_to_confident_irs_variant_ids_mapping)
    else:
        irs_sample_confident_variants = {}

    all_confident_variants = {}
    for sample in set(irs_sample_confident_variants.keys()).union(vapor_confident_variants.keys()):
        sample_vapor_good = set(vapor_confident_variants[sample].__dict__['good_variant_ids']) \
            if sample in vapor_confident_variants else set()
        sample_vapor_bad = set(vapor_confident_variants[sample].__dict__['bad_variant_ids']) \
            if sample in vapor_confident_variants else set()
        sample_irs_good = set(irs_sample_confident_variants[sample].__dict__['good_variant_ids']) \
            if sample in irs_sample_confident_variants else set()
        sample_irs_bad = set(irs_sample_confident_variants[sample].__dict__['bad_variant_ids']) \
            if sample in irs_sample_confident_variants else set()
        all_good = sample_vapor_good.union(sample_irs_good)
        all_bad = sample_vapor_bad.union(sample_irs_bad)
        all_confident_variants[sample] = get_truth_overlap.SampleConfidentVariants(good_variant_ids=all_good,
                                                                                   bad_variant_ids=all_bad)

    get_truth_overlap.output_confident_variants(all_confident_variants, output_file=arguments.output)
    return all_confident_variants


if __name__ == "__main__":
    main()
