#!/usr/bin/env python

import sys
import argparse
import logging

import pandas as pd

from pysam import VariantFile, VariantRecord

from sv_utils import genomics_io, get_truth_overlap
from typing import List, Text, Optional, Dict, Tuple, Iterable, Set, Mapping


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
    matched_record_ids = 0
    matched_record_ids_good = 0
    matched_record_ids_bad = 0
    matched_samples_good = 0
    matched_samples_bad = 0
    total_calls = 0
    with VariantFile(vcf) as vcf:
        for record in vcf:
            if record.id in valid_irs_variant_ids:
                matched_record_ids += 1
                called_samples = get_called_samples(record)
                total_calls += len(called_samples)
                for sample_list in samples_list_to_report_mapping:
                    if record.id in samples_list_to_report_mapping[sample_list][0]:
                        matched_record_ids_good += 1
                        for sample in called_samples:
                            if sample in sample_list:
                                matched_samples_good += 1
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
                        matched_record_ids_bad += 1
                        for sample in called_samples:
                            if sample in sample_list:
                                matched_samples_bad += 1
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
    logging.info(f"Valid vcf IRS record ids: {matched_record_ids}")
    logging.info(f"Total calls in valid IRS variants: {total_calls}")
    logging.info(f"Number times a good variant was matched in an IRS batch: {matched_record_ids_good}")
    logging.info(f"Number times a bad variant was matched in an IRS batch: {matched_record_ids_bad}")
    logging.info(f"Matched good IRS variant sample calls: {matched_samples_good}")
    logging.info(f"Matched bad IRS variant sample calls: {matched_samples_bad}")
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
                        help="Max number of supporting vapor reads required for neg example")
    parser.add_argument("--vapor-read-support-neg-cov-thresh", type=int, default=5,
                        help="Max number of covering vapor reads required for neg example")
    parser.add_argument("--irs-output", type=str, default="-",
                        help="File to output IRS results to. If omitted or set to '-', print to stdout")
    parser.add_argument("--vapor-output", type=str, default="-",
                        help="File to output Vapor results to. If omitted or set to '-', print to stdout")
    parser.add_argument("--vapor-max-cnv-size", type=int, default=5000,
                        help="Maximum size CNV to trust vapor results for")
    parser.add_argument("--irs-sample-batch-lists", type=str,
                        help="list of lists of samples used in each IRS test batch")
    parser.add_argument("--irs-contigs-file", type=str,
                        help="list of contigs to restrict IRS variants to")
    parser.add_argument("--irs-test-report-list", type=str,
                        help="list of IRS results files")
    parser.add_argument("--irs-good-pvalue-threshold", type=float, default=0.000001,
                        help="Maximum pvalue to choose a good record from the IRS report")
    parser.add_argument("--irs-bad-pvalue-threshold", type=float, default=0.2,
                        help="Minimum pvalue to choose a bad record from the IRS report")
    parser.add_argument("--irs-min-probes", type=int, default=5,
                        help="IRS results file")
    parser.add_argument("--irs-min-cnv-size", type=int, default=10000,
                        help="Minimum size CNV to trust IRS results for")
    parser.add_argument("-l", "--log-level", required=False, default="INFO",
                        help="Specify level of logging information, ie. info, warning, error (not case-sensitive). "
                             "Default: INFO")

    parsed_arguments = parser.parse_args(argv[1:] if len(argv) > 1 else ["--help"])
    return parsed_arguments


def main(argv: Optional[List[Text]] = None) -> get_truth_overlap.ConfidentVariants:
    arguments = __parse_arguments(sys.argv if argv is None else argv)

    logging.basicConfig(format='[%(levelname)s:%(asctime)s] %(message)s', level=getattr(logging, arguments.log_level))

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

    if arguments.irs_sample_batch_lists is not None:
        batch_list_paths = read_list_file(arguments.irs_sample_batch_lists)
        irs_report_paths = read_list_file(arguments.irs_test_report_list)

        samples_list_to_confident_irs_variant_ids_mapping = {
            frozenset(read_list_file(sample_list)):
                get_confident_variant_ids_from_irs_report(report_list,
                                                          arguments.irs_good_pvalue_threshold,
                                                          arguments.irs_bad_pvalue_threshold,
                                                          arguments.irs_min_probes,
                                                          valid_irs_variant_ids)
            for sample_list, report_list in zip(batch_list_paths, irs_report_paths)
        }
        logging.debug("IRS dump:")
        logging.debug(str(samples_list_to_confident_irs_variant_ids_mapping))
        logging.debug("IRS sources:")
        for x, y in zip(batch_list_paths, irs_report_paths):
            logging.debug(f"{x}\t{y}")

        # for each variant in the IRS table that passes filters as good,
        # find all non ref samples and add variant ID to good list
        logging.info(f"Sample sets: {len(samples_list_to_confident_irs_variant_ids_mapping)}")
        logging.info(f"Valid IRS variant ids: {len(valid_irs_variant_ids)}")
        irs_sample_confident_variants = get_irs_sample_confident_variants(arguments.vcf,
                                                                          valid_irs_variant_ids,
                                                                          samples_list_to_confident_irs_variant_ids_mapping)
    else:
        logging.info("No IRS batches were provided.")
        irs_sample_confident_variants = {}
    logging.info(f"Samples with confident IRS variants: {len(irs_sample_confident_variants)}")

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
    logging.info(f"Samples with confident Vapor variants: {len(vapor_confident_variants)}")

    logging.info("Writing IRS labels...")
    get_truth_overlap.output_confident_variants(irs_sample_confident_variants, output_file=arguments.irs_output)
    logging.info("Writing Vapor labels...")
    get_truth_overlap.output_confident_variants(vapor_confident_variants, output_file=arguments.vapor_output)
    return irs_sample_confident_variants, vapor_confident_variants


if __name__ == "__main__":
    main()
