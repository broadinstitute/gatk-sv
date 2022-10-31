#!/usr/bin/env python

import sys
import argparse
import gzip
import numpy
import json
import warnings

import pandas as pd

from pysam import VariantFile, VariantRecord

from sv_utils import common, genomics_io, interval_overlaps, get_truth_overlap
from typing import List, Text, Optional, Dict, Tuple, TypeVar, Union, Collection, Callable, Iterable, Set, Mapping


def get_confident_variants_vapor(vapor_files: Optional[Dict[str, str]],
                                 precision: float,
                                 valid_variant_ids: set) -> get_truth_overlap.ConfidentVariants:
    vapor_info = {
        sample_id:
            get_truth_overlap.select_confident_vapor_variants(vapor_file=vapor_file,
                                                              valid_variant_ids=valid_variant_ids,
                                                              precision=precision)
        for sample_id, vapor_file in vapor_files.items()
    }
    return vapor_info


# non veriant genotypes (taken from svtk.utils)
NULL_GT = [(0, 0), (None, None), (0, ), (None, )]


def get_called_samples(record: VariantRecord) -> set:
    samples = set()
    for sample in record.samples.keys():
        if record.samples[sample]['GT'] not in NULL_GT and \
           ('FT' not in record.samples[sample] or record.samples[sample]['FT'] == 'PASS'):
            samples.add(sample)
    return samples


def load_irs_test_report(irs_test_report):
    irs_results = genomics_io.tsv_to_pandas(data_file=irs_test_report, require_header=True, header_start='')
    irs_results.set_index("ID", inplace=True)
    return irs_results


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
    parser.add_argument("--output-summary", type=str, required=True,
                        help="File to output summary results to")
    parser.add_argument("--output-detail", type=str, required=True,
                        help="File (gzipped) to output detail results to")
    parser.add_argument("--vapor-max-cnv-size", type=int, default="5000",
                        help="Maximum size CNV to trust vapor results for")
    parser.add_argument("--irs-sample-batch-lists", type=str,
                        help="list of lists of samples used in each IRS test batch", required=True)
    parser.add_argument("--irs-contigs-file", type=str,
                        help="list of contigs to restrict IRS variants to")
    parser.add_argument("--irs-test-report-list", type=str,
                        help="list of IRS results files", required=True)
    parser.add_argument("--irs-good-pvalue-threshold", type=float, default=0.001,
                        help="Maximum pvalue to choose a good record from the IRS report")
    parser.add_argument("--irs-min-probes", type=int, default=4,
                        help="IRS results file")
    parser.add_argument("--irs-min-cnv-size", type=int, default="50000",
                        help="Minimum size CNV to trust IRS results for")

    parsed_arguments = parser.parse_args(argv[1:] if len(argv) > 1 else ["--help"])
    return parsed_arguments


def main(argv: Optional[List[Text]] = None) -> get_truth_overlap.ConfidentVariants:
    arguments = __parse_arguments(sys.argv if argv is None else argv)

    if arguments.irs_contigs_file is not None:
        irs_contigs = frozenset([line.split("\t")[0] for line in read_list_file(arguments.irs_contigs_file)])
    else:
        irs_contigs = None

    vapor_files = get_truth_overlap.get_vapor_files(arguments.vapor_json)
    vapor_variants = {sample:  genomics_io.vapor_to_pandas(vapor_files[sample]) for sample in vapor_files}
    for sample in vapor_variants:
        variants_df = vapor_variants[sample]
        p_gt_bad = 10 ** -variants_df[get_truth_overlap.Keys.vapor_gq]
        variants_df['p_non_ref'] = numpy.where(variants_df[get_truth_overlap.Keys.vapor_gt] == "0/0", p_gt_bad, 1.0 - p_gt_bad)

    sample_list_file_to_report_file_mapping = zip(read_list_file(arguments.irs_sample_batch_lists),
                                                  read_list_file(arguments.irs_test_report_list))

    samples_list_to_irs_report_data = {
        frozenset(read_list_file(sample_list)):
            load_irs_test_report(report)
        for sample_list, report in sample_list_file_to_report_file_mapping
    }

    vapor_valid_variants = 0
    irs_valid_variants = 0
    vapor_tested_var_gts = 0
    irs_tested_var_gts = 0
    vaport_supported_var_gts = 0
    irs_supported_var_gts = 0
    with VariantFile(arguments.vcf) as vcf, \
            open(arguments.output_summary, 'w') as out_summary, \
            gzip.open(arguments.output_detail, 'wb') as out_detail:
        out_summary.write("VAPOR_TESTABLE_VARIANTS\tIRS_TESTABLE_VARIANTS\tVAPOR_TESTED_VAR_GTS" +
                          "\tIRS_BATCH_TESTED_VAR_GTS\tVAPOR_SUPPORTED_GTS\tIRS_SUPPORTED_VAR_GTS\n")
        out_detail.write("SVID\tSAMPLE\t\SVTYPE\tSVLEN\tCALLED_GT\tCALLED_GQ\tVAPOR_SUPPORT_GT\tVAPOR_SUPPORT_GQ\tIRS_PVALUE\n")

        for record in vcf:
            svtype = record.info['SVTYPE']
            vapor_valid = False
            irs_valid = False
            if (svtype == 'DEL' or svtype == 'DUP') and record.info['SVLEN'] < arguments.vapor_max_cnv_size:
                vapor_valid = True
                vapor_valid_variants = vapor_valid_variants + 1
            if (svtype == 'DEL' or svtype == 'DUP') and record.info['SVLEN'] >= arguments.irs_min_cnv_size \
                    and (irs_contigs is None or record.contig in irs_contigs):
                irs_valid = True
                irs_valid_variants = irs_valid_variants + 1
            if vapor_valid or irs_valid:
                for sample in get_called_samples(record):
                    vapor_support = False
                    vapor_support_gt = 'NA'
                    vapor_support_gq = 'NA'
                    if vapor_valid and sample in vapor_variants:
                        vapor_rec = vapor_variants[sample].loc[record.id]
                        vapor_gt = vapor_rec['VaPoR_GT']
                        vapor_tested_var_gts = vapor_tested_var_gts + 1
                        vapor_support = vapor_gt in {"0/1", "1/1"} and \
                                        vapor_rec['p_non_ref'] > 1 - arguments.vapor_min_precision
                        if vapor_support:
                            vapor_support_gt = vapor_rec['VaPoR_GT']
                            vapor_support_gq = vapor_rec['VaPoR_GQ']
                            vaport_supported_var_gts = vaport_supported_var_gts + 1
                    irs_support = False
                    irs_support_pval = 'NA'
                    if irs_valid:
                        for sample_list in samples_list_to_irs_report_data:
                            if sample in sample_list:
                                irs_data = samples_list_to_irs_report_data[sample_list]
                                #print(irs_data)
                                if record.id in irs_data.index:
                                    irs_tested_var_gts = irs_tested_var_gts + 1
                                    irs_rec = irs_data.loc[record.id]
                                    #print(irs_rec)
                                    irs_support = irs_rec['NPROBES'] >= arguments.irs_min_probes and \
                                                irs_rec["PVALUE"] <= arguments.irs_good_pvalue_threshold
                                    if irs_support:
                                        irs_support_pval = irs_rec['PVALUE']
                                        irs_supported_var_gts = irs_supported_var_gts + 1

                    if vapor_support or irs_support:
                        out_detail.write("{}\t{}\t{}\t{}\t{}/{}\t{}\t{}\t{}\t{}\n".format(record.id,
                                                                                  sample,
                                                                                  record.info['SVTYPE'],
                                                                                  record.info['SVLEN'],
                                                                                  record.samples[sample]['GT'][0],
                                                                                  record.samples[sample]['GT'][1],
                                                                                  record.samples[sample]['GQ'],
                                                                                  vapor_support_gt,
                                                                                  vapor_support_gq,
                                                                                  irs_support_pval))
        out_summary.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(vapor_valid_variants,
                                                            irs_valid_variants,
                                                            vapor_tested_var_gts,
                                                            irs_tested_var_gts,
                                                            vaport_supported_var_gts,
                                                            irs_supported_var_gts))


if __name__ == "__main__":
    main()
