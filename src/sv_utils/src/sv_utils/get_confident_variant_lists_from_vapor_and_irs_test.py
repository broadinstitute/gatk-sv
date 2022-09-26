#!/usr/bin/env python

import sys
import argparse
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
            get_truth_overlap.select_confident_vapor_variants(vapor_file=vapor_file, valid_variant_ids=valid_variant_ids,
                                                                                    precision=precision)
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
    parser.add_argument("--min-vapor-precision", type=float, default=get_truth_overlap.Default.min_vapor_precision,
                        help="Minimum allowed precision for selecting good or bad variants from VaPoR")
    parser.add_argument("--output", "-O", type=str, default="-",
                        help="File to output results to. If omitted or set to '-', print to stdout")
    parser.add_argument("--vapor-max-cnv-size", type=int, default="5000",
                        help="Maximum size CNV to trust vapor results for")
    parser.add_argument("--irs-test-report", type=str,
                        help="IRS results file", required=True)
    parser.add_argument("--irs-pvalue-threshold", type=float, default=0.001,
                        help="IRS results file")
    parser.add_argument("--irs-min-probes", type=int, default=4,
                        help="IRS results file")
    parser.add_argument("--irs-min-cnv-size", type=int, default="50000",
                        help="Minimum size CNV to trust IRS results for")

    parsed_arguments = parser.parse_args(argv[1:] if len(argv) > 1 else ["--help"])
    if parsed_arguments.vapor_json is None:
        raise ValueError("Must supply one or more --vapor-json")
    return parsed_arguments


def main(argv: Optional[List[Text]] = None) -> get_truth_overlap.ConfidentVariants:
    arguments = __parse_arguments(sys.argv if argv is None else argv)

    valid_vapor_variant_ids = set()
    valid_irs_variant_ids = set()

    # scan the vcf to get a list of valid variants
    with VariantFile(arguments.vcf) as vcf:
        for record in vcf:
            svtype = record.info['SVTYPE']
            if not (svtype == 'DEL' or svtype == 'DUP') or record.info['SVLEN'] <= arguments.vapor_max_cnv_size:
                valid_vapor_variant_ids.add(record.id)
            if (svtype == 'DEL' or svtype == 'DUP') and record.info['SVLEN'] >= arguments.irs_min_cnv_size:
                valid_irs_variant_ids.add(record.id)

    vapor_files = get_truth_overlap.get_vapor_files(arguments.vapor_json)

    vapor_confident_variants = get_confident_variants_vapor(
        vapor_files=vapor_files,
        precision=arguments.min_vapor_precision,
        valid_variant_ids=valid_vapor_variant_ids
    )

    # read the IRS table
    irs_good_variant_ids = get_good_variant_ids_from_irs_report(arguments.irs_test_report,
                                                                arguments.irs_pvalue_threshold,
                                                                arguments.irs_min_probes,
                                                                arguments.irs_min_cnv_size,
                                                                valid_irs_variant_ids)

    # for each variant in the IRS table that passes filters as good,
    # find all non ref samples and add variant ID to good list
    irs_confident_variants = get_irs_sample_confident_variants(arguments.vcf,
                                                               irs_good_variant_ids,
                                                               )

    all_confident_variants = {}
    for sample in set(irs_confident_variants.keys()).union(vapor_confident_variants.keys()):
        sample_vapor_good = set(vapor_confident_variants[sample].__dict__['good_variant_ids']) if sample in vapor_confident_variants else set()
        sample_vapor_bad = set(vapor_confident_variants[sample].__dict__['bad_variant_ids']) if sample in vapor_confident_variants else set()
        sample_irs_good = set(irs_confident_variants[sample].__dict__['good_variant_ids']) if sample in irs_confident_variants else set()
        sample_irs_bad = set(irs_confident_variants[sample].__dict__['bad_variant_ids']) if sample in irs_confident_variants else set()
        all_good = sample_vapor_good.union(sample_irs_good)
        all_bad = sample_vapor_bad.union(sample_irs_bad)
        all_confident_variants[sample] = get_truth_overlap.SampleConfidentVariants(good_variant_ids=all_good,
                                                                                   bad_variant_ids=all_bad)

    get_truth_overlap.output_confident_variants(all_confident_variants, output_file=arguments.output)
    return all_confident_variants


def get_irs_sample_confident_variants(vcf: str,
                                      irs_good_variant_ids: set) -> get_truth_overlap.ConfidentVariants:
    irs_confident_variants = {}
    with VariantFile(vcf) as vcf:
        for record in vcf:
            if record.id in irs_good_variant_ids:
                samples = get_called_samples(record)
                for sample in samples:
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
                                bad_variant_ids=set())
    return irs_confident_variants


def get_good_variant_ids_from_irs_report(irs_test_report: str,
                                         irs_pvalue_threshold: float,
                                         min_probes: int,
                                         min_size: int,
                                         valid_irs_variant_ids: set):
    irs_results = genomics_io.tsv_to_pandas(data_file=irs_test_report, require_header=True, header_start='')
    irs_results.set_index("ID", inplace=True)
    irs_good_variant_ids = valid_irs_variant_ids.intersection(
        irs_results.loc[((irs_results["END"] - irs_results["START"]) >= min_size) &
                        (irs_results["NPROBES"] >= min_probes) &
                        (pd.notna(irs_results["PVALUE"])) &
                        (irs_results["PVALUE"] <= irs_pvalue_threshold)].index
    )
    return irs_good_variant_ids


if __name__ == "__main__":
    main()
