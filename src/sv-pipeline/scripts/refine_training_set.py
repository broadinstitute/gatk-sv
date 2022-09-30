#!/bin/env python

from collections import defaultdict
import argparse
import json
import pysam
import sys
from typing import Any, List, Text, Set, Dict, Optional, Tuple

"""
Creates a GQRecalibrator training sites json file with SV calls from PacBio data.

Looks for calls clustered together in the input vcfs and comes up with consensus positive/negative labels 
for small/med dup/del/ins calls. Note large del/dup/ins and all inv from the input json are passed through 
automatically. Other SV types are not included.
"""

GOOD_VARIANTS_KEY = 'good_variant_ids'
BAD_VARIANTS_KEY = 'bad_variant_ids'
VAPOR_LABEL = 'vapor'


def parse_vapor_labels(json_path: Text,
                       sample_id: Text) -> Dict:
    # Load good/bad training site lists for the sample
    true_labels = set()
    false_labels = set()
    with open(json_path) as f:
        training_json = json.load(f)
        if sample_id not in training_json:
            raise ValueError(f"Could not find sample {sample_id} in training json {json_path}")
        for vid in training_json[sample_id][GOOD_VARIANTS_KEY]:
            true_labels.add(vid)
        for vid in training_json[sample_id][BAD_VARIANTS_KEY]:
            false_labels.add(vid)
    return true_labels, false_labels


def get_vids_and_bypass(vcf: pysam.VariantFile,
                        cnv_size_cutoff: int,
                        non_cnv_size_cutoff: int) -> Tuple[Set[Text], Set[Text]]:
    vids = set()
    large_cnv_vids = set()
    bypass_vids = set()
    cnv_types = ['DEL', 'DUP']
    for record in vcf:
        vids.add(record.id)
        svtype = record.info['SVTYPE']
        svlen = record.info['SVLEN'] if 'SVLEN' in record.info else record.stop - record.pos
        if svtype in cnv_types and svlen >= cnv_size_cutoff:
            large_cnv_vids.add(record.id)
        if svtype in ['INS', 'INV'] and svlen >= non_cnv_size_cutoff:
            bypass_vids.add(record.id)
    return vids, large_cnv_vids, bypass_vids


def parse_truth_support(labels: Dict[Text, List[Text]],
                        vcf: pysam.VariantFile,
                        main_vids: Set[Text],
                        truth_alg: Text) -> Dict:
    for record in vcf:
        if record.id in main_vids and record.info.get("STATUS", "") == "TP":
            labels[record.id].append(truth_alg)


def refine_labels(vids: Set[Text],
                  sample_id: Text,
                  large_cnv_vids: Set[Text],
                  bypass_vids: Set[Text],
                  vapor_true: Set[Text],
                  vapor_false: Set[Text],
                  strict_labels: Dict[Text, List[Text]],
                  loose_labels: Dict[Text, List[Text]],
                  strict_support_threshold: int,
                  loose_support_threshold: int) -> Dict:
    labels = {GOOD_VARIANTS_KEY: list(), BAD_VARIANTS_KEY: list()}
    for vid in vids:
        in_vapor_true = vid in vapor_true
        in_vapor_false = vid in vapor_false
        if vid in large_cnv_vids:
            # Calls where we are going to use array-based label from the vapor file
            if in_vapor_true:
                labels[GOOD_VARIANTS_KEY].append(vid)
            elif in_vapor_false:
                labels[BAD_VARIANTS_KEY].append(vid)
        elif (vid in bypass_vids) or not (in_vapor_true or in_vapor_false):
            # Vapor no-calls and unlabelable variants (large INV/INS)
            continue
        else:
            vapor_score = 1 if in_vapor_true else 0
            support_strict = len(strict_labels[vid]) + vapor_score
            support_loose = len(loose_labels[vid]) + vapor_score
            if support_strict >= strict_support_threshold:
                labels[GOOD_VARIANTS_KEY].append(vid)
            elif support_loose <= loose_support_threshold:
                labels[BAD_VARIANTS_KEY].append(vid)
    return {sample_id: labels}


def _parse_arg_list(arg: Text) -> List[Text]:
    if arg is None:
        return set()
    else:
        return arg.split(',')


def _parse_arguments(argv: List[Text]) -> argparse.Namespace:
    # noinspection PyTypeChecker
    parser = argparse.ArgumentParser(
        description="Refine truth set labels using a VCF containing clustered variants.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--loose-concordance-vcfs", nargs='+', default=[], required=True,
                        help="Space-separated list of vcfs, this will be the cleaned vcf evaluated with  "
                             "SVConcordance using loose concordance parameters against each tool")
    parser.add_argument("--strict-concordance-vcfs", nargs='+', default=[], required=True,
                        help="Space-separated list of vcfs, this will be the cleaned vcf evaluated with "
                             " SVConcordance using strict concordance parameters against each tool")
    parser.add_argument("--main-vcf", type=str, required=True,  help="Usually this will be the cleaned vcf")
    parser.add_argument("--vapor-json", type=str, required=True, help="Vapor truth set json")
    parser.add_argument("--json-out", type=str, required=True, help="Output json path")
    parser.add_argument("--table-out", type=str, required=True, help="Output tsv path")
    parser.add_argument("--sample-id", type=str, required=True, help="Sample id")
    parser.add_argument("--min-strict-algorithm-count", type=int, default=1,
                        help="Minimum number of strictly concordant truth algorithms required for call to be true")
    parser.add_argument("--max-loose-algorithm-count", type=int, default=0,
                        help="Maximum number of loosely concordant truth algorithms allowed for call to be false")
    parser.add_argument("--cnv-size-cutoff", type=int, default=5000,
                        help="Retain DEL and DUP variants in the input truth json that are above this size")
    parser.add_argument("--non-cnv-size-cutoff", type=int, default=5000,
                        help="Retain INS/INV variants in the input truth json that are above this size")
    parser.add_argument("--truth-algorithms", type=str, default="pbsv,sniffles,pav",
                        help="Comma-delimited list of truth ALGORITHMS values")

    if len(argv) <= 1:
        parser.parse_args(["--help"])
        sys.exit(0)
    parsed_arguments = parser.parse_args(argv[1:])
    if len(parsed_arguments.loose_concordance_vcfs) != len(parsed_arguments.strict_concordance_vcfs):
        raise ValueError(f"Must have same number of strict concordance vcfs as loose concordance vcfs")
    return parsed_arguments


def main(argv: Optional[List[Text]] = None):
    if argv is None:
        argv = sys.argv
    arguments = _parse_arguments(argv)

    # Training json labels
    vapor_true, vapor_false = parse_vapor_labels(json_path=arguments.vapor_json, sample_id=arguments.sample_id)

    # Get vids we're looking for
    with pysam.VariantFile(arguments.main_vcf) as vcf:
        main_vids, large_cnv_vids, bypass_vids = \
            get_vids_and_bypass(vcf, cnv_size_cutoff=arguments.cnv_size_cutoff,
                                non_cnv_size_cutoff=arguments.non_cnv_size_cutoff)

    # Parse clustered vcfs and generate labels
    truth_algs = _parse_arg_list(arguments.truth_algorithms)
    strict_labels = defaultdict(list)
    if len(arguments.loose_concordance_vcfs) != len(truth_algs):
        raise ValueError(f"Must enter same number of strict and loose concordance vcfs as truth algorithms")
    for i in range(len(arguments.strict_concordance_vcfs)):
        with pysam.VariantFile(arguments.strict_concordance_vcfs[i]) as vcf:
            parse_truth_support(labels=strict_labels, vcf=vcf, main_vids=main_vids, truth_alg=truth_algs[i])
    loose_labels = defaultdict(list)
    for i in range(len(arguments.loose_concordance_vcfs)):
        with pysam.VariantFile(arguments.loose_concordance_vcfs[i]) as vcf:
            parse_truth_support(labels=loose_labels, vcf=vcf, main_vids=main_vids, truth_alg=truth_algs[i])

    # Get consensus
    refined_labels = refine_labels(vids=main_vids, sample_id=arguments.sample_id, large_cnv_vids=large_cnv_vids,
                                   bypass_vids=bypass_vids,
                                   vapor_true=vapor_true, vapor_false=vapor_false,
                                   strict_labels=strict_labels, loose_labels=loose_labels,
                                   strict_support_threshold=arguments.min_strict_algorithm_count,
                                   loose_support_threshold=arguments.max_loose_algorithm_count)
    with open(arguments.json_out, 'w') as f:
        f.write(json.dumps(refined_labels))
    refined_label_sets = {key: set(val) for key, val in refined_labels[arguments.sample_id].items()}
    with open(arguments.table_out, 'w') as f:
        f.write("vid\tsample\tlabel\tvapor\t" + "\t".join([f"{a}_strict" for a in truth_algs])
                + "\t" + "\t".join([f"{a}_loose" for a in truth_algs]) + "\n")
        for vid in main_vids:
            label_col = "1" if vid in refined_label_sets[GOOD_VARIANTS_KEY] \
                else "0" if vid in refined_label_sets[BAD_VARIANTS_KEY] else "NA"
            vapor_col = "1" if vid in vapor_true else "0" if vid in vapor_false else "NA"
            strict_cols = "\t".join(["1" if alg in strict_labels[vid] else "0" for alg in truth_algs])
            loose_cols = "\t".join(["1" if alg in loose_labels[vid] else "0" for alg in truth_algs])
            f.write(f"{vid}\t{arguments.sample_id}\t{label_col}\t{vapor_col}\t{strict_cols}\t{loose_cols}\n")


if __name__ == "__main__":
    main()
