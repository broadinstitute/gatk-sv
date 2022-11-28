#!/bin/env python

from collections import defaultdict
import argparse
import json
import pysam
import sys
from typing import Any, List, Text, Set, Dict, Optional, Tuple


GOOD_VARIANTS_KEY = 'good_variant_ids'
BAD_VARIANTS_KEY = 'bad_variant_ids'


def _unlabeled():
    return 'unlabeled'


def parse_input_labels(json_path: Text,
                       sample_id: Text) -> Dict:
    # Load good/bad training site lists for the sample
    labels = defaultdict(_unlabeled)
    with open(json_path) as f:
        training_json = json.load(f)
        if sample_id not in training_json:
            raise ValueError(f"Could not find sample {sample_id} in training json {json_path}")
        for vid in training_json[sample_id][GOOD_VARIANTS_KEY]:
            labels[vid] = 'True'
        for vid in training_json[sample_id][BAD_VARIANTS_KEY]:
            labels[vid] = 'False'
    return labels


def get_vids_and_large_cnvs(vcf: pysam.VariantFile,
                            size_cutoff: int) -> Tuple[Set[Text], Set[Text]]:
    vids = set()
    large_cnv_vids = set()
    cnv_types = ['DEL', 'DUP']
    for record in vcf:
        vids.add(record.id)
        if record.info['SVTYPE'] in cnv_types:
            svlen = record.info['SVLEN'] if 'SVLEN' in record.info else record.stop - record.pos
            if svlen > size_cutoff:
                large_cnv_vids.add(record.id)
    return vids, large_cnv_vids


def parse_truth_support(vcf: pysam.VariantFile,
                         main_vids: Set[Text],
                         truth_algs: Set[Text],
                         min_alg_count: int,
                         max_alg_count: int) -> Dict:
    truth_dict = defaultdict(_unlabeled)
    for record in vcf:
        vids = [v for v in record.info['MEMBERS'] if v in main_vids]
        if len(vids) == 0:
            continue
        algs = record.info['ALGORITHMS']
        if isinstance(algs, str):
            algs = tuple(algs)
        algs = [a for a in algs if a in truth_algs]
        if len(algs) >= min_alg_count:
            label = 'True'
        elif len(algs) <= max_alg_count:
            label = 'False'
        else:
            continue
        for v in vids:
            truth_dict[v] = label
    return truth_dict


def refine_labels(vids: Set[Text],
                  large_cnv_vids: Set[Text],
                  labels1: Dict,
                  labels2: Dict) -> Dict:
    return {key: labels1[key] for key in vids
            if (labels1[key] == labels2[key] and (labels1[key] == 'True' or labels1[key] == 'False'))
            or key in large_cnv_vids}


def write_json(path: Text,
               sample_id: Text,
               labels: Dict) -> None:
    out = {sample_id: {GOOD_VARIANTS_KEY: [vid for vid, label in labels.items() if label == 'True'],
                       BAD_VARIANTS_KEY: [vid for vid, label in labels.items() if label == 'False']}}
    with open(path, 'w') as f:
        f.write(json.dumps(out))


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
    parser.add_argument("--clustered-vcf", type=str, required=True,
                        help="Usually this will be the cleaned vcf clustered with truth calls")
    parser.add_argument("--main-vcf", type=str, required=True,  help="Usually this will be the cleaned vcf")
    parser.add_argument("--truth-json", type=str, required=True, help="GQRecalibrator truth set json")
    parser.add_argument("--out", type=str, required=True, help="Output json path")
    parser.add_argument("--sample-id", type=str, required=True, help="Sample id")
    parser.add_argument("--min-algorithm-count", type=int, default=1,
                        help="Minimum number of truth algorithms required for call to be true")
    parser.add_argument("--max-algorithm-count", type=int, default=0,
                        help="Maximum number of truth algorithms for call to be false")
    parser.add_argument("--cnv-size-cutoff", type=int, default=10000,
                        help="Retain DEL and DUP variants in the input truth json that are above this size")
    parser.add_argument("--truth-algorithms", type=str, default="pbsv,sniffles,pav",
                        help="Comma-delimited list of truth ALGORITHMS values")

    if len(argv) <= 1:
        parser.parse_args(["--help"])
        sys.exit(0)
    parsed_arguments = parser.parse_args(argv[1:])
    if parsed_arguments.min_algorithm_count <= parsed_arguments.max_algorithm_count:
        raise ValueError(f"Min algorithm count ({parsed_arguments.min_algorithm_count} must be "
                         f"strictly greater than max algorithm count ({parsed_arguments.max_algorithm_count})")
    return parsed_arguments


def main(argv: Optional[List[Text]] = None):
    if argv is None:
        argv = sys.argv
    arguments = _parse_arguments(argv)

    # Training json labels
    original_labels = parse_input_labels(json_path=arguments.truth_json, sample_id=arguments.sample_id)

    # Get vids we're looking for
    with pysam.VariantFile(arguments.main_vcf) as vcf:
        main_vids, large_cnv_vids = get_vids_and_large_cnvs(vcf, size_cutoff=arguments.cnv_size_cutoff)

    # Parse clustered vcf and generate labels
    with pysam.VariantFile(arguments.clustered_vcf) as vcf:
        truth_algs = _parse_arg_list(arguments.truth_algorithms)
        clustering_labels = parse_truth_support(vcf=vcf, main_vids=main_vids,
                                                truth_algs=truth_algs,
                                                min_alg_count=arguments.min_algorithm_count,
                                                max_alg_count=arguments.max_algorithm_count)

    # Get consensus
    refined_labels = refine_labels(vids=main_vids, large_cnv_vids=large_cnv_vids,
                                   labels1=original_labels, labels2=clustering_labels)
    write_json(path=arguments.out, sample_id=arguments.sample_id, labels=refined_labels)


if __name__ == "__main__":
    main()
