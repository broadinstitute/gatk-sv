#!/bin/env python

import argparse
from collections import Counter
import json
import gzip
import pysam
import sys
from typing import Any, List, Text, Set, Dict, Optional

import numpy as np


def load_json(path):
    with open(path) as f:
        training_sites = json.load(f)
    return {
        s: {'good_variant_ids': set(training_sites[s]['good_variant_ids']),
            'bad_variant_ids': set(training_sites[s]['bad_variant_ids'])} for s in training_sites
    }


def count_samples(vcf_path):
    with pysam.VariantFile(vcf_path) as vcf:
        return len(vcf.header.samples)


def create_vcf_tsv(out_path, truth_json_path, vcf_path, min_size_medium,
                   min_size_large, downsample_factor=None, labeled_only=False):
    # Info fields to load (double-check these if you get an invalid header error)
    info_fields = [
        'END2', 'CHR2', 'ALGORITHMS', 'EVIDENCE', 'SVLEN', 'SVTYPE', 'AF', 'STATUS', 'NON_REF_GENOTYPE_CONCORDANCE',
        'ALGORITHMS'
    ]

    # Format fields to load (double-check these if you get an invalid header error)
    format_fields = [
        'SAMPLE', 'GT', 'EV', 'GQ', 'SL', 'RD_GQ', 'SR_GQ', 'PE_GQ', 'OGQ'
    ]

    # Types to load
    valid_types = set(['DEL', 'DUP', 'INS', 'INV', 'CPX', 'CTX'])

    # Load truth labels
    training_sites = load_json(truth_json_path)

    _gt_non_ref_or_no_call_map = dict()

    def _is_non_ref_or_no_call(gt):
        s = _gt_non_ref_or_no_call_map.get(gt, None)
        if s is None:
            s = any([a is not None and a > 0 for a in gt]) or all([a is None for a in gt])
            _gt_non_ref_or_no_call_map[gt] = s
        return s

    def _reformat_field(val, key):
        if isinstance(val, tuple):
            if key == 'GT':
                return "/".join([str(e) for e in val])
            else:
                return ",".join([str(e) for e in val])
        else:
            return val

    # Parse vcf
    type_counter = Counter()
    zipped_out_path = out_path if out_path.endswith(".gz") else out_path + ".gz"
    with pysam.VariantFile(vcf_path) as vcf, gzip.open(zipped_out_path, 'w') as fout:
        # Write column names
        fout.write(bytes("\t".join(['CHROM', 'POS', 'END', 'VID', 'FILTER_CLASS', 'LABEL'] +
                                   format_fields + info_fields) + "\n", 'utf-8'))
        samples = list(set(vcf.header.samples).intersection(set(training_sites.keys())))
        i = 0
        for r in vcf:
            i += 1
            if downsample_factor is not None and i % downsample_factor > 0:
                continue
            svtype = r.info['SVTYPE']
            if svtype not in valid_types:
                continue
            svlen = r.info.get('SVLEN', None)
            if svlen is None:
                svlen = np.nan
            r_class = svtype
            if svtype in ['DEL', 'DUP']:
                r_class += '_' + ('s' if svlen < min_size_medium else 'm' if svlen < min_size_large else 'l')
            r_data = [[r.chrom, r.pos, r.stop, r.id, r_class, 1 if r.id in training_sites[s]['good_variant_ids']
                      else 0 if r.id in training_sites[s]['bad_variant_ids'] else -1]
                      + [s if k == 'SAMPLE' else _reformat_field(r.samples[s].get(k, None), k) for k in format_fields]
                      + [_reformat_field(r.info.get(k, None), k) for k in info_fields]
                      for s in samples if _is_non_ref_or_no_call(r.samples[s]['GT']) and
                      ((not labeled_only) or r.id in training_sites[s]['good_variant_ids']
                       or r.id in training_sites[s]['bad_variant_ids'])]
            type_counter[r_class] += len(r_data)
            if i % 10000 == 0:
                print(f"Processed {i} records; position {r.chrom}:{r.pos}")
                print("\t" + ", ".join(sorted([f"{key}: {val}" for key, val in type_counter.items()])))
            if len(r_data) > 0:
                fout.write(bytes("\n".join(["\t".join([str(x) for x in entry]) for entry in r_data]) + "\n", 'utf-8'))


def _parse_arguments(argv: List[Text]) -> argparse.Namespace:
    # noinspection PyTypeChecker
    parser = argparse.ArgumentParser(
        description="Writes VCF fields to tsv for consumption by optimize_sl_cutoffs.py",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--vcf', type=str, help='Input vcf with SL annotations')
    parser.add_argument('--truth-json', type=str, help='Input truth json')
    parser.add_argument('--out', type=str, help='Output table path, will be gzipped (.gz)')
    parser.add_argument("--downsample-factor", type=int,
                        help="If provided, load only every DOWNSAMPLE_FACTOR variants from the vcf")
    parser.add_argument("--labeled-only", action='store_true',
                        help="Limit output to labeled sites; intended for very large call sets with "
                             "sparsely labeled genotypes.")
    parser.add_argument("--medium-size", type=float, default=500,
                        help="Min size for medium DEL/DUP")
    parser.add_argument("--large-size", type=float, default=10000,
                        help="Min size for large DEL/DUP")
    if len(argv) <= 1:
        parser.parse_args(["--help"])
        sys.exit(0)
    parsed_arguments = parser.parse_args(argv[1:])
    return parsed_arguments


def main(argv: Optional[List[Text]] = None):
    if argv is None:
        argv = sys.argv
    args = _parse_arguments(argv)

    # Write out genotypes table to disk
    print("Reading VCF and writing to data table...")
    create_vcf_tsv(out_path=args.out, truth_json_path=args.truth_json, vcf_path=args.vcf,
                   min_size_medium=args.medium_size, min_size_large=args.large_size,
                   downsample_factor=args.downsample_factor, labeled_only=args.labeled_only)
    print("Done!")


if __name__ == "__main__":
    main()
