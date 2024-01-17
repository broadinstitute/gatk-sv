#!/bin/env python

import argparse
import sys
import pysam
from typing import List, Text, Set, Optional
import logging


BOTHSIDES_SUPPORT_KEY = "BOTHSIDES_SUPPORT"
HIGH_SR_BACKGROUND_KEY = "HIGH_SR_BACKGROUND"

BOTHSIDES_SUPPORT_HEADER = f"##INFO=<ID={BOTHSIDES_SUPPORT_KEY},Number=0,Type=Flag,Description=\"Variant has " \
                           f"read-level support for both sides of breakpoint\">"
HIGH_SR_BACKGROUND_HEADER = f"##INFO=<ID={HIGH_SR_BACKGROUND_KEY},Number=0,Type=Flag,Description=\"High number of " \
                            f"SR splits in background samples indicating messy region\">"


def process(vcf, fout, bothsides_support_variants, high_sr_background_variants):

    def _log_unvisited(vid_list, n_not_found, key):
        if len(vid_list) > 0:
            logging.info(f"{len(vid_list)} {key} records were not found")
            if len(vid_list) > 100:
                logging.info("Logging only the first 50 records")
                vid_list = vid_list[:50]
            sys.stderr.write("\n".join(vid_list) + "\n")
        if len(vid_list) != n_not_found:
            logging.warning(f"Number of {key} CPX/INS variants not found {len(vid_list)} in the cleaned VCF does not "
                            f"match the number of unmatched variants with the flag already set ({n_not_found})")
        else:
            logging.info(f"Number of {key} CPX/INS variants not found in the cleaned VCF exactly matches the number of "
                         f"unmatched variants with the flag already set ({n_not_found})")

    visited = set()
    n_duplicates = 0
    n_found_bothsides = 0
    n_found_high_background = 0
    n_not_found_with_bothsides = 0
    n_not_found_with_high_background = 0
    for record in vcf:
        record_key = get_record_key(record)
        if record_key in visited:
            if n_duplicates < 50:
                logging.warning(f"Duplicate key {record_key}")
                n_duplicates += 1
                if n_duplicates == 50:
                    logging.warning("Suppressing further duplicate warnings")
        visited.add(record_key)
        if record_key in bothsides_support_variants:
            record.info[BOTHSIDES_SUPPORT_KEY] = 1
            n_found_bothsides += 1
        elif BOTHSIDES_SUPPORT_KEY in record.info:
            if record.info['SVTYPE'] in ['INS', 'CPX']:
                n_not_found_with_bothsides += 1
        if record_key in high_sr_background_variants:
            record.info[HIGH_SR_BACKGROUND_KEY] = 1
            n_found_high_background += 1
        elif HIGH_SR_BACKGROUND_KEY in record.info:
            if record.info['SVTYPE'] in ['INS', 'CPX']:
                n_not_found_with_high_background += 1
        fout.write(record)
    unvisited_bothsides = sorted(list(bothsides_support_variants - visited))
    unvisited_high_background = sorted(list(high_sr_background_variants - visited))
    logging.info(f"{n_found_bothsides} INS/CPX records newly flagged {BOTHSIDES_SUPPORT_KEY}")
    logging.info(f"{n_not_found_with_bothsides} INS/CPX records were not in the list but were already flagged "
                 f"{BOTHSIDES_SUPPORT_KEY}")
    logging.info(f"{n_found_high_background} INS/CPX records newly flagged {HIGH_SR_BACKGROUND_KEY}")
    logging.info(f"{n_not_found_with_high_background} INS/CPX records were not in the list but were already flagged "
                 f"{HIGH_SR_BACKGROUND_KEY}")
    _log_unvisited(unvisited_bothsides, n_not_found_with_bothsides, BOTHSIDES_SUPPORT_KEY)
    _log_unvisited(unvisited_high_background, n_not_found_with_high_background, HIGH_SR_BACKGROUND_KEY)


def get_record_key(record):
    return ":".join([str(x) for x in [record.chrom, record.pos,
                                      '' if record.info['SVTYPE'] == 'INS' else record.stop,
                                      record.info['SVTYPE'], record.info.get('SVLEN', ''),
                                      record.info.get('CHR2', record.chrom), record.info.get('END2', ''),
                                      record.info.get('CPX_TYPE', '')]])


def process_cpx_vcf(vcf, bothsides_support_set, high_sr_background_set):

    def _check_record(r, vid_set, variant_set):
        if r.id in vid_set:
            variant_set.add(get_record_key(r))

    bothsides_support_variants = set()
    high_sr_background_variants = set()
    for record in vcf:
        if record.info['SVTYPE'] not in ['INS', 'CPX']:
            continue
        _check_record(record, bothsides_support_set, bothsides_support_variants)
        _check_record(record, high_sr_background_set, high_sr_background_variants)
    logging.info(f"Found {len(bothsides_support_variants)} {BOTHSIDES_SUPPORT_KEY} INS/CPX variants")
    logging.info(f"Found {len(high_sr_background_variants)} {HIGH_SR_BACKGROUND_KEY} INS/CPX variants")
    return bothsides_support_variants, high_sr_background_variants


def _parse_variant_id_file(path: Text) -> Set[Text]:
    with open(path) as f:
        return {line.strip().split('\t')[-1] for line in f}


def _parse_arguments(argv: List[Text]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Apply bothsides SR support and high SR background annotations to CPX/INS records",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--cleaned-vcf', type=str, help='CleanVcf vcf (defaults to stdin)')
    parser.add_argument('--out', type=str, help='Output file (defaults to stdout)')
    parser.add_argument('--cpx-resolve-vcf', type=str, required=True, help='ResolveComplexVariants vcf')
    parser.add_argument('--bothsides-list', type=str, required=True, help='Bothsides SR pass variant list')
    parser.add_argument('--high-background-list', type=str, required=True, help='High SR background variant list')
    parser.add_argument("--log-level",
                        help="Specify level of logging information, ie. info, warning, error (not case-sensitive)",
                        required=False, default="INFO")
    if len(argv) <= 1:
        parser.parse_args(["--help"])
        sys.exit(0)
    parsed_arguments = parser.parse_args(argv[1:])
    return parsed_arguments


def main(argv: Optional[List[Text]] = None):
    if argv is None:
        argv = sys.argv
    args = _parse_arguments(argv)

    log_level = args.log_level
    numeric_level = getattr(logging, log_level.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % log_level)
    logging.basicConfig(level=numeric_level,
                        format='%(levelname)s: %(message)s')

    if args.cleaned_vcf is None:
        vcf = pysam.VariantFile(sys.stdin)
    else:
        vcf = pysam.VariantFile(args.cleaned_vcf)

    header = vcf.header
    if BOTHSIDES_SUPPORT_KEY in header.filters:
        raise ValueError(f"Found {BOTHSIDES_SUPPORT_KEY} is defined as a FILTER in the header. "
                         f"First remove all SR status filters with bcftools annotate -x "
                         f"\"FILTER/{BOTHSIDES_SUPPORT_KEY},FILTER/{HIGH_SR_BACKGROUND_KEY}\"")
    if HIGH_SR_BACKGROUND_KEY in header.filters:
        raise ValueError(f"Found {HIGH_SR_BACKGROUND_KEY} is defined as a FILTER in the header. "
                         f"First remove all SR status filters with bcftools annotate -x "
                         f"\"FILTER/{BOTHSIDES_SUPPORT_KEY},FILTER/{HIGH_SR_BACKGROUND_KEY}\"")

    header.add_line(BOTHSIDES_SUPPORT_HEADER)
    header.add_line(HIGH_SR_BACKGROUND_HEADER)
    if args.out is None:
        fout = pysam.VariantFile(sys.stdout, 'w', header=header)
    else:
        fout = pysam.VariantFile(args.out, 'w', header=header)

    bothsides_support_set = _parse_variant_id_file(args.bothsides_list)
    high_sr_background_set = _parse_variant_id_file(args.high_background_list)
    logging.info(f"Loaded {len(bothsides_support_set)} {BOTHSIDES_SUPPORT_KEY} variant ids")
    logging.info(f"Loaded {len(high_sr_background_set)} {HIGH_SR_BACKGROUND_KEY} variant ids")
    with pysam.VariantFile(args.cpx_resolve_vcf) as cpx_vcf:
        bothsides_support_variants, high_sr_background_variants = process_cpx_vcf(cpx_vcf, bothsides_support_set=bothsides_support_set, high_sr_background_set=high_sr_background_set)
    process(vcf, fout, bothsides_support_variants=bothsides_support_variants, high_sr_background_variants=high_sr_background_variants)
    vcf.close()
    fout.close()


if __name__ == "__main__":
    main()
