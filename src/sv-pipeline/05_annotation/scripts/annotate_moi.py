#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Annotate mode of inheritance (MOI) for trio de novo structural variant callsets.

For every record, compare the case sample's genotype with the provided parent
genotypes and add two INFO fields:

    MOI (String):
        DE_NOVO               case is non-ref, every provided parent is ref
        INHERITED_FROM_MOTHER     case is non-ref, mother non-ref only
        INHERITED_FROM_FATHER     case is non-ref, father non-ref only
        INHERITED_FROM_BOTH       case is non-ref, both provided parents non-ref
        PARENT_ONLY             case is ref, at least one provided parent non-ref
        UNASSESSABLE            case genotype missing/unknown
    MOI_CONFIDENCE (String):
        CONFIRMED             determination does not depend on a parent that
                                was not assayed (or whose genotype was missing)
        UNCONFIRMED           a parent was not assayed (or its genotype was
                                missing), so the call cannot be fully confirmed

The case sample must be present in the VCF. Parents are optional, but a parent
provided on the command line must be present in the VCF sample columns (the
program exits with an error otherwise).

A sample counts as non-reference at a record if its genotype carries a
non-reference allele of that record. This is a per-record allele-presence
test, not allele-index specific: on the rare records retained with multiple
ALT alleles (MULTIALLELIC), a case carrying ALT1 and a parent carrying ALT2
will be considered matching.

Outputs:
    <prefix>.vcf.gz  annotated, bgzipped VCF
    <prefix>.vcf.gz.tbi  tabix index
    <prefix>.moi_summary.tsv  per-MOI record counts
"""

import argparse
import sys

import pysam

MOI_INFO_LINE = (
    '##INFO=<ID=MOI,Number=1,Type=String,'
    'Description="Mode of inheritance relative to the case sample: '
    'DE_NOVO|INHERITED_FROM_MOTHER|INHERITED_FROM_FATHER|INHERITED_FROM_BOTH|'
    'PARENT_ONLY|UNASSESSABLE">'
)
MOI_CONFIDENCE_INFO_LINE = (
    '##INFO=<ID=MOI_CONFIDENCE,Number=1,Type=String,'
    'Description="Confidence in the MOI call: CONFIRMED when all relevant '
    'parents were assayed, UNCONFIRMED when a parent was not assayed">'
)

# sentinel states for a sample's genotype at a record
_MISSING = "missing"
_REF = "ref"
_NONREF = "nonref"

MOI_VALUES = (
    "DE_NOVO",
    "INHERITED_FROM_MOTHER",
    "INHERITED_FROM_FATHER",
    "INHERITED_FROM_BOTH",
    "PARENT_ONLY",
    "UNASSESSABLE",
)


def allele_presence(gt):
    """Classify a genotype as ref, nonref, or missing.

    ref: at least one present allele and all are 0.
    nonref: at least one present allele is non-zero.
    missing: no present (non-None) alleles at all, e.g. "./." or absent GT.
    """
    if not gt:
        return _MISSING
    present = [allele for allele in gt if allele is not None]
    if not present:
        return _MISSING
    if all(allele == 0 for allele in present):
        return _REF
    return _NONREF


def get_gt(samples_view, sample_id):
    """Return the GT array for a sample, or None if absent."""
    try:
        return samples_view[sample_id]["GT"]
    except (KeyError, ValueError, TypeError):
        return None


def classify(case, mother, father, mother_provided, father_provided):
    """Return (moi, confidence) for one record.

    mother/father are the observed states (ref/nonref/missing); *_provided
    say whether that parent was requested as part of the trio.
    """
    mother_incomplete = (not mother_provided) or (mother == _MISSING)
    father_incomplete = (not father_provided) or (father == _MISSING)
    incomplete = mother_incomplete or father_incomplete

    def confidence():
        return "CONFIRMED" if not incomplete else "UNCONFIRMED"

    if case == _MISSING:
        return "UNASSESSABLE", confidence()

    mother_nonref = mother == _NONREF
    father_nonref = father == _NONREF

    if case == _REF:
        if mother_nonref or father_nonref:
            return "PARENT_ONLY", confidence()
        # ref in every assayed sample (record retained for another reason)
        return "UNASSESSABLE", confidence()

    # case is non-ref
    if mother_nonref and father_nonref:
        return "INHERITED_FROM_BOTH", confidence()
    if mother_nonref:
        return "INHERITED_FROM_MOTHER", confidence()
    if father_nonref:
        return "INHERITED_FROM_FATHER", confidence()
    return "DE_NOVO", confidence()


def index_tabix(path):
    """Index a bgzipped VCF (same call used elsewhere in the pysam 0.15 image)."""
    pysam.tabix_index(path, preset="vcf", force=True)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf', help='Input vcf.gz (bgzipped, tabix index optional).')
    parser.add_argument('prefix', help='Output prefix; <prefix>.vcf.gz, '
                                       '<prefix>.vcf.gz.tbi and '
                                       '<prefix>.moi_summary.tsv are written.')
    parser.add_argument('--case', required=True, help='Case sample name.')
    parser.add_argument('--mother', help='Mother sample name (optional).')
    parser.add_argument('--father', help='Father sample name (optional).')

    args = parser.parse_args()

    vcf = pysam.VariantFile(args.vcf)

    if args.case not in vcf.header.samples:
        sys.exit(f"Error: case sample {args.case} not found in {args.vcf} "
                 f"(samples: {vcf.header.samples})")
    for parent in (args.mother, args.father):
        if parent is not None and parent not in vcf.header.samples:
            sys.exit(f"Error: parent sample {parent} not found in {args.vcf} "
                     f"(samples: {vcf.header.samples})")

    # add INFO lines to the header (idempotent; older pysam lacks header.info)
    for info_id, info_line in (("MOI", MOI_INFO_LINE),
                               ("MOI_CONFIDENCE", MOI_CONFIDENCE_INFO_LINE)):
        try:
            already_present = info_id in vcf.header.info
        except AttributeError:
            already_present = False
        if not already_present:
            vcf.header.add_line(info_line)

    out_path = f"{args.prefix}.vcf.gz"
    counts = {moi: 0 for moi in MOI_VALUES}

    with pysam.VariantFile(out_path, "w", header=vcf.header) as fout:
        for record in vcf:
            case = allele_presence(get_gt(record.samples, args.case))
            if args.mother:
                mother = allele_presence(get_gt(record.samples, args.mother))
            else:
                mother = _MISSING
            if args.father:
                father = allele_presence(get_gt(record.samples, args.father))
            else:
                father = _MISSING

            moi, confidence = classify(
                case, mother, father,
                mother_provided=(args.mother is not None),
                father_provided=(args.father is not None))

            record.info["MOI"] = moi
            record.info["MOI_CONFIDENCE"] = confidence
            counts[moi] += 1
            fout.write(record)

    index_tabix(out_path)

    with open(f"{args.prefix}.moi_summary.tsv", "w") as summary:
        summary.write("MOI\tCOUNT\n")
        for moi in MOI_VALUES:
            summary.write(f"{moi}\t{counts[moi]}\n")


if __name__ == '__main__':
    main()
