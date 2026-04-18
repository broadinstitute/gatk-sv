#!/usr/bin/env python3
"""
Analyze trio genotype inheritance patterns.

This version uses conservative Mendelian logic for biallelic genotypes with
possible missing alleles ('.').

Rows are removed from the summary when there is no reliable alternative
genotype in any family member. A reliable alt genotype is a genotype without '.'
that contains at least one alt allele, e.g. 0/1 or 1/1. Genotypes such as 1/.
are not treated as reliable evidence of an alt transmission state.

Categories:
    - de_novo                    : child has more alt copies than parents can supply
    - paternal_inheritance_error : father must transmit alt, child has too few alts
    - maternal_inheritance_error : mother must transmit alt, child has too few alts
    - error_from_either_side     : both parents must transmit alt, child has too few alts
    - corrected_inherited        : both child alleles known and traceable one-to-father, one-to-mother
    - missing_allele_uninformative: no Mendelian error detected but cannot fully trace child alleles
    - filtered_no_reliable_alt   : excluded from summary counts

Usage:
    python3 analyze_trio_inheritance.py <input_file> [summary_tsv] [annotated_tsv]

    summary_tsv   — per-family count table (stdout if omitted)
    annotated_tsv — input rows with an added 'inheritance_category' column
                    (stdout if omitted and summary_tsv is given; skipped to avoid
                    mixing outputs when both are omitted)

Input columns (tab-separated, with header):
    father_gt  mo_gt  child_gt  count_variants  family_ID
"""

import sys
from collections import defaultdict


def normalize_gt(gt):
    gt = gt.strip().replace('|', '/')
    if gt == '.':
        return './.'
    return gt


def parse_alleles(gt):
    gt = normalize_gt(gt)
    alleles = gt.split('/')
    if len(alleles) == 1:
        alleles = [alleles[0], alleles[0]]
    if len(alleles) != 2:
        raise ValueError(f"Unexpected genotype format: {gt}")
    return alleles


def has_reliable_alt(gt):
    alleles = parse_alleles(gt)
    return '.' not in alleles and any(allele != '0' for allele in alleles)


def has_any_missing(gt):
    return '.' in parse_alleles(gt)


def known_alleles(gt):
    """Return list of non-'.' alleles in a genotype."""
    return [a for a in parse_alleles(gt) if a != '.']


def _can_trace_both_alleles(father_gt, mo_gt, child_gt):
    """Return True iff both child alleles are known (no '.') and can be assigned
    one-to-father and one-to-mother using only known (non-'.') parental alleles."""
    child_alleles = parse_alleles(child_gt)
    if '.' in child_alleles:
        return False
    father_known = known_alleles(father_gt)
    mother_known = known_alleles(mo_gt)
    if not father_known or not mother_known:
        return False
    a, b = child_alleles
    if a in father_known and b in mother_known:
        return True
    if b in father_known and a in mother_known:
        return True
    return False


def alt_copy_range(gt):
    """Return (min_alt_copies, max_alt_copies) for a genotype with possible '.'."""
    min_alt = 0
    max_alt = 0
    for allele in parse_alleles(gt):
        if allele == '.':
            max_alt += 1
        elif allele == '0':
            pass
        else:
            min_alt += 1
            max_alt += 1
    return min_alt, max_alt


def transmitted_alt_range(gt):
    """Return (min_alt_transmitted, max_alt_transmitted) for one parent.

    Only a definite 1/1 forces alt transmission.
    """
    min_alt, max_alt = alt_copy_range(gt)
    if min_alt == 2 and max_alt == 2:
        return 1, 1
    if max_alt == 0:
        return 0, 0
    return 0, 1


def classify(father_gt, mo_gt, child_gt):
    """Classify one trio GT pattern using parental transmission constraints."""
    if not any(has_reliable_alt(gt) for gt in (father_gt, mo_gt, child_gt)):
        return 'filtered_no_reliable_alt'

    f_min, f_max = transmitted_alt_range(father_gt)
    m_min, m_max = transmitted_alt_range(mo_gt)
    c_min, c_max = alt_copy_range(child_gt)

    parent_min = f_min + m_min
    parent_max = f_max + m_max

    if c_min > parent_max:
        return 'de_novo'

    if c_max < parent_min:
        father_forced = f_min == 1
        mother_forced = m_min == 1
        if father_forced and mother_forced:
            return 'error_from_either_side'
        if father_forced:
            return 'paternal_inheritance_error'
        if mother_forced:
            return 'maternal_inheritance_error'
        return 'error_from_either_side'

    if _can_trace_both_alleles(father_gt, mo_gt, child_gt):
        return 'corrected_inherited'
    return 'missing_allele_uninformative'


SUMMARY_CATEGORIES = [
    'de_novo',
    'paternal_inheritance_error',
    'maternal_inheritance_error',
    'error_from_either_side',
    'corrected_inherited',
    'missing_allele_uninformative',
]


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    infile        = sys.argv[1]
    summary_out   = sys.argv[2] if len(sys.argv) > 2 else None
    annotated_out = sys.argv[3] if len(sys.argv) > 3 else None

    # Original summary (current logic)
    counts_orig         = defaultdict(lambda: defaultdict(int))
    filtered_orig       = defaultdict(int)

    # Strict summary: remove any row containing '.' in father/mother/child first
    counts_strict       = defaultdict(lambda: defaultdict(int))
    filtered_strict_rel = defaultdict(int)
    filtered_strict_mis = defaultdict(int)

    families       = []
    annotated_rows = []

    with open(infile) as fh:
        orig_header = fh.readline().rstrip('\n')
        header_cols = orig_header.split('\t')

        # Support both schemas:
        # 1) father_gt mo_gt child_gt count_variants family_ID
        # 2) family father mother child count
        header_map = {name: i for i, name in enumerate(header_cols)}

        def pick_idx(options):
            for name in options:
                if name in header_map:
                    return header_map[name]
            return None

        idx_father = pick_idx(['father_gt', 'father'])
        idx_mother = pick_idx(['mo_gt', 'mother'])
        idx_child = pick_idx(['child_gt', 'child'])
        idx_count = pick_idx(['count_variants', 'count'])
        idx_family = pick_idx(['family_ID', 'family'])

        if None in (idx_father, idx_mother, idx_child, idx_count, idx_family):
            raise ValueError(
                'Input header not recognized. Expected columns like '\
                'father_gt/mo_gt/child_gt/count_variants/family_ID '\
                'or family/father/mother/child/count.'
            )

        for line in fh:
            line = line.rstrip('\n')
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < len(header_cols):
                continue

            father_gt = parts[idx_father]
            mo_gt = parts[idx_mother]
            child_gt = parts[idx_child]
            n = int(parts[idx_count])
            fam = parts[idx_family]

            if fam not in families:
                families.append(fam)

            cat = classify(father_gt, mo_gt, child_gt)
            annotated_rows.append((parts, cat))

            # Original summary
            if cat == 'filtered_no_reliable_alt':
                filtered_orig[fam] += n
            else:
                counts_orig[fam][cat] += n

            # Strict summary: drop rows with any missing allele in trio first
            if has_any_missing(father_gt) or has_any_missing(mo_gt) or has_any_missing(child_gt):
                filtered_strict_mis[fam] += n
                continue

            if cat == 'filtered_no_reliable_alt':
                filtered_strict_rel[fam] += n
            else:
                counts_strict[fam][cat] += n

    # ── Annotated input file ──────────────────────────────────────────────────
    annotated_lines = [orig_header + '\tinheritance_category']
    for parts, cat in annotated_rows:
        annotated_lines.append('\t'.join(parts) + '\t' + cat)
    annotated_text = '\n'.join(annotated_lines) + '\n'

    if annotated_out:
        with open(annotated_out, 'w') as fh:
            fh.write(annotated_text)
        print(f"Annotated file written to: {annotated_out}")
    elif not summary_out:
        # No file args: print annotated rows to stdout (summary goes to stdout too — print both)
        print(annotated_text, end='')

    # ── Per-family summary counts (two summary types) ────────────────────────
    summary_lines = ['\t'.join([
        'family_ID',
        'summary_type',
        *SUMMARY_CATEGORIES,
        'total_informative',
        'filtered_no_reliable_alt',
        'filtered_any_missing',
    ])]

    for fam in families:
        # Original summary
        row_counts_orig = [counts_orig[fam].get(cat, 0) for cat in SUMMARY_CATEGORIES]
        total_orig = sum(row_counts_orig)
        summary_lines.append('\t'.join([
            fam,
            'original',
            *[str(c) for c in row_counts_orig],
            str(total_orig),
            str(filtered_orig[fam]),
            '0',
        ]))

        # Strict summary (rows with any '.' removed first)
        row_counts_strict = [counts_strict[fam].get(cat, 0) for cat in SUMMARY_CATEGORIES]
        total_strict = sum(row_counts_strict)
        summary_lines.append('\t'.join([
            fam,
            'no_missing_members',
            *[str(c) for c in row_counts_strict],
            str(total_strict),
            str(filtered_strict_rel[fam]),
            str(filtered_strict_mis[fam]),
        ]))

    summary_text = '\n'.join(summary_lines) + '\n'

    if summary_out:
        with open(summary_out, 'w') as fh:
            fh.write(summary_text)
        print(f"Summary written to: {summary_out}")
    else:
        print(summary_text, end='')


if __name__ == '__main__':
    main()
