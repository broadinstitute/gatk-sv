#!/usr/bin/env python3

import argparse
import gzip
import re
import sys
from collections import defaultdict


VARIANT_CLASSES = [
    "snv",
    "del_1_49bp",
    "ins_1_49bp",
    "del_gt49bp",
    "ins_gt49bp",
]

FUNCTION_CLASSES = [
    "lof",
    "missense",
    "rest_coding",
    "intronic",
    "intergenic",
    "unclassified",
]

# VEP consequence priority (high -> low).
SEVERITY_ORDER = [
    "transcript_ablation",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "stop_gained",
    "frameshift_variant",
    "stop_lost",
    "start_lost",
    "transcript_amplification",
    "feature_elongation",
    "feature_truncation",
    "inframe_insertion",
    "inframe_deletion",
    "missense_variant",
    "protein_altering_variant",
    "splice_region_variant",
    "incomplete_terminal_codon_variant",
    "start_retained_variant",
    "stop_retained_variant",
    "synonymous_variant",
    "coding_sequence_variant",
    "mature_miRNA_variant",
    "5_prime_UTR_variant",
    "3_prime_UTR_variant",
    "non_coding_transcript_exon_variant",
    "intron_variant",
    "NMD_transcript_variant",
    "non_coding_transcript_variant",
    "upstream_gene_variant",
    "downstream_gene_variant",
    "TFBS_ablation",
    "TFBS_amplification",
    "TF_binding_site_variant",
    "regulatory_region_ablation",
    "regulatory_region_amplification",
    "regulatory_region_variant",
    "intergenic_variant",
]

SEVERITY_RANK = {term: i for i, term in enumerate(SEVERITY_ORDER)}

LOF_TERMS = {
    "transcript_ablation",
    "stop_gained",
    "frameshift_variant",
    "splice_donor_variant",
    "splice_acceptor_variant",
    "stop_lost",
    "start_lost",
    "transcript_amplification",
}

MISSENSE_TERMS = {
    "missense_variant",
    "inframe_insertion",
    "inframe_deletion",
    "protein_altering_variant",
}

REST_CODING_TERMS = {
    "synonymous_variant",
    "stop_retained_variant",
    "start_retained_variant",
    "incomplete_terminal_codon_variant",
    "coding_sequence_variant",
}

INTRONIC_TERMS = {
    "intron_variant",
    "splice_region_variant",
    "splice_donor_region_variant",
    "splice_polypyrimidine_tract_variant",
    "splice_donor_5th_base_variant",
    "NMD_transcript_variant",
    "non_coding_transcript_exon_variant",
    "non_coding_transcript_variant",
}

INTERGENIC_TERMS = {
    "intergenic_variant",
    "3_prime_UTR_variant",
    "5_prime_UTR_variant",
    "upstream_gene_variant",
    "downstream_gene_variant",
    "regulatory_region_variant",
    "regulatory_region_ablation",
    "regulatory_region_amplification",
    "TF_binding_site_variant",
    "TFBS_ablation",
    "TFBS_amplification",
}


def open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def parse_info(info_str):
    info = {}
    for item in info_str.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            info[k] = v
        elif item:
            info[item] = True
    return info


def classify_variant(ref, alt):
    if not alt or alt == ".":
        return None
    if alt.startswith("<") and alt.endswith(">"):
        return None
    if "[" in alt or "]" in alt:
        return None
    if "*" in alt:
        return None

    lr = len(ref)
    la = len(alt)

    if lr == 1 and la == 1:
        return "snv"

    if lr > la:
        delta = lr - la
        if 1 <= delta <= 49:
            return "del_1_49bp"
        if delta > 49:
            return "del_gt49bp"
        return None

    if la > lr:
        delta = la - lr
        if 1 <= delta <= 49:
            return "ins_1_49bp"
        if delta > 49:
            return "ins_gt49bp"
        return None

    return None


def parse_vep_format_line(line):
    m = re.search(r"Format:\\s*([^\">]+)", line)
    if not m:
        return None
    fmt = [x.strip() for x in m.group(1).split("|")]
    return fmt if fmt else None


def build_vep_maps(vep_raw, allele_idx, consequence_idx):
    by_allele = defaultdict(list)
    all_terms = []

    if not vep_raw:
        return by_allele, all_terms

    for ann in vep_raw.split(","):
        fields = ann.split("|")
        allele = ""
        consequences = ""

        if allele_idx is not None and allele_idx < len(fields):
            allele = fields[allele_idx]

        if consequence_idx is not None and consequence_idx < len(fields):
            consequences = fields[consequence_idx]

        terms = [x for x in consequences.split("&") if x]
        if not terms:
            continue

        all_terms.extend(terms)
        if allele:
            by_allele[allele].extend(terms)

    return by_allele, all_terms


def most_severe_term(terms):
    if not terms:
        return None

    best = None
    best_rank = 10**9
    for term in terms:
        rank = SEVERITY_RANK.get(term, len(SEVERITY_ORDER) + 100)
        if rank < best_rank:
            best_rank = rank
            best = term
    return best


def functional_bucket(term):
    if term is None:
        return "unclassified"
    if term in LOF_TERMS:
        return "lof"
    if term in MISSENSE_TERMS:
        return "missense"
    if term in REST_CODING_TERMS:
        return "rest_coding"
    if term in INTRONIC_TERMS:
        return "intronic"
    if term in INTERGENIC_TERMS:
        return "intergenic"
    return "unclassified"


def write_summary(out_path, counts, totals, meta):
    with open(out_path, "wt") as out:
        out.write("# input_vcf\t%s\n" % meta["input_vcf"])
        out.write("# vep_key\t%s\n" % meta["vep_key"])
        out.write("# excluded_trv_records\t%d\n" % meta["excluded_trv_records"])
        out.write("# skipped_non_target_alleles\t%d\n" % meta["skipped_non_target_alleles"])
        out.write("# processed_variant_alleles\t%d\n" % meta["processed_variant_alleles"])
        out.write("variant_class\tlof\tmissense\trest_coding\tintronic\tintergenic\tunclassified\ttotal\n")

        for vclass in VARIANT_CLASSES:
            row_total = 0
            row = [vclass]
            for fclass in FUNCTION_CLASSES:
                n = counts[vclass][fclass]
                row.append(str(n))
                row_total += n
            row.append(str(row_total))
            out.write("\t".join(row) + "\n")

        grand = ["all_classes"]
        grand_total = 0
        for fclass in FUNCTION_CLASSES:
            n = totals[fclass]
            grand.append(str(n))
            grand_total += n
        grand.append(str(grand_total))
        out.write("\t".join(grand) + "\n")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Count SNV/INS/DEL classes and VEP functional classes from a VCF. "
            "Excludes records with 'TRV' in variant ID and uses REF/ALT lengths for type/size."
        )
    )
    parser.add_argument("--input", required=True, help="Input VCF/VCF.GZ")
    parser.add_argument("--output", required=True, help="Output TSV summary")
    parser.add_argument("--vep-key", default="vep", help="INFO key for VEP annotations (default: vep)")
    parser.add_argument(
        "--max-records",
        type=int,
        default=0,
        help="Optional limit for variant records (0 means process all)",
    )
    args = parser.parse_args()

    counts = {k: {f: 0 for f in FUNCTION_CLASSES} for k in VARIANT_CLASSES}
    totals = {f: 0 for f in FUNCTION_CLASSES}

    excluded_trv_records = 0
    skipped_non_target_alleles = 0
    processed_variant_alleles = 0
    record_seen = 0

    vep_fields = None
    allele_idx = None
    consequence_idx = None

    with open_text(args.input) as fh:
        for line in fh:
            if line.startswith("##INFO=<ID=%s," % args.vep_key):
                parsed = parse_vep_format_line(line)
                if parsed:
                    vep_fields = parsed
            if line.startswith("#CHROM"):
                break

        if vep_fields:
            lower_fields = [x.lower() for x in vep_fields]
            allele_idx = lower_fields.index("allele") if "allele" in lower_fields else None
            consequence_idx = (
                lower_fields.index("consequence") if "consequence" in lower_fields else None
            )
        else:
            allele_idx = 0
            consequence_idx = 1

        for line in fh:
            if not line or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue

            record_seen += 1
            if args.max_records and record_seen > args.max_records:
                break

            var_id = fields[2]
            if "TRV" in var_id:
                excluded_trv_records += 1
                continue

            ref = fields[3]
            alts = fields[4].split(",")
            info = parse_info(fields[7])

            by_allele, all_terms = build_vep_maps(info.get(args.vep_key, ""), allele_idx, consequence_idx)

            for alt in alts:
                vclass = classify_variant(ref, alt)
                if vclass is None:
                    skipped_non_target_alleles += 1
                    continue

                allele_terms = by_allele.get(alt, [])
                terms = allele_terms if allele_terms else all_terms
                sev = most_severe_term(terms)
                fclass = functional_bucket(sev)

                counts[vclass][fclass] += 1
                totals[fclass] += 1
                processed_variant_alleles += 1

    meta = {
        "input_vcf": args.input,
        "vep_key": args.vep_key,
        "excluded_trv_records": excluded_trv_records,
        "skipped_non_target_alleles": skipped_non_target_alleles,
        "processed_variant_alleles": processed_variant_alleles,
    }
    write_summary(args.output, counts, totals, meta)

    sys.stderr.write("Done. Wrote summary: %s\n" % args.output)


if __name__ == "__main__":
    main()
