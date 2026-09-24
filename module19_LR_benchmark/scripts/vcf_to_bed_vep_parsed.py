#!/usr/bin/env python3
"""Convert annotated VCF(.gz) to BED, splitting VEP into Consequence/IMPACT/SYMBOL.

For each variant the most severe VEP annotation is selected
(HIGH > MODERATE > LOW > MODIFIER).  Multiple transcripts with the same top
IMPACT are collapsed: Consequence values are joined with "&", unique SYMBOL
values with ",".

Output columns (40 total, plus any --extra-info-fields appended at the end):
  #CHROM  START  END  ID  REF  ALT  QUAL  FILTER
  allele_type  allele_length  SOURCE  REGION  TRID  dbGaP_ID
  gnomAD_V4_match_type  gnomAD_V4_match_ID  gnomAD_V4_match_source
  AF  AC  AN
  PREDICTED_BREAKEND_EXONIC  PREDICTED_COPY_GAIN  PREDICTED_DUP_PARTIAL
  PREDICTED_INTERGENIC  PREDICTED_INTRAGENIC_EXON_DUP  PREDICTED_INTRONIC
  PREDICTED_INV_SPAN  PREDICTED_LOF  PREDICTED_MSV_EXON_OVERLAP
  PREDICTED_NEAREST_TSS  PREDICTED_NONCODING_BREAKPOINT  PREDICTED_NONCODING_SPAN
  PREDICTED_PARTIAL_DISPERSED_DUP  PREDICTED_PARTIAL_EXON_DUP  PREDICTED_PROMOTER
  PREDICTED_TSS_DUP  PREDICTED_UTR
  vep_Consequence  vep_IMPACT  vep_SYMBOL
  [extra INFO fields, in the order given to --extra-info-fields]

Usage:
    python3 vcf_to_bed_vep_parsed.py <input.vcf[.gz]> [output.bed]
        [--extra-info-fields FIELD1,FIELD2,...]

--extra-info-fields lets you pull additional site-level INFO fields straight
through into extra BED columns (e.g. nhomref,nhet,nhomalt), without having to
hardcode them in FIELDS. Each name must match an INFO ID already present in
the VCF header/records; missing values are written as ".".

If output is omitted, replaces .vcf.gz/.vcf with .vep_parsed.bed.
"""

import argparse
import gzip
import os

# INFO fields to extract (in order); these become columns 9-37
FIELDS = [
    "allele_type", "allele_length", "SOURCE", "REGION", "TRID", "dbGaP_ID",
    "gnomAD_V4_match_type", "gnomAD_V4_match_ID", "gnomAD_V4_match_source",
    "AF", "AC", "AN",
    "PREDICTED_BREAKEND_EXONIC", "PREDICTED_COPY_GAIN", "PREDICTED_DUP_PARTIAL",
    "PREDICTED_INTERGENIC", "PREDICTED_INTRAGENIC_EXON_DUP", "PREDICTED_INTRONIC",
    "PREDICTED_INV_SPAN", "PREDICTED_LOF", "PREDICTED_MSV_EXON_OVERLAP",
    "PREDICTED_NEAREST_TSS", "PREDICTED_NONCODING_BREAKPOINT",
    "PREDICTED_NONCODING_SPAN", "PREDICTED_PARTIAL_DISPERSED_DUP",
    "PREDICTED_PARTIAL_EXON_DUP", "PREDICTED_PROMOTER", "PREDICTED_TSS_DUP",
    "PREDICTED_UTR",
]

IMPACT_RANK = {"HIGH": 4, "MODERATE": 3, "LOW": 2, "MODIFIER": 1}


def parse_info(info_str):
    d = {}
    for token in info_str.split(";"):
        if "=" in token:
            k, v = token.split("=", 1)
            d[k] = v
        else:
            d[token] = "TRUE"
    return d


def best_vep(vep_str):
    """Return (Consequence, IMPACT, SYMBOL) for the most severe VEP annotation."""
    if not vep_str or vep_str == ".":
        return ".", ".", "."

    # Each transcript annotation is comma-separated
    transcripts = vep_str.split(",")
    best_rank = -1
    best_csq = []
    best_sym = set()
    best_impact = "."

    for t in transcripts:
        parts = t.split("|")
        if len(parts) < 4:
            continue
        csq    = parts[1] if len(parts) > 1 else "."
        impact = parts[2] if len(parts) > 2 else "."
        symbol = parts[3] if len(parts) > 3 else "."

        rank = IMPACT_RANK.get(impact, 0)
        if rank > best_rank:
            best_rank   = rank
            best_csq    = [csq]
            best_sym    = {symbol} if symbol else set()
            best_impact = impact
        elif rank == best_rank:
            if csq not in best_csq:
                best_csq.append(csq)
            if symbol:
                best_sym.add(symbol)

    if best_rank < 0:
        return ".", ".", "."

    consequence = "&".join(best_csq) if best_csq else "."
    symbol      = ",".join(sorted(best_sym)) if best_sym else "."
    return consequence, best_impact, symbol


def vcf_to_bed(infile, outfile, extra_info_fields=None):
    extra_info_fields = extra_info_fields or []
    header_cols = (
        ["#CHROM", "START", "END", "ID", "REF", "ALT", "QUAL", "FILTER"]
        + FIELDS
        + ["vep_Consequence", "vep_IMPACT", "vep_SYMBOL"]
        + extra_info_fields
    )

    opener = gzip.open if infile.endswith(".gz") else open
    with opener(infile, "rt") as fin, open(outfile, "w") as fout:
        fout.write("\t".join(header_cols) + "\n")
        for line in fin:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            chrom, pos, vid, ref, alt, qual, filt, info_str = cols[:8]
            start = int(pos) - 1  # 0-based BED start

            info = parse_info(info_str)
            end  = int(info.get("END", start + len(ref)))

            info_vals = [info.get(f, ".") for f in FIELDS]
            csq, impact, symbol = best_vep(info.get("vep", "."))
            extra_vals = [info.get(f, ".") for f in extra_info_fields]

            fout.write(
                "\t".join(
                    [chrom, str(start), str(end), vid, ref, alt, qual, filt]
                    + info_vals
                    + [csq, impact, symbol]
                    + extra_vals
                )
                + "\n"
            )

    print(f"Saved: {outfile}", flush=True)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert annotated VCF(.gz) to BED with parsed VEP columns.",
    )
    parser.add_argument("input_vcf", help="Annotated input VCF(.gz)")
    parser.add_argument("output_bed", nargs="?", default=None,
                         help="Output BED path (default: replace .vcf.gz/.vcf with .vep_parsed.bed)")
    parser.add_argument("--extra-info-fields", default="",
                         help="Comma-separated list of additional INFO field names to append as "
                              "extra BED columns, e.g. nhomref,nhet,nhomalt")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    infile = args.input_vcf
    if args.output_bed:
        outfile = args.output_bed
    else:
        base = infile
        for ext in (".vcf.gz", ".vcf.bgz", ".vcf"):
            if base.endswith(ext):
                base = base[: -len(ext)]
                break
        outfile = base + ".vep_parsed.bed"

    extra_fields = [f.strip() for f in args.extra_info_fields.split(",") if f.strip()]
    vcf_to_bed(infile, outfile, extra_info_fields=extra_fields)
