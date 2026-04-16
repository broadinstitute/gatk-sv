#!/usr/bin/env python3

import sys
import gzip

def open_vcf(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)

def parse_info(info_str):
    info = {}
    for item in info_str.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            info[k] = v
        else:
            info[item] = True
    return info

if len(sys.argv) != 3:
    sys.stderr.write(
        "Usage: parse_vcf_info.py <input.vcf[.gz]> <output.tsv>\n"
    )
    sys.exit(1)

vcf_path = sys.argv[1]
out_path = sys.argv[2]

with open_vcf(vcf_path) as vcf, open(out_path, "w") as out:
    out.write(
        "CHROM\tPOS\tREF\tALT\t"
        "allele_length\tallele_type\t"
        "gnomAD_V4_match_ID\tgnomAD_V4_match_source\n"
    )

    for line in vcf:
        if line.startswith("#"):
            continue

        fields = line.rstrip().split("\t")
        chrom, pos, _, ref, alts, _, _, info_str = fields[:8]
        info = parse_info(info_str)

        allele_length = info.get("allele_length", "")
        allele_type = info.get("allele_type", "")
        gnomad_id = info.get("gnomAD_V4_match_ID", "")
        gnomad_source = info.get("gnomAD_V4_match_source", "")

        for alt in alts.split(","):
            out.write(
                f"{chrom}\t{pos}\t{ref}\t{alt}\t"
                f"{allele_length}\t{allele_type}\t"
                f"{gnomad_id}\t{gnomad_source}\n"
            )
