#!/usr/bin/env python3
"""Extract selected VCF columns and INFO fields into a tab-delimited text table.

Output columns:
- chr, pos, ID, filter
- INFO fields: allele_type, allele_length, AC, AN, AF, VEP Consequence, VEP Gene
- All INFO fields whose key starts with PREDICTED_
"""

import argparse
import gzip
from typing import Dict, Iterable, List, Set, TextIO


def open_text(path: str, mode: str) -> TextIO:
    if path.endswith(".gz"):
        return gzip.open(path, mode)  # type: ignore[return-value]
    return open(path, mode)


def parse_info(info_str: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    if not info_str or info_str == ".":
        return out

    for item in info_str.split(";"):
        if not item:
            continue
        if "=" in item:
            k, v = item.split("=", 1)
            out[k] = v
        else:
            # Flag-style INFO fields
            out[item] = "1"
    return out


def first_present(info: Dict[str, str], keys: Iterable[str], default: str = ".") -> str:
    for k in keys:
        if k in info and info[k] != "":
            return info[k]
    return default


def parse_vep_consequence_gene(vep_value: str) -> tuple[str, str]:
    """Parse INFO/VEP value where each annotation is pipe-delimited.

    For each annotation entry:
    - consequence is the 2nd item (index 1)
    - gene is the 5th item (index 4)
    """
    if not vep_value or vep_value == ".":
        return ".", "."

    consequences: List[str] = []
    genes: List[str] = []

    for ann in vep_value.split(","):
        parts = ann.split("|")

        consequence = parts[1].strip() if len(parts) > 1 else ""
        gene = parts[4].strip() if len(parts) > 4 else ""

        if consequence and consequence not in consequences:
            consequences.append(consequence)
        if gene and gene not in genes:
            genes.append(gene)

    consequence_str = ",".join(consequences) if consequences else "."
    gene_str = ",".join(genes) if genes else "."
    return consequence_str, gene_str


def collect_predicted_keys(vcf_path: str) -> List[str]:
    keys: Set[str] = set()

    with open_text(vcf_path, "rt") as fin:
        for line in fin:
            if line.startswith("##INFO=<ID=PREDICTED_"):
                # Example: ##INFO=<ID=PREDICTED_XYZ,Number=...>
                try:
                    left = line.split("##INFO=<ID=", 1)[1]
                    key = left.split(",", 1)[0].strip()
                    if key.startswith("PREDICTED_"):
                        keys.add(key)
                except Exception:
                    pass
                continue

            if line.startswith("#CHROM"):
                break

            if line.startswith("#"):
                continue

            # Fallback scan in case predicted fields are not declared in header.
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            info = parse_info(fields[7])
            for k in info:
                if k.startswith("PREDICTED_"):
                    keys.add(k)

    return sorted(keys)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Extract selected VCF columns and INFO values into tab-delimited text, "
            "including all PREDICTED_* INFO keys."
        )
    )
    parser.add_argument("--input-vcf", required=True, help="Input VCF(.gz)")
    parser.add_argument("--output-txt", required=True, help="Output TSV/TXT path")
    args = parser.parse_args()

    predicted_keys = collect_predicted_keys(args.input_vcf)

    # Synonym lists to tolerate different key naming styles.
    allele_type_keys = ["allele_type", "ALLELE_TYPE", "SVTYPE"]
    allele_length_keys = ["allele_length", "ALLELE_LENGTH", "SVLEN"]
    ac_keys = ["AC", "ac"]
    an_keys = ["AN", "an"]
    af_keys = ["AF", "af"]

    base_header = [
        "chr",
        "pos",
        "ID",
        "filter",
        "allele_type",
        "allele_length",
        "AC",
        "AN",
        "AF",
        "VEP_Consequence",
        "VEP_Gene",
    ]
    header = base_header + predicted_keys

    n_records = 0
    with open_text(args.input_vcf, "rt") as fin, open_text(args.output_txt, "wt") as fout:
        fout.write("\t".join(header) + "\n")

        for line in fin:
            if not line or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue

            chrom = fields[0]
            pos = fields[1]
            vid = fields[2]
            flt = fields[6]
            info = parse_info(fields[7])
            vep_consequence, vep_gene = parse_vep_consequence_gene(info.get("VEP", "."))

            row = [
                chrom,
                pos,
                vid,
                flt,
                first_present(info, allele_type_keys),
                first_present(info, allele_length_keys),
                first_present(info, ac_keys),
                first_present(info, an_keys),
                first_present(info, af_keys),
                vep_consequence,
                vep_gene,
            ]

            for k in predicted_keys:
                row.append(info.get(k, "."))

            fout.write("\t".join(row) + "\n")
            n_records += 1

    print(
        {
            "input_vcf": args.input_vcf,
            "output_txt": args.output_txt,
            "records_written": n_records,
            "predicted_columns": len(predicted_keys),
        }
    )


if __name__ == "__main__":
    main()
