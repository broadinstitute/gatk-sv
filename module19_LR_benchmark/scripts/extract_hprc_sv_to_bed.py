#!/usr/bin/env python3
"""Filter SVs by SOURCE and export a BED-like table with selected INFO/VEP annotations."""

import argparse
import gzip
import re
from collections import OrderedDict
from typing import Dict, Iterable, List, Optional, Tuple


def open_text(path: str, mode: str):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def parse_info(info_str: str) -> OrderedDict:
    info = OrderedDict()
    if not info_str:
        return info
    for item in info_str.split(";"):
        if item == "":
            continue
        if "=" in item:
            k, v = item.split("=", 1)
            info[k] = v
        else:
            info[item] = True
    return info


def get_first_value(raw: Optional[str], default: str = ".") -> str:
    if raw is None:
        return default
    token = raw.split(",", 1)[0].strip()
    return token if token else default


def parse_vep_header_format(info_header_line: str) -> List[str]:
    # Example: ... Format: Allele|Consequence|IMPACT|SYMBOL|...">
    m = re.search(r"Format:\s*([^\"]+)", info_header_line)
    if not m:
        return []
    return [x.strip() for x in m.group(1).split("|")]


def is_sv_record(alts: List[str], info: Dict[str, str]) -> bool:
    svtype = info.get("SVTYPE")
    if svtype:
        return svtype not in {"SNV", "MNV"}
    # Fallback heuristic for symbolic SV alleles.
    return any(a.startswith("<") and a.endswith(">") for a in alts)


def parse_vep_annotations(
    vep_value: str,
    field_index: Dict[str, int],
) -> Tuple[str, str]:
    if not vep_value:
        return ".", "."

    consequences = set()
    elements = set()

    idx_consequence = field_index.get("Consequence")
    idx_symbol = field_index.get("SYMBOL")
    idx_feature_type = field_index.get("Feature_type")
    idx_feature = field_index.get("Feature")

    for ann in vep_value.split(","):
        cols = ann.split("|")

        if idx_consequence is not None and idx_consequence < len(cols):
            cval = cols[idx_consequence].strip()
            if cval:
                for term in cval.split("&"):
                    term = term.strip()
                    if term:
                        consequences.add(term)

        if idx_symbol is not None and idx_symbol < len(cols):
            sval = cols[idx_symbol].strip()
            if sval and sval != "-":
                elements.add(sval)

        feature_text = ""
        if idx_feature is not None and idx_feature < len(cols):
            feature = cols[idx_feature].strip()
            if feature and feature != "-":
                if idx_feature_type is not None and idx_feature_type < len(cols):
                    ftype = cols[idx_feature_type].strip()
                    if ftype and ftype != "-":
                        feature_text = f"{ftype}:{feature}"
                    else:
                        feature_text = feature
                else:
                    feature_text = feature
        if feature_text:
            elements.add(feature_text)

    consequence_str = ",".join(sorted(consequences)) if consequences else "."
    element_str = ",".join(sorted(elements)) if elements else "."
    return consequence_str, element_str


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Filter an annotated VCF to SVs with SOURCE=HPRC_SV_Integration and "
            "export BED-like annotations."
        )
    )
    parser.add_argument("--input-vcf", required=True, help="Input annotated VCF(.gz)")
    parser.add_argument(
        "--output-vcf",
        required=True,
        help="Filtered output VCF.gz containing only selected SV records",
    )
    parser.add_argument(
        "--output-bed",
        required=True,
        help="Output BED.gz table with selected columns",
    )
    parser.add_argument(
        "--source-value",
        default="HPRC_SV_Integration",
        help="Required value for INFO/SOURCE",
    )
    parser.add_argument(
        "--vep-key",
        default="vep",
        help="INFO key containing VEP annotations (default: vep)",
    )
    args = parser.parse_args()

    vep_fields: List[str] = []
    predicted_keys: List[str] = []
    vep_header_prefix = f"##INFO=<ID={args.vep_key},"

    total = 0
    source_matched = 0
    sv_matched = 0

    with open_text(args.input_vcf, "rt") as fin, gzip.open(args.output_vcf, "wt") as fvcf, gzip.open(
        args.output_bed, "wt"
    ) as fbed:
        for line in fin:
            if line.startswith("##"):
                if line.startswith("##INFO=<ID=PREDICTED_"):
                    m = re.search(r"##INFO=<ID=([^,>]+)", line)
                    if m:
                        predicted_keys.append(m.group(1))
                if line.startswith(vep_header_prefix):
                    vep_fields = parse_vep_header_format(line)
                fvcf.write(line)
                continue

            if line.startswith("#CHROM"):
                fvcf.write(line)
                base_cols = [
                    "chrom",
                    "pos",
                    "end",
                    "id",
                    "SVTYPE",
                    "SVLEN",
                    "vep_consequences",
                    "vep_disrupted_genes_elements",
                ]
                fbed.write("\t".join(base_cols + predicted_keys) + "\n")
                continue

            total += 1
            fields = line.rstrip("\n").split("\t")
            chrom, pos, vid, ref, alt, _qual, _filter, info_raw = fields[:8]
            alts = alt.split(",")
            info = parse_info(info_raw)

            if info.get("SOURCE") != args.source_value:
                continue
            source_matched += 1

            if not is_sv_record(alts, info):
                continue
            sv_matched += 1

            fvcf.write(line)

            end = get_first_value(info.get("END"), default=pos)
            svtype = get_first_value(info.get("SVTYPE"), default=".")
            svlen = get_first_value(info.get("SVLEN"), default=".")

            field_index = {name: idx for idx, name in enumerate(vep_fields)} if vep_fields else {}
            vep_cons, vep_elements = parse_vep_annotations(
                str(info.get(args.vep_key, "")),
                field_index,
            )

            pred_values = []
            for k in predicted_keys:
                v = info.get(k)
                if v is True:
                    pred_values.append("1")
                elif v is None:
                    pred_values.append(".")
                else:
                    pred_values.append(str(v))

            row = [chrom, pos, end, vid, svtype, svlen, vep_cons, vep_elements] + pred_values
            fbed.write("\t".join(row) + "\n")

    print(
        {
            "total_records": total,
            "source_matched": source_matched,
            "sv_matched": sv_matched,
            "output_vcf": args.output_vcf,
            "output_bed": args.output_bed,
            "predicted_columns": len(predicted_keys),
        }
    )


if __name__ == "__main__":
    main()
