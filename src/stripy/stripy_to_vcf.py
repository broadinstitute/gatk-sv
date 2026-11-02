#!/usr/bin/env python3

import argparse
import json
import re
from pathlib import Path

import pysam


def parse_coords(coord):
    m = re.match(r"(chr)?(?P<chrom>[^:]+):(?P<start>\d+)-(?P<end>\d+)", str(coord))
    if not m:
        return None
    chrom = m.group("chrom")
    if not chrom.startswith("chr"):
        chrom = "chr" + chrom
    return chrom, int(m.group("start")), int(m.group("end"))


def chrom_key(c):
    m = re.match(r"chr(\d+)$", c)
    if m:
        return (0, int(m.group(1)))
    order = {"chrX": (1, 23), "chrY": (1, 24), "chrM": (2, 25), "chrMT": (2, 25)}
    return order.get(c, (9, c))


def normalize_filter_label(raw_filter):
    if raw_filter is None:
        return None
    filter_label = str(raw_filter).strip()
    if not filter_label or filter_label == "." or filter_label.upper() == "PASS":
        return None
    return filter_label


def to_vcf_filter_id(filter_label):
    filter_id = re.sub(r"[^0-9A-Za-z_]+", "_", filter_label).strip("_")
    if not filter_id:
        filter_id = "STRIPY_FILTER"
    if filter_id[0].isdigit():
        filter_id = "STRIPY_" + filter_id
    return filter_id


def build_filter_id_map(loci):
    raw_filters = set()
    for locus in loci:
        filter_label = normalize_filter_label(locus["filter"])
        if filter_label is not None:
            raw_filters.add(filter_label)

    filter_id_map = {}
    used_filter_ids = set()
    for filter_label in sorted(raw_filters):
        base_filter_id = to_vcf_filter_id(filter_label)
        filter_id = base_filter_id
        suffix = 2
        while filter_id in used_filter_ids:
            filter_id = f"{base_filter_id}_{suffix}"
            suffix += 1
        filter_id_map[filter_label] = filter_id
        used_filter_ids.add(filter_id)
    return filter_id_map


def load_inputs(json_path):
    with open(json_path, "r") as f:
        J = json.load(f)
    loci = []

    def to_float(x):
        if x is None or x == "":
            return None
        try:
            return float(x)
        except Exception:
            return float("nan")

    def to_int(x):
        if x is None or x == "":
            return None
        try:
            return int(x)
        except Exception:
            try:
                return int(float(x))
            except Exception:
                return None

    def ci_tuple(a):
        if not isinstance(a, dict):
            return None
        ci = a.get("CI") or {}
        ci_min = to_int(ci.get("Min"))
        ci_max = to_int(ci.get("Max"))
        if ci_min is None or ci_max is None:
            return None
        return ci_min, ci_max

    def outlier(a):
        if not isinstance(a, dict):
            return None
        return to_int(a.get("IsPopulationOutlier"))

    for entry in J.get("GenotypingResults", []):
        for locus_name, locus in entry.items():
            tl = locus["TargetedLocus"]
            coords = tl["Coordinates"]
            parsed_coords = parse_coords(coords)
            if parsed_coords is None:
                raise ValueError(f"Unable to parse STRipy coordinates for {locus_name}: {coords!r}")
            (chrom, start, end) = parsed_coords
            motif = tl.get("Motif")
            alleles = locus.get("Alleles") or []
            a1 = alleles[0] if len(alleles) > 0 else None
            a2 = alleles[1] if len(alleles) > 1 else None

            diseases = []
            corr = tl.get("CorrespondingDisease") or {}
            for meta in corr.values():
                disease_symbol = meta.get("DiseaseSymbol")
                if disease_symbol:
                    diseases.append(disease_symbol)
            diseases_s = "|".join(sorted(set(diseases)))
            meta = locus.get("Metadata") or {}
            filt = locus.get("Filter")
            coverage = to_int(meta.get("Coverage"))
            locus_id = tl.get("LocusID") or locus_name
            loci.append({
                "chrom": chrom,
                "pos": start,
                "end": end,
                "id": str(locus_id),
                "motif": motif,
                "period": len(motif) if motif else None,
                "a1_rep": to_float(a1.get("Repeats")) if isinstance(a1, dict) else None,
                "a2_rep": to_float(a2.get("Repeats")) if isinstance(a2, dict) else None,
                "a1_ci": ci_tuple(a1),
                "a2_ci": ci_tuple(a2),
                "a1_out": outlier(a1),
                "a2_out": outlier(a2) if a2 is not None else None,
                "a1_z": to_float(a1.get("PopulationZscore")) if isinstance(a1, dict) else None,
                "a2_z": to_float(a2.get("PopulationZscore")) if isinstance(a2, dict) else None,
                "coverage": coverage,
                "filter": filt,
                "diseases": diseases_s,
            })
    loci.sort(key=lambda x: (chrom_key(x["chrom"]), x["pos"]))
    return loci


VCF_HEADER = {
    "INFO": [
        {
            "ID": "END",
            "Number": "1",
            "Type": "Integer",
            "Description": "Stop position of the interval"
        },
        {
            "ID": "SVTYPE",
            "Number": "1",
            "Type": "String",
            "Description": "Type of structural variant"
        },
        {
            "ID": "RU",
            "Number": "1",
            "Type": "String",
            "Description": "Repeat unit in the reference orientation"
        },
        {
            "ID": "PERIOD",
            "Number": "1",
            "Type": "Integer",
            "Description": "Length of the repeat unit"
        },
        {
            "ID": "DISEASES",
            "Number": ".",
            "Type": "String",
            "Description": "Associated disease symbols for this locus (| separated)"
        },
        {
            "ID": "LOCUS",
            "Number": "1",
            "Type": "String",
            "Description": "Gene/locus identifier from STRipy"
        },
    ],
    "FORMAT": [
        {
            "ID": "GT",
            "Number": "1",
            "Type": "String",
            "Description": "Unphased genotype"
        },
        {
            "ID": "REPCN",
            "Number": "2",
            "Type": "Float",
            "Description": "Number of repeat units spanned by each allele"
        },
        {
            "ID": "REPCI1",
            "Number": "2",
            "Type": "Integer",
            "Description": "95% CI min,max on repeat counts of first allele"
        },
        {
            "ID": "REPCI2",
            "Number": "2",
            "Type": "Integer",
            "Description": "95% CI min,max on repeat counts of second allele"
        },
        {
            "ID": "OUTLIER",
            "Number": "2",
            "Type": "Integer",
            "Description": "Allelic population outlier flags (0/1) assigned by STRipy"
        },
        {
            "ID": "ZSCORE",
            "Number": "2",
            "Type": "Float",
            "Description": "Allelic population Z-scores assigned by STRipy"
        },
        {
            "ID": "DP",
            "Number": "1",
            "Type": "Integer",
            "Description": "Total Depth"
        },
        {
            "ID": "STR_FILTER",
            "Number": ".",
            "Type": "String",
            "Description": "Filter status assigned by STRipy"
        },
    ]

}


def write_with_pysam(loci, out_path, sample_name, contigs=None):
    header = pysam.VariantHeader()
    filter_id_map = build_filter_id_map(loci)

    header.add_meta("source", value="STRipy2VCF")
    header.add_meta("ALT", items=[("ID", "STR"), ("Description", "Short tandem repeat")])
    for filter_label, filter_id in filter_id_map.items():
        header.add_meta(
            "FILTER",
            items=[("ID", filter_id), ("Description", f"Filter status assigned by STRipy: {filter_label}")],
        )
    for info in VCF_HEADER["INFO"]:
        items = list(info.items())
        header.add_meta("INFO", items=items)

    for format in VCF_HEADER["FORMAT"]:
        items = list(format.items())
        header.add_meta("FORMAT", items=items)

    if contigs is not None:
        for name, length in contigs:
            header.contigs.add(name, length=length)
    else:
        chrom_list = set(locus['chrom'] for locus in loci)
        for chrom_name in chrom_list:
            header.contigs.add(chrom_name)

    # This is a single sample tool
    header.add_sample(sample_name)
    vf = pysam.VariantFile(out_path, mode="w", header=header)
    for loc in loci:
        rec = header.new_record()
        rec.contig = loc["chrom"]
        rec.start = loc["pos"] - 1
        rec.stop = loc["end"]
        rec.id = str(loc["id"])
        rec.ref = "N"
        rec.alts = ("<STR>",)
        filter_label = normalize_filter_label(loc["filter"])
        rec.filter.add(filter_id_map[filter_label] if filter_label else "PASS")
        rec.info["SVTYPE"] = "STR"
        if loc["motif"]:
            rec.info["RU"] = loc["motif"]
            if loc["period"]:
                rec.info["PERIOD"] = int(loc["period"])
        rec.info["DISEASES"] = loc["diseases"]
        rec.info["LOCUS"] = loc["id"]
        sample = rec.samples[sample_name]
        sample["GT"] = (None, None)
        if loc["a1_rep"] is not None or loc["a2_rep"] is not None:
            sample["REPCN"] = (loc["a1_rep"], loc["a2_rep"])
        if loc["a1_ci"] is not None:
            sample["REPCI1"] = loc["a1_ci"]
        if loc["a2_ci"] is not None:
            sample["REPCI2"] = loc["a2_ci"]
        if loc["a1_out"] is not None or loc["a2_out"] is not None:
            sample["OUTLIER"] = (loc["a1_out"], loc["a2_out"])
        if loc["a1_z"] is not None or loc["a2_z"] is not None:
            sample["ZSCORE"] = (loc["a1_z"], loc["a2_z"])
        if loc["coverage"] is not None:
            sample["DP"] = int(loc["coverage"])
        if loc["filter"]:
            sample["STR_FILTER"] = [str(loc["filter"])]
        else:
            sample["STR_FILTER"] = ["PASS"]
        vf.write(rec)
    vf.close()


def main():
    ap = argparse.ArgumentParser(description="Convert STRipy JSON output into a VCF with SV-style STR annotations.")
    ap.add_argument("--json", required=True, help="STRipy JSON report")
    ap.add_argument("-o", "--out", required=True, help="Output VCF (compressed and tabix indexed)")
    ap.add_argument("--sample-name", default=None, help="Sample name for VCF")
    ap.add_argument("--reference", default=None, help="Optional reference fasta")
    args = ap.parse_args()

    if args.sample_name is None:
        args.sample_name = Path(args.out).stem

    loci = load_inputs(args.json)
    contigs = None
    if args.reference:
        reference = pysam.FastaFile(args.reference)
        contigs = [(name, reference.get_reference_length(name)) for name in reference.references]
        reference.close()
    write_with_pysam(loci, args.out, args.sample_name, contigs=contigs)


if __name__ == "__main__":
    main()
