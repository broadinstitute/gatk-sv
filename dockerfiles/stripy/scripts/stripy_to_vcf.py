#!/usr/bin/env python3

import argparse
import json
import os
import re
import sys
from datetime import datetime
from pathlib import Path

import pysam
from pyfaidx import Fasta

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

def load_inputs(json_path):
    with open(json_path, "r") as f:
        J = json.load(f)
    loci = []
    for entry in J.get("GenotypingResults", []):
        for locus_name, locus in entry.items():
            tl = locus["TargetedLocus"]
            coords = tl["Coordinates"]
            (chrom, start, end) = parse_coords(coords)
            motif = tl["Motif"]
            (a1, a2) = locus["Alleles"]
            def to_float(x):
                try:
                    return float(x)
                except:
                    return float("nan")
            def ci_str(a):
                ci = a["CI"]
                return f"{ci.get('Min','')}-{ci.get('Max','')}"
            def outlier(a):
                return int(a["IsPopulationOutlier"])
            diseases = []
            corr = tl["CorrespondingDisease"]
            for sym, meta in corr.items():
                diseases.append(meta["DiseaseSymbol"])
            diseases_s = "|".join(sorted(set(diseases)))
            meta = locus["Metadata"]
            filt = locus["Filter"]
            coverage = meta["Coverage"]
            locus_id = tl.get("LocusID") or locus_name
            loci.append({
                "chrom": chrom,
                "pos": start,
                "end": end,
                "id": str(locus_id),
                "motif": motif,
                "period": len(motif),
                "a1_rep": to_float(a1["Repeats"]),
                "a2_rep": to_float(a2["Repeats"]),
                "a1_ci": ci_str(a1),
                "a2_ci": ci_str(a2),
                "a1_out": outlier(a1),
                "a2_out": outlier(a2),
                "a1_z": to_float(a1["PopulationZscore"]),
                "a2_z": to_float(a2["PopulationZscore"]),
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

    header.add_meta("source", value="STRipy2VCF")
    header.add_meta("ALT", items=[("ID", "STR"), ("Description", "Short tandem repeat")])
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
        rec.filter.add(loc["filter"] if loc["filter"] else "PASS")
        rec.info["SVTYPE"] = "STR"
        if loc["motif"]:
            rec.info["RU"] = loc["motif"]
            if loc["period"]:
                rec.info["PERIOD"] = int(loc["period"])
        rec.info["DISEASES"] = loc["diseases"]
        rec.info["LOCUS"] = loc["id"]
        rec.samples[sample_name]["GT"] = (0, 0)
        rec.samples[sample_name]["REPCN"] = (loc["a1_rep"], loc["a2_rep"])
        def parse_ci_tuple(ci_s):
            m = re.match(r"^\s*([+-]?[0-9]+(?:\.[0-9]+)?)\s*-\s*([+-]?[0-9]+(?:\.[0-9]+)?)\s*$", str(ci_s))
            if not m:
                return (0, 0)
            a = int(float(m.group(1)))
            b = int(float(m.group(2)))
            return (a, b)
        rec.samples[sample_name]["REPCI1"] = parse_ci_tuple(loc["a1_ci"])
        rec.samples[sample_name]["REPCI2"] = parse_ci_tuple(loc["a2_ci"])
        rec.samples[sample_name]["OUTLIER"] = (int(loc["a1_out"]), int(loc["a2_out"]))
        rec.samples[sample_name]["ZSCORE"] = (loc["a1_z"], loc["a2_z"]) 
        rec.samples[sample_name]["DP"] = int(loc["coverage"]) 
        if loc["filter"]:
            rec.samples[sample_name]["STR_FILTER"] = [str(loc["filter"])]
        else:
            rec.samples[sample_name]["STR_FILTER"] = ["PASS"]
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
        reference = Fasta(args.reference, as_raw=True, read_ahead=1000000)
        contigs = [(name, len(reference[name])) for name in reference.keys()]
        reference.close()
    write_with_pysam(loci, args.out, args.sample_name, contigs=contigs)

if __name__ == "__main__":
    main()