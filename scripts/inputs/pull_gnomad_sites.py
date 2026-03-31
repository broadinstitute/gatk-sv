#!/usr/bin/env python3

import argparse
from pathlib import Path

import hail as hl


DEFAULT_INPUT = (
    "gs://gcp-public-data--gnomad/release/4.1.1/ht/genomes/"
    "gnomad.genomes.v4.1.1.sites.ht"
)
DEFAULT_OUTPUT = "gnomad.genomes.v4.1.1.autosomal.snp.biallelic.pass.af_gt_0.5.tsv.bgz"

def get_adj_freq_index(ht: hl.Table) -> int:
    freq_meta = hl.eval(ht.index_globals().freq_meta)
    for i, meta in enumerate(freq_meta):
        if meta.get("group") == "adj" and len(meta) == 1:
            return i
    raise ValueError("Could not find the plain adj frequency entry in ht.freq_meta")


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Download a simple subset of gnomAD v4.1.1 genome sites to the current "
            "working directory."
        )
    )
    parser.add_argument("--input-ht", default=DEFAULT_INPUT, help="Input Hail Table URI")
    parser.add_argument(
        "--output",
        default=DEFAULT_OUTPUT,
        help="Output TSV name written in the current working directory",
    )
    args = parser.parse_args()

    output_path = Path.cwd() / args.output

    hl.init(default_reference="GRCh38")

    ht = hl.read_table(args.input_ht)
    adj_index = get_adj_freq_index(ht)

    filter_expr = hl.all(
        lambda condition: condition,
        [
            ht.locus.in_autosome(),
            hl.len(ht.alleles) == 2,
            hl.is_snp(ht.alleles[0], ht.alleles[1]),
            hl.len(ht.filters) == 0,
            ht.freq[adj_index].AF > 0.5,
        ],
    )
    ht = ht.filter(filter_expr)

    ht = ht.select(
        contig=ht.locus.contig,
        position=ht.locus.position,
        ref=ht.alleles[0],
        alt=ht.alleles[1],
        af=ht.freq[adj_index].AF,
    )

    ht.export(str(output_path))
    print(f"Wrote {output_path}")

    hl.stop()


if __name__ == "__main__":
    main()
