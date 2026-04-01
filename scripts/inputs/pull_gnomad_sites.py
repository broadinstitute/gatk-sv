#!/usr/bin/env python3

import argparse
import os
import shutil
from pathlib import Path
from urllib.parse import urlparse

import hail as hl


DEFAULT_INPUT = (
    "gs://gcp-public-data--gnomad/release/4.1.1/ht/genomes/"
    "gnomad.genomes.v4.1.1.sites.ht"
)
DEFAULT_OUTPUT = "gnomad.genomes.v4.1.1.autosomal.snp.biallelic.pass.af_gt_0.5.vcf.bgz"


def require_java() -> None:
    java_home = os.environ.get("JAVA_HOME")
    java_on_path = shutil.which("java")
    if java_home or java_on_path:
        return

    raise SystemExit(
        "Java is required to run Hail, but neither JAVA_HOME is set nor `java` was "
        "found on PATH. Install a supported Java runtime and set JAVA_HOME before "
        "running this script."
    )


def require_supported_filesystem(input_ht: str) -> None:
    scheme = urlparse(input_ht).scheme
    if scheme != "gs":
        return

    if hl.utils.hadoop_scheme_supported("gs"):
        return

    raise SystemExit(
        "This Hail/Spark install cannot read `gs://` paths because the Google Cloud "
        "Storage Hadoop connector is not available. Install the GCS connector for "
        "your Spark/Hadoop setup, then rerun the script. A common fix is to add the "
        "shaded GCS connector jar to Spark's jars/classpath on the machine running "
        "Hail."
    )


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
        help="Output VCF name written in the current working directory",
    )
    args = parser.parse_args()

    output_path = Path.cwd() / args.output
    metadata = {
        "info": {
            "AF": {
                "Description": "gnomAD adj allele frequency",
                "Number": "A",
                "Type": "Float",
            }
        }
    }

    require_java()
    hl.init()
    hl.default_reference("GRCh38")
    require_supported_filesystem(args.input_ht)

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
        filters=ht.filters,
        info=hl.struct(AF=[hl.float64(ht.freq[adj_index].AF)]),
    )

    hl.export_vcf(ht, str(output_path), metadata=metadata)
    print(f"Wrote {output_path}")

    hl.stop()


if __name__ == "__main__":
    main()
