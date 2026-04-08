"""
Pull common-SNP site lists from a gnomAD Hail Table for ploidy estimation.

This is a **local-only** utility that requires `hail` (and a working Java
runtime).  Hail is intentionally *not* included in the standard package
dependencies because it is large, requires Spark/Java, and is not needed
inside the Docker images used by WDL workflows.  Install it separately
before running this subcommand::

    pip install hail

The command applies the following filters to the input Hail Table:

* Biallelic SNPs only
* PASS filter status
* Autosome / chrX non-PAR loci: adj AF > ``--af-cutoff`` (default 0.5)
* chrY non-PAR loci:            adj AF > ``--af-cutoff-y`` (default 0.01)

The result is exported as a block-gzipped VCF.
"""

import argparse
import os
import shutil
import tempfile
from pathlib import Path
from urllib.parse import urlparse

DEFAULT_INPUT_HT = (
    "gs://gcp-public-data--gnomad/release/4.1.1/ht/genomes/"
    "gnomad.genomes.v4.1.1.sites.ht"
)

DEFAULT_AF_CUTOFF = 0.5
DEFAULT_AF_CUTOFF_Y = 0.01


# ---------------------------------------------------------------------------
# Pre-flight checks
# ---------------------------------------------------------------------------

def _require_hail():
    """Import hail lazily so the rest of the package never depends on it."""
    try:
        import hail as hl  # type: ignore[import-not-found]  # noqa: F811
        return hl
    except ImportError:
        raise SystemExit(
            "The 'pull-snps' subcommand requires Hail, which is not installed.\n"
            "Hail is intentionally excluded from the standard package dependencies\n"
            "because it requires Spark/Java and is only used on local machines,\n"
            "not inside workflow Docker containers.\n\n"
            "Install it with:\n"
            "    pip install hail\n"
        )


def _require_pysam():
    """Import pysam lazily for BGZF compression and tabix indexing."""
    try:
        import pysam  # type: ignore[import-not-found]
        return pysam
    except ImportError:
        raise SystemExit(
            "The 'pull-snps' subcommand requires `pysam` to produce a real BGZF "
            "VCF and tabix index.\n"
            "Install it with:\n"
            "    pip install pysam\n"
        )


def _require_java() -> None:
    java_home = os.environ.get("JAVA_HOME")
    java_on_path = shutil.which("java")
    if java_home or java_on_path:
        return
    raise SystemExit(
        "Java is required to run Hail, but neither JAVA_HOME is set nor "
        "`java` was found on PATH.  Install a supported Java runtime and "
        "set JAVA_HOME before running this subcommand."
    )


def _require_supported_filesystem(hl, input_ht: str) -> None:
    scheme = urlparse(input_ht).scheme
    if scheme != "gs":
        return
    if hl.utils.hadoop_scheme_supported("gs"):
        return
    raise SystemExit(
        "This Hail/Spark install cannot read `gs://` paths because the "
        "Google Cloud Storage Hadoop connector is not available.  Install "
        "the GCS connector for your Spark/Hadoop setup, then rerun."
    )


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _get_adj_freq_index(hl, ht) -> int:
    """Return the index of the plain 'adj' entry in ``ht.freq_meta``."""
    freq_meta = hl.eval(ht.index_globals().freq_meta)
    for i, meta in enumerate(freq_meta):
        if meta.get("group") == "adj" and len(meta) == 1:
            return i
    raise ValueError(
        "Could not find the plain adj frequency entry in ht.freq_meta"
    )


def _bgzip_and_index_vcf(input_vcf: str, output_vcf: str) -> None:
    """Compress a plain-text VCF with BGZF and create a tabix index."""
    pysam = _require_pysam()
    pysam.tabix_compress(input_vcf, output_vcf, force=True)
    pysam.tabix_index(output_vcf, preset="vcf", force=True)


# ---------------------------------------------------------------------------
# Core logic
# ---------------------------------------------------------------------------

def pull_snps(
    hl,
    input_ht: str,
    output_path: str,
    af_cutoff: float = DEFAULT_AF_CUTOFF,
    af_cutoff_y: float = DEFAULT_AF_CUTOFF_Y,
) -> None:
    """Filter a gnomAD Hail Table and export a sites-only VCF.

    Parameters
    ----------
    hl : module
        The ``hail`` module (imported lazily by the caller).
    input_ht : str
        URI of the input Hail Table (e.g. a ``gs://`` path).
    output_path : str
        Local path for the output block-gzipped VCF.
    af_cutoff : float
        Minimum adj AF for autosomal and chrX non-PAR loci.
    af_cutoff_y : float
        Minimum adj AF for chrY non-PAR loci.
    """
    ht = hl.read_table(input_ht)
    adj_index = _get_adj_freq_index(hl, ht)
    af = ht.freq[adj_index].AF

    # Shared filters: biallelic, SNP, PASS
    common_filters = [
        hl.len(ht.alleles) == 2,
        hl.is_snp(ht.alleles[0], ht.alleles[1]),
        hl.len(ht.filters) == 0,
    ]

    # Contig-specific AF filters
    is_autosome_or_x = ht.locus.in_autosome() | ht.locus.in_x_nonpar()
    is_y = ht.locus.in_y_nonpar()

    af_filter = hl.if_else(
        is_autosome_or_x,
        af > af_cutoff,
        hl.if_else(is_y, af > af_cutoff_y, False),
    )

    filter_expr = hl.all(lambda c: c, [*common_filters, af_filter])
    ht = ht.filter(filter_expr)

    # Re-derive AF from the filtered table (Hail requires all
    # expressions in a select to come from the same table object).
    filtered_af = ht.freq[adj_index].AF

    # Keep only AF in the INFO field
    ht = ht.select(
        filters=ht.filters,
        info=hl.struct(AF=[hl.float64(filtered_af)]),
    )

    metadata = {
        "info": {
            "AF": {
                "Description": "gnomAD adj allele frequency",
                "Number": "A",
                "Type": "Float",
            }
        }
    }
    output_path = str(Path(output_path))
    with tempfile.TemporaryDirectory(prefix="pull_snps_") as tmpdir:
        temp_vcf = str(Path(tmpdir) / "sites.vcf")
        hl.export_vcf(ht, temp_vcf, metadata=metadata)
        _bgzip_and_index_vcf(temp_vcf, output_path)

    print(f"Wrote {output_path}")
    print(f"Wrote {output_path}.tbi")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Pull common-SNP site lists from a gnomAD Hail Table.\n\n"
            "NOTE: This subcommand requires 'hail' which is NOT included in\n"
            "the standard package dependencies (it needs Spark/Java and is\n"
            "only intended for local machines, not workflow Docker containers).\n"
            "Install local-only extras separately:  pip install hail pysam"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--input-ht",
        default=DEFAULT_INPUT_HT,
        help=(
            "Input Hail Table URI "
            f"(default: {DEFAULT_INPUT_HT})"
        ),
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output BGZF-compressed VCF file name; a tabix index is written automatically",
    )
    parser.add_argument(
        "--af-cutoff",
        type=float,
        default=DEFAULT_AF_CUTOFF,
        help=(
            "Minimum adj allele-frequency threshold for autosomal and "
            f"chrX non-PAR loci (default: {DEFAULT_AF_CUTOFF})"
        ),
    )
    parser.add_argument(
        "--af-cutoff-y",
        type=float,
        default=DEFAULT_AF_CUTOFF_Y,
        help=(
            "Minimum adj allele-frequency threshold for chrY non-PAR loci "
            f"(default: {DEFAULT_AF_CUTOFF_Y})"
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Entry point for the ``pull-snps`` subcommand."""
    args = _parse_args()

    # Lazy import — keeps hail out of the critical path for every other
    # subcommand and ensures a clear error when it is not installed.
    hl = _require_hail()
    _require_java()

    hl.init()
    hl.default_reference("GRCh38")
    _require_supported_filesystem(hl, args.input_ht)

    output_path = str(Path.cwd() / args.output)
    pull_snps(
        hl,
        input_ht=args.input_ht,
        output_path=output_path,
        af_cutoff=args.af_cutoff,
        af_cutoff_y=args.af_cutoff_y,
    )

    hl.stop()


if __name__ == "__main__":
    main()
