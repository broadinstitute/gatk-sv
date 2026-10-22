#!/usr/bin/env python3
"""
Estimate BAF (Biallelic Allele Fraction) from a single sample's SD (SiteDepth) file.

Reads a single-sample SD compressed table (.sd.txt.gz) produced by GATK's
CollectSVEvidence, and a dbSNP VCF defining REF/ALT alleles. For each site
present in both files, computes:
  baf = ALT_depth / (REF_depth + ALT_depth)
Omitting sites where the fraction is 0 or 1 (uninformative homozygous calls).

Both files are chromosome/position-sorted, so we traverse them simultaneously
(merge join) without loading the VCF into memory.

Output is written to stdout as uncompressed TSV (4 columns):
  contig  pos(0-based)  baf  sample

The caller should pipe through bgzip + tabix to produce the BGZF-compressed,
tabix-indexed output expected by PrintSVEvidence.

Usage:
  baf_from_sd.py --sd-file sample.sd.txt.gz \\
                 --sd-loci-vcf loci.vcf.gz \\
                 --sample-name SAMPLE \\
                 | bgzip > output.baf.txt.gz
"""

import gzip
import logging
import sys

import pysam

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    stream=sys.stderr,
)
logger = logging.getLogger(__name__)

# Base name lookup for depth columns (A=col3, C=col4, G=col5, T=col6)
BASES = ("A", "C", "G", "T")


def chrom_sort_key(chrom):
    """Natural sort key for chromosome names."""
    if chrom.startswith("chr"):
        chrom = chrom[3:]
    if chrom == "X":
        return 23
    if chrom == "Y":
        return 24
    if chrom in ("M", "MT"):
        return 25
    try:
        return int(chrom)
    except ValueError:
        return 99


def vcf_iterator(vcf_path):
    """
    Iterate over biallelic SNPs from a VCF using pysam.
    Yields (chrom, pos_1based, ref, alt).
    """
    vcf = pysam.VariantFile(vcf_path)
    for record in vcf:
        ref = record.ref.upper()
        alt = str(record.alts[0]).upper() if record.alts else None
        # Only biallelic SNPs
        if (
            ref and alt and len(ref) == 1 and len(alt) == 1 and record.alts
        ):
            yield (record.contig, record.start + 1, ref, alt)  # 1-based


def sd_iterator(sd_path):
    """
    Iterate over SD file entries.
    Yields (chrom, pos_0based, {A: depth, C: depth, G: depth, T: depth}).
    """
    fh = (
        gzip.open(sd_path, "rt")
        if sd_path.endswith(".gz")
        else open(sd_path, "r")
    )
    with fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 7:
                continue
            chrom = parts[0]
            pos = int(parts[1])
            depths = {
                BASES[i]: int(parts[i + 3]) for i in range(4)
            }
            yield (chrom, pos, depths)


def compute_baf_merge_join(sd_path, vcf_path, sample_name):
    """
    Merge-join the sorted SD file and sorted VCF to compute BAF.

    SD positions are 0-based (GATK convention).
    VCF positions are 1-based (VCF convention).
    We normalize to 1-based for comparison, then output 0-based for BAF.
    """
    baf_records = []
    sites_matched = 0
    sites_kept = 0
    sites_filtered = 0

    vcf_iter = vcf_iterator(vcf_path)
    sd_iter = sd_iterator(sd_path)

    vcf_entry = next(vcf_iter, None)
    sd_entry = next(sd_iter, None)

    while vcf_entry is not None and sd_entry is not None:
        vcf_chrom, vcf_pos_1 = vcf_entry[0], vcf_entry[1]  # 1-based
        sd_chrom, sd_pos_0, sd_depths = sd_entry[0], sd_entry[1], sd_entry[2]  # 0-based

        # Normalize to 1-based for comparison
        vcf_key = (chrom_sort_key(vcf_chrom), vcf_pos_1)
        sd_key = (chrom_sort_key(sd_chrom), sd_pos_0 + 1)

        if vcf_key == sd_key:
            # Match — compute BAF
            ref, alt = vcf_entry[2], vcf_entry[3]
            ref_depth = sd_depths.get(ref, 0)
            alt_depth = sd_depths.get(alt, 0)
            total = ref_depth + alt_depth

            if total > 0:
                baf = alt_depth / total
                if 0.0 < baf < 1.0:
                    baf_records.append((sd_chrom, sd_pos_0, baf, sample_name))
                    sites_kept += 1
                else:
                    sites_filtered += 1
            else:
                sites_filtered += 1
            sites_matched += 1

            vcf_entry = next(vcf_iter, None)
            sd_entry = next(sd_iter, None)

        elif vcf_key < sd_key:
            # VCF ahead, advance VCF
            vcf_entry = next(vcf_iter, None)
        else:
            # SD ahead, advance SD
            sd_entry = next(sd_iter, None)

    logger.info(
        f"Merge-join: {sites_matched} loci matched, "
        f"{sites_kept} BAF values kept, {sites_filtered} filtered"
    )

    # Records are already sorted (merge join preserves order)
    return baf_records


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Estimate BAF from a single sample's SD (SiteDepth) file."
    )
    parser.add_argument(
        "--sd-file",
        required=True,
        help="Path to the SD file (.sd.txt.gz) from CollectSVEvidence",
    )
    parser.add_argument(
        "--sd-loci-vcf",
        required=True,
        help="Path to the loci VCF defining SNP positions (vcf or vcf.gz)",
    )
    parser.add_argument(
        "--sample-name",
        required=True,
        help="Sample name for the case (used in BAF output)",
    )

    args = parser.parse_args()

    logger.info(f"Merge-joining SD: {args.sd_file} with VCF: {args.sd_loci_vcf}")
    baf_records = compute_baf_merge_join(
        args.sd_file, args.sd_loci_vcf, args.sample_name
    )

    if not baf_records:
        logger.warning("No BAF records produced (all sites were 0, 1, or no depth)")

    # Write to stdout (uncompressed TSV, caller pipes to bgzip)
    for chrom, pos, baf, sample in baf_records:
        print(f"{chrom}\t{pos}\t{baf:.6f}\t{sample}")

    logger.info(f"Wrote {len(baf_records)} BAF records to stdout")


if __name__ == "__main__":
    main()
