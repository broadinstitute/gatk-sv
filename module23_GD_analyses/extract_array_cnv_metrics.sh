#!/usr/bin/env bash
# Extract per-sample copy-number-relevant metrics (GT, BAF, LRR by default)
# from a SNP-array VCF for a given genomic region.
#
# Usage:
#   extract_array_cnv_metrics.sh <in.vcf[.gz]> <region> <out.tsv> [--extra]
#
#   <in.vcf[.gz]>  Array VCF. Plain, gzipped, or bgzipped all work -- the
#                  script bgzips + indexes a working copy if needed.
#   <region>       e.g. chr1:115000-120000
#   <out.tsv>      Output path (long format: one row per variant x sample)
#   --extra        Also include R, THETA, X, Y, NORMX, NORMY, IGC
#
# Output columns (default):
#   CHROM POS ID REF ALT SAMPLE GT BAF LRR
set -euo pipefail

if [[ $# -lt 3 ]]; then
  echo "Usage: $0 <in.vcf[.gz]> <region> <out.tsv> [--extra]" >&2
  exit 1
fi

IN_VCF="$1"
REGION="$2"
OUT="$3"
EXTRA="${4:-}"

WORKDIR=$(mktemp -d)
trap 'rm -rf "$WORKDIR"' EXIT

# Detect compression and make sure we have a bgzipped, tabix-indexed copy
# to query by region (bcftools -r requires an index).
magic=$(head -c2 "$IN_VCF" | od -An -tx1 | tr -d ' ')
if [[ "$magic" == "1f8b" ]]; then
  if bgzip -t "$IN_VCF" 2>/dev/null; then
    VCF_BGZ="$IN_VCF"
  else
    echo "Input is gzip (not bgzip) -- recompressing..." >&2
    VCF_BGZ="$WORKDIR/input.vcf.gz"
    gzip -dc "$IN_VCF" | bgzip -c > "$VCF_BGZ"
  fi
else
  echo "Input is uncompressed -- bgzipping..." >&2
  VCF_BGZ="$WORKDIR/input.vcf.gz"
  bgzip -c "$IN_VCF" > "$VCF_BGZ"
fi

if [[ ! -f "${VCF_BGZ}.tbi" && ! -f "${VCF_BGZ}.csi" ]]; then
  echo "No index found -- indexing..." >&2
  if [[ "$VCF_BGZ" == "$IN_VCF" ]]; then
    # don't write next to a possibly read-only/shared input; index a symlink
    # in the workdir instead
    ln -s "$(realpath "$VCF_BGZ")" "$WORKDIR/input.vcf.gz"
    VCF_BGZ="$WORKDIR/input.vcf.gz"
  fi
  tabix -p vcf "$VCF_BGZ"
fi

if [[ -n "$EXTRA" && "$EXTRA" != "--extra" ]]; then
  echo "Unrecognized argument: $EXTRA" >&2
  exit 1
fi

if [[ "$EXTRA" == "--extra" ]]; then
  header="CHROM\tPOS\tID\tREF\tALT\tSAMPLE\tGT\tBAF\tLRR\tR\tTHETA\tX\tY\tNORMX\tNORMY\tIGC"
  fmt='[%CHROM\t%POS\t%ID\t%REF\t%ALT\t%SAMPLE\t%GT\t%BAF\t%LRR\t%R\t%THETA\t%X\t%Y\t%NORMX\t%NORMY\t%IGC\n]'
else
  header="CHROM\tPOS\tID\tREF\tALT\tSAMPLE\tGT\tBAF\tLRR"
  fmt='[%CHROM\t%POS\t%ID\t%REF\t%ALT\t%SAMPLE\t%GT\t%BAF\t%LRR\n]'
fi

printf "%b\n" "$header" > "$OUT"
bcftools query -r "$REGION" -f "$fmt" "$VCF_BGZ" >> "$OUT"

n=$(($(wc -l < "$OUT") - 1))
echo "Wrote $n (variant x sample) rows to $OUT" >&2
