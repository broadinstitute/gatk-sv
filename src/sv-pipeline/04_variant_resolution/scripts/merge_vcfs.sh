#!/bin/bash
#
# merge_vcfs.sh
#
# Merge VCFs across batches for genotyping in a single batch
#
# Copies all variants from other batches into a single VCF for the base batch
#

set -e

vcflist=$1
prefix=$2

# Copy variants from each VCF from other batches, replacing genotypes
# with appropriate number of null genotypes for this batch
while read vcf; do
  gsutil cp $vcf .
done < $vcflist

while read vcf; do
  basename $vcf
done < $vcflist > vcfs.list

/opt/sv-pipeline/04_variant_resolution/scripts/merge_vcfs.py vcfs.list ${prefix}.vcf

vcf-sort -c ${prefix}.vcf | bgzip -c > ${prefix}.vcf.gz
