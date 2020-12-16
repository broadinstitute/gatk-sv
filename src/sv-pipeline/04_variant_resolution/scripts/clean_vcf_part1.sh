#!/bin/bash
#
# clean_vcf_part1.sh
#

set -euxo pipefail

# inputs
vcf=$1			# bgzipped vcf
backgroundlist=$2	# something about ugly regions
famfile=$3		# ped file must include all samples in vcf -- we only look at cols 2 and 5: sample name and sex
allosome_fai=$4		# names of allosomes in fai format -- we only look at column 1: the name of the allosome
bothsides=$5

# outputs
# includelist.txt: the names of all the samples in the input vcf
# sexchr.revise.txt: the names of the events where genotypes got tweaked on allosomes
# int.vcr.gz: a revised vcf, bgzipped

bgzip -dc@3 "$vcf" \
    | awk -v allosomeFile="$allosome_fai" -v pedFile="$famfile" -v bgdFile="$backgroundlist" -v bothFile="$bothsides" -f /opt/sv-pipeline/04_variant_resolution/scripts/clean_vcf_part1.awk \
    | bgzip -c@3 > int.vcf.gz
