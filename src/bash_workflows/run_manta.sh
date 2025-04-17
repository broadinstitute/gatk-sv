#!/bin/bash

set -Eeuo pipefail

sample_id=$1
bam_or_cram_file=$2
bam_or_cram_index=$3
reference_fasta=$4
region_bed=$5
region_bed_index=$6
num_cpu=${7:-8} # default 8 cores
jobs_per_cpu=${8:-1.3}

num_jobs=$(awk -v a="$num_cpu" -v b="$jobs_per_cpu" 'BEGIN {printf "%d", int(a * b + 0.5)}')

# prepare the analysis job
/usr/local/bin/manta/bin/configManta.py \
  --bam "$bam_or_cram_file" \
  --referenceFasta "$reference_fasta" \
  --runDir . \
  --callRegions "$region_bed"

# always tell manta there are 2 GiB per job, otherwise it will
# scale back the requested number of jobs, even if they won't
# need that much memory
./runWorkflow.py \
  --mode local \
  --jobs "$num_jobs" \
  --memGb $(( num_jobs * 2 ))

# inversion conversion, then compression and index
python2 /usr/local/bin/manta/libexec/convertInversion.py \
  /usr/local/bin/samtools \
  "$reference_fasta" \
  results/variants/diploidSV.vcf.gz \
  | bcftools reheader -s <(echo "$sample_id") \
  > diploidSV.vcf

bgzip -c diploidSV.vcf > "$sample_id.manta.vcf.gz"
tabix -p vcf "${sample_id}.manta.vcf.gz"
