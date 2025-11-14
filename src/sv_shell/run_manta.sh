#!/bin/bash

# This script is a bash implementation of the following workflow/task in WDL:
# Filename: wdl/Manta.wdl
# Workflow: Manta

set -Exeuo pipefail

sample_id=$1
bam_or_cram_file=$2
bam_or_cram_index=$3
reference_fasta=$4
region_bed=$5
region_bed_index=$6
outputs_json_filename=$7
num_cpu=${8:-8} # default 8 cores
jobs_per_cpu=${9:-1.3}

# You may automatically adjust the job count and memory size as the following.
# We currently use the hardcoded values for simplicity.
#num_jobs=$(awk -v a="$num_cpu" -v b="$jobs_per_cpu" 'BEGIN {printf "%d", int(a * b + 0.5)}')
#mem_size=$(( num_jobs * 2 ))
num_jobs=8
mem_size=16

echo "=============== Running manta"

# you need to define a separate directory for each manta run
# (e.g., when multiple runs invoked in a single docker image).
# Otherwise, you will get the following error:
# > configManta.py: error: Run directory already contains workflow script file '/runWorkflow.py'.
# > Each analysis must be configured in a separate directory.
TMPDIR=`mktemp -d -p .` || exit 1

working_dir=$(mktemp -d /wd_manta_XXXXXXXX)
working_dir="$(realpath ${working_dir})"
output_dir=$(mktemp -d /output_manta_XXXXXXXX)
output_dir="$(realpath ${output_dir})"
cd "${working_dir}"

# prepare the analysis job
/usr/local/bin/manta/bin/configManta.py \
  --bam "$bam_or_cram_file" \
  --referenceFasta "$reference_fasta" \
  --runDir ${working_dir} \
  --callRegions "$region_bed"

# always tell manta there are 2 GiB per job, otherwise it will
# scale back the requested number of jobs, even if they won't
# need that much memory
${working_dir}/runWorkflow.py \
  --mode local \
  --jobs "$num_jobs" \
  --memGb ${mem_size}

# inversion conversion, then compression and index
python2 /usr/local/bin/manta/libexec/convertInversion.py \
  /opt/samtools/bin/samtools \
  "$reference_fasta" \
  ${working_dir}/results/variants/diploidSV.vcf.gz \
  | bcftools reheader -s <(echo "$sample_id") \
  > ${working_dir}/diploidSV.vcf

bgzip -c ${working_dir}/diploidSV.vcf > "${working_dir}/$sample_id.manta.vcf.gz"
tabix -p vcf "${working_dir}/${sample_id}.manta.vcf.gz"

output_vcf_filename="$(realpath ${output_dir}/$sample_id.manta.vcf.gz)"
output_vcf_index_filename="$(realpath ${output_dir}/$sample_id.manta.vcf.gz.tbi)"
mv "${working_dir}/$sample_id.manta.vcf.gz" "${output_vcf_filename}"
mv "${working_dir}/${sample_id}.manta.vcf.gz.tbi" "${output_vcf_index_filename}"

outputs_filename="${output_dir}/outputs.json"
outputs_json=$(jq -n \
  --arg vcf "${output_vcf_filename}" \
  --arg vcf_idx "${output_vcf_index_filename}" \
  '{vcf: $vcf, index: $vcf_idx}' )
echo "${outputs_json}" > "${outputs_filename}"
cp "${outputs_filename}" "${outputs_json_filename}"
