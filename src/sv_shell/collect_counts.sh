#!/bin/bash

# This script is a bash implementation of the following workflow/task in WDL:
# Filename: wdl/CollectCoverage.wdl
# Task: CollectCounts

set -Exeuo pipefail

intervals=$1
cram_or_bam=$2
cram_or_bam_idx=$3
sample_id=$4
ref_fasta=$5
ref_fasta_fai=$6
ref_fasta_dict=$7
outputs_json_filename=$8
gatk4_jar_override=${9:-/root/gatk.jar}
disabled_read_filters=${10:-""}

echo "=============== Running collect_counts.sh"
echo "intervals:             " "${intervals}"
echo "cram_or_bam:           " "${cram_or_bam}"
echo "cram_or_bam_idx:       " "${cram_or_bam_idx}"
echo "sample_id:             " "${sample_id}"
echo "ref_fasta:             " "${ref_fasta}"
echo "ref_fasta_fai:         " "${ref_fasta_fai}"
echo "ref_fasta_dict:        " "${ref_fasta_dict}"
echo "gatk4_jar_override:    " "${gatk4_jar_override}"
echo "disabled_read_filters: " "${disabled_read_filters}"

working_dir=$(mktemp -d /wd_collect_counts_XXXXXXXX)
working_dir="$(realpath ${working_dir})"
output_dir=$(mktemp -d /output_collect_counts_XXXXXXXX)
output_dir="$(realpath ${output_dir})"
cd "${working_dir}"

# In the WDL implementation, this is computed based on a few factors,
# the following is the result of the computation using the default values.
command_mem_mb=10240

declare -a disabled_read_filters_arr
if [ -n "$disabled_read_filters" ]; then
  disabled_read_filters_arr=( "--disable-read-filter" "${disabled_read_filters}" )
else
  disabled_read_filters_arr=()
fi

export GATK_LOCAL_JAR="$gatk4_jar_override"

java -Xmx${command_mem_mb}m -jar /opt/gatk.jar CollectReadCounts \
  -L "${intervals}" \
  --input "${cram_or_bam}" \
  --read-index "${cram_or_bam_idx}" \
  --reference "${ref_fasta}" \
  --format TSV \
  --interval-merging-rule OVERLAPPING_ONLY \
  --output "${sample_id}".counts.tsv \
  "${disabled_read_filters_arr[@]}"

sed -ri "s/@RG\tID:GATKCopyNumber\tSM:.+/@RG\tID:GATKCopyNumber\tSM:${sample_id}/g" "${sample_id}.counts.tsv"
bgzip --force "${sample_id}.counts.tsv"

counts_filename="${output_dir}/${sample_id}.counts.tsv.gz"
mv "${sample_id}.counts.tsv.gz" "${counts_filename}"

outputs_filename="${output_dir}/outputs.json"
outputs_json=$(jq -n \
  --arg c "${counts_filename}" \
  '{counts: $c}' )
echo "${outputs_json}" > "${outputs_filename}"
cp "${outputs_filename}" "${outputs_json_filename}"
