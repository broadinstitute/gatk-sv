#!/bin/bash

# This script is a bash implementation of the following workflow/task in WDL:
# Filename: wdl/CollectSVEvidence.wdl
# Workflow: CollectSVEvidence

set -Exeuo pipefail

sample_id=${1}
bam_or_cram_file=${2}
bam_or_cram_index=${3}
reference_fasta=${4}
reference_index=${5}
reference_dict=${6}
sd_locs_vcf=${7}
preprocessed_intervals=${8}
outputs_json_filename=${9}
site_depth_min_mapq=${10:-6}
site_depth_min_baseq=${11:-10}
primary_contigs_list="${12:-}"
gatk_jar_override="${13:-/root/gatk.jar}"
command_mem_mb=${14:-3250}

working_dir=$(mktemp -d wd_collect_sv_evidence_XXXXXXXX)
working_dir="$(realpath ${working_dir})"
output_dir=$(mktemp -d output_collect_sv_evidence_XXXXXXXX)
output_dir="$(realpath ${output_dir})"
cd "${working_dir}"

export GATK_LOCAL_JAR="$gatk_jar_override"

primary_contigs_arg=""
if [[ -n "$primary_contigs_list" ]]; then
  primary_contigs_arg="-L $primary_contigs_list"
fi

java -Xmx${command_mem_mb}m -jar /opt/gatk.jar CollectSVEvidence \
    -I "${bam_or_cram_file}" \
    --sample-name "${sample_id}" \
    -F "${sd_locs_vcf}" \
    -SR "${sample_id}.sr.txt.gz" \
    -PE "${sample_id}.pe.txt.gz" \
    -SD "${sample_id}.sd.txt.gz" \
    --site-depth-min-mapq "${site_depth_min_mapq}" \
    --site-depth-min-baseq "${site_depth_min_baseq}" \
    -R "${reference_fasta}" \
    ${primary_contigs_arg} \
    --read-filter NonZeroReferenceLengthAlignmentReadFilter

split_out_filename="${output_dir}/${sample_id}.sr.txt.gz"
split_out_index_filename="${output_dir}/${sample_id}.sr.txt.gz.tbi"
disc_out_filename="${output_dir}/${sample_id}.pe.txt.gz"
disc_out_index_filename="${output_dir}/${sample_id}.pe.txt.gz.tbi"
sd_out_filename="${output_dir}/${sample_id}.sd.txt.gz"
sd_out_index_filename="${output_dir}/${sample_id}.sd.txt.gz.tbi"

outputs_filename="${output_dir}/outputs.json"
outputs_json=$(jq -n \
  --arg split_out "${split_out_filename}" \
  --arg split_out_index "${split_out_index_filename}" \
  --arg disc_out "${disc_out_filename}" \
  --arg disc_out_index "${disc_out_index_filename}" \
  --arg sd_out "${sd_out_filename}" \
  --arg sd_out_index "${sd_out_index_filename}" \
  '{split_out: $split_out, split_out_index: $split_out_index, disc_out: $disc_out, disc_out_index: $disc_out_index, sd_out: $sd_out, sd_out_index: $sd_out_index}')
echo "${outputs_json}" > "${outputs_filename}"
cp "${outputs_filename}" "${outputs_json_filename}"