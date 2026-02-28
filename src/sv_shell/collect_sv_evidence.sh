#!/bin/bash

# This script is a bash implementation of the following workflow/task in WDL:
# Filename: wdl/CollectSVEvidence.wdl
# Workflow: CollectSVEvidence

set -Exeuo pipefail

# -------------------------------------------------------
# ==================== Input & Setup ====================
# -------------------------------------------------------

input_json=${1}
output_json_filename=${2-""}
output_dir=${3:-""}

input_json="$(realpath ${input_json})"

if [ -z "${output_dir}" ]; then
  output_dir=$(mktemp -d ${SV_SHELL_BASE_DIR}/output_collect_sv_evidence_XXXXXXXX)
else
  mkdir -p "${output_dir}"
fi
output_dir="$(realpath ${output_dir})"

if [ -z "${output_json_filename}" ]; then
  output_json_filename="${output_dir}/output.json"
else
  output_json_filename="$(realpath ${output_json_filename})"
fi

working_dir=$(mktemp -d ${SV_SHELL_BASE_DIR}/wd_collect_sv_evidence_XXXXXXXX)
working_dir="$(realpath ${working_dir})"
cd "${working_dir}"
echo "collect_sv_evidence working directory: ${working_dir}"


sample_id=$(jq -r ".sample_id" "$input_json")
bam_or_cram_file=$(jq -r ".bam_or_cram_file" "$input_json")
bam_or_cram_index=$(jq -r ".bam_or_cram_index" "$input_json")
reference_fasta=$(jq -r ".reference_fasta" "$input_json")
reference_index=$(jq -r ".reference_index" "$input_json")
reference_dict=$(jq -r ".reference_dict" "$input_json")
sd_locs_vcf=$(jq -r ".sd_locs_vcf" "$input_json")
preprocessed_intervals=$(jq -r ".preprocessed_intervals" "$input_json")
site_depth_min_mapq=$(jq -r ".site_depth_min_mapq // 6" "$input_json")
site_depth_min_baseq=$(jq -r ".site_depth_min_baseq // 10" "$input_json")
primary_contigs_list=$(jq -r '.primary_contigs_list // ""' "$input_json")
command_mem_mb=$(jq -r '.command_mem_mb // 3250' "$input_json")


# -------------------------------------------------------
# ======================= Command =======================
# -------------------------------------------------------


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

tabix -f -0 -s1 -b2 -e2 "${sample_id}.sr.txt.gz"
tabix -f -0 -s1 -b2 -e2 "${sample_id}.pe.txt.gz"
tabix -f -0 -s1 -b2 -e2 "${sample_id}.sd.txt.gz"

# -------------------------------------------------------
# ======================= Output ========================
# -------------------------------------------------------

split_out_filename="${output_dir}/${sample_id}.sr.txt.gz"
mv "${sample_id}.sr.txt.gz" "${split_out_filename}"
mv "${sample_id}.sr.txt.gz.tbi" "${split_out_filename}.tbi"

disc_out_filename="${output_dir}/${sample_id}.pe.txt.gz"
mv "${sample_id}.pe.txt.gz" "${disc_out_filename}"
mv "${sample_id}.pe.txt.gz.tbi" "${disc_out_filename}.tbi"

sd_out_filename="${output_dir}/${sample_id}.sd.txt.gz"
mv "${sample_id}.sd.txt.gz" "${sd_out_filename}"
mv "${sample_id}.sd.txt.gz.tbi" "${sd_out_filename}.tbi"

jq -n \
  --arg split_out "${split_out_filename}" \
  --arg split_out_index "${split_out_filename}.tbi" \
  --arg disc_out "${disc_out_filename}" \
  --arg disc_out_index "${disc_out_filename}.tbi" \
  --arg sd_out "${sd_out_filename}" \
  --arg sd_out_index "${sd_out_filename}.tbi" \
  '{
      "split_out": $split_out,
      "split_out_index": $split_out_index,
      "disc_out": $disc_out,
      "disc_out_index": $disc_out_index,
      "sd_out": $sd_out,
      "sd_out_index": $sd_out_index
  }' > "${output_json_filename}"

echo "Successfully finished collect sv evidence, output json filename: ${output_json_filename}"