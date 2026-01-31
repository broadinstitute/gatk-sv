#!/bin/bash

set -Exeuo pipefail

function getJavaMem() {
  # get JVM memory in MiB by getting total memory from /proc/meminfo
  # and multiplying by java_mem_fraction

  local mem_fraction=${java_mem_fraction:=0.85}
  cat /proc/meminfo | \
    awk -v MEM_FIELD="$1" -v frac="${mem_fraction}" '{
      f[substr($1, 1, length($1)-1)] = $2
    } END {
      printf "%dM", f[MEM_FIELD] * frac / 1024
    }'
}
JVM_MAX_MEM=$(getJavaMem MemTotal)
echo "JVM memory: $JVM_MAX_MEM"

MergeEvidence() {
  ulimit -n 100000

  local -n _files=$1
  local -n _samples=$2
  local _primary_contigs_fai=$3
  local _reference_dict=$4
  local _output_filename=$5

  evidence_file=$(mktemp evidence_XXXXX.list)
  printf "%s\n" "${_files[@]}" > "${evidence_file}"

  samples_file=$(mktemp samples_XXXXX.list)
  printf "%s\n" "${_samples[@]}" > "${samples_file}"

  # Note that this piece has a few differences with the version in WDL.
  #  - It does not support renaming or sub-setting contings as they are not needed for the single-sample pipeline.
  #  - It only indexes files if the file is missing an index.
  awk '/txt\.gz$/' "${evidence_file}" | while read fil; do
    if [ ! -f "${fil}.tbi" ]; then
      echo "Tabix indexing ${fil} ..."
      tabix -f -0 -s1 -b2 -e2 $fil
      echo "Done."
    fi
  done

  java "-Xmx${JVM_MAX_MEM}" -jar /opt/gatk.jar PrintSVEvidence \
    -F "${evidence_file}" \
    --sample-names "${samples_file}" \
    --sequence-dictionary "${_reference_dict}" \
    -O "${_output_filename}"
}


SDtoBAF() {
  ulimit -n 100000

  local -n _sd_files=$1
  local _n _samples=$2
  local _sd_locs_vcf=$3
  local _reference_dict=$4
  local _min_het_probability=$5
  local _output_filename=$6

  inputs_file=$(mktemp inputs_XXXXX.list)
  printf "%s\n" "${_sd_files[@]}" > "${inputs_file}"

  samples_file=$(mktemp samples_XXXXX.list)
  printf "%s\n" "${_samples[@]}" > "${samples_file}"

  awk '/txt\.gz$/' "${inputs_file}" | while read fil; do
    tabix -f -s1 -b2 -e2 $fil
  done

  java "-Xmx${JVM_MAX_MEM}" -jar /opt/gatk.jar SiteDepthtoBAF \
      -F "${inputs_file}" \
      --sample-names "${samples_file}" \
      --sequence-dictionary "${_reference_dict}" \
      --baf-sites-vcf "${_sd_locs_vcf}" \
      --min-het-probability "${_min_het_probability}" \
      -O "${_output_filename}"
}


# -------------------------------------------------------
# ==================== Input & Setup ====================
# -------------------------------------------------------


input_json=${1}
output_json_filename=${2-""}
output_dir=${3:-""}

input_json="$(realpath ${input_json})"

if [ -z "${output_dir}" ]; then
  output_dir=$(mktemp -d ${SV_SHELL_BASE_DIR}/output_batch_evidence_merging_XXXXXXXX)
else
  mkdir -p "${output_dir}"
fi
output_dir="$(realpath ${output_dir})"

if [ -z "${output_json_filename}" ]; then
  output_json_filename="${output_dir}/output.json"
else
  output_json_filename="$(realpath ${output_json_filename})"
fi

working_dir=$(mktemp -d ${SV_SHELL_BASE_DIR}/wd_batch_evidence_merging_XXXXXXXX)
working_dir="$(realpath ${working_dir})"
cd "${working_dir}"


batch=($(jq -r '.batch' "$input_json"))
samples=($(jq -r '.samples[]' "$input_json"))
PE_files=($(jq -r '.PE_files[]' "$input_json"))
SR_files=($(jq -r '.SR_files[]' "$input_json"))
SD_files=($(jq -r '.SD_files[]' "$input_json"))
sd_locs_vcf=$(jq -r ".sd_locs_vcf" "${input_json}")
reference_dict=$(jq -r ".reference_dict" "${input_json}")
primary_contigs_fai=$(jq -r ".primary_contigs_fai" "${input_json}")
min_het_probability=$(jq -r ".min_het_probability" "${input_json}")


# -------------------------------------------------------
# ======================= Command =======================
# -------------------------------------------------------


# ---- MergeSREvidence
echo "Started running Merge SR Evidence."
merged_sr_output_filename="${batch}.sr.txt.gz"
merged_sr_output_filename="$(realpath ${merged_sr_output_filename})"
# Note that "SR_files" and "samples" in the following is not passing a string, it is rather passing the samples variable.
MergeEvidence "SR_files" "samples" "${primary_contigs_fai}" "${reference_dict}" "${merged_sr_output_filename}"
echo "Finished running Merge SR Evidence."


# ---- MergePEEvidence
echo "Started running Merge PE Evidence."
merged_pe_output_filename="${batch}.pe.txt.gz"
merged_pe_output_filename="$(realpath ${merged_pe_output_filename})"
# Note that "PE_files" and "samples" in the following is not passing a string, it is rather passing the samples variable.
MergeEvidence "PE_files" "samples" "${primary_contigs_fai}" "${reference_dict}" "${merged_pe_output_filename}"
echo "Finished running Merge SR Evidence."


echo "Started running SD to BAF."
sd_to_baf_output_filename="${batch}.baf.txt.gz"
sd_to_baf_output_filename="$(realpath ${sd_to_baf_output_filename})"
# Note that "SD_files" and "samples" in the following is not passing a string, it is rather passing the samples variable.
SDtoBAF "SD_files" "samples" "${sd_locs_vcf}" "${reference_dict}" "${min_het_probability}" "${sd_to_baf_output_filename}"
echo "Finished running SD to BAF."


# -------------------------------------------------------
# ======================= Output ========================
# -------------------------------------------------------


sd_to_baf_output_filename_outdir="${output_dir}/$(basename "${sd_to_baf_output_filename}")"
mv "${sd_to_baf_output_filename}" "${sd_to_baf_output_filename_outdir}"
mv "${sd_to_baf_output_filename}.tbi" "${sd_to_baf_output_filename_outdir}.tbi"

merged_sr_output_filename_outdir="${output_dir}/$(basename "${merged_sr_output_filename}")"
mv "${merged_sr_output_filename}" "${merged_sr_output_filename_outdir}"
mv "${merged_sr_output_filename}.tbi" "${merged_sr_output_filename_outdir}.tbi"

merged_pe_output_filename_outdir="${output_dir}/$(basename "${merged_pe_output_filename}")"
mv "${merged_pe_output_filename}" "${merged_pe_output_filename_outdir}"
mv "${merged_pe_output_filename}.tbi" "${merged_pe_output_filename_outdir}.tbi"

outputs_json=$(jq -n \
  --arg merged_BAF "${sd_to_baf_output_filename_outdir}" \
  --arg merged_BAF_index "${sd_to_baf_output_filename_outdir}.tbi" \
  --arg merged_SR "${merged_sr_output_filename_outdir}" \
  --arg merged_SR_index "${merged_sr_output_filename_outdir}.tbi" \
  --arg merged_PE "${merged_pe_output_filename_outdir}" \
  --arg merged_PE_index "${merged_pe_output_filename_outdir}.tbi" \
  '{
     "merged_BAF": $merged_BAF,
     "merged_BAF_index": $merged_BAF_index,
     "merged_SR": $merged_SR,
     "merged_SR_index": $merged_SR_index,
     "merged_PE": $merged_PE,
     "merged_PE_index": $merged_PE_index
   }' \
)
echo "${outputs_json}" > "${output_json_filename}"

echo "Output JSON filename: ${output_json_filename}"
