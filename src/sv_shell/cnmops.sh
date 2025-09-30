#!/bin/bash

set -Eeuo pipefail


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


function CNSampleNormal() {
  # TODO: check if these variables are correctly set
  local _chr=$1
  local _mode=$2
  local _r=$3
  local _wd=$4

  echo "----------- Starting CN Sample Normal -------------"
  echo "chr: ${_chr}"
  echo "mode: ${_mode}"
  echo "r: ${_r}"
  echo "Working dir: ${_wd}"
  echo "---------------------------------------------------"

  cd "${_wd}"

  java "-Xmx${JVM_MAX_MEM}" -jar /opt/gatk.jar PrintSVEvidence \
    --sequence-dictionary "${ref_dict}" \
    --evidence-file "${bincov_matrix}" \
    -L "${_chr}" \
    -O "${_chr}.RD.txt"

  if [ "${_mode}" == "normal" ]; then
    mv "${_chr}.RD.txt" "${_chr}.${_mode}.RD.txt"
  else
    awk -v sex="${_mode}" '$5==sex' "${ped_file}" | cut -f2 > ids.to.include
    col=$(head -n 1 "${_chr}.RD.txt" | tr '\t' '\n'|cat -n| grep -wf ids.to.include | awk -v ORS="," '{print $1}' | sed 's/,$//g' | sed 's:\([0-9]\+\):$&:g')
    col_a="{print \$1,\$2,\$3,$col}"
    awk -f <(echo "$col_a") "${_chr}.RD.txt" | tr ' ' '\t' > "${_chr}.${_mode}.RD.txt"
  fi

  # redirect stdout and stderr to cnmops.out so that EMPTY_OUTPUT_ERROR can be detected, but use tee to also output them to
  # terminal so that errors can be debugged
  EMPTY_OUTPUT_ERROR="No CNV regions in result object. Rerun cn.mops with different parameters!"
  set +e
  echo "Starting to run cnMOPS_workflow"
  bash /opt/WGD/bin/cnMOPS_workflow.sh -S "${exclude_list}" -x "${exclude_list}" -r "${_r}" -o . -M "${_chr}.${_mode}.RD.txt" </dev/null wor2>&1 | tee cnmops.out
  echo "Finished running cnMOPS_workflow"
  RC=$?
  set -e
  if [ ! $RC -eq 0 ]; then
    if grep -q "$EMPTY_OUTPUT_ERROR" "cnmops.out"; then
      touch calls/cnMOPS.cnMOPS.gff
    else
      echo "cnMOPS_workflow.sh returned a non-zero code that was not due to an empty call file."
      exit $RC
    fi
  fi

  echo "----------- Finished CN Sample Normal -------------"
}






# -------------------------------------------------------
# ==================== Input & Setup ====================
# -------------------------------------------------------


input_json=${1}
output_json_filename=${2-""}
output_dir=${3:-""}

input_json="$(realpath ${input_json})"

if [ -z "${output_dir}" ]; then
  output_dir=$(mktemp -d /output_batch_evidence_merging_XXXXXXXX)
else
  mkdir -p "${output_dir}"
fi
output_dir="$(realpath ${output_dir})"

if [ -z "${output_json_filename}" ]; then
  output_json_filename="${output_dir}/output.json"
else
  output_json_filename="$(realpath ${output_json_filename})"
fi

working_dir=$(mktemp -d /wd_batch_evidence_merging_XXXXXXXX)
working_dir="$(realpath ${working_dir})"
cd "${working_dir}"


allo_file=($(jq -r '.allo_file' "$input_json"))
chrom_file=($(jq -r '.chrom_file' "$input_json"))
exclude_list=($(jq -r '.exclude_list' "$input_json"))
ped_file=($(jq -r '.ped_file' "$input_json"))
r1=($(jq -r '.r1' "$input_json"))
r2=($(jq -r '.r2' "$input_json"))
ref_dict=($(jq -r '.ref_dict' "$input_json"))
bincov_matrix=($(jq -r '.bincov_matrix' "$input_json"))


# -------------------------------------------------------
# ======================= Command =======================
# -------------------------------------------------------


allos=($(awk '{print $1}' "${allo_file}"))

for allo in "${allos[@]}"; do
  # Male R2
  working_dir=$(mktemp -d /wd_cn_sample_normal_${allo}_${r2}_XXXXXXXX)
  working_dir="$(realpath ${working_dir})"
  CNSampleNormal "${allo}" "1" "${r2}" "${working_dir}"

  # Male R1
  working_dir=$(mktemp -d /wd_cn_sample_normal_${allo}_${r1}_XXXXXXXX)
  working_dir="$(realpath ${working_dir})"
  CNSampleNormal "${allo}" "1" "${r1}" "${working_dir}"
done

chroms=($(awk '{print $1}' "${chrom_file}"))
for chrom in "${chroms[@]}"; do
  # Normal R2
  working_dir=$(mktemp -d /wd_cn_sample_normal_${chrom}_${r2}_XXXXXXXX)
  working_dir="$(realpath ${working_dir})"
  CNSampleNormal "${chrom}" "normal" "${r2}" "${working_dir}"

  # Normal R1
  working_dir=$(mktemp -d /wd_cn_sample_normal_${chrom}_${r1}_XXXXXXXX)
  working_dir="$(realpath ${working_dir})"
  CNSampleNormal "${chrom}" "normal" "${r1}" "${working_dir}"
done



first_row_string="${Allos[0]}"

echo "${first_row_string}"

# Split that string into a new array called 'fields'
# read -ra splits on whitespace (tabs and spaces)
read -ra fields <<< "$first_row_string"

# Now access an element from the 'fields' array (e.g., the 2nd field is at index 1)
echo "${fields[1]}"