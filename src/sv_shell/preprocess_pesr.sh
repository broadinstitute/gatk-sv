#!/bin/bash

set -Eeuo pipefail


function StandardizeVCFs() {
  declare -n _raw_vcfs=$1
  declare -n _samples=$2
  local _caller=$3
  local _contigs=$4
  local _min_svsize=$5
  local _prefix=$6


  echo "----------- Starting StandardizeVCFs -------------"
  echo "_raw_vcfs: ${_raw_vcfs[@]}"
  echo "_samples: ${_samples[@]}"
  echo "_caller: ${_caller}"
  echo "_contigs: ${_contigs}"
  echo "_min_svsize: ${_min_svsize}"
  echo "_prefix: ${_prefix}"
  echo "-------------------------------------------------"

  #  vcfs=(~{sep=" " raw_vcfs})
  #  sample_ids=(~{sep=" " samples})
  num_samples="${#_samples[@]}"

  mkdir out
  for (( i=0; i<"${num_samples}"; i++ ));
  do
    vcf="${_raw_vcfs[$i]}"
    sample_id="${_samples[$i]}"
    svtk standardize \
      --sample-names "${sample_id}" \
      --prefix "${_caller}_${sample_id}" \
      --contigs "${_contigs}" \
      --min-size "${min_svsize}" \
      "${vcf}" \
      tmp.vcf \
      "${_caller}"
    sample_no=$(printf %03d $i)
    bcftools sort tmp.vcf -Oz -o "out/std_${sample_no}.${_caller}.${sample_id}.vcf.gz"
  done
  tar czf "${_prefix}.tar.gz" -C out/ .

  echo "----------- Finished StandardizeVCFs -------------"
}

# -------------------------------------------------------
# ==================== Input & Setup ====================
# -------------------------------------------------------


input_json=${1}
output_json_filename=${2-""}
output_dir=${3:-""}

input_json="$(realpath ${input_json})"

if [ -z "${output_dir}" ]; then
  output_dir=$(mktemp -d /output_preprocess_pesr_XXXXXXXX)
else
  mkdir -p "${output_dir}"
fi
output_dir="$(realpath ${output_dir})"

if [ -z "${output_json_filename}" ]; then
  output_json_filename="${output_dir}/output.json"
else
  output_json_filename="$(realpath ${output_json_filename})"
fi

working_dir=$(mktemp -d /wd_preprocess_pesr_XXXXXXXX)
working_dir="$(realpath ${working_dir})"
cd "${working_dir}"
echo "Preprocess PE-SR Working directory: ${working_dir}"

samples=($(jq -r '.samples[]' "$input_json"))
dragen_vcfs=($(jq -r '.dragen_vcfs[] // ""' "$input_json"))
manta_vcfs=($(jq -r '.manta_vcfs[] // ""' "$input_json"))
scramble_vcfs=($(jq -r '.scramble_vcfs[] // ""' "$input_json"))
wham_vcfs=($(jq -r '.wham_vcfs[] // ""' "$input_json"))
contigs=$(jq -r '.contigs // ""' "$input_json")
min_svsize=$(jq -r '.min_svsize // ""' "$input_json")
batch=$(jq -r '.batch // ""' "$input_json")


# -------------------------------------------------------
# ======================= Command =======================
# -------------------------------------------------------


# note that in the following methods, "ARRAY_NAME" is passing the array, not a string.
# future implementation improvement: the following should be implemented as a for loop, check wdl implementation


dragen_out_output=""
if (( "${#dragen_vcfs[@]}" > 0 )); then
  algorithm="dragen"
  prefix="${batch}.${algorithm}_std"
  working_dir=$(mktemp -d /wd_preprocess_pesr_dragen_XXXXXXXX)
  working_dir="$(realpath ${working_dir})"
  cd "${working_dir}"

  StandardizeVCFs "dragen_vcfs" "samples" "${algorithm}" "${contigs}" "${min_svsize}" "${prefix}"

  dragen_out="${working_dir}/${prefix}.tar.gz"
  dragen_out_output="${output_dir}/$(basename "${dragen_out}")"
  mv "${dragen_out}" "${dragen_out_output}"
fi


manta_out_output=""
if (( "${#manta_vcfs[@]}" > 0 )); then
  algorithm="manta"
  prefix="${batch}.${algorithm}_std"
  working_dir=$(mktemp -d /wd_preprocess_pesr_manta_XXXXXXXX)
  working_dir="$(realpath ${working_dir})"
  cd "${working_dir}"

  StandardizeVCFs "manta_vcfs" "samples" "${algorithm}" "${contigs}" "${min_svsize}" "${prefix}"

  manta_out="${working_dir}/${prefix}.tar.gz"
  manta_out_output="${output_dir}/$(basename "${manta_out}")"
  mv "${manta_out}" "${manta_out_output}"
fi


scramble_out_output=""
if (( "${#manta_vcfs[@]}" > 0 )); then
  algorithm="scramble"
  prefix="${batch}.${algorithm}_std"
  working_dir=$(mktemp -d /wd_preprocess_pesr_scramble_XXXXXXXX)
  working_dir="$(realpath ${working_dir})"
  cd "${working_dir}"

  StandardizeVCFs "scramble_vcfs" "samples" "${algorithm}" "${contigs}" "${min_svsize}" "${prefix}"

  scramble_out="${working_dir}/${prefix}.tar.gz"
  scramble_out_output="${output_dir}/$(basename "${scramble_out}")"
  mv "${scramble_out}" "${scramble_out_output}"
fi


wham_out_output=""
if (( "${#wham_vcfs[@]}" > 0 )); then
  algorithm="wham"
  prefix="${batch}.${algorithm}_std"
  working_dir=$(mktemp -d /wd_preprocess_pesr_wham_XXXXXXXX)
  working_dir="$(realpath ${working_dir})"
  cd "${working_dir}"

  StandardizeVCFs "wham_vcfs" "samples" "${algorithm}" "${contigs}" "${min_svsize}" "${prefix}"

  wham_out="${working_dir}/${prefix}.tar.gz"
  wham_out_output="${output_dir}/$(basename "${wham_out}")"
  mv "${wham_out}" "${wham_out_output}"
fi


# -------------------------------------------------------
# ======================= Output ========================
# -------------------------------------------------------


outputs_json=$(jq -n \
  --arg dragen "${dragen_out_output}" \
  --arg manta "${manta_out_output}" \
  --arg scramble "${scramble_out_output}" \
  --arg wham "${wham_out_output}" \
  '{std_dragen_vcf_tar: $dragen, std_manta_vcf_tar: $manta, std_scramble_vcf_tar: $scramble, std_wham_vcf_tar: $wham}' )
echo "${outputs_json}" > "${output_json_filename}"

echo "Finished preprocess PESR successfully, output json filename: ${output_json_filename}"
