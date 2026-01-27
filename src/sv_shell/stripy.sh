#!/bin/bash

set -Exeuo pipefail

# -------------------------------------------------------
# ==================== Input & Setup ====================
# -------------------------------------------------------


input_json=${1}
output_json_filename=${2-""}
output_dir=${3:-""}

input_json="$(realpath ${input_json})"

if [ -z "${output_dir}" ]; then
  output_dir=$(mktemp -d /output_stripy_XXXXXXXX)
else
  mkdir -p "${output_dir}"
fi
output_dir="$(realpath ${output_dir})"

if [ -z "${output_json_filename}" ]; then
  output_json_filename="${output_dir}/output.json"
else
  output_json_filename="$(realpath ${output_json_filename})"
fi

working_dir=$(realpath $(mktemp -d "/wd_stripy_XXXXXXXX"))
cd "${working_dir}"
echo "stripy Working directory: ${working_dir}"


ped_file=$(jq -r ".ped_file" "$input_json")
sample_name=$(jq -r ".sample_name" "$input_json")
bam_or_cram_file=$(jq -r ".bam_or_cram_file" "$input_json")
genome_build=$(jq -r '.genome_build // "hg38"' "$input_json")
reference_fasta=$(jq -r ".reference_fasta" "$input_json")
analysis=$(jq -r '.analysis // "standard"' "$input_json")
locus=$(jq -r '.locus // "AFF2,AR,ARX_1,ARX_2,ATN1,ATXN1,ATXN10,ATXN2,ATXN3,ATXN7,ATXN8OS,BEAN1,C9ORF72,CACNA1A,CBL,CNBP,COMP,DAB1,DIP2B,DMD,DMPK,FGF14,FMR1,FOXL2,FXN,GIPC1,GLS,HOXA13_1,HOXA13_2,HOXA13_3,HOXD13,HTT,JPH3,LRP12,MARCHF6,NIPA1,NOP56,NOTCH2NLC,NUTM2B-AS1,PABPN1,PHOX2B,PPP2R2B,PRDM12,RAPGEF2,RFC1,RILPL1,RUNX2,SAMD12,SOX3,STARD7,TBP,TBX1,TCF4,TNRC6A,XYLT1,YEATS2,ZIC2,ZIC3"' "$input_json")


# -------------------------------------------------------
# ======================= Command =======================
# -------------------------------------------------------

# GetSampleSex
# ---------------------------------------------------------------------------------------------------------------------
unknown_sex="female"

awk -F '\t' -v smp="${sample_name}" -v unk="${unknown_sex}" '{if ($2==smp) { if ($5 == 1) {print "male"} else if ($5 == 2) {print "female"} else {print unk}}}' "${ped_file}" > "${sample_name}.sex.txt"

# Fail if the sample id wasn't found
if ! [ -s "${sample_name}.sex.txt" ]; then
  echo "ERROR: Sample ${sample_name} not found in ped file ${ped_file}"
  exit 1
fi

GetSampleSex_out_file="$(realpath "${sample_name}.sex.txt")"
GetSampleSex_out_string=$(cat "${GetSampleSex_out_file}")


# RunStripy
# ---------------------------------------------------------------------------------------------------------------------
stripy_output_dir="STRipy_output"
mkdir -p "${stripy_output_dir}"
stripy \
  --input "${bam_or_cram_file}" \
  --genome "${genome_build}" \
  --reference "${reference_fasta}" \
  --output "${stripy_output_dir}" \
  --analysis "${analysis}" \
  --output-json true \
  --output-tsv true \
  --output-html true \
  --output-vcf true \
  --verbose false \
  --num-threads $(nproc) \
  --locus "${locus}" \
  --sex "${GetSampleSex_out_string}"

ACTUAL_FILENAME=$(basename "${bam_or_cram_file}")
ACTUAL_BASE=$(echo "${ACTUAL_FILENAME}" | sed 's/\.[^.]*$//')
TARGET_BASE="${sample_name}"
echo "ACTUAL_FILENAME: ${ACTUAL_FILENAME}"
echo "ACTUAL_BASE: ${ACTUAL_BASE}"
echo "TARGET_BASE: ${TARGET_BASE}"
ls "${stripy_output_dir}"
if [ -f "${stripy_output_dir}/${ACTUAL_FILENAME}.json" ]; then
  mv "${stripy_output_dir}/${ACTUAL_FILENAME}.json" "${stripy_output_dir}/${TARGET_BASE}.json"
fi
if [ -f "${stripy_output_dir}/${ACTUAL_FILENAME}.tsv" ]; then
  mv "${stripy_output_dir}/${ACTUAL_FILENAME}.tsv" "${stripy_output_dir}/${TARGET_BASE}.tsv"
fi
if [ -f "${stripy_output_dir}/${ACTUAL_FILENAME}.html" ]; then
  mv "${stripy_output_dir}/${ACTUAL_FILENAME}.html" "${stripy_output_dir}/${TARGET_BASE}.html"
fi
if [ -f "${stripy_output_dir}/${ACTUAL_BASE}.vcf" ]; then
  mv "${stripy_output_dir}/${ACTUAL_BASE}.vcf" "${stripy_output_dir}/${TARGET_BASE}.vcf"
fi


# -------------------------------------------------------
# ======================= Output ========================
# -------------------------------------------------------

json_path="${stripy_output_dir}/${sample_name}.json"
tsv_path="${stripy_output_dir}/${sample_name}.tsv"
html_path="${stripy_output_dir}/${sample_name}.html"
vcf_path="${stripy_output_dir}/${sample_name}.vcf"

json_path_output_dir="$(realpath "${output_dir}/$(basename "${json_path}")")"
mv "${json_path}" "${json_path_output_dir}"

tsv_path_output_dir="$(realpath "${output_dir}/$(basename "${tsv_path}")")"
mv "${tsv_path}" "${tsv_path_output_dir}"

html_path_output_dir="$(realpath "${output_dir}/$(basename "${html_path}")")"
mv "${html_path}" "${html_path_output_dir}"

vcf_path_output_dir="$(realpath "${output_dir}/$(basename "${vcf_path}")")"
mv "${vcf_path}" "${vcf_path_output_dir}"


jq -n \
  --arg json "${json_path_output_dir}" \
  --arg tsv "${tsv_path_output_dir}" \
  --arg html "${html_path_output_dir}" \
  --arg vcf "${vcf_path_output_dir}" \
  '{
      json_output: $json,
      tsv_output: $tsv,
      html_output: $html,
      vcf_output: $vcf
  }' > "${output_json_filename}"

echo "Finished stripy, output json filename: ${output_json_filename}"
