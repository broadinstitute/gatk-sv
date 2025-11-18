#!/bin/bash

set -Exeuo pipefail

# -------------------------------------------------------
# ==================== Input & Setup ====================
# -------------------------------------------------------

input_json=${1}
output_json_filename=${2:-""}
output_dir=${3:-""}

input_json="$(realpath ${input_json})"

if [ -z "${output_dir}" ]; then
  output_dir=$(mktemp -d /output_sv_cluster_XXXXXXXX)
else
  mkdir -p "${output_dir}"
fi
output_dir="$(realpath ${output_dir})"

if [ -z "${output_json_filename}" ]; then
  output_json_filename="${output_dir}/output.json"
else
  output_json_filename="$(realpath ${output_json_filename})"
fi

working_dir="$(mktemp -d /wd_sv_cluster_XXXXXXXX)"
working_dir="$(realpath ${working_dir})"
cd "${working_dir}"

vcfs=($(jq -r '(.vcfs // [])[]' "$input_json"))
vcfs_tar=$(jq -r '.vcfs_tar // ""' "$input_json")


output_prefix=$(jq -r '.output_prefix' "$input_json")
ploidy_table=$(jq -r '.ploidy_table' "$input_json")
reference_fasta=$(jq -r '.reference_fasta' "$input_json")

contig=$(jq -r '.contig // ""' "$input_json")
variant_prefix=$(jq -r '.variant_prefix // ""' "$input_json")
default_no_call=$(jq -r '.default_no_call // ""' "$input_json")
omit_members=$(jq -r '.omit_members // ""' "$input_json")
enable_cnv=$(jq -r '.enable_cnv // ""' "$input_json")
fast_mode=$(jq -r '.fast_mode // ""' "$input_json")
algorithm=$(jq -r '.algorithm // ""' "$input_json")
defrag_padding_fraction=$(jq -r '.defrag_padding_fraction // ""' "$input_json")
depth_sample_overlap=$(jq -r '.depth_sample_overlap // ""' "$input_json")
depth_interval_overlap=$(jq -r '.depth_interval_overlap // ""' "$input_json")
depth_size_similarity=$(jq -r '.depth_size_similarity // ""' "$input_json")
depth_breakend_window=$(jq -r '.depth_breakend_window // ""' "$input_json")
mixed_sample_overlap=$(jq -r '.mixed_sample_overlap // ""' "$input_json")
mixed_interval_overlap=$(jq -r '.mixed_interval_overlap // ""' "$input_json")
mixed_size_similarity=$(jq -r '.mixed_size_similarity // ""' "$input_json")
mixed_breakend_window=$(jq -r '.mixed_breakend_window // ""' "$input_json")
pesr_sample_overlap=$(jq -r '.pesr_sample_overlap // ""' "$input_json")
pesr_interval_overlap=$(jq -r '.pesr_interval_overlap // ""' "$input_json")
pesr_size_similarity=$(jq -r '.pesr_size_similarity // ""' "$input_json")
pesr_breakend_window=$(jq -r '.pesr_breakend_window // ""' "$input_json")
insertion_length_summary_strategy=$(jq -r '.insertion_length_summary_strategy // ""' "$input_json")
breakpoint_summary_strategy=$(jq -r '.breakpoint_summary_strategy // ""' "$input_json")
alt_allele_summary_strategy=$(jq -r '.alt_allele_summary_strategy // ""' "$input_json")
additional_args=$(jq -r '.additional_args // ""' "$input_json")
java_mem_fraction=$(jq -r '.java_mem_fraction // 0.85' "$input_json")


# -------------------------------------------------------
# ======================= Command =======================
# -------------------------------------------------------


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



if (( ${#vcfs[@]} > 0 )); then
  printf "%s\n" "${vcfs[@]}" | awk '{print "-V "$0}' > arguments.txt
elif [[ -n "${vcfs_tar}" ]]; then
  mkdir vcfs
  tar xzf "${vcfs_tar}" -C vcfs/
  ls vcfs/*.vcf.gz | awk '{print "-V "$0}' > arguments.txt
else
  echo "ERROR: neither vcfs nor vcfs_tar was provided"
  exit 1
fi


default_no_call_arg=""
if [[ "${default_no_call}" == "true" ]]; then
  default_no_call_arg="--default-no-call"
fi

omit_members_arg=""
if [[ "${omit_members}" == "true" ]]; then
  omit_members_arg="--omit-members"
fi

enable_cnv_arg=""
if [[ "${enable_cnv}" == "true" ]]; then
  enable_cnv_arg="--enable-cnv"
fi

fast_mode_arg=""
if [[ "${fast_mode}" == "true" ]]; then
  fast_mode_arg="--fast-mode"
fi

java "-Xmx${JVM_MAX_MEM}" -jar /opt/gatk.jar SVCluster \
  --arguments_file arguments.txt \
  --output "${output_prefix}.vcf.gz" \
  --ploidy-table "${ploidy_table}" \
  --reference "${reference_fasta}" \
  ${contig:+-L $contig} \
  ${fast_mode_arg} \
  ${enable_cnv_arg} \
  ${omit_members_arg} \
  ${default_no_call_arg} \
  ${variant_prefix:+--variant-prefix $variant_prefix} \
  ${algorithm:+--algorithm $algorithm} \
  ${defrag_padding_fraction:+--defrag-padding-fraction $defrag_padding_fraction} \
  ${defrag_sample_overlap:+--defrag-sample-overlap $defrag_sample_overlap} \
  ${depth_sample_overlap:+--depth-sample-overlap $depth_sample_overlap} \
  ${depth_interval_overlap:+--depth-interval-overlap $depth_interval_overlap} \
  ${depth_size_similarity:+--depth-size-similarity $depth_size_similarity} \
  ${depth_breakend_window:+--depth-breakend-window $depth_breakend_window} \
  ${mixed_sample_overlap:+--mixed-sample-overlap $mixed_sample_overlap} \
  ${mixed_interval_overlap:+--mixed-interval-overlap $mixed_interval_overlap} \
  ${mixed_size_similarity:+--mixed-size-similarity $mixed_size_similarity} \
  ${mixed_breakend_window:+--mixed-breakend-window $mixed_breakend_window} \
  ${pesr_sample_overlap:+--pesr-sample-overlap $pesr_sample_overlap} \
  ${pesr_interval_overlap:+--pesr-interval-overlap $pesr_interval_overlap} \
  ${pesr_size_similarity:+--pesr-size-similarity $pesr_size_similarity} \
  ${pesr_breakend_window:+--pesr-breakend-window $pesr_breakend_window} \
  ${insertion_length_summary_strategy:+--insertion-length-summary-strategy $insertion_length_summary_strategy} \
  ${breakpoint_summary_strategy:+--breakpoint-summary-strategy $breakpoint_summary_strategy} \
  ${alt_allele_summary_strategy:+--alt-allele-summary-strategy $alt_allele_summary_strategy} \
  ${additional_args}

cluster_out_in_wd="$(realpath ${output_prefix}.vcf.gz)"
cluster_out_index_in_wd="$(realpath ${output_prefix}.vcf.gz.tbi)"


# -------------------------------------------------------
# ======================= Output ========================
# -------------------------------------------------------


cluster_out_in_output_dir="${output_dir}/$(basename "${cluster_out_in_wd}")"
mv "${cluster_out_in_wd}" "${cluster_out_in_output_dir}"
mv "${cluster_out_in_wd}.tbi" "${cluster_out_in_output_dir}.tbi"

outputs_json=$(jq -n \
  --arg out "${cluster_out_in_output_dir}" \
  --arg out_index "${cluster_out_in_output_dir}.tbi" \
  '{
     "out": $out,
     "out_index": $out_index
   }' \
)
echo "${outputs_json}" > "${output_json_filename}"

echo "Successfully finished SVCluster, output json filename: ${output_json_filename}"
