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
  output_dir=$(mktemp -d /output_join_raw_calls_XXXXXXXX)
else
  mkdir -p "${output_dir}"
fi
output_dir="$(realpath ${output_dir})"

if [ -z "${output_json_filename}" ]; then
  output_json_filename="${output_dir}/output.json"
else
  output_json_filename="$(realpath ${output_json_filename})"
fi

working_dir=$(mktemp -d /wd_join_raw_calls_XXXXXXXX)
working_dir="$(realpath ${working_dir})"
cd "${working_dir}"

prefix=$(jq -r ".prefix" "${input_json}")

clustered_depth_vcf=$(jq -r '.clustered_depth_vcf // ""' "${input_json}")
clustered_depth_vcf_index=$(jq -r '.clustered_depth_vcf_index // ""' "${input_json}")

clustered_dragen_vcf=$(jq -r '.clustered_dragen_vcf // ""' "${input_json}")
clustered_dragen_vcf_index=$(jq -r '.clustered_dragen_vcf_index // ""' "${input_json}")

clustered_manta_vcf=$(jq -r '.clustered_manta_vcf // ""' "${input_json}")
clustered_manta_vcf_index=$(jq -r '.clustered_manta_vcf_index // ""' "${input_json}")

clustered_scramble_vcf=$(jq -r '.clustered_scramble_vcf // ""' "${input_json}")
clustered_scramble_vcf_index=$(jq -r '.clustered_scramble_vcf_index // ""' "${input_json}")

clustered_wham_vcf=$(jq -r '.clustered_wham_vcf // ""' "${input_json}")
clustered_wham_vcf_index=$(jq -r '.clustered_wham_vcf_index // ""' "${input_json}")

ped_file=$(jq -r ".ped_file" "${input_json}")
contig_list=$(jq -r ".contig_list" "${input_json}")

reference_fasta=$(jq -r ".reference_fasta" "${input_json}")
reference_fasta_fai=$(jq -r ".reference_fasta_fai" "${input_json}")
reference_dict=$(jq -r ".reference_dict" "${input_json}")
java_mem_fraction=$(jq -r '.java_mem_fraction // "null"' "${input_json}")

# -------------------------------------------------------
# ======================= Command =======================
# -------------------------------------------------------


# ConcatInputVcfs
# ---------------------------------------------------------------------------------------------------------------------

vcf_files=(
  "${clustered_depth_vcf}"
  "${clustered_dragen_vcf}"
  "${clustered_manta_vcf}"
  "${clustered_scramble_vcf}"
  "${clustered_wham_vcf}"
)

vcf_files_list="vcfs_list.txt"
touch "${vcf_files_list}"
for vcf_file_path in "${vcf_files[@]}"; do
  if [[ -n "${vcf_file_path}" ]]; then
    echo "${vcf_file_path}" >> "${vcf_files_list}"
  fi
done

concat_vcf="${prefix}.join_raw_calls.concat.vcf.gz"
concat_vcf_idx="${prefix}.join_raw_calls.concat.vcf.gz.tbi"
VCFS="~{write_lines(vcfs)}"

bcftools concat --no-version --allow-overlaps -Oz --file-list "${vcf_files_list}" > "${concat_vcf}"

tabix "${concat_vcf}"


# CreatePloidyTableFromPed
# ---------------------------------------------------------------------------------------------------------------------

CreatePloidyTableFromPed_out="$(realpath "${prefix}.ploidy.tsv")"
python /opt/sv-pipeline/scripts/ploidy_table_from_ped.py \
  --ped "${ped_file}" \
  --out "${CreatePloidyTableFromPed_out}" \
  --contigs "${contig_list}"


# FormatVcfForGatk
# ---------------------------------------------------------------------------------------------------------------------

FormatVcfForGatk_gatk_formatted_vcf="$(realpath "${prefix}.join_raw_calls.gatk_formatted.vcf.gz")"
FormatVcfForGatk_gatk_formatted_vcf_index="$(realpath "${prefix}.join_raw_calls.gatk_formatted.vcf.gz.tbi")"

python /opt/sv-pipeline/scripts/format_svtk_vcf_for_gatk.py \
  --vcf "${concat_vcf}" \
  --out "${FormatVcfForGatk_gatk_formatted_vcf}" \
  --ploidy-table "${CreatePloidyTableFromPed_out}"

tabix "${FormatVcfForGatk_gatk_formatted_vcf}"


# SVCluster
# ---------------------------------------------------------------------------------------------------------------------

sv_cluster_output_dir=$(mktemp -d "/output_sv_cluster_XXXXXXXX")
sv_cluster_output_dir="$(realpath ${sv_cluster_output_dir})"
sv_cluster_inputs_json="$(realpath "${sv_cluster_output_dir}/sv_cluster_inputs.json")"
sv_cluster_output_json="$(realpath "${sv_cluster_output_dir}/sv_cluster_output.json")"

sv_cluster_wd_dir=$(mktemp -d "/wd_sv_cluster_XXXXXXXX")
sv_cluster_wd_dir="$(realpath ${sv_cluster_wd_dir})"

jq -n \
  --arg vcfs "${FormatVcfForGatk_gatk_formatted_vcf}" \
  --arg ploidy_table "${CreatePloidyTableFromPed_out}" \
  --arg output_prefix "${prefix}.join_raw_calls" \
  --argjson fast_mode true \
  --arg algorithm "SINGLE_LINKAGE" \
  --argjson pesr_sample_overlap 0 \
  --argjson mixed_sample_overlap 0 \
  --argjson depth_sample_overlap 0 \
  --arg reference_fasta "${reference_fasta}" \
  --arg reference_fasta_fai "${reference_fasta_fai}" \
  --arg reference_dict "${reference_dict}" \
  --argjson java_mem_fraction "${java_mem_fraction}" \
  --arg variant_prefix "${prefix}_" \
  '{
      "vcfs": [$vcfs],
      "ploidy_table": $ploidy_table,
      "output_prefix": $output_prefix,
      "fast_mode": $fast_mode,
      "algorithm": $algorithm,
      "pesr_sample_overlap": $pesr_sample_overlap,
      "mixed_sample_overlap": $mixed_sample_overlap,
      "depth_sample_overlap": $depth_sample_overlap,
      "reference_fasta": $reference_fasta,
      "reference_fasta_fai": $reference_fasta_fai,
      "reference_dict": $reference_dict,
      "java_mem_fraction": $java_mem_fraction,
      "variant_prefix": $variant_prefix
  }' > "${sv_cluster_inputs_json}"

  bash /opt/sv_shell/sv_cluster.sh "${sv_cluster_inputs_json}" "${sv_cluster_output_json}" "${sv_cluster_output_dir}"

  echo "Finished SV clustering; output json: ${sv_cluster_output_dir}"

  sv_cluster_vcf_out=$(jq -r ".out" "${sv_cluster_output_json}")
