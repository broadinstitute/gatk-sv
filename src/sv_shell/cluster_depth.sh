#!/bin/bash

set -Exeuo pipefail

function ShardBedThenConvertToVcf()
{
  local _del_or_dup=$1
  local _bed=$2
  local -n vcfs_array=$3
  local -n vcfs_idx_array=$4

  cd "${working_dir}"

  # Task: ScatterCompressedBedOmitHeaders
  shard_bed_convert_vcf_working_dir=$(mktemp -d /wd_shard_bed_convert_vcf_XXXXXXXX)
  shard_bed_convert_vcf_working_dir="$(realpath ${shard_bed_convert_vcf_working_dir})"
  cd "${shard_bed_convert_vcf_working_dir}"
  n_digits=6
  mkdir out
  _prefix="${batch}.cluster_batch.depth.${_del_or_dup}.shard_"
  gunzip -c "${_bed}" \
    | awk '$0!~"#"' \
    | split -d -a "${n_digits}" -l "${records_per_bed_shard}" - "out/${_prefix}"
  for file in out/${_prefix}*; do
    mv $file $file.bed
    bgzip $file.bed
  done

  # Note that the following generates an empty array if it does not find a file matching the pattern,
  # similar to WDL equivalent, rather than failing.
  beds_array=($(ls out/${_prefix}*.bed.gz 2>/dev/null))

  vcfs_array=()
  vcfs_idx_array=()

  # Task: CNVBedToGatkVcf
  counter=0
  for _bed_file in "${beds_array[@]}"; do
    _out_prefix="${batch}.cluster_batch.depth.gatk_formatted.${_del_or_dup}.shard_${counter}"
    python /opt/sv-pipeline/scripts/convert_bed_to_gatk_vcf.py \
      --bed "${_bed_file}" \
      --out "${_out_prefix}.vcf.gz" \
      --samples "${sample_list}" \
      --contigs "${contig_list}" \
      --vid-prefix "${batch}_raw_depth_${_del_or_dup}_shard_${counter}_" \
      --ploidy-table "${ploidy_table}" \
      --fai "${reference_fasta_fai}"

    tabix "${_out_prefix}.vcf.gz"

    vcfs_array+=("$(realpath "${_out_prefix}.vcf.gz")")
    vcfs_idx_array+=("$(realpath "${_out_prefix}.vcf.gz.tbi")")

    counter=$((counter + 1))
  done
}


# -------------------------------------------------------
# ==================== Input & Setup ====================
# -------------------------------------------------------


input_json=${1}
output_json_filename=${2:-""}
output_dir=${3:-""}

input_json="$(realpath ${input_json})"

if [ -z "${output_dir}" ]; then
  output_dir=$(mktemp -d /output_cluster_depth_XXXXXXXX)
else
  mkdir -p "${output_dir}"
fi
output_dir="$(realpath ${output_dir})"

if [ -z "${output_json_filename}" ]; then
  output_json_filename="${output_dir}/output.json"
else
  output_json_filename="$(realpath ${output_json_filename})"
fi

working_dir=$(mktemp -d /wd_cluster_depth_XXXXXXXX)
working_dir="$(realpath ${working_dir})"
cd "${working_dir}"

batch=$(jq -r ".batch" "${input_json}")
del_bed=$(jq -r ".del_bed" "${input_json}")
dup_bed=$(jq -r ".dup_bed" "${input_json}")
records_per_bed_shard=$(jq -r ".records_per_bed_shard" "${input_json}")
sample_list=$(jq -r ".sample_list" "${input_json}")
contig_list=$(jq -r ".contig_list" "${input_json}")
ploidy_table=$(jq -r ".ploidy_table" "${input_json}")
contig_list=$(jq -r '.contig_list // ""' "${input_json}")
contig_subset_list=$(jq -r '.contig_subset_list // ""' "${input_json}")
depth_interval_overlap=$(jq -r ".depth_interval_overlap" "${input_json}")
reference_fasta=$(jq -r ".reference_fasta" "${input_json}")
reference_fasta_fai=$(jq -r ".reference_fasta_fai" "${input_json}")
reference_dict=$(jq -r ".reference_dict" "${input_json}")
java_mem_fraction=$(jq -r '.java_mem_fraction // ""' "${input_json}")

# -------------------------------------------------------
# ======================= Command =======================
# -------------------------------------------------------


# ScatterCompressedBedOmitHeaders and CNVBedToGatkVcf
# ---------------------------------------------------------------------------------------------------------------------

del_vcfs_array=()
del_vcfs_idx_array=()
ShardBedThenConvertToVcf "DEL" "${del_bed}" "del_vcfs_array" "del_vcfs_idx_array"

dup_vcfs_array=()
dup_vcfs_idx_array=()
ShardBedThenConvertToVcf "DUP" "${dup_bed}" "dup_vcfs_array" "dup_vcfs_idx_array"


# Per contig clustering
# ---------------------------------------------------------------------------------------------------------------------

vcfs=("${del_vcfs_array[@]}" "${dup_vcfs_array[@]}")
vcfs_string="$(printf '%s\n' "${vcfs[@]}" | jq -R . | jq -s . | jq -c .)"

contigs_list_file="${contig_subset_list:-$contig_list}"
contigs=($(cat "${contigs_list_file}"))

for contig in "${contigs[@]}"; do
  echo "Starting to cluster ${contig}."

  cd "${working_dir}"
  cluster_contig_output_dir=$(mktemp -d /output_cluster_contig_XXXXXXXX)
  cluster_contig_output_dir="$(realpath ${cluster_contig_output_dir})"
  contig_cluster_inputs_json="$(realpath "${cluster_contig_output_dir}/contig_cluster_inputs.json")"
  contig_cluster_output_json="$(realpath "${cluster_contig_output_dir}/contig_cluster_output.json")"

  jq -n \
    --argjson vcfs "${vcfs_string}" \
    --arg ploidy_table "${ploidy_table}" \
    --arg output_prefix "${batch}.cluster_batch.depth.${contig}.clustered" \
    --arg contig "${contig}" \
    --argjson fast_mode true \
    --arg algorithm "SINGLE_LINKAGE" \
    --argjson depth_sample_overlap 0 \
    --argjson depth_interval_overlap "${depth_interval_overlap}" \
    --argjson depth_breakend_window 10000000 \
    --arg reference_fasta "${reference_fasta}" \
    --arg reference_fasta_fai "${reference_fasta_fai}" \
    --arg reference_dict "${reference_dict}" \
    --arg java_mem_fraction "${java_mem_fraction}" \
    --arg variant_prefix "${batch}_depth_${contig}_" \
    '{
        "vcfs": $vcfs,
        "ploidy_table": $ploidy_table,
        "output_prefix": $output_prefix,
        "contig": $contig,
        "fast_mode": $fast_mode,
        "algorithm": $algorithm,
        "depth_sample_overlap": $depth_sample_overlap,
        "depth_interval_overlap": $depth_interval_overlap,
        "depth_breakend_window": $depth_breakend_window,
        "reference_fasta": $reference_fasta,
        "reference_fasta_fai": $reference_fasta_fai,
        "reference_dict": $reference_dict,
        "java_mem_fraction": $java_mem_fraction,
        "variant_prefix": $variant_prefix
    }' > "${contig_cluster_inputs_json}"

    bash /opt/sv_shell/sv_cluster.sh "${contig_cluster_inputs_json}" "${contig_cluster_output_json}" "${cluster_contig_output_dir}"

    echo "Finished clustering ${contig}; output json: ${cluster_contig_output_dir}"
done
