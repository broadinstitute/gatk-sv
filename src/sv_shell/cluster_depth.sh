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
reference_fasta_fai=$(jq -r ".reference_fasta_fai" "${input_json}")

# -------------------------------------------------------
# ======================= Command =======================
# -------------------------------------------------------


# ScatterCompressedBedOmitHeaders
# ---------------------------------------------------------------------------------------------------------------------

del_vcfs_array=()
del_vcfs_idx_array=()
ShardBedThenConvertToVcf "DEL" "${del_bed}" "del_vcfs_array" "del_vcfs_idx_array"

dup_vcfs_array=()
dup_vcfs_idx_array=()
ShardBedThenConvertToVcf "DUP" "${dup_bed}" "dup_vcfs_array" "dup_vcfs_idx_array"
