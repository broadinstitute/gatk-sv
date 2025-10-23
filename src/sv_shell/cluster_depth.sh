#!/bin/bash

set -Exeuo pipefail

function ScatterCompressedBedOmitHeaders()
{
  local _bed=$1
  local _prefix=$2
  local -n _out_array=$3

  scatter_compressed_bed_working_dir=$(mktemp -d /wd_scatter_compressed_bed_XXXXXXXX)
  scatter_compressed_bed_working_dir="$(realpath ${scatter_compressed_bed_working_dir})"
  cd "${scatter_compressed_bed_working_dir}"
  n_digits=6
  mkdir out
  gunzip -c "${_bed}" \
    | awk '$0!~"#"' \
    | split -d -a "${n_digits}" -l "${records_per_bed_shard}" - "out/${_prefix}"
  for file in out/${_prefix}*; do
    mv $file $file.bed
    bgzip $file.bed
  done

  # Note that the following generates an empty array if it does not find a file matching the pattern,
  # similar to WDL equivalent, rather than failing.
  _out_array=($(ls out/${_prefix}*.bed.gz 2>/dev/null))
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


# -------------------------------------------------------
# ======================= Command =======================
# -------------------------------------------------------


# ScatterCompressedBedOmitHeaders
# ---------------------------------------------------------------------------------------------------------------------

del_depth_bed_list=()
ScatterCompressedBedOmitHeaders "${del_bed}" "${batch}.cluster_batch.depth.del.shard_" "del_depth_bed_list"

dup_depth_bed_list=()
ScatterCompressedBedOmitHeaders "${dup_bed}" "${batch}.cluster_batch.depth.dup.shard_" "dup_depth_bed_list"


