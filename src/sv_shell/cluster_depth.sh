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
records_per_bed_shard=$(jq -r ".records_per_bed_shard" "${input_json}")


# -------------------------------------------------------
# ======================= Command =======================
# -------------------------------------------------------


# ScatterCompressedBedOmitHeaders
# ---------------------------------------------------------------------------------------------------------------------

scatter_compressed_bed_working_dir=$(mktemp -d /wd_scatter_compressed_bed_XXXXXXXX)
scatter_compressed_bed_working_dir="$(realpath ${scatter_compressed_bed_working_dir})"
cd "${scatter_compressed_bed_working_dir}"
n_digits=6
del_depth_prefix="${batch}.cluster_batch.depth.del.shard_"
mkdir out
gunzip -c "${del_bed}" \
    | awk '$0!~"#"' \
    | split -d -a "${n_digits}" -l "${records_per_bed_shard}" - "out/${del_depth_prefix}"
for file in out/${del_depth_prefix}*; do
    mv $file $file.bed
    bgzip $file.bed
done
