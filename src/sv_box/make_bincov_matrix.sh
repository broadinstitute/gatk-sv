#!/bin/bash

set -Eeuo pipefail

inputs_json=${1}
output_dir=${2:-""}


working_dir=$(mktemp -d wd_make_bincov_matrix_XXXXXXXX)
working_dir="$(realpath ${working_dir})"
output_dir=$(mktemp -d output_make_bincov_matrix_XXXXXXXX)
output_dir="$(realpath ${output_dir})"
cd "${working_dir}"

bin_file_name="${working_dir}/locs.bed.gz"


count_file=$(jq -r ".count_file" "${inputs_json}")
binsize=$(jq -r ".binsize" "${inputs_json}")
bincov_matrix_samples=$(jq -r ".bincov_matrix_samples" "${inputs_json}")


echo "${count_file}"


# ---------------------------
# --------- SetBins ---------
# ---------------------------

# kill the dictionary | kill the header | adjust to bed format: 0-based half-open intervals
zcat "${count_file}" \
  | sed '/^@/d' \
  | sed '/^CONTIG	START	END	COUNT$/d' \
  | sed '/^#/d' \
  | awk -v x="1" 'BEGIN{OFS="\t"}{$2=$2-x; print $1,$2,$3}' > tmp_locs

# determine bin size, and drop all bins not exactly equal to this size
if [[ "${binsize}" == "null" ]]; then
  # use the most common bin size from the bins
  binsize=$(
    sed -n '1,1000p' tmp_locs | awk '{ print $3-$2 }' \
    | sort | uniq -c | sort -nrk1,1 \
    | sed -n '1p' | awk '{ print $2 }'
  )
fi

# store binsize where cromwell can read it
# TODO: the folloiwng is not needed, just output the defined bin size.
#echo "${binsize}" > ~{binsize_output_file_name}

# write final bed file with header, and compress it
awk -v FS="\t" -v b="${binsize}" 'BEGIN{ print "#Chr\tStart\tEnd" } { if ($3-$2==b) print $0 }' tmp_locs \
    | bgzip -c \
    > "${bin_file_name}"

# if bincov_matrix_samples was passed, convert to tab-separated string
bincov_header_file_name="bincov_header_file.tsv"
if [[ "${bincov_matrix_samples}" != "null" ]]; then
  echo "${bincov_matrix_samples[@]}" | jq -r '. | @tsv' > "${bincov_header_file_name}"
else
  touch "${bincov_header_file_name}"
fi







#bin_locs="${output_dir}/locs.bed.gz"
#mv bin_file_name bin_locs
#bincov_matrix_header_file="${output_dir}/bincov_header_file.tsv"




outputs_filename="${output_dir}/outputs.json"
#outputs_json=$(jq -n \
#  --arg bin_locs "${bin_locs}" \
#
#  '{counts: $c}' )
#echo "${outputs_json}" > "${outputs_filename}"
#cp "${outputs_filename}" "${outputs_json_filename}"

