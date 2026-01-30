#!/bin/bash

# note that the WDL implementation of this script has gcnv support,
# this script has temporarily dropped any gcnv related logic, as
# we're postponing implementing gcnv.

set -Exeuo pipefail

function MergeSample() {
  local _sample_id=$1
  local _std_cnmops=$2
  local _large_cnmops=$3
  local _max_dist="${4:-"0.25"}"

  echo "----------- Starting MergeSample -------------"
  echo "_std_cnmops: ${_std_cnmops}"
  echo "_large_cnmops: ${_large_cnmops}"
  echo "_max_dist: ${_max_dist}"
  echo "_sample_id: ${_sample_id}"
  echo "----------------------------------------------"

  zcat "${_std_cnmops}" "${_large_cnmops}" | awk -F "\t" -v OFS="\t" -v sample="${_sample_id}" '{if ($5==sample) print}' > cnmops.cnv
  sort -k1,1V -k2,2n cnmops.cnv > "${_sample_id}.bed"
  bedtools merge -i "${_sample_id}.bed" -d 0 -c 4,5,6,7 -o distinct > "${_sample_id}.merged.bed"
  /opt/sv-pipeline/00_preprocessing/scripts/defragment_cnvs.py \
    --max-dist "${_max_dist}" \
    "${_sample_id}.merged.bed" \
    "${_sample_id}.merged.defrag.bed"

  sort -k1,1V -k2,2n "${_sample_id}.merged.defrag.bed" > "${_sample_id}.merged.defrag.sorted.bed"

  echo "----------- Finished MergeSample -------------"
}


function MergeSet() {
  declare -n _beds=$1
  local _svtype=$2
  local _batch=$3

  echo "----------- Starting MergeSet ----------------"
  echo "_beds: ${_beds[@]}"
  echo "_svtype: ${_svtype}"
  echo "_batch: ${_batch}"
  echo "----------------------------------------------"

  zcat -f "${_beds[@]}" \
    | sort -k1,1V -k2,2n \
    | awk -v OFS="\t" -v svtype="${_svtype}" -v batch="${_batch}" '{$4=batch"_"svtype"_"NR; print}' \
    | cat <(echo -e "#chr\\tstart\\tend\\tname\\tsample\\tsvtype\\tsources") - \
    | bgzip -c > "${_batch}.${_svtype}.bed.gz";
  tabix -p bed "${_batch}.${_svtype}.bed.gz"

echo "----------- Finished MergeSet ----------------"
}

# -------------------------------------------------------
# ==================== Input & Setup ====================
# -------------------------------------------------------


input_json=${1}
output_json_filename=${2-""}
output_dir=${3:-""}

input_json="$(realpath ${input_json})"

if [ -z "${output_dir}" ]; then
  output_dir=$(mktemp -d ${SV_SHELL_BASE_DIR}/output_merge_depth_XXXXXXXX)
else
  mkdir -p "${output_dir}"
fi
output_dir="$(realpath ${output_dir})"

if [ -z "${output_json_filename}" ]; then
  output_json_filename="${output_dir}/output.json"
else
  output_json_filename="$(realpath ${output_json_filename})"
fi

working_dir=$(mktemp -d ${SV_SHELL_BASE_DIR}/wd_merge_depth_XXXXXXXX)
working_dir="$(realpath ${working_dir})"
cd "${working_dir}"
echo "Merge Depth Working directory: ${working_dir}"

batch=($(jq -r '.batch' "$input_json"))
samples=($(jq -r '.samples[]' "$input_json"))
std_cnmops_del=$(jq -r '.std_cnmops_del' "$input_json")
std_cnmops_dup=$(jq -r '.std_cnmops_dup' "$input_json")
large_cnmops_del=$(jq -r '.large_cnmops_del' "$input_json")
large_cnmops_dup=$(jq -r '.large_cnmops_dup' "$input_json")
defragment_max_dist=$(jq -r '.defragment_max_dist // ""' "$input_json")


# -------------------------------------------------------
# ======================= Command =======================
# -------------------------------------------------------

# Merge Sample DEL
merged_del_beds=()
for sample in "${samples[@]}"; do
  working_dir=$(mktemp -d ${SV_SHELL_BASE_DIR}/wd_merge_sample_del_XXXXXXXX)
  working_dir="$(realpath ${working_dir})"
  cd "${working_dir}"

  MergeSample "${sample}" "${std_cnmops_del}" "${large_cnmops_del}" "${defragment_max_dist}"
  merged_del_beds+=("${working_dir}/${sample}.merged.defrag.sorted.bed")
done


# Merge Sample DUP
merged_dup_bed=()
for sample in "${samples[@]}"; do
  working_dir=$(mktemp -d ${SV_SHELL_BASE_DIR}/wd_merge_sample_dup_XXXXXXXX)
  working_dir="$(realpath ${working_dir})"
  cd "${working_dir}"

  MergeSample "${sample}" "${std_cnmops_dup}" "${large_cnmops_dup}" "${defragment_max_dist}"
  merged_dup_beds+=("${working_dir}/${sample}.merged.defrag.sorted.bed")
done


# Merge DEL
working_dir=$(mktemp -d ${SV_SHELL_BASE_DIR}/wd_merge_del_XXXXXXXX)
working_dir="$(realpath ${working_dir})"
cd "${working_dir}"
# note that "merged_del_beds" in the following is passing the array, not a string.
MergeSet "merged_del_beds" "DEL" "${batch}"
merge_set_del="${working_dir}/${batch}.DEL.bed.gz"
merge_set_del_idx="${working_dir}/${batch}.DEL.bed.gz.tbi"


# Merge DUP
working_dir=$(mktemp -d ${SV_SHELL_BASE_DIR}/wd_merge_dup_XXXXXXXX)
working_dir="$(realpath ${working_dir})"
cd "${working_dir}"
# note that "merged_del_beds" in the following is passing the array, not a string.
MergeSet "merged_dup_beds" "DUP" "${batch}"
merge_set_dup="${working_dir}/${batch}.DUP.bed.gz"
merge_set_dup_idx="${working_dir}/${batch}.DUP.bed.gz.tbi"


# -------------------------------------------------------
# ======================= Output ========================
# -------------------------------------------------------

del_out="${output_dir}/$(basename "${merge_set_del}")"
del_index_out="${output_dir}/$(basename "${merge_set_del_idx}")"

dup_out="${output_dir}/$(basename "${merge_set_dup}")"
dup_index_out="${output_dir}/$(basename "${merge_set_dup_idx}")"

mv "${merge_set_del}" "${del_out}"
mv "${merge_set_del_idx}" "${del_index_out}"
mv "${merge_set_dup}" "${dup_out}"
mv "${merge_set_dup_idx}" "${dup_index_out}"

outputs_json=$(jq -n \
  --arg del "${del_out}" \
  --arg del_index "${del_index_out}" \
  --arg dup "${dup_out}" \
  --arg dup_index "${dup_index_out}" \
  '{Del: $del, Del_idx: $del_index, Dup: $dup, Dup_idx: $dup_index}' )
echo "${outputs_json}" > "${output_json_filename}"

echo "Finished merge depth successfully, output json filename: ${output_json_filename}"
