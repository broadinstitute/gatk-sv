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
  output_dir=$(mktemp -d ${SV_SHELL_BASE_DIR}/output_cluster_batch_XXXXXXXX)
else
  mkdir -p "${output_dir}"
fi
output_dir="$(realpath ${output_dir})"

if [ -z "${output_json_filename}" ]; then
  output_json_filename="${output_dir}/output.json"
else
  output_json_filename="$(realpath ${output_json_filename})"
fi

working_dir=$(mktemp -d ${SV_SHELL_BASE_DIR}/wd_cluster_batch_XXXXXXXX)
working_dir="$(realpath ${working_dir})"
cd "${working_dir}"

batch=$(jq -r ".batch" "${input_json}")
contig_list=$(jq -r ".contig_list" "${input_json}")
ped_file=$(jq -r ".ped_file" "${input_json}")
del_bed=$(jq -r ".del_bed" "${input_json}")
dup_bed=$(jq -r ".dup_bed" "${input_json}")
dragen_vcf_tar=$(jq -r '.dragen_vcf_tar // ""' "${input_json}")
manta_vcf_tar=$(jq -r '.manta_vcf_tar // ""' "${input_json}")
melt_vcf_tar=$(jq -r '.melt_vcf_tar // ""' "${input_json}")
wham_vcf_tar=$(jq -r '.wham_vcf_tar // ""' "${input_json}")
scramble_vcf_tar=$(jq -r '.scramble_vcf_tar // ""' "${input_json}")
chr_x=$(jq -r '.chr_x // ""' "${input_json}")
chr_y=$(jq -r '.chr_y // ""' "${input_json}")
retain_female_chr_y=$(jq -r '.retain_female_chr_y' "${input_json}")
pesr_min_size=$(jq -r '.pesr_min_size // "50"' "${input_json}")
pesr_exclude_intervals=$(jq -r '.pesr_exclude_intervals' "${input_json}")
contig_subset_list=$(jq -r '.contig_subset_list // ""' "${input_json}")
pesr_interval_overlap=$(jq -r '.pesr_interval_overlap' "${input_json}")
pesr_breakend_window=$(jq -r '.pesr_breakend_window' "${input_json}")
pesr_clustering_algorithm=$(jq -r '.pesr_clustering_algorithm // ""' "${input_json}")
reference_fasta=$(jq -r '.reference_fasta' "${input_json}")
reference_fasta_fai=$(jq -r '.reference_fasta_fai' "${input_json}")
reference_dict=$(jq -r '.reference_dict' "${input_json}")
depth_records_per_bed_shard=$(jq -r '.depth_records_per_bed_shard // "1000000"' "${input_json}")
depth_exclude_intervals=$(jq -r '.depth_exclude_intervals' "${input_json}")
depth_exclude_overlap_fraction=$(jq -r '.depth_exclude_overlap_fraction' "${input_json}")
depth_clustering_algorithm=$(jq -r '.depth_clustering_algorithm' "${input_json}")
depth_interval_overlap=$(jq -r '.depth_interval_overlap' "${input_json}")


# -------------------------------------------------------
# ======================= Command =======================
# -------------------------------------------------------


# GetSampleIdsFromVcfTar
# ---------------------------------------------------------------------------------------------------------------------
if [ -n "${dragen_vcf_tar}" ]; then
  vcf_tar="${dragen_vcf_tar}"
elif [ -n "${manta_vcf_tar}" ]; then
  vcf_tar="${manta_vcf_tar}"
elif [ -n "${wham_vcf_tar}" ]; then
  vcf_tar="${wham_vcf_tar}"
else
  echo "Error: No VCF tar file provided" >&2
  exit 1
fi
mkdir vcfs
tar xzf "${vcf_tar}" -C vcfs/
ls vcfs/*.vcf.gz | xargs -n1 bcftools query -l | sort -u > "${batch}.samples.txt"

GetSampleIdsFromVcfTar_out_file="$(realpath "${batch}.samples.txt")"



# CreatePloidyTableFromPed
# ---------------------------------------------------------------------------------------------------------------------

if [[ "${retain_female_chr_y}" == "true" ]]; then
  ploidy_table_output_file="${batch}.ploidy.FEMALE_chrY_1.tsv"
else
  ploidy_table_output_file="${batch}.ploidy.tsv"
fi
ploidy_table_output_file="$(realpath ${ploidy_table_output_file})"

# parameter expansion: use value after :+ if chr_x provided, empty string otherwise
python /opt/sv-pipeline/scripts/ploidy_table_from_ped.py \
  --ped "${ped_file}" \
  --out tmp.tsv \
  --contigs "${contig_list}" \
  ${chr_x:+--chr-x $chr_x} \
  ${chr_y:+--chr-y $chr_y}

# TODO : For now we retain female Y genotypes for metric generation
if [ "${retain_female_chr_y}" = true ]; then
    sed -e 's/\t0/\t1/g' tmp.tsv > "${ploidy_table_output_file}"
else
    mv tmp.tsv "${ploidy_table_output_file}"
fi


# Cluster PESR Dragen
# ---------------------------------------------------------------------------------------------------------------------

dragen_pesr_clustered_vcf=""
dragen_pesr_clustered_vcf_idx=""
if [ -n "${dragen_vcf_tar}" ]; then
  echo "Running PE/SR clustering on Dragen VCF."

  cd "${working_dir}"
  dragen_pesr_output_dir=$(mktemp -d "${SV_SHELL_BASE_DIR}/output_dragen_pesr_XXXXXXXX")
  dragen_pesr_output_dir="$(realpath ${dragen_pesr_output_dir})"
  dragen_pesr_inputs_json="$(realpath "${dragen_pesr_output_dir}/dragen_pesr_inputs.json")"
  dragen_pesr_output_json="$(realpath "${dragen_pesr_output_dir}/dragen_pesr_output.json")"

  jq -n \
    --arg vcf_tar "${dragen_vcf_tar}" \
    --arg ploidy_table "${ploidy_table_output_file}" \
    --arg batch "${batch}" \
    --arg caller "dragen" \
    --argjson min_size "${pesr_min_size}" \
    --arg exclude_intervals "${pesr_exclude_intervals}" \
    --arg contig_list "${contig_list}" \
    --argjson contig_subset_list "${contig_subset_list:-null}" \
    --argjson pesr_interval_overlap ${pesr_interval_overlap} \
    --argjson pesr_breakend_window ${pesr_breakend_window} \
    --arg clustering_algorithm "${pesr_clustering_algorithm}" \
    --arg reference_fasta "${reference_fasta}" \
    --arg reference_fasta_fai "${reference_fasta_fai}" \
    --arg reference_dict "${reference_dict}" \
    --argjson sort_vcf_list false \
    --argjson vcfs '[]' \
    '{
          "batch": $batch,
          "caller": $caller,
          "clustering_algorithm": $clustering_algorithm,
          "contig_list": $contig_list,
          "contig_subset_list": $contig_subset_list,
          "exclude_intervals": $exclude_intervals,
          "min_size": $min_size,
          "pesr_breakend_window": $pesr_breakend_window,
          "pesr_interval_overlap": $pesr_interval_overlap,
          "ploidy_table": $ploidy_table,
          "reference_dict": $reference_dict,
          "reference_fasta": $reference_fasta,
          "reference_fasta_fai": $reference_fasta_fai,
          "sort_vcf_list": $sort_vcf_list,
          "vcf_tar": $vcf_tar,
          "vcfs": $vcfs
      }' > "${dragen_pesr_inputs_json}"

    bash /opt/sv_shell/cluster_pesr.sh "${dragen_pesr_inputs_json}" "${dragen_pesr_output_json}"

    dragen_pesr_clustered_vcf=$(jq -r '.clustered_vcf' "${dragen_pesr_output_json}")
    dragen_pesr_clustered_vcf_idx=$(jq -r '.clustered_vcf_index' "${dragen_pesr_output_json}")

    echo "Finished running PE/SR clustering on Dragen VCF."
fi


# Cluster PESR Manta
# ---------------------------------------------------------------------------------------------------------------------

manta_pesr_clustered_vcf=""
manta_pesr_clustered_vcf_idx=""
if [ -n "${manta_vcf_tar}" ]; then
  echo "Running PE/SR clustering on Manta VCF."

  cd "${working_dir}"
  manta_pesr_output_dir=$(mktemp -d "${SV_SHELL_BASE_DIR}/output_manta_pesr_XXXXXXXX")
  manta_pesr_output_dir="$(realpath ${manta_pesr_output_dir})"
  manta_pesr_inputs_json="$(realpath "${manta_pesr_output_dir}/manta_pesr_inputs.json")"
  manta_pesr_output_json="$(realpath "${manta_pesr_output_dir}/manta_pesr_output.json")"

  jq -n \
    --arg vcf_tar "${manta_vcf_tar}" \
    --arg ploidy_table "${ploidy_table_output_file}" \
    --arg batch "${batch}" \
    --arg caller "manta" \
    --argjson min_size "${pesr_min_size}" \
    --arg exclude_intervals "${pesr_exclude_intervals}" \
    --arg contig_list "${contig_list}" \
    --argjson contig_subset_list "${contig_subset_list:-null}" \
    --argjson pesr_interval_overlap ${pesr_interval_overlap} \
    --argjson pesr_breakend_window ${pesr_breakend_window} \
    --arg clustering_algorithm "${pesr_clustering_algorithm}" \
    --arg reference_fasta "${reference_fasta}" \
    --arg reference_fasta_fai "${reference_fasta_fai}" \
    --arg reference_dict "${reference_dict}" \
    --argjson sort_vcf_list false \
    --argjson vcfs '[]' \
    '{
          "batch": $batch,
          "caller": $caller,
          "clustering_algorithm": $clustering_algorithm,
          "contig_list": $contig_list,
          "contig_subset_list": $contig_subset_list,
          "exclude_intervals": $exclude_intervals,
          "min_size": $min_size,
          "pesr_breakend_window": $pesr_breakend_window,
          "pesr_interval_overlap": $pesr_interval_overlap,
          "ploidy_table": $ploidy_table,
          "reference_dict": $reference_dict,
          "reference_fasta": $reference_fasta,
          "reference_fasta_fai": $reference_fasta_fai,
          "sort_vcf_list": $sort_vcf_list,
          "vcf_tar": $vcf_tar,
          "vcfs": $vcfs
      }' > "${manta_pesr_inputs_json}"

    bash /opt/sv_shell/cluster_pesr.sh "${manta_pesr_inputs_json}" "${manta_pesr_output_json}"

    manta_pesr_clustered_vcf=$(jq -r '.clustered_vcf' "${manta_pesr_output_json}")
    manta_pesr_clustered_vcf_idx=$(jq -r '.clustered_vcf_index' "${manta_pesr_output_json}")

    echo "Finished running PE/SR clustering on Manta VCF."
fi


# Cluster PESR Wham
# ---------------------------------------------------------------------------------------------------------------------

wham_pesr_clustered_vcf=""
wham_pesr_clustered_vcf_idx=""
if [ -n "${wham_vcf_tar}" ]; then
  echo "Running PE/SR clustering on Wham VCF."

  cd "${working_dir}"
  wham_pesr_output_dir=$(mktemp -d "${SV_SHELL_BASE_DIR}/output_wham_pesr_XXXXXXXX")
  wham_pesr_output_dir="$(realpath ${wham_pesr_output_dir})"
  wham_pesr_inputs_json="$(realpath "${wham_pesr_output_dir}/wham_pesr_inputs.json")"
  wham_pesr_output_json="$(realpath "${wham_pesr_output_dir}/wham_pesr_output.json")"

  jq -n \
    --arg vcf_tar "${wham_vcf_tar}" \
    --arg ploidy_table "${ploidy_table_output_file}" \
    --arg batch "${batch}" \
    --arg caller "wham" \
    --argjson min_size "${pesr_min_size}" \
    --arg exclude_intervals "${pesr_exclude_intervals}" \
    --arg contig_list "${contig_list}" \
    --argjson contig_subset_list "${contig_subset_list:-null}" \
    --argjson pesr_interval_overlap ${pesr_interval_overlap} \
    --argjson pesr_breakend_window ${pesr_breakend_window} \
    --arg clustering_algorithm "${pesr_clustering_algorithm}" \
    --arg reference_fasta "${reference_fasta}" \
    --arg reference_fasta_fai "${reference_fasta_fai}" \
    --arg reference_dict "${reference_dict}" \
    --argjson sort_vcf_list false \
    --argjson vcfs '[]' \
    '{
          "batch": $batch,
          "caller": $caller,
          "clustering_algorithm": $clustering_algorithm,
          "contig_list": $contig_list,
          "contig_subset_list": $contig_subset_list,
          "exclude_intervals": $exclude_intervals,
          "min_size": $min_size,
          "pesr_breakend_window": $pesr_breakend_window,
          "pesr_interval_overlap": $pesr_interval_overlap,
          "ploidy_table": $ploidy_table,
          "reference_dict": $reference_dict,
          "reference_fasta": $reference_fasta,
          "reference_fasta_fai": $reference_fasta_fai,
          "sort_vcf_list": $sort_vcf_list,
          "vcf_tar": $vcf_tar,
          "vcfs": $vcfs
      }' > "${wham_pesr_inputs_json}"

    bash /opt/sv_shell/cluster_pesr.sh "${wham_pesr_inputs_json}" "${wham_pesr_output_json}"

    wham_pesr_clustered_vcf=$(jq -r '.clustered_vcf' "${wham_pesr_output_json}")
    wham_pesr_clustered_vcf_idx=$(jq -r '.clustered_vcf_index' "${wham_pesr_output_json}")

    echo "Finished running PE/SR clustering on Wham VCF."
fi


# Cluster PESR Scramble
# ---------------------------------------------------------------------------------------------------------------------

scramble_pesr_clustered_vcf=""
scramble_pesr_clustered_vcf_idx=""
if [ -n "${scramble_vcf_tar}" ]; then
  echo "Running PE/SR clustering on Scramble VCF."

  cd "${working_dir}"
  scramble_pesr_output_dir=$(mktemp -d "${SV_SHELL_BASE_DIR}/output_scramble_pesr_XXXXXXXX")
  scramble_pesr_output_dir="$(realpath ${scramble_pesr_output_dir})"
  scramble_pesr_inputs_json="$(realpath "${scramble_pesr_output_dir}/scramble_pesr_inputs.json")"
  scramble_pesr_output_json="$(realpath "${scramble_pesr_output_dir}/scramble_pesr_output.json")"

  jq -n \
    --arg vcf_tar "${scramble_vcf_tar}" \
    --arg ploidy_table "${ploidy_table_output_file}" \
    --arg batch "${batch}" \
    --arg caller "scramble" \
    --argjson min_size "${pesr_min_size}" \
    --arg exclude_intervals "${pesr_exclude_intervals}" \
    --arg contig_list "${contig_list}" \
    --argjson contig_subset_list "${contig_subset_list:-null}" \
    --argjson pesr_interval_overlap ${pesr_interval_overlap} \
    --argjson pesr_breakend_window ${pesr_breakend_window} \
    --arg clustering_algorithm "${pesr_clustering_algorithm}" \
    --arg reference_fasta "${reference_fasta}" \
    --arg reference_fasta_fai "${reference_fasta_fai}" \
    --arg reference_dict "${reference_dict}" \
    --argjson sort_vcf_list false \
    --argjson vcfs '[]' \
    '{
          "batch": $batch,
          "caller": $caller,
          "clustering_algorithm": $clustering_algorithm,
          "contig_list": $contig_list,
          "contig_subset_list": $contig_subset_list,
          "exclude_intervals": $exclude_intervals,
          "min_size": $min_size,
          "pesr_breakend_window": $pesr_breakend_window,
          "pesr_interval_overlap": $pesr_interval_overlap,
          "ploidy_table": $ploidy_table,
          "reference_dict": $reference_dict,
          "reference_fasta": $reference_fasta,
          "reference_fasta_fai": $reference_fasta_fai,
          "sort_vcf_list": $sort_vcf_list,
          "vcf_tar": $vcf_tar,
          "vcfs": $vcfs
      }' > "${scramble_pesr_inputs_json}"

    bash /opt/sv_shell/cluster_pesr.sh "${scramble_pesr_inputs_json}" "${scramble_pesr_output_json}"

    scramble_pesr_clustered_vcf=$(jq -r '.clustered_vcf' "${scramble_pesr_output_json}")
    scramble_pesr_clustered_vcf_idx=$(jq -r '.clustered_vcf_index' "${scramble_pesr_output_json}")

    echo "Finished running PE/SR clustering on Scramble VCF."
fi


# Cluster Depth
# ---------------------------------------------------------------------------------------------------------------------

echo "Running cluster depth"

cd "${working_dir}"
cluster_depth_output_dir=$(mktemp -d "${SV_SHELL_BASE_DIR}/output_cluster_depth_XXXXXXXX")
cluster_depth_output_dir="$(realpath ${cluster_depth_output_dir})"
cluster_depth_inputs_json="$(realpath "${cluster_depth_output_dir}/cluster_depth_inputs.json")"
cluster_depth_output_json="$(realpath "${cluster_depth_output_dir}/cluster_depth_output.json")"

jq -n \
  --arg del_bed "${del_bed}" \
  --arg dup_bed "${dup_bed}" \
  --arg batch "${batch}" \
  --arg ploidy_table "${ploidy_table_output_file}" \
  --arg contig_list "${contig_list}" \
  --argjson contig_subset_list "${contig_subset_list:-null}" \
  --arg sample_list "${GetSampleIdsFromVcfTar_out_file}" \
  --arg records_per_bed_shard "${depth_records_per_bed_shard}" \
  --arg exclude_intervals "${depth_exclude_intervals}" \
  --arg exclude_overlap_fraction "${depth_exclude_overlap_fraction}" \
  --arg clustering_algorithm "${depth_clustering_algorithm}" \
  --arg depth_interval_overlap "${depth_interval_overlap}" \
  --arg reference_fasta "${reference_fasta}" \
  --arg reference_fasta_fai "${reference_fasta_fai}" \
  --arg reference_dict "${reference_dict}" \
  '{
        "del_bed": $del_bed,
        "dup_bed": $dup_bed,
        "batch": $batch,
        "ploidy_table": $ploidy_table,
        "contig_list": $contig_list,
        "contig_subset_list": $contig_subset_list,
        "sample_list": $sample_list,
        "records_per_bed_shard": $records_per_bed_shard,
        "exclude_intervals": $exclude_intervals,
        "exclude_overlap_fraction": $exclude_overlap_fraction,
        "clustering_algorithm": $clustering_algorithm,
        "depth_interval_overlap": $depth_interval_overlap,
        "reference_fasta": $reference_fasta,
        "reference_fasta_fai": $reference_fasta_fai,
        "reference_dict": $reference_dict
    }' > "${cluster_depth_inputs_json}"

bash /opt/sv_shell/cluster_depth.sh "${cluster_depth_inputs_json}" "${cluster_depth_output_json}"


# -------------------------------------------------------
# ======================= Output ========================
# -------------------------------------------------------


clustered_depth_vcf_wd=$(jq -r ".clustered_vcf" "${cluster_depth_output_json}")
clustered_depth_vcf_out="${output_dir}/$(basename "${clustered_depth_vcf_wd}")"
mv "${clustered_depth_vcf_wd}" "${clustered_depth_vcf_out}"

clustered_depth_vcf_index_wd=$(jq -r ".clustered_vcf_index" "${cluster_depth_output_json}")
clustered_depth_vcf_index_out="${output_dir}/$(basename "${clustered_depth_vcf_index_wd}")"
mv "${clustered_depth_vcf_index_wd}" "${clustered_depth_vcf_index_out}"

dragen_pesr_clustered_vcf_out=""
if [[ -n "${dragen_pesr_clustered_vcf}" ]]; then
  dragen_pesr_clustered_vcf_out="${output_dir}/$(basename "${dragen_pesr_clustered_vcf}")"
  mv "${dragen_pesr_clustered_vcf}" "${dragen_pesr_clustered_vcf_out}"
fi

dragen_pesr_clustered_vcf_idx_out=""
if [[ -n "${dragen_pesr_clustered_vcf_idx}" ]]; then
  dragen_pesr_clustered_vcf_idx_out="${output_dir}/$(basename "${dragen_pesr_clustered_vcf_idx}")"
  mv "${dragen_pesr_clustered_vcf_idx}" "${dragen_pesr_clustered_vcf_idx_out}"
fi

manta_pesr_clustered_vcf_out=""
if [[ -n "${manta_pesr_clustered_vcf}" ]]; then
  manta_pesr_clustered_vcf_out="${output_dir}/$(basename "${manta_pesr_clustered_vcf}")"
  mv "${manta_pesr_clustered_vcf}" "${manta_pesr_clustered_vcf_out}"
fi

manta_pesr_clustered_vcf_idx_out=""
if [[ -n "${manta_pesr_clustered_vcf_idx}" ]]; then
  manta_pesr_clustered_vcf_idx_out="${output_dir}/$(basename "${manta_pesr_clustered_vcf_idx}")"
  mv "${manta_pesr_clustered_vcf_idx}" "${manta_pesr_clustered_vcf_idx_out}"
fi

wham_pesr_clustered_vcf_out=""
if [[ -n "${wham_pesr_clustered_vcf}" ]]; then
  wham_pesr_clustered_vcf_out="${output_dir}/$(basename "${wham_pesr_clustered_vcf}")"
  mv "${wham_pesr_clustered_vcf}" "${wham_pesr_clustered_vcf_out}"
fi

wham_pesr_clustered_vcf_idx_out=""
if [[ -n "${wham_pesr_clustered_vcf_idx}" ]]; then
  wham_pesr_clustered_vcf_idx_out="${output_dir}/$(basename "${wham_pesr_clustered_vcf_idx}")"
  mv "${wham_pesr_clustered_vcf_idx}" "${wham_pesr_clustered_vcf_idx_out}"
fi

scramble_pesr_clustered_vcf_out=""
if [[ -n "${scramble_pesr_clustered_vcf}" ]]; then
  scramble_pesr_clustered_vcf_out="${output_dir}/$(basename "${scramble_pesr_clustered_vcf}")"
  mv "${scramble_pesr_clustered_vcf}" "${scramble_pesr_clustered_vcf_out}"
fi

scramble_pesr_clustered_vcf_idx_out=""
if [[ -n "${scramble_pesr_clustered_vcf_idx}" ]]; then
  scramble_pesr_clustered_vcf_idx_out="${output_dir}/$(basename "${scramble_pesr_clustered_vcf_idx}")"
  mv "${scramble_pesr_clustered_vcf_idx}" "${scramble_pesr_clustered_vcf_idx_out}"
fi

outputs_json=$(jq -n \
  --arg clustered_depth_vcf "${clustered_depth_vcf_out}" \
  --arg clustered_depth_vcf_index "${clustered_depth_vcf_index_out}" \
  --arg clustered_dragen_vcf "${dragen_pesr_clustered_vcf_out}" \
  --arg clustered_dragen_vcf_index "${dragen_pesr_clustered_vcf_idx_out}" \
  --arg clustered_manta_vcf "${manta_pesr_clustered_vcf_out}" \
  --arg clustered_manta_vcf_index "${manta_pesr_clustered_vcf_idx_out}" \
  --arg clustered_wham_vcf "${wham_pesr_clustered_vcf_out}" \
  --arg clustered_wham_vcf_index "${wham_pesr_clustered_vcf_idx_out}" \
  --arg clustered_scramble_vcf "${scramble_pesr_clustered_vcf_out}" \
  --arg clustered_scramble_vcf_index "${scramble_pesr_clustered_vcf_idx_out}" \
  '{
      "clustered_depth_vcf": $clustered_depth_vcf,
      "clustered_depth_vcf_index": $clustered_depth_vcf_index,
      "clustered_dragen_vcf": $clustered_dragen_vcf,
      "clustered_dragen_vcf_index": $clustered_dragen_vcf_index,
      "clustered_manta_vcf": $clustered_manta_vcf,
      "clustered_manta_vcf_index": $clustered_manta_vcf_index,
      "clustered_wham_vcf": $clustered_wham_vcf,
      "clustered_wham_vcf_index": $clustered_wham_vcf_index,
      "clustered_scramble_vcf": $clustered_scramble_vcf,
      "clustered_scramble_vcf_index": $clustered_scramble_vcf_index
  }' > "${output_json_filename}"
)

echo "Successfully finished Gather Batch Evidence, output json filename: ${output_json_filename}"