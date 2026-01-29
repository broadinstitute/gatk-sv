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
  output_dir=$(mktemp -d ${SV_SHELL_BASE_DIR}/output_combine_batches_XXXXXXXX)
else
  mkdir -p "${output_dir}"
fi
output_dir="$(realpath ${output_dir})"

if [ -z "${output_json_filename}" ]; then
  output_json_filename="${output_dir}/output.json"
else
  output_json_filename="$(realpath ${output_json_filename})"
fi

working_dir=$(mktemp -d ${SV_SHELL_BASE_DIR}/wd_combine_batches_XXXXXXXX)
working_dir="$(realpath ${working_dir})"
cd "${working_dir}"

batches=$(jq -r ".batches" "${input_json}")
cohort_name=$(jq -r ".cohort_name" "${input_json}")
pesr_vcfs=$(jq -r ".pesr_vcfs" "${input_json}")
depth_vcfs=$(jq -r ".depth_vcfs" "${input_json}")
raw_sr_bothside_pass_files=($(jq -r ".raw_sr_bothside_pass_files[]" "${input_json}"))
raw_sr_background_fail_files=($(jq -r ".raw_sr_background_fail_files[]" "${input_json}"))
min_sr_background_fail_batches=$(jq -r ".min_sr_background_fail_batches" "${input_json}")
ped_file=$(jq -r ".ped_file" "${input_json}")
contig_list=$(jq -r ".contig_list" "${input_json}")
chr_x=$(jq -r ".chr_x" "${input_json}")
chr_y=$(jq -r ".chr_y" "${input_json}")
reference_dict=$(jq -r '.reference_dict' "$input_json")
reference_fasta=$(jq -r '.reference_fasta' "$input_json")
reference_fasta_fai=$(jq -r '.reference_fasta_fai' "$input_json")
java_mem_fraction=$(jq -r '.java_mem_fraction // ""' "$input_json")
clustering_config_part1=$(jq -r ".clustering_config_part1" "${input_json}")
clustering_config_part2=$(jq -r ".clustering_config_part2" "${input_json}")
stratification_config_part1=$(jq -r ".stratification_config_part1" "${input_json}")
stratification_config_part2=$(jq -r ".stratification_config_part2" "${input_json}")
track_bed_files=($(jq -r ".track_bed_files[]" "${input_json}"))
track_names=($(jq -r ".track_names[]" "${input_json}"))


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


# -------------------------------------------------------
# ======================= Command =======================
# -------------------------------------------------------


# CombineSRBothsidePass
# ---------------------------------------------------------------------------------------------------------------------

combine_sr_both_side_pass_prefix="${cohort_name}.sr_bothside_pass"


# GetNonRefVariantLists
non_ref_variants_list="$(realpath "${combine_sr_both_side_pass_prefix}.non_ref_vids.shard_0.list")"
bcftools view -G -i 'SUM(AC)>0||SUM(FORMAT/SR_GT)>0' "${pesr_vcfs}" |\
 bcftools query -f '%ID\n' > "${non_ref_variants_list}"


# CalculateBothsideSupportFraction
calculate_both_side_support_fraction_output="$(realpath "${combine_sr_both_side_pass_prefix}.sr_bothside_support.txt")"
echo "${non_ref_variants_list}" > non_ref_vid_lists.txt
printf "%s\n" "${raw_sr_bothside_pass_files[@]}" > raw_sr_bothside_pass_files.txt
python /opt/sv-pipeline/04_variant_resolution/scripts/calculate_sr_bothside_support.py \
    non_ref_vid_lists.txt \
    raw_sr_bothside_pass_files.txt \
    > "${calculate_both_side_support_fraction_output}"
CombineSRBothsidePass_output="${calculate_both_side_support_fraction_output}"


# CombineBackgroundFail
CombineBackgroundFail_out="$(realpath "${cohort_name}.background_fail.txt")"
min_background_fail_first_col=$(awk "BEGIN {print ${min_sr_background_fail_batches} * ${#raw_sr_background_fail_files[@]}}")

printf "%s\n" "${raw_sr_background_fail_files[@]}" > raw_sr_background_fail_files.txt
while read SHARD; do
  if [ -n "${SHARD}" ]; then
    cat "${SHARD}"
  fi
done < raw_sr_background_fail_files.txt \
  | sort | uniq -c | \
  awk -v OFS='\t' -v min_val="${min_background_fail_first_col}" '{if($1 >= min_val) print $2}' \
  > "${CombineBackgroundFail_out}"


# CreatePloidyTableFromPed
CreatePloidyTableFromPed_out="$(realpath "${cohort_name}.ploidy.FEMALE_chrY_1.tsv")"
python /opt/sv-pipeline/scripts/ploidy_table_from_ped.py \
  --ped "${ped_file}" \
  --out tmp.tsv \
  --contigs "${contig_list}" \
  --chr-x chr_x \
  --chr-y chr_y

#retain_female_chr_y
sed -e 's/\t0/\t1/g' tmp.tsv > "${CreatePloidyTableFromPed_out}"


# GatkFormatting.FormatVcf
all_vcfs=("$pesr_vcfs" "$depth_vcfs")
reformatted_vcfs_array=()
reformatted_vcfs_idx_array=()
for vcf in "${all_vcfs[@]}"; do
  cd "${working_dir}"

  format_vcf_working_dir=$(mktemp -d ${SV_SHELL_BASE_DIR}/wd_format_vcf_XXXXXXXX)
  format_vcf_working_dir="$(realpath ${format_vcf_working_dir})"
  cd "${format_vcf_working_dir}"

  # Convert format
  formatted_vcf_output="$(basename "$vcf" ".vcf.gz").reformat_gatk.vcf.gz"
  formatted_vcf_output="$(realpath "${formatted_vcf_output}")"
  python /opt/sv-pipeline/scripts/format_svtk_vcf_for_gatk.py \
    --vcf "${vcf}" \
    --out "${formatted_vcf_output}" \
    --ploidy-table "${CreatePloidyTableFromPed_out}" \
    --bothside-pass-list "${CombineSRBothsidePass_output}" \
    --background-fail-list "${CombineBackgroundFail_out}" \
    --add-sr-pos --scale-down-gq

  tabix "${formatted_vcf_output}"

  reformatted_vcfs_array+=("${formatted_vcf_output}")
  reformatted_vcfs_idx_array+=("${formatted_vcf_output}.tbi")
done



# JoinVcfs
# ---------------------------------------------------------------------------------------------------------------------

cd "${working_dir}"
join_vcfs_output_dir=$(mktemp -d ${SV_SHELL_BASE_DIR}/output_join_vcfs_XXXXXXXX)
join_vcfs_output_dir="$(realpath ${join_vcfs_output_dir})"
join_vcfs_inputs_json="$(realpath "${join_vcfs_output_dir}/sv_cluster_inputs.json")"
join_vcfs_outputs_json="$(realpath "${join_vcfs_output_dir}/sv_cluster_outputs.json")"

vcf_json_array=$(printf "%s\n" "${reformatted_vcfs_array[@]}" | jq -nR '[inputs]')

jq -n \
  --argjson vcfs "${vcf_json_array}" \
  --arg ploidy_table "${CreatePloidyTableFromPed_out}" \
  --arg output_prefix "${cohort_name}.combine_batches.join_vcfs" \
  --arg fast_mode false \
  --argjson pesr_sample_overlap 0 \
  --argjson pesr_interval_overlap 1 \
  --argjson pesr_breakend_window 0 \
  --argjson depth_sample_overlap 0 \
  --argjson depth_interval_overlap 1 \
  --argjson depth_breakend_window 0 \
  --argjson mixed_sample_overlap 0 \
  --argjson mixed_interval_overlap 1 \
  --argjson mixed_breakend_window 0 \
  --arg reference_fasta "${reference_fasta}" \
  --arg reference_fasta_fai "${reference_fasta_fai}" \
  --arg reference_dict "${reference_dict}" \
  --arg java_mem_fraction "${java_mem_fraction}" \
  '{
      "vcfs": $vcfs,
      "ploidy_table": $ploidy_table,
      "output_prefix": $output_prefix,
      "fast_mode": $fast_mode,
      "pesr_sample_overlap": $pesr_sample_overlap,
      "pesr_interval_overlap": $pesr_interval_overlap,
      "pesr_breakend_window": $pesr_breakend_window,
      "depth_sample_overlap": $depth_sample_overlap,
      "depth_interval_overlap": $depth_interval_overlap,
      "depth_breakend_window": $depth_breakend_window,
      "mixed_sample_overlap": $mixed_sample_overlap,
      "mixed_interval_overlap": $mixed_interval_overlap,
      "mixed_breakend_window": $mixed_breakend_window,
      "reference_fasta": $reference_fasta,
      "reference_fasta_fai": $reference_fasta_fai,
      "reference_dict": $reference_dict,
      "java_mem_fraction": $java_mem_fraction
  }' > "${join_vcfs_inputs_json}"

bash /opt/sv_shell/sv_cluster.sh "${join_vcfs_inputs_json}" "${join_vcfs_outputs_json}" "${join_vcfs_output_dir}"

join_vcfs_out=$(jq -r ".out" "${join_vcfs_outputs_json}")

echo "Successfully finished Join VCFs clustering."


# ClusterSites
# First round of clustering
# ---------------------------------------------------------------------------------------------------------------------

cd "${working_dir}"
cluster_sites_output_dir=$(mktemp -d ${SV_SHELL_BASE_DIR}/output_cluster_sites_XXXXXXXX)
cluster_sites_output_dir="$(realpath ${cluster_sites_output_dir})"
cluster_sites_inputs_json="$(realpath "${cluster_sites_output_dir}/sv_cluster_inputs.json")"
cluster_sites_outputs_json="$(realpath "${cluster_sites_output_dir}/sv_cluster_outputs.json")"

join_vcfs_out_array=("${join_vcfs_out}")
cluster_sites_input_vcf_array=$(printf "%s\n" "${join_vcfs_out_array[@]}" | jq -nR '[inputs]')

jq -n \
  --argjson vcfs "${cluster_sites_input_vcf_array}" \
  --arg ploidy_table "${CreatePloidyTableFromPed_out}" \
  --arg output_prefix "${cohort_name}.combine_batches.cluster_sites" \
  --arg fast_mode false \
  --arg breakpoint_summary_strategy "REPRESENTATIVE" \
  --argjson pesr_sample_overlap 0.5 \
  --argjson pesr_interval_overlap 0.1 \
  --argjson pesr_breakend_window 300 \
  --argjson depth_sample_overlap 0.5 \
  --argjson depth_interval_overlap 0.5 \
  --argjson depth_breakend_window 500000 \
  --argjson mixed_sample_overlap 0.5 \
  --argjson mixed_interval_overlap 0.5 \
  --argjson mixed_breakend_window 1000000 \
  --arg reference_fasta "${reference_fasta}" \
  --arg reference_fasta_fai "${reference_fasta_fai}" \
  --arg reference_dict "${reference_dict}" \
  --arg java_mem_fraction "${java_mem_fraction}" \
  --arg variant_prefix "${cohort_name}_" \
  '{
      "vcfs": $vcfs,
      "ploidy_table": $ploidy_table,
      "output_prefix": $output_prefix,
      "fast_mode": $fast_mode,
      "breakpoint_summary_strategy": $breakpoint_summary_strategy,
      "pesr_sample_overlap": $pesr_sample_overlap,
      "pesr_interval_overlap": $pesr_interval_overlap,
      "pesr_breakend_window": $pesr_breakend_window,
      "depth_sample_overlap": $depth_sample_overlap,
      "depth_interval_overlap": $depth_interval_overlap,
      "depth_breakend_window": $depth_breakend_window,
      "mixed_sample_overlap": $mixed_sample_overlap,
      "mixed_interval_overlap": $mixed_interval_overlap,
      "mixed_breakend_window": $mixed_breakend_window,
      "reference_fasta": $reference_fasta,
      "reference_fasta_fai": $reference_fasta_fai,
      "reference_dict": $reference_dict,
      "java_mem_fraction": $java_mem_fraction,
      "variant_prefix": $variant_prefix
  }' > "${cluster_sites_inputs_json}"

bash /opt/sv_shell/sv_cluster.sh "${cluster_sites_inputs_json}" "${cluster_sites_outputs_json}" "${cluster_sites_output_dir}"

ClusterSites_out=$(jq -r ".out" "${cluster_sites_outputs_json}")

echo "Successfully finished Cluster Sites clustering."


track_intervals_array=()
for x in "${track_bed_files[@]}"; do
  track_intervals_array+=(--track-intervals "${x}")
done

track_names_array=()
for x in "${track_names[@]}"; do
  track_names_array+=(--track-name "${x}")
done


# GroupedSVClusterPart1
# Second round of clustering
# ---------------------------------------------------------------------------------------------------------------------

cd "${working_dir}"
grouped_sv_cluster_p1_wd=$(mktemp -d ${SV_SHELL_BASE_DIR}/wd_grouped_sv_cluster_p1_XXXXXXXX)
grouped_sv_cluster_p1_wd="$(realpath ${grouped_sv_cluster_p1_wd})"
cd "${grouped_sv_cluster_p1_wd}"

GroupedSVClusterPart1_out_prefix="${cohort_name}.combine_batches.recluster_part_1"

java "-Xmx${JVM_MAX_MEM}" -jar /opt/gatk.jar GroupedSVCluster \
  --reference "${reference_fasta}" \
  --ploidy-table "${CreatePloidyTableFromPed_out}" \
  -V "${ClusterSites_out}" \
  -O "${GroupedSVClusterPart1_out_prefix}.vcf.gz" \
  --clustering-config "${clustering_config_part1}" \
  --stratify-config "${stratification_config_part1}" \
  "${track_intervals_array[@]}" \
  "${track_names_array[@]}" \
  --stratify-overlap-fraction 0 \
  --stratify-num-breakpoint-overlaps 1 \
  --stratify-num-breakpoint-overlaps-interchromosomal 1 \
  --breakpoint-summary-strategy "REPRESENTATIVE"

GroupedSVClusterPart1_vcf_out=$(realpath "${GroupedSVClusterPart1_out_prefix}.vcf.gz")
GroupedSVClusterPart1_vcf_out_index=$(realpath "${GroupedSVClusterPart1_out_prefix}.vcf.gz.tbi")

echo "Successfully finished grouped SV cluster part 1."


# GroupedSVClusterPart2
# Third round of clustering
# ---------------------------------------------------------------------------------------------------------------------

cd "${working_dir}"
grouped_sv_cluster_p2_wd=$(mktemp -d ${SV_SHELL_BASE_DIR}/wd_grouped_sv_cluster_p2_XXXXXXXX)
grouped_sv_cluster_p2_wd="$(realpath ${grouped_sv_cluster_p2_wd})"
cd "${grouped_sv_cluster_p2_wd}"

GroupedSVClusterPart2_out_prefix="${cohort_name}.combine_batches.recluster_part_2"

java "-Xmx${JVM_MAX_MEM}" -jar /opt/gatk.jar GroupedSVCluster \
  --reference "${reference_fasta}" \
  --ploidy-table "${CreatePloidyTableFromPed_out}" \
  -V "${GroupedSVClusterPart1_vcf_out}" \
  -O "${GroupedSVClusterPart2_out_prefix}.vcf.gz" \
  --clustering-config "${clustering_config_part2}" \
  --stratify-config "${stratification_config_part2}" \
  "${track_intervals_array[@]}" \
  "${track_names_array[@]}" \
  --stratify-overlap-fraction 0 \
  --stratify-num-breakpoint-overlaps 1 \
  --stratify-num-breakpoint-overlaps-interchromosomal 1 \
  --breakpoint-summary-strategy "REPRESENTATIVE"

GroupedSVClusterPart2_vcf_out="$(realpath "${GroupedSVClusterPart2_out_prefix}.vcf.gz")"
GroupedSVClusterPart2_vcf_out_index="$(realpath "${GroupedSVClusterPart2_out_prefix}.vcf.gz.tbi")"

echo "Successfully finished grouped SV cluster part 2."


# GatkToSvtkVcf
# ---------------------------------------------------------------------------------------------------------------------
# Use "depth" as source to match legacy headers
# AC/AF cause errors due to being lists instead of single values

cd "${working_dir}"
gatk_to_svtk_vcf_prefix="${cohort_name}.combine_batches.svtk_formatted"

python /opt/sv-pipeline/scripts/format_gatk_vcf_for_svtk.py \
  --vcf "${GroupedSVClusterPart2_vcf_out}" \
  --out "${gatk_to_svtk_vcf_prefix}.vcf.gz" \
  --source "depth" \
  --contigs "${contig_list}" \
  --remove-infos "AC,AF,AN,HIGH_SR_BACKGROUND,BOTHSIDES_SUPPORT,SR1POS,SR2POS" \
  --remove-formats "CN" \
  --set-pass
tabix "${gatk_to_svtk_vcf_prefix}.vcf.gz"

GatkToSvtkVcf_out="$(realpath "${gatk_to_svtk_vcf_prefix}.vcf.gz")"
GatkToSvtkVcf_out_index="$(realpath "${gatk_to_svtk_vcf_prefix}.vcf.gz.tbi")"


# ExtractSRVariantLists
# ---------------------------------------------------------------------------------------------------------------------

ExtractSRVariantLists_out_prefix="${cohort_name}.combine_batches"
high_sr_background_list="$(realpath "${ExtractSRVariantLists_out_prefix}.high_sr_background.txt")"
bothsides_sr_support="$(realpath "${ExtractSRVariantLists_out_prefix}.bothsides_sr_support.txt")"

bcftools query -f '%ID\t%HIGH_SR_BACKGROUND\t%BOTHSIDES_SUPPORT\n' "${GroupedSVClusterPart2_vcf_out}" > flags.txt
awk -F'\t' '($2 != "."){print $1}' flags.txt > "${high_sr_background_list}"
awk -F'\t' '($3 != "."){print $1}' flags.txt > "${bothsides_sr_support}"



# -------------------------------------------------------
# ======================= Output ========================
# -------------------------------------------------------

combined_vcfs_out="${output_dir}/$(basename "${GatkToSvtkVcf_out}")"
mv "${GatkToSvtkVcf_out}" "${combined_vcfs_out}"

combined_vcf_indexes_out="${output_dir}/$(basename "${GatkToSvtkVcf_out_index}")"
mv "${GatkToSvtkVcf_out_index}" "${combined_vcf_indexes_out}"

cluster_background_fail_lists_out="${output_dir}/$(basename "${high_sr_background_list}")"
mv "${high_sr_background_list}" "${cluster_background_fail_lists_out}"

cluster_bothside_pass_lists_out="${output_dir}/$(basename "${bothsides_sr_support}")"
mv "${bothsides_sr_support}" "${cluster_bothside_pass_lists_out}"

outputs_json=$(jq -n \
  --arg combined_vcfs "${combined_vcfs_out}" \
  --arg combined_vcf_indexes "${combined_vcf_indexes_out}" \
  --arg cluster_background_fail_lists "${cluster_background_fail_lists_out}" \
  --arg cluster_bothside_pass_lists "${cluster_bothside_pass_lists_out}" \
  '{
      "combined_vcfs": $combined_vcfs,
      "combined_vcf_indexes": $combined_vcf_indexes,
      "cluster_background_fail_lists": $cluster_background_fail_lists,
      "cluster_bothside_pass_lists": $cluster_bothside_pass_lists
  }' > "${output_json_filename}"
)

echo "Successfully finished Combine Batches, output json filename: ${output_json_filename}"
