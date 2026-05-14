#!/bin/bash

# This script is a bash implementation of the following workflow/task in WDL:
# Filename: wdl/Scramble.wdl
# Workflow: Scramble

set -Exeuo pipefail

if [ -z "${SV_SHELL_CLEAN_UP_WORKING_DIR:-}" ]; then
  SV_SHELL_CLEAN_UP_WORKING_DIR=true
fi

# -------------------------------------------------------
# ==================== Input & Setup ====================
# -------------------------------------------------------

input_json=${1}
outputs_json_filename=${2:-""}
output_dir=${3:-""}

input_json="$(realpath "${input_json}")"

if [ -z "${output_dir}" ]; then
  output_dir=$(mktemp -d ${SV_SHELL_BASE_DIR}/output_scramble_XXXXXXXX)
else
  mkdir -p "${output_dir}"
fi
output_dir="$(realpath "${output_dir}")"

if [ -z "${outputs_json_filename}" ]; then
  outputs_json_filename="${output_dir}/output.json"
else
  outputs_json_filename="$(realpath "${outputs_json_filename}")"
fi

sample_name=$(jq -r ".sample_name" "${input_json}")
bam_or_cram_file=$(jq -r ".bam_or_cram_file" "${input_json}")
bam_or_cram_index=$(jq -r ".bam_or_cram_index" "${input_json}")
original_bam_or_cram_file=$(jq -r ".original_bam_or_cram_file" "${input_json}")
original_bam_or_cram_index=$(jq -r ".original_bam_or_cram_index" "${input_json}")
counts_file=$(jq -r ".counts_file" "${input_json}")
input_vcf=$(jq -r ".input_vcf" "${input_json}")
reference_fasta=$(jq -r ".reference_fasta" "${input_json}")
reference_index=$(jq -r ".reference_index" "${input_json}")
regions_list=$(jq -r ".regions_list" "${input_json}")

# Critical parameter for sensitivity/specificity
# Recommended values for aligners:
#   BWA-MEM: 90
#   DRAGEN-3.7.8: 60
alignment_score_cutoff=$(jq -r ".alignment_score_cutoff" "${input_json}")

mei_bed=$(jq -r ".mei_bed" "${input_json}")
min_clipped_reads_fraction=$(jq -r '.min_clipped_reads_fraction // 0.22' "${input_json}")
percent_align_cutoff=$(jq -r '.percent_align_cutoff // 70' "${input_json}")
part2_threads=$(jq -r '.part2_threads // 7' "${input_json}")
make_scramble_vcf_args=$(jq -r '.make_scramble_vcf_args // ""' "${input_json}")

echo "=============== Running scramble.sh"

initial_wd=$PWD
output_dir=$(mktemp -d ${SV_SHELL_BASE_DIR}/output_scramble_XXXXXXXX)
output_dir="$(realpath ${output_dir})"

scramble_dir="/app/scramble-gatk-sv"


# -------------------------------------------------------
# ======================= Command =======================
# -------------------------------------------------------


# We bypass version detection and conservatively assume Dragen 3.7.8.
# This adds a few extra steps but is safer logic in case of issues with other Dragen aligner versions.
# Check aligner
#gatk PrintReadsHeader \
#  -I "${bam_or_cram_file}" \
#  --read-index "${bam_or_cram_index}" \
#  -O "${sample_name}.header.sam" \
#  -R "${reference_fasta}"
#
#count=$(awk '$0~"@PG" && $0~"ID: DRAGEN SW build" && $0~"VN: 05.021.604.3.7.8"' "${sample_id}.header.sam" | wc -l)
#if [ "${count}" -gt 0 ]; then
#  is_dragen_3_7_8="true"
#else
#  is_dragen_3_7_8="false"
#fi

# -------------
# ScramblePart1
# -------------

working_dir_p1=$(mktemp -d ${SV_SHELL_BASE_DIR}/wd_scramble_p1_XXXXXXXX)
working_dir_p1="$(realpath ${working_dir_p1})"
cd "${working_dir_p1}"

zcat "${counts_file}" \
  | awk '$0!~"@"' \
  | sed 1d \
  | awk 'NR % 100 == 0' \
  | cut -f4 \
  | Rscript -e "cat(round(${min_clipped_reads_fraction}*median(data.matrix(read.csv(file(\"stdin\"))))))" \
  > cutoff.txt
MIN_CLIPPED_READS=$(cat cutoff.txt)
echo "MIN_CLIPPED_READS: ${MIN_CLIPPED_READS}"

gzipped_clusters_file="${working_dir_p1}/${sample_name}_scramble_clusters.tsv.gz"
gzipped_clusters_file="$(realpath ${gzipped_clusters_file})"

# Identify clusters of split reads
while read region; do
  time "${scramble_dir}"/cluster_identifier/src/build/cluster_identifier -l -s ${MIN_CLIPPED_READS} -r "${region}" -t "${reference_fasta}" "${bam_or_cram_file}" \
    | gzip >> "${gzipped_clusters_file}"
done < "${regions_list}"

# -------------
# ScramblePart2
# -------------

cd "${initial_wd}"
working_dir_p2=$(mktemp -d ${SV_SHELL_BASE_DIR}/wd_scramble_p2_XXXXXXXX)
working_dir_p2="$(realpath ${working_dir_p2})"
cd "${working_dir_p2}"

clusters_file="${working_dir_p2}/${sample_name}_scramble_clusters.tsv"
meiRef="${scramble_dir}/cluster_analysis/resources/MEI_consensus_seqs.fa"

# create a blast db from the reference
cat "${reference_fasta}" | makeblastdb -in - -parse_seqids -title ref -dbtype nucl -out ref

gunzip -c "${gzipped_clusters_file}" > "${clusters_file}"

# Produce ${clusterFile}_MEIs.txt
Rscript --vanilla "${scramble_dir}"/cluster_analysis/bin/SCRAMble.R \
  --out-name "${clusters_file}" \
  --cluster-file "${clusters_file}" \
  --install-dir "${scramble_dir}"/cluster_analysis/bin \
  --mei-refs "${meiRef}" \
  --ref "${working_dir_p2}/ref" \
  --no-vcf --eval-meis \
  --cores "${part2_threads}" \
  --pct-align "${percent_align_cutoff}" \
  -n "${MIN_CLIPPED_READS}" \
  --mei-score "${alignment_score_cutoff}"

mv ${clusters_file}_MEIs.txt "${sample_name}".scramble.tsv
gzip "${sample_name}".scramble.tsv
scramble_table="$(realpath ${sample_name}.scramble.tsv.gz)"


# ---------------
# MakeScrambleVcf
# ---------------

cd "${initial_wd}"
working_dir_make_vcf=$(mktemp -d ${SV_SHELL_BASE_DIR}/wd_scramble_make_vcf_XXXXXXXX)
working_dir_make_vcf="$(realpath ${working_dir_make_vcf})"
cd "${working_dir_make_vcf}"

python /opt/sv-pipeline/scripts/make_scramble_vcf.py \
  --table "${scramble_table}" \
  --input-vcf "${input_vcf}" \
  --alignments-file "${original_bam_or_cram_file}" \
  --sample "${sample_name}" \
  --reference "${reference_fasta}" \
  --mei-bed "${mei_bed}" \
  --out unsorted.vcf.gz \
  ${make_scramble_vcf_args}
bcftools sort unsorted.vcf.gz -Oz -o "${sample_name}".scramble.vcf.gz
tabix "${sample_name}".scramble.vcf.gz

vcf_filename="${output_dir}/${sample_name}.scramble.vcf.gz"
vcf_idx_filename="${output_dir}/${sample_name}.scramble.vcf.gz.tbi"
clusters_filename="${output_dir}/${sample_name}.scramble_clusters.tsv.gz"
table_filename="${output_dir}/${sample_name}.scramble.tsv.gz"

mv "${sample_name}.scramble.vcf.gz" "${vcf_filename}"
mv "${sample_name}.scramble.vcf.gz.tbi" "${vcf_idx_filename}"
mv "${gzipped_clusters_file}" "${clusters_filename}"
mv "${scramble_table}" "${table_filename}"

if [ "${SV_SHELL_CLEAN_UP_WORKING_DIR}" == "true" ]; then
  rm -rf "${working_dir_p1}"
  rm -rf "${working_dir_p2}"
  rm -rf "${working_dir_make_vcf}"
fi


# -------------------------------------------------------
# ======================= Output ========================
# -------------------------------------------------------


jq -n \
  --arg vcf "${vcf_filename}" \
  --arg vcf_idx "${vcf_idx_filename}" \
  --arg clusters "${clusters_filename}" \
  --arg table "${table_filename}" \
  '{vcf: $vcf, index: $vcf_idx, clusters: $clusters, table: $table}' > "${outputs_json_filename}"

echo "Successfully finished Scramble."
