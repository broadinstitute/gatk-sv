#!/bin/bash

# This script is a bash implementation of the following workflow/task in WDL:
# Filename: wdl/Scramble.wdl
# Workflow: Scramble

set -Exeuo pipefail

sample_name=${1}
bam_or_cram_file=${2}
bam_or_cram_index=${3}
original_bam_or_cram_file=${4}
original_bam_or_cram_index=${5}
counts_file=${6}
input_vcf=${7}
reference_fasta=${8}
reference_index=${9}
regions_list=${10}

# Critical parameter for sensitivity/specificity
# Recommended values for aligners:
#   BWA-MEM: 90
#   DRAGEN-3.7.8: 60
alignment_score_cutoff=${11}

mei_bed=${12}

outputs_json_filename=${13}

min_clipped_reads_fraction=${14:-0.22}
percent_align_cutoff=${15:-70}

part2_threads=${16:-7}

scramble_vcf_script=${17:-"/opt/sv-pipeline/scripts/make_scramble_vcf.py"}
make_scramble_vcf_args=${18:-""}

min_clipped_reads=${19:-20}
# TODO: make this settable from the caller script.
#  If set to "" it will result in calculating this based on the read depth.
#  The default value is set to 20 to match the downsampled data for testing purpose.

echo "=============== Running scramble.sh"
echo "sample_name:                " "${sample_name}"
echo "bam_or_cram_file:           " "${bam_or_cram_file}"
echo "bam_or_cram_index:          " "${bam_or_cram_index}"
echo "original_bam_or_cram_file:  " "${original_bam_or_cram_file}"
echo "original_bam_or_cram_index: " "${original_bam_or_cram_index}"
echo "counts_file:                " "${counts_file}"
echo "input_vcf:                  " "${input_vcf}"
echo "reference_fasta:            " "${reference_fasta}"
echo "reference_index:            " "${reference_index}"
echo "regions_list:               " "${regions_list}"
echo "alignment_score_cutoff:     " "${alignment_score_cutoff}"
echo "mei_bed:                    " "${mei_bed}"
echo "min_clipped_reads_fraction: " "${min_clipped_reads_fraction}"
echo "percent_align_cutoff:       " "${percent_align_cutoff}"
echo "part2_threads:              " "${part2_threads}"
echo "scramble_vcf_script:        " "${scramble_vcf_script}"
echo "make_scramble_vcf_args:     " "${make_scramble_vcf_args}"


initial_wd=$PWD
output_dir=$(mktemp -d output_scramble_XXXXXXXX)
output_dir="$(realpath ${output_dir})"

scramble_dir="/app/scramble-gatk-sv"

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

working_dir_p1=$(mktemp -d wd_scramble_p1_XXXXXXXX)
working_dir_p1="$(realpath ${working_dir_p1})"
cd "${working_dir_p1}"

# Calibrate clipped reads cutoff based on median coverage
if [[ "${min_clipped_reads}" == "" ]]; then
  zcat "${counts_file}" \
    | awk '$0!~"@"' \
    | sed 1d \
    | awk 'NR % 100 == 0' \
    | cut -f4 \
    | Rscript -e "cat(round(${min_clipped_reads_fraction}*median(data.matrix(read.csv(file(\"stdin\"))))))" \
    > cutoff.txt
  MIN_CLIPPED_READS=$(cat cutoff.txt)
  echo "MIN_CLIPPED_READS: ${MIN_CLIPPED_READS}"
else
  MIN_CLIPPED_READS="${min_clipped_reads}"
fi

gzipped_clusters_file="${working_dir_p1}/${sample_name}_scramble_clusters.tsv.gz)"
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
working_dir_p2=$(mktemp -d wd_scramble_p2_XXXXXXXX)
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
working_dir_make_vcf=$(mktemp -d wd_scramble_make_vcf_XXXXXXXX)
working_dir_make_vcf="$(realpath ${working_dir_make_vcf})"
cd "${working_dir_make_vcf}"

python "${scramble_vcf_script}" \
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

outputs_filename="${output_dir}/outputs.json"
outputs_json=$(jq -n \
  --arg vcf "${vcf_filename}" \
  --arg vcf_idx "${vcf_idx_filename}" \
  --arg clusters "${clusters_filename}" \
  --arg table "${table_filename}" \
  '{vcf: $vcf, index: $vcf_idx, clusters: $clusters, table: $table}')
echo "${outputs_json}" > "${outputs_filename}"
cp "${outputs_filename}" "${outputs_json_filename}"
