#!/bin/bash

set -Eeuo pipefail

# ./scramble.sh NA12878.final.cram NA12878.final.cram.crai NA12878.final.cram NA12878.final.cram.crai NA12878.counts.tsv.gz NA12878.manta.vcf.gz test Homo_sapiens_assembly38.fasta Homo_sapiens_assembly38.fasta.fai primary_contigs.list 90 hg38.repeatmasker.mei.with_SVA.pad_50_merged.bed.gz

bam_or_cram_file=${1}
bam_or_cram_index=${2}
original_bam_or_cram_file=${3}
original_bam_or_cram_index=${4}
counts_file=${5}
manta_vcf=${6}
sample_name=${7}
reference_fasta=${8}
reference_index=${9}
regions_list=${10}

# Critical parameter for sensitivity/specificity
# Recommended values for aligners:
#   BWA-MEM: 90
#   DRAGEN-3.7.8: 60
alignment_score_cutoff=${11}

mei_bed=${12}

min_clipped_reads_fraction=${13:-0.22}
percent_align_cutoff=${14:-70}

part2_threads=${15:-7}

scramble_vcf_script=${16:-"/opt/sv-pipeline/scripts/make_scramble_vcf.py"}
make_scramble_vcf_args=${17:-""}

# In case it is re-run, the script will wait for a response
# on override the existing file, that may not work in a pipeline.
rm -f test.scramble.tsv.gz

# ScramblePart1
# -------------

# Calibrate clipped reads cutoff based on median coverage
zcat "${counts_file}" \
  | awk '$0!~"@"' \
  | sed 1d \
  | awk 'NR % 100 == 0' \
  | cut -f4 \
  | Rscript -e "cat(round(${min_clipped_reads_fraction}*median(data.matrix(read.csv(file(\"stdin\"))))))" \
  > cutoff.txt
export MIN_CLIPPED_READS=$(cat cutoff.txt)
echo "MIN_CLIPPED_READS: ${MIN_CLIPPED_READS}"

# Identify clusters of split reads
while read region; do
  time /app/scramble-gatk-sv/cluster_identifier/src/build/cluster_identifier -l -s ${MIN_CLIPPED_READS} -r "${region}" -t "${reference_fasta}" "${bam_or_cram_file}" \
    | gzip >> "${sample_name}".scramble_clusters.tsv.gz
done < "${regions_list}"

export clusters_file="${sample_name}.scramble_clusters.tsv.gz"

# ScramblePart2
# -------------

export xDir=$PWD
export clusterFile=$xDir/clusters
export scrambleDir="/app/scramble-gatk-sv"
export meiRef=$scrambleDir/cluster_analysis/resources/MEI_consensus_seqs.fa

# create a blast db from the reference
cat "${reference_fasta}" | makeblastdb -in - -parse_seqids -title ref -dbtype nucl -out ref

gunzip -c "${clusters_file}" > $clusterFile

# Produce ${clusterFile}_MEIs.txt
Rscript --vanilla "${scrambleDir}"/cluster_analysis/bin/SCRAMble.R --out-name "${clusterFile}" \
        --cluster-file "${clusterFile}" --install-dir "${scrambleDir}"/cluster_analysis/bin \
        --mei-refs "${meiRef}" --ref "${xDir}"/ref --no-vcf --eval-meis --cores "${part2_threads}" \
        --pct-align "${percent_align_cutoff}" -n "${MIN_CLIPPED_READS}" --mei-score "${alignment_score_cutoff}"

# Save raw outputs
mv ${clusterFile}_MEIs.txt "${sample_name}".scramble.tsv
gzip "${sample_name}".scramble.tsv


# MakeScrambleVcf
# --------------

export scramble_table="${sample_name}.scramble.tsv.gz"

python "${scramble_vcf_script}" \
  --table "${scramble_table}" \
  --manta-vcf "${manta_vcf}" \
  --alignments-file "${original_bam_or_cram_file}" \
  --sample "${sample_name}" \
  --reference "${reference_fasta}" \
  --mei-bed "${mei_bed}" \
  --out unsorted.vcf.gz \
  ${make_scramble_vcf_args}
bcftools sort unsorted.vcf.gz -Oz -o "${sample_name}".scramble.vcf.gz
tabix "${sample_name}".scramble.vcf.gz
