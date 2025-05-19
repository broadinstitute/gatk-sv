#!/bin/bash

set -Eeuo pipefail

# ./scramble.sh NA12878.final.cram NA12878.final.cram.crai NA12878.final.cram NA12878.final.cram.crai NA12878.counts.tsv.gz NA12878.manta.vcf.gz test Homo_sapiens_assembly38.fasta Homo_sapiens_assembly38.fasta.fai primary_contigs.list 90 hg38.repeatmasker.mei.with_SVA.pad_50_merged.bed.gz

sample_name=${1}
bam_or_cram_file=${2}
bam_or_cram_index=${3}
original_bam_or_cram_file=${4}
original_bam_or_cram_index=${5}
counts_file=${6}
manta_vcf=${7}
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

min_clipped_reads=${18:-20}
# TODO: let this to be settable from the caller,
#  if you set it to "" it will result in calculating this
#  based on the read depth, override it mainly for testing using downsampled files.
#  Note that the default value is set to 20 to match the downsampled data for testing purpose.


echo "=============== Running scramble.sh"
echo "sample_name:                " "${sample_name}"
echo "bam_or_cram_file:           " "${bam_or_cram_file}"
echo "bam_or_cram_index:          " "${bam_or_cram_index}"
echo "original_bam_or_cram_file:  " "${original_bam_or_cram_file}"
echo "original_bam_or_cram_index: " "${original_bam_or_cram_index}"
echo "counts_file:                " "${counts_file}"
echo "manta_vcf:                  " "${manta_vcf}"
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


# In case it is re-run, the script will wait for a response
# on override the existing file, that may not work in a pipeline.
rm -f test.scramble.tsv.gz

# TODO: this seems an overkill
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

# ScramblePart1
# -------------

# Calibrate clipped reads cutoff based on median coverage
if [[ "${min_clipped_reads}" == "" ]]; then
  zcat "${counts_file}" \
    | awk '$0!~"@"' \
    | sed 1d \
    | awk 'NR % 100 == 0' \
    | cut -f4 \
    | Rscript -e "cat(round(${min_clipped_reads_fraction}*median(data.matrix(read.csv(file(\"stdin\"))))))" \
    > cutoff.txt
  export MIN_CLIPPED_READS=$(cat cutoff.txt)
  echo "MIN_CLIPPED_READS: ${MIN_CLIPPED_READS}"
else
  export MIN_CLIPPED_READS="${min_clipped_reads}"
fi

# Identify clusters of split reads
while read region; do
  time /app/scramble-gatk-sv/cluster_identifier/src/build/cluster_identifier -l -s ${MIN_CLIPPED_READS} -r "${region}" -t "${reference_fasta}" "${bam_or_cram_file}" \
    | gzip >> "${sample_name}".scramble_clusters.tsv.gz
done < "${regions_list}"

export clusters_file="${sample_name}.scramble_clusters.tsv.gz"

# ScramblePart2
# -------------

export xDir="" # $PWD
export clusterFile="${xDir}/clusters"
export scrambleDir="/app/scramble-gatk-sv"
export meiRef=$scrambleDir/cluster_analysis/resources/MEI_consensus_seqs.fa

# create a blast db from the reference
cat "${reference_fasta}" | makeblastdb -in - -parse_seqids -title ref -dbtype nucl -out ref

gunzip -c "${clusters_file}" > "${clusterFile}"

# Produce ${clusterFile}_MEIs.txt
Rscript --vanilla "${scrambleDir}"/cluster_analysis/bin/SCRAMble.R \
  --out-name "${clusterFile}" \
  --cluster-file "${clusterFile}" \
  --install-dir "${scrambleDir}"/cluster_analysis/bin \
  --mei-refs "${meiRef}" \
  --ref "${xDir}/ref" \
  --no-vcf --eval-meis \
  --cores "${part2_threads}" \
  --pct-align "${percent_align_cutoff}" \
  -n "${MIN_CLIPPED_READS}" \
  --mei-score "${alignment_score_cutoff}"

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
