#!/bin/bash

set -Eeuo pipefail

intervals=$1
cram_or_bam=$2
cram_or_bam_idx=$3
sample_id=$4
ref_fasta=$5
ref_fasta_fai=$6
ref_fasta_dict=$7  # TODO: this is not used in commands, not sure if it is needed or not?
gatk4_jar_override=${8:-/root/gatk.jar}
disabled_read_filters=${9:-""}

rm -rf "${sample_id}".counts.tsv

# TODO: in the original code, this is computed based on a few factors,
#  the following is the result of the computation using the default values.
command_mem_mb=10240

declare -a disabled_read_filters_arr
if [ -n "$disabled_read_filters" ]; then
  disabled_read_filters_arr=( "--disable-read-filter" "${disabled_read_filters}" )
else
  disabled_read_filters_arr=()
fi

export GATK_LOCAL_JAR="$gatk4_jar_override"

java -Xmx${command_mem_mb}m -jar /opt/gatk.jar CollectReadCounts \
  -L "${intervals}" \
  --input "${cram_or_bam}" \
  --read-index "${cram_or_bam_idx}" \
  --reference "${ref_fasta}" \
  --format TSV \
  --interval-merging-rule OVERLAPPING_ONLY \
  --output "${sample_id}".counts.tsv \
  "${disabled_read_filters_arr[@]}"

sed -ri "s/@RG\tID:GATKCopyNumber\tSM:.+/@RG\tID:GATKCopyNumber\tSM:${sample_id}/g" "${sample_id}.counts.tsv"
bgzip "${sample_id}.counts.tsv"
