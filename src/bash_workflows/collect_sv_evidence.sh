#!/bin/bash

set -Eeuo pipefail

sample_id=$1
bam_or_cram_file=$2
bam_or_cram_index=$3
reference_fasta=$4
reference_index=$5
reference_dict=$6
sd_locs_vcf=$7
preprocessed_intervals=$8
site_depth_min_mapq=${9:-6}
site_depth_min_baseq=${10:-10}
primary_contigs_list="${11:-}"
gatk_jar_override="${12:-/root/gatk.jar}"
command_mem_mb=${13:-3250}

export GATK_LOCAL_JAR="$gatk_jar_override"

primary_contigs_arg=""
if [[ -n "$primary_contigs_list" ]]; then
  primary_contigs_arg="-L $primary_contigs_list"
fi

java -Xmx${command_mem_mb}m -jar /opt/gatk.jar CollectSVEvidence \
    -I "${bam_or_cram_file}" \
    --sample-name "${sample_id}" \
    -F "${sd_locs_vcf}" \
    -SR "${sample_id}.sr.txt.gz" \
    -PE "${sample_id}.pe.txt.gz" \
    -SD "${sample_id}.sd.txt.gz" \
    --site-depth-min-mapq "${site_depth_min_mapq}" \
    --site-depth-min-baseq "${site_depth_min_baseq}" \
    -R "${reference_fasta}" \
    ${primary_contigs_arg} \
    --read-filter NonZeroReferenceLengthAlignmentReadFilter
