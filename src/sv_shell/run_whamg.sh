#!/bin/bash

# This script is a bash implementation of the following workflow/task in WDL:
# Filename: wdl/Whamg.wdl
# Workflow: Whamg

set -Exeuo pipefail

sample_id=$1
cram_file=$2
cram_index=$3
reference_fasta=$4
reference_index=$5
include_bed_file=$6
primary_contigs_list=$7
outputs_json_filename=$8
cpu_cores=${9:-4}

include_bed_file_abs_path=$(realpath $include_bed_file)
reference_fasta_abs_path=$(realpath $reference_fasta)

working_dir=$(mktemp -d /wd_collect_counts_XXXXXXXX)
working_dir="$(realpath ${working_dir})"
output_dir=$(mktemp -d /output_collect_counts_XXXXXXXX)
output_dir="$(realpath ${output_dir})"
cd "${working_dir}"

chr_list=$(paste -s -d ',' "${primary_contigs_list}")

# print some info that may be useful for debugging
df -h
echo "whamg $(whamg 2>&1 | grep Version)"

echo "Converting cram to bam ..."
# covert cram to bam
samtools view -b1@ ${cpu_cores} -T "${reference_fasta}" "${cram_file}" > sample.bam
echo "Finished converting cram to bam."

# index bam file
echo "Indexing the bam file ..."
samtools index -@ ${cpu_cores} sample.bam
echo "Finished indexing the bam file."

# ensure that index files are present in appropriate locations
ln -s sample.bam.bai sample.bai
if [ ! -e "${reference_fasta}.fai" ]; then
  ln -s "${reference_index}" "${reference_fasta}.fai"
fi

# run whamg on all specified intervals
echo "Running whamg on specified intervals ..."
mkdir -p tmpVcfs
cd tmpVcfs

awk 'BEGIN{FS=OFS="\t"}{printf("%07d\t%s\n",NR,$1":"$2"-"$3)}' "${include_bed_file_abs_path}" |\
  while read -r line interval; do
    vcfFile="$line.wham.vcf.gz"
    whamg \
        -c "${chr_list}" \
        -x ${cpu_cores} \
        -a "${reference_fasta_abs_path}" \
        -f ../sample.bam \
        -r $interval \
      | bgzip -c > $vcfFile
    bcftools index -t $vcfFile
  done
echo "Finished running whamg on specified intervals."

# We need to update both the VCF sample ID and the TAGS INFO field in the WHAM output VCFs.
# WHAM uses both to store the sample identifier, and by default uses the SM identifier from the BAM file.
# We need to update both to reflect the potentially-renamed sample identifier used by the pipeline (sample_id) --
# svtk standardize_vcf uses the TAGS field to identify the sample for WHAM VCFs.

ls -1 *.wham.vcf.gz > vcf.list
bcftools concat -a -Ov -f vcf.list | \
  sed -e 's/^#CHROM\t.*/#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t${sample_id}/' -e 's/;TAGS=[^;]*;/;TAGS=${sample_id};/' | \
  bgzip -c > "../${sample_id}.wham.vcf.gz"
cd ..
bcftools index -t "${sample_id}.wham.vcf.gz"

vcf_filename="${output_dir}/${sample_id}.wham.vcf.gz"
mv "${sample_id}.wham.vcf.gz" "${vcf_filename}"
index_filename="${output_dir}/${sample_id}.wham.vcf.gz.tbi"
mv "${sample_id}.wham.vcf.gz.tbi" "${index_filename}"

outputs_filename="${output_dir}/outputs.json"
outputs_json=$(jq -n \
  --arg vcf "${vcf_filename}" \
  --arg index "${index_filename}" \
  '{vcf: $vcf, index: $index}' )
echo "${outputs_json}" > "${outputs_filename}"
cp "${outputs_filename}" "${outputs_json_filename}"
