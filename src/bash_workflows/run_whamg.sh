#!/bin/bash

set -Eeuo pipefail

sample_id=$1
cram_file=$2
cram_index=$3
reference_fasta=$4
reference_index=$5
include_bed_file=$6
primary_contigs_list=$7
cpu_cores=${8:-4}

include_bed_file_abs_path=$(realpath $include_bed_file)
reference_fasta_abs_path=$(realpath $reference_fasta)

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
