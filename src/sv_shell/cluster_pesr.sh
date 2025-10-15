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
  output_dir=$(mktemp -d /output_cluster_pesr_XXXXXXXX)
else
  mkdir -p "${output_dir}"
fi
output_dir="$(realpath ${output_dir})"

if [ -z "${output_json_filename}" ]; then
  output_json_filename="${output_dir}/output.json"
else
  output_json_filename="$(realpath ${output_json_filename})"
fi

working_dir=$(mktemp -d /wd_cluster_pesr_XXXXXXXX)
working_dir="$(realpath ${working_dir})"
cd "${working_dir}"

batch=$(jq -r ".batch" "${input_json}")
caller=$(jq -r ".caller" "${input_json}")
reference_fasta_fai=$(jq -r ".reference_fasta_fai" "${input_json}")
vcf_tar=$(jq -r ".vcf_tar" "${input_json}")
ploidy_table=$(jq -r ".ploidy_table" "${input_json}")
exclude_intervals=$(jq -r ".exclude_intervals" "${input_json}")
min_size=$(jq -r ".min_size" "${input_json}")


# -------------------------------------------------------
# ======================= Command =======================
# -------------------------------------------------------

# PreparePESRVcfs
# ---------------------------------------------------------------------------------------------------------------------

output_prefix="${batch}.cluster_batch.${caller}.prep_vcfs"

cut -f1,2 "${reference_fasta_fai}" > genome.file
mkdir in/ out/
tar xzf "${vcf_tar}" -C in/
i=0
for VCF in in/*.vcf.gz; do
  NAME=$(basename $VCF .vcf.gz)
  SAMPLE_NUM=`printf %05d $i`
  # Convert format
  python /opt/sv-pipeline/scripts/format_svtk_vcf_for_gatk.py \
    --vcf $VCF \
    --out tmp.vcf.gz \
    --ploidy-table "${ploidy_table}"
  # Interval, contig, and size filtering
  bcftools query -f '%CHROM\t%POS\t%POS\t%ID\t%SVTYPE\n%CHROM\t%END\t%END\t%ID\t%SVTYPE\n%CHR2\t%END2\t%END2\t%ID\t%SVTYPE\n' tmp.vcf.gz \
    | awk '$1!="." && $2!="."' \
    | sort -k1,1V -k2,2n -k3,3n \
    > ends.bed
  bedtools intersect -sorted -u -wa -g genome.file -wa -a ends.bed -b "${exclude_intervals}" | cut -f4 | sort | uniq \
    > excluded_vids.list
  bcftools view -i "ID!=@excluded_vids.list && (INFO/SVLEN='.' || INFO/SVLEN=-1 || abs(INFO/SVLEN)>=${min_size})" tmp.vcf.gz \
    -Oz -o "out/${SAMPLE_NUM}.${NAME}.vcf.gz"
  tabix out/$SAMPLE_NUM.$NAME.vcf.gz
  i=$((i+1))
done
tar czf "${output_prefix}.tar.gz" -C out/ .
