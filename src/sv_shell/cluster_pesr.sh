#!/bin/bash

set -Exeuo pipefail


function ExcludeIntervalsByEndpoints() {
  local _vcf=$1
  local _reference_fasta_fai=$2
  local _intervals=$3
  local _output_prefix=$4

  cut -f1,2 "${_reference_fasta_fai}" > genome.file
  bcftools query -f \
    '%CHROM\t%POS\t%POS\t%ID\t%SVTYPE\n%CHROM\t%END\t%END\t%ID\t%SVTYPE\n%CHR2\t%END2\t%END2\t%ID\t%SVTYPE\n' \
    "${_vcf}" \
    | awk '$1!="."' \
    | sort -k1,1V -k2,2n -k3,3n \
    > ends.bed
  bedtools intersect -sorted -u -wa -g genome.file -wa -a ends.bed -b "${_intervals}" \
    | cut -f4 \
    | sort \
    | uniq \
    > excluded_vids.list
  bcftools view -i '%ID!=@excluded_vids.list' "${_vcf}" -Oz -o "${_output_prefix}.vcf.gz"
  tabix "${_output_prefix}.vcf.gz"
}

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
reference_fasta=$(jq -r ".reference_fasta" "${input_json}")
reference_fasta_fai=$(jq -r ".reference_fasta_fai" "${input_json}")
reference_dict=$(jq -r ".reference_dict" "${input_json}")
vcf_tar=$(jq -r ".vcf_tar" "${input_json}")
ploidy_table=$(jq -r ".ploidy_table" "${input_json}")
exclude_intervals=$(jq -r ".exclude_intervals" "${input_json}")
min_size=$(jq -r ".min_size" "${input_json}")
contig_list=$(jq -r '.contig_list // ""' "${input_json}")
contig_subset_list=$(jq -r '.contig_subset_list // ""' "${input_json}")
clustering_algorithm=$(jq -r '.clustering_algorithm // ""' "${input_json}")
pesr_interval_overlap=$(jq -r '.pesr_interval_overlap' "${input_json}")
pesr_breakend_window=$(jq -r '.pesr_breakend_window' "${input_json}")
java_mem_fraction=$(jq -r '.java_mem_fraction // 0.85' "$input_json")


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
  bcftools view -i "ID!=@excluded_vids.list && (INFO/SVLEN='.' || INFO/SVLEN=-1 || INFO/SVLEN>=${min_size})" tmp.vcf.gz \
    -Oz -o "out/${SAMPLE_NUM}.${NAME}.vcf.gz"
  tabix out/$SAMPLE_NUM.$NAME.vcf.gz
  i=$((i+1))
done
tar czf "${output_prefix}.tar.gz" -C out/ .

prepare_pesr_vcf_out="$(realpath ${output_prefix}.tar.gz)"

# ---------------------------------------------------------------------------------------------------------------------

svtk_format_vcfs=()
svtk_format_vcf_indexes=()

contigs_list="${contig_subset_list:-$contig_list}"
contigs=()
while IFS= read -r line; do
    contigs+=("$line")
done < "${contigs_list}"

for contig in "${contigs[@]}"; do

  # SVCluster
  # -------------------------------------------------------------------------------------------------------------------
  sv_cluster_output_dir=$(mktemp -d /output_sv_cluster_XXXXXXXX)
  sv_cluster_output_dir="$(realpath ${sv_cluster_output_dir})"
  sv_cluster_inputs_json="$(realpath "${sv_cluster_output_dir}/sv_cluster_inputs.json")"
  sv_cluster_outputs_json="$(realpath "${sv_cluster_output_dir}/sv_cluster_outputs.json")"

  jq -n \
    --arg vcfs_tar "${prepare_pesr_vcf_out}" \
    --arg ploidy_table "${ploidy_table}" \
    --arg output_prefix "${batch}.cluster_batch.${caller}.${contig}.clustered" \
    --arg contig "${contig}" \
    --arg fast_mode true \
    --arg algorithm "${clustering_algorithm}" \
    --argjson pesr_sample_overlap 0 \
    --argjson pesr_interval_overlap ${pesr_interval_overlap} \
    --argjson pesr_breakend_window ${pesr_breakend_window} \
    --arg reference_fasta "${reference_fasta}" \
    --arg reference_fasta_fai "${reference_fasta_fai}" \
    --arg reference_dict "${reference_dict}" \
    --arg java_mem_fraction "${java_mem_fraction}" \
    --arg variant_prefix "${batch}_${caller}_${contig}_" \
    '{
        "vcfs_tar": $vcfs_tar,
        "ploidy_table": $ploidy_table,
        "output_prefix": $output_prefix,
        "contig": $contig,
        "fast_mode": $fast_mode,
        "algorithm": $algorithm,
        "pesr_sample_overlap": $pesr_sample_overlap,
        "pesr_interval_overlap": $pesr_interval_overlap,
        "pesr_breakend_window": $pesr_breakend_window,
        "reference_fasta": $reference_fasta,
        "reference_fasta_fai": $reference_fasta_fai,
        "reference_dict": $reference_dict,
        "java_mem_fraction": $java_mem_fraction,
    }' > "${sv_cluster_inputs_json}"

  bash /opt/sv_shell/sv_cluster.sh "${sv_cluster_inputs_json}" "${sv_cluster_outputs_json}" "${sv_cluster_output_dir}"

  sv_cluster_out=$(jq -r ".out" "${sv_cluster_outputs_json}")


  # ExcludeIntervalsByEndpoints
  # -------------------------------------------------------------------------------------------------------------------
  # Remove variants from VCF that overlap with regions in the exclude_intervals file.
  working_dir=$(mktemp -d /wd_ExcludeIntervalsByEndpoints_XXXXXXXX)
  working_dir="$(realpath ${working_dir})"
  cd "${working_dir}"
  output_prefix="${batch}.cluster_batch.${caller}.${contig}.exclude_intervals"
  ExcludeIntervalsByEndpoints "${sv_cluster_out}" "${reference_fasta_fai}" "${exclude_intervals}" "${output_prefix}"

  out_exclud_intervas_by_endpoint="${working_dir}/${output_prefix}.vcf.gz"
  out_index_exclud_intervas_by_endpoint="${working_dir}/${output_prefix}.vcf.gz.tbi"

  echo "${out_exclud_intervas_by_endpoint}"


  # GatkToSvtkVcf
  # -------------------------------------------------------------------------------------------------------------------
  working_dir=$(mktemp -d /wd_GatkToSvtkVcf_XXXXXXXX)
  working_dir="$(realpath ${working_dir})"
  cd "${working_dir}"
  output_prefix="${batch}.cluster_batch.${caller}.${contig}.svtk_formatted"
  python /opt/sv-pipeline/scripts/format_gatk_vcf_for_svtk.py \
    --vcf "${out_exclud_intervas_by_endpoint}" \
    --out "${output_prefix}.vcf.gz" \
    --source "${caller}" \
    --contigs "${contig_list}" \
    --remove-format CN

  svtk_format_vcfs+=("${working_dir}/${output_prefix}.vcf.gz")
  svtk_format_vcf_indexes+=("${working_dir}/${output_prefix}.vcf.gz.tbi")
done




