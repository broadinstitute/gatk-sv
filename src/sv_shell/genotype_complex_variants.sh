#!/bin/bash

set -Exeuo pipefail

# -------------------------------------------------------
# ==================== Input & Setup ====================
# -------------------------------------------------------


input_json=${1}
output_json_filename=${2-""}
output_dir=${3:-""}

input_json="$(realpath ${input_json})"

if [ -z "${output_dir}" ]; then
  output_dir=$(mktemp -d /output_genotype_complex_variants_XXXXXXXX)
else
  mkdir -p "${output_dir}"
fi
output_dir="$(realpath ${output_dir})"

if [ -z "${output_json_filename}" ]; then
  output_json_filename="${output_dir}/output.json"
else
  output_json_filename="$(realpath ${output_json_filename})"
fi

working_dir=$(mktemp -d /wd_genotype_complex_variants_XXXXXXXX)
working_dir="$(realpath ${working_dir})"
cd "${working_dir}"
echo "Genotype complex variants working directory: ${working_dir}"


depth_vcfs=$(jq -r ".depth_vcfs" "$input_json")
ped_file=$(jq -r ".ped_file" "$input_json")
batches=$(jq -r ".batches" "$input_json")
complex_resolve_vcfs=$(jq -r ".complex_resolve_vcfs" "$input_json")
cohort_name=$(jq -r ".cohort_name" "$input_json")
median_coverage_files=$(jq -r ".median_coverage_files" "$input_json")
bin_exclude=$(jq -r ".bin_exclude" "$input_json")
ref_dict=$(jq -r ".ref_dict" "$input_json")
bincov_files=$(jq -r ".bincov_files" "$input_json")
depth_gt_rd_sep_files=$(jq -r ".depth_gt_rd_sep_files" "$input_json")

n_rd_test_bins=100000

function getJavaMem() {
  # get JVM memory in MiB by getting total memory from /proc/meminfo
  # and multiplying by java_mem_fraction

  local mem_fraction=${java_mem_fraction:=0.6}
  cat /proc/meminfo | \
    awk -v MEM_FIELD="$1" -v frac="${mem_fraction}" '{
      f[substr($1, 1, length($1)-1)] = $2
    } END {
      printf "%dM", f[MEM_FIELD] * frac / 1024
    }'
}
JVM_MAX_MEM=$(getJavaMem MemTotal)
echo "JVM memory: $JVM_MAX_MEM"

# -------------------------------------------------------
# ======================= Command =======================
# -------------------------------------------------------

# GetSampleIdsFromVcf
# ---------------------------------------------------------------------------------------------------------------------
sample_list="$(basename "${depth_vcfs}" .vcf.gz).samples.txt"
bcftools query -l "${depth_vcfs}" > "${sample_list}"


# SubsetPedFile
# ---------------------------------------------------------------------------------------------------------------------
ped_subset_filename="$(basename "${ped_file}" .ped).${batches}.ped"
awk 'FNR==NR {a[$1]; next}; $2 in a' "${sample_list}" "${ped_file}" > "${ped_subset_filename}"


# ScatterCpxGenotyping
# ---------------------------------------------------------------------------------------------------------------------

# GenotypeCpxCnvs
# Workflow to perform depth-based genotyping for a single vcf shard scattered
# across batches on predicted CPX CNVs

# Convert VCF to bed of CPX CNV intervals
# Get CNV intervals from complex SV for depth genotyping
GetCpxCnvIntervals_cpx_cnv_bed="${cohort_name}.complex_CNV_intervals_to_test.bed.gz"
bash /opt/sv-pipeline/04_variant_resolution/scripts/gather_cpx_intervals_for_rd_gt.sh \
  "${complex_resolve_vcfs}" \
  "${GetCpxCnvIntervals_cpx_cnv_bed}"


# GenotypeCpxCnvsPerBatch
# ---------------------------------------------------------------------------------------------------------------------
# GetSampleIdsFromMedianCoverageFile
GetSampleIdsFromMedianCoverageFile_out_file="$(realpath "${batches}.samples.txt")"
head -1 "${median_coverage_files}" | sed -e 's/\t/\n/g' > "${GetSampleIdsFromMedianCoverageFile_out_file}"


# RdTestGenotype
# ---------------------------------------------------------------------------------------------------------------------
bed="${cohort_name}.complex_CNV_intervals_to_test.bed"
gunzip -c "${GetCpxCnvIntervals_cpx_cnv_bed}" > "${bed}"


#bincov_files="/inputs/tmp_rdgenotype/NA12878.RD.txt.gz"


RdTestGenotype_prefix=$(basename "${bed}" .bed)

grep -v "^#" "${bed}" | sort -k1,1V -k2,2n | bedtools merge -i stdin -d 1000000 > merged.bed

if [ -s merged.bed ]; then
  java "-Xmx${JVM_MAX_MEM}" -jar /opt/gatk.jar PrintSVEvidence \
    --sequence-dictionary "${ref_dict}" \
    --evidence-file "${bincov_files}" \
    -L merged.bed \
    -O local.RD.txt.gz

  tabix -f -0 -s1 -b2 -e3 local.RD.txt.gz
else
  touch local.RD.txt
  bgzip local.RD.txt
  tabix -p bed local.RD.txt.gz
fi

# just make sure the file is indexed, no need to index it in-place
# tabix -f -p bed "${bin_exclude}"

Rscript /opt/RdTest/RdTest.R \
  -b "${bed}" \
  -c local.RD.txt.gz \
  -m "${median_coverage_files}" \
  -f "${ped_file}" \
  -n "${RdTestGenotype_prefix}" \
  -w "${GetSampleIdsFromMedianCoverageFile_out_file}" \
  -i ${n_rd_test_bins} \
  -r "${depth_gt_rd_sep_files}" \
  -y "${bin_exclude}" \
  -g TRUE

if [ -f "${RdTestGenotype_prefix}.geno" ] && [ -f "${RdTestGenotype_prefix}.gq" ] ; then
  python /opt/sv-pipeline/04_variant_resolution/scripts/merge_RdTest_genotypes.py \
    "${RdTestGenotype_prefix}.geno" "${RdTestGenotype_prefix}.gq" rd.geno.cnv.bed
  sort -k1,1V -k2,2n rd.geno.cnv.bed | uniq | bgzip -c > rd.geno.cnv.bed.gz
else
  # In case RdTest does not produce output because the input is empty
  echo "" | bgzip -c > rd.geno.cnv.bed.gz
fi

RdTestGenotype_melted_genotypes="$(realpath "rd.geno.cnv.bed.gz")"

zcat "${RdTestGenotype_melted_genotypes}" \
  | sort -Vk1,1 -k2,2n -k3,3n \
  | bgzip -c \
  > "${cohort_name}.rd_genos.bed.gz"

# TODO: this is temp, after you resolved the above scattering, decide where you get this, maybe `RdTestGenotype_melted_genotypes`
#GenotypeBatch_genotypes="/inputs/tmp_rdgenotype/NA12878.rd_genos.bed.gz"
GenotypeBatch_genotypes="$(realpath "${cohort_name}.rd_genos.bed.gz")"


# ZcatCompressedFiles
# ---------------------------------------------------------------------------------------------------------------------
# the following is a far simplified version of the WDL implementation for a case where you have only one file
ZcatCompressedFiles_outfile="${cohort_name}.CPX_intervals.merged_rd_genos.bed.gz"
zcat "${GenotypeBatch_genotypes}" \
  | sort -Vk1,1 -k2,2n -k3,3n -k4,4V -k5,5V \
  | bgzip -c \
  > "${ZcatCompressedFiles_outfile}"


# ParseGenotypes
# ---------------------------------------------------------------------------------------------------------------------

# TODO TEMP
#GetCpxCnvIntervals_cpx_cnv_bed="/inputs/tmp_parsegenotypes/NA12878.chr1.chr1.complex_CNV_intervals_to_test.bed.gz"

python /opt/sv-pipeline/04_variant_resolution/scripts/process_posthoc_cpx_depth_regenotyping.py \
  --vcf "${complex_resolve_vcfs}" \
  --intervals "${GetCpxCnvIntervals_cpx_cnv_bed}" \
  --genotypes "${ZcatCompressedFiles_outfile}" \
  --ped "${ped_file}" \
  --out out.vcf.gz \
  --reclassification-table "${cohort_name}.CPXregenotyping_reclassification_table.txt"

# Compress for storage
gzip "${cohort_name}.CPXregenotyping_reclassification_table.txt"

# Re-sort and index since coordinates may have changed
mkdir temp
bcftools sort \
  --temp-dir temp \
  --output-type z \
  --output-file "${cohort_name}.postCPXregenotyping.vcf.gz" \
  out.vcf.gz
tabix "${cohort_name}.postCPXregenotyping.vcf.gz"

ParseGenotypes_cpx_depth_gt_resolved_vcf="$(realpath "${cohort_name}.postCPXregenotyping.vcf.gz")"
ParseGenotypes_cpx_depth_gt_resolved_vcf_idx="${ParseGenotypes_cpx_depth_gt_resolved_vcf}.tbi"
ParseGenotypes_reclassification_table="$(realpath "${cohort_name}.CPXregenotyping_reclassification_table.txt.gz")"



# TODO: if you end up not sharding above, basically this workflow is done here,
#  you can produce output, because all other steps remaining are merging outputs of various shards
