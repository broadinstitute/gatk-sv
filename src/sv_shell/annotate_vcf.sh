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
  output_dir=$(mktemp -d /output_annotate_vcf_XXXXXXXX)
else
  mkdir -p "${output_dir}"
fi
output_dir="$(realpath ${output_dir})"

if [ -z "${output_json_filename}" ]; then
  output_json_filename="${output_dir}/output.json"
else
  output_json_filename="$(realpath ${output_json_filename})"
fi

working_dir=$(mktemp -d /wd_annotate_vcf_XXXXXXXX)
working_dir="$(realpath ${working_dir})"
cd "${working_dir}"
echo "Median coverage Working directory: ${working_dir}"


vcf=$(jq -r '.vcf' "$input_json")
prefix=$(jq -r '.prefix' "$input_json")
external_af_ref_bed=$(jq -r '.external_af_ref_bed' "$input_json")
par_bed=$(jq -r '.par_bed' "$input_json")

# Note that the following files are optional in the WDL version,
# however, since they are used in the single-sample pipeline,
# they are made required for simplicity.
protein_coding_gtf=$(jq -r '.protein_coding_gtf' "$input_json")
noncoding_bed=$(jq -r '.noncoding_bed' "$input_json")
# in the same context, `promoter_window` and `max_breakend_as_cnv_length`
# are excluded as they are not provided by default in the single-sample mode.

java_mem_fraction=$(jq -r '.java_mem_fraction // ""' "$input_json")


function getJavaMem() {
  # get JVM memory in MiB by getting total memory from /proc/meminfo
  # and multiplying by java_mem_fraction

  local mem_fraction=${java_mem_fraction:=0.85}
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

# SplitRefBed
# ---------------------------------------------------------------------------------------------------------------------
# TODO: assuming you always provide external_af_ref_bed

split_ref_bed_working_dir=$(mktemp -d /wd_split_ref_bed_XXXXXXXX)
split_ref_bed_working_dir="$(realpath ${split_ref_bed_working_dir})"
cd "${split_ref_bed_working_dir}"

prefix="$(basename "${external_af_ref_bed}" .bed.gz)"
ref_bed_del="$(realpath "${prefix}.DEL.bed")"
ref_bed_dup="$(realpath "${prefix}.DUP.bed")"
ref_bed_ins="$(realpath "${prefix}.INS.bed")"
ref_bed_inv="$(realpath "${prefix}.INV_CPX.bed")"
ref_bed_bnd="$(realpath "${prefix}.BND_CTX.bed")"

zcat "${external_af_ref_bed}" | head -1 > header || true

cat header <(zcat "${external_af_ref_bed}" | awk '{if ($6=="DEL") print}') > "${ref_bed_del}"
cat header <(zcat "${external_af_ref_bed}" | awk '{if ($6=="DUP") print}') > "${ref_bed_dup}"
cat header <(zcat "${external_af_ref_bed}" | awk '{if ($6=="INS" || $6=="INS:ME" || $6=="INS:ME:ALU" || $6=="INS:ME:LINE1" || $6=="INS:ME:SVA" || $6=="ALU" || $6=="LINE1" || $6=="SVA" || $6=="HERVK" ) print}') > "${ref_bed_ins}"
cat header <(zcat "${external_af_ref_bed}" | awk '{if ($6=="INV" || $6=="CPX") print}' ) > "${ref_bed_inv}"
cat header <(zcat "${external_af_ref_bed}" | awk '{if ($6=="BND" || $6=="CTX") print}' ) > "${ref_bed_bnd}"


# AnnotateFunctionalConsequences
# ---------------------------------------------------------------------------------------------------------------------

AnnotateFunctionalConsequences_annotated_vcf="$(realpath "${prefix}.annotated.vcf.gz")"
AnnotateFunctionalConsequences_annotated_vcf_index="$(realpath "${prefix}.annotated.vcf.gz.tbi")"

tabix -f -p vcf "${vcf}"

java "-Xmx${JVM_MAX_MEM}" -jar /opt/gatk.jar SVAnnotate \
  -V "${vcf}" \
  -O "${AnnotateFunctionalConsequences_annotated_vcf}" \
  --protein-coding-gtf "${protein_coding_gtf}" \
  --non-coding-bed "${noncoding_bed}"


# ComputeAFs
# ---------------------------------------------------------------------------------------------------------------------

ComputeAFs_af_vcf="${prefix}.wAFs.vcf.gz"
ComputeAFs_af_vcf_idx="${prefix}.wAFs.vcf.gz.tbi"

/opt/sv-pipeline/05_annotation/scripts/compute_AFs.py "${AnnotateFunctionalConsequences_annotated_vcf}" stdout \
  --par "${par_bed}" \
| bgzip -c \
> "${ComputeAFs_af_vcf}"

tabix -p vcf "${ComputeAFs_af_vcf}"
