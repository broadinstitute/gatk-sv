#!/bin/bash

set -Exeuo pipefail

function BedtoolsClosest()
{
  local _bed_a=$1
  local _bed_b=$2
  local _out=$3

  paste <(head -1 "${_bed_a}") <(head -1 "${_bed_b}") | sed -e "s/#//g" > "${_out}"
  bedtools closest -wo -a <(sort -k1,1 -k2,2n "${_bed_a}") -b <(sort -k1,1 -k2,2n "${_bed_b}") >> "${_out}"
}

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
external_af_population=($(jq -r '.external_af_population[]' "$input_json"))
external_af_ref_prefix=$(jq -r '.external_af_ref_prefix' "$input_json")

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

# AnnotateExternalAFPerShard
# ---------------------------------------------------------------------------------------------------------------------

# SplitQueryVcf
# ---------------------------------
# splits a VCF into smaller BED files each containing a set type of variants.

SplitQueryVcf_prefix="$(basename "${ComputeAFs_af_vcf}" .vcf.gz)"
SplitQueryVcf_bed="$(realpath "${SplitQueryVcf_prefix}.bed")"
SplitQueryVcf_del="$(realpath "${SplitQueryVcf_prefix}.DEL.bed")"
SplitQueryVcf_dup="$(realpath "${SplitQueryVcf_prefix}.DUP.bed")"
SplitQueryVcf_ins="$(realpath "${SplitQueryVcf_prefix}.INS.bed")"
SplitQueryVcf_inv="$(realpath "${SplitQueryVcf_prefix}.INV_CPX.bed")"
SplitQueryVcf_bnd="$(realpath "${SplitQueryVcf_prefix}.BND_CTX.bed")"

svtk vcf2bed -i SVTYPE -i SVLEN "${ComputeAFs_af_vcf}" tmp.bed
cut -f1-4,7-8 tmp.bed > "${SplitQueryVcf_bed}"

head -1 "${SplitQueryVcf_bed}" > header

cat header <(awk '{if ($5=="DEL") print}' "${SplitQueryVcf_bed}" )> "${SplitQueryVcf_del}"
cat header <(awk '{if ($5=="DUP") print}' "${SplitQueryVcf_bed}" )> "${SplitQueryVcf_dup}"
cat header <(awk '{if ($5=="INS" || $5=="INS:ME" || $5=="INS:ME:ALU" || $5=="INS:ME:LINE1" || $5=="INS:ME:SVA" || $5=="ALU" || $5=="LINE1" || $5=="SVA" || $5=="HERVK" ) print}' "${SplitQueryVcf_bed}" )> "${SplitQueryVcf_ins}"
cat header <(awk '{if ($5=="INV" || $5=="CPX") print}' "${SplitQueryVcf_bed}" )> "${SplitQueryVcf_inv}"
cat header <(awk '{if ($5=="BND" || $5=="CTX") print}' "${SplitQueryVcf_bed}" )> "${SplitQueryVcf_bnd}"


# compare_[del, dup, ins, inv, bnd]
# ---------------------------------

compare_del_bed="$(realpath bedtools_closest_out_del.bed)"
compare_dup_bed="$(realpath bedtools_closest_out_dup.bed)"
compare_ins_bed="$(realpath bedtools_closest_out_ins.bed)"
compare_inv_bed="$(realpath bedtools_closest_out_inv.bed)"
compare_bnd_bed="$(realpath bedtools_closest_out_bnd.bed)"

BedtoolsClosest "${SplitQueryVcf_del}" "${ref_bed_del}" "${compare_del_bed}"
BedtoolsClosest "${SplitQueryVcf_dup}" "${ref_bed_dup}" "${compare_dup_bed}"
BedtoolsClosest "${SplitQueryVcf_ins}" "${ref_bed_ins}" "${compare_ins_bed}"
BedtoolsClosest "${SplitQueryVcf_inv}" "${ref_bed_inv}" "${compare_inv_bed}"
BedtoolsClosest "${SplitQueryVcf_bnd}" "${ref_bed_bnd}" "${compare_bnd_bed}"


printf "%s\n" "${external_af_population[@]}" > pop_list.txt


# SelectMatchedSVs
# ---------------------------------

# calcu_del
calcu_del_output_comp="${compare_del_bed%.bed}.comparison"
Rscript /opt/sv-pipeline/05_annotation/scripts/R1.bedtools_closest_CNV.R \
  -i "${compare_del_bed}" \
  -o "${calcu_del_output_comp}" \
  -p pop_list.txt

# calcu_dup
calcu_dup_output_comp="${compare_dup_bed%.bed}.comparison"
Rscript /opt/sv-pipeline/05_annotation/scripts/R1.bedtools_closest_CNV.R \
  -i "${compare_dup_bed}" \
  -o "${calcu_dup_output_comp}" \
  -p pop_list.txt

# calcu_inv
calcu_inv_output_comp="${compare_inv_bed%.bed}.comparison"
Rscript /opt/sv-pipeline/05_annotation/scripts/R1.bedtools_closest_CNV.R \
  -i "${compare_inv_bed}" \
  -o "${calcu_inv_output_comp}" \
  -p pop_list.txt


# SelectMatchedINSs
# ---------------------------------

# calcu_ins
calcu_ins_output_comp="${compare_ins_bed%.bed}.comparison"
Rscript /opt/sv-pipeline/05_annotation/scripts/R2.bedtools_closest_INS.R \
  -i "${compare_ins_bed}" \
  -o "${calcu_ins_output_comp}" \
  -p pop_list.txt

# calcu_ins
calcu_bnd_output_comp="${compare_bnd_bed%.bed}.comparison"
Rscript /opt/sv-pipeline/05_annotation/scripts/R2.bedtools_closest_INS.R \
  -i "${compare_bnd_bed}" \
  -o "${calcu_bnd_output_comp}" \
  -p pop_list.txt


# ModifyVcf
# ---------------------------------------------------------------------------------------------------------------------

cat "${calcu_del_output_comp}" > labeled.bed
cat "${calcu_dup_output_comp}" >> labeled.bed
cat "${calcu_ins_output_comp}" >> labeled.bed
cat "${calcu_inv_output_comp}" >> labeled.bed
cat "${calcu_bnd_output_comp}" >> labeled.bed

ModifyVcf_annotated_vcf="$(realpath "${prefix}.annotated.vcf")"

python /opt/sv_shell/annotate_external_af_modify_vcf.py \
  --vcf "${ComputeAFs_af_vcf}" \
  --ref-prefix "${external_af_ref_prefix}" \
  --labeled-bed labeled.bed \
  --output-filename "${ModifyVcf_annotated_vcf}"

bgzip "${ModifyVcf_annotated_vcf}"
ModifyVcf_annotated_vcf="${ModifyVcf_annotated_vcf}.gz"
tabix "${ModifyVcf_annotated_vcf}"



# -------------------------------------------------------
# ======================= Output ========================
# -------------------------------------------------------

ModifyVcf_annotated_vcf_output_dir="${output_dir}/$(basename "${ModifyVcf_annotated_vcf}")"
mv "${ModifyVcf_annotated_vcf}" "${ModifyVcf_annotated_vcf_output_dir}"
mv "${ModifyVcf_annotated_vcf}.tbi" "${ModifyVcf_annotated_vcf_output_dir}.tbi"

jq -n \
  --arg annotated_vcf "${ModifyVcf_annotated_vcf_output_dir}" \
  --arg annotated_vcf_tbi "${ModifyVcf_annotated_vcf_output_dir}.tbi" \
  '{
      annotated_vcf: $annotated_vcf,
      annotated_vcf_tbi: $annotated_vcf_tbi
  }' > "${output_json_filename}"

echo "Finished AnnotateVcf successfully, output json filename: ${output_json_filename}"
