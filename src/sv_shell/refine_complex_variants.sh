#!/bin/bash

set -Exeuo pipefail

function SeekDepthSuppForCpx()
{
  local _cpx_lg_cnv=$1
  local _raw_depth_bed=$2
  local _prefix=$3

  zcat "${_cpx_lg_cnv}" | awk '{print $6}' | sort | uniq > sample_list.tsv

  echo -e '#chr\tpos\tend\tSVTYPE\tSVID\tsample\tbatch\tdepth_cov' > "${_prefix}.depth_supp.bed"

  while read sample_name; do
    zcat "${_cpx_lg_cnv}"    | awk -v sample="$sample_name" '{if ($6==sample) print}' > query.bed
    zcat "${_raw_depth_bed}" | awk -v sample="$sample_name" '{if ($5==sample) print}' > ref.bed
    bedtools coverage -a query.bed -b ref.bed \
      | awk '{print $1,$2,$3,$4,$5,$6,$7,$NF}' \
      | sed -e 's/ /\t/g' \
      >> "${_prefix}.depth_supp.bed"
  done < sample_list.tsv

  bgzip "${_prefix}.depth_supp.bed"
}

# -------------------------------------------------------
# ==================== Input & Setup ====================
# -------------------------------------------------------


input_json=${1}
output_json_filename=${2-""}
output_dir=${3:-""}

input_json="$(realpath ${input_json})"

if [ -z "${output_dir}" ]; then
  output_dir=$(mktemp -d /output_refine_complex_variants_XXXXXXXX)
else
  mkdir -p "${output_dir}"
fi
output_dir="$(realpath ${output_dir})"

if [ -z "${output_json_filename}" ]; then
  output_json_filename="${output_dir}/output.json"
else
  output_json_filename="$(realpath ${output_json_filename})"
fi

working_dir=$(mktemp -d /wd_refine_complex_variants_XXXXXXXX)
working_dir="$(realpath ${working_dir})"
cd "${working_dir}"
echo "Refine complex variants working directory: ${working_dir}"

prefix=$(jq -r ".prefix" "$input_json")
batch_name=$(jq -r ".batch_name_list" "$input_json")
batch_sample_lists=$(jq -r ".batch_sample_lists" "$input_json")
PE_metrics=$(jq -r ".PE_metrics" "$input_json")
vcf=$(jq -r ".vcf" "$input_json")
Depth_DEL_beds=$(jq -r ".Depth_DEL_beds" "$input_json")
Depth_DUP_beds=$(jq -r ".Depth_DUP_beds" "$input_json")
min_pe_cpx=$(jq -r ".min_pe_cpx" "$input_json")
min_pe_ctx=$(jq -r ".min_pe_ctx" "$input_json")


# -------------------------------------------------------
# ======================= Command =======================
# -------------------------------------------------------

# GetSampleBatchPEMap
# ---------------------------------------------------------------------------------------------------------------------
echo "${PE_metrics}" > batch_pe_files.txt
python /opt/sv-pipeline/scripts/get_sample_batch_pe_map.py \
  --batch-name-list "${batch_name}" \
  --batch-sample-lists "${batch_sample_lists}" \
  --batch-pe-files batch_pe_files.txt \
  --prefix "${prefix}"

GetSampleBatchPEMap_sample_batch_pe_map="$(realpath "${prefix}.sample_batch_pe_map.tsv")"


# VcfToBed
# ---------------------------------------------------------------------------------------------------------------------
vcf_basename=$(basename "${vcf}" .vcf.gz)
svtk vcf2bed "${vcf}" -i ALL --include-filters "${vcf_basename}.bed"
bgzip "${vcf_basename}.bed"

VcfToBed_bed_output="$(realpath "${vcf_basename}.bed.gz")"


# SplitCpxCtx
# ---------------------------------------------------------------------------------------------------------------------
SplitCpxCtx_prefix=$(basename "${VcfToBed_bed_output}" .bed.gz)

## The commented out lines in the following are from WDL, however, running them needs disabling pipefail.
## The alternatives in the following use substitution that do not need disabling pipefail.
# zcat ~{bed} | head -1 > ~{prefix}.cpx_ctx.bed
# filterColumn=$(zcat "${VcfToBed_bed_output}" | head -1 | tr "\t" "\n" | awk '$1=="FILTER" {print NR}')
head -1 < <(zcat "${VcfToBed_bed_output}") > ${SplitCpxCtx_prefix}.cpx_ctx.bed
filterColumn=$(head -1 < <(zcat "${VcfToBed_bed_output}") | tr "\t" "\n" | awk '$1=="FILTER" {print NR}')

zcat "${VcfToBed_bed_output}" | awk 'NR > 1' | { grep CPX || true; } | awk -v filter_column="${filterColumn}" '$filter_column !~ /UNRESOLVED/' >> "${SplitCpxCtx_prefix}.cpx_ctx.bed"

zcat "${VcfToBed_bed_output}" | awk 'NR > 1' | { grep CTX || true; } >> "${SplitCpxCtx_prefix}.cpx_ctx.bed"

# INS with INV in SOURCE - will be converted to CPX later so need to evaluate evidence
zcat "${VcfToBed_bed_output}" | awk 'NR > 1' | { grep INS || true; } | { grep INV || true; } >> "${SplitCpxCtx_prefix}.cpx_ctx.bed"

bgzip "${SplitCpxCtx_prefix}.cpx_ctx.bed"

SplitCpxCtx_cpx_ctx_bed="$(realpath "${SplitCpxCtx_prefix}.cpx_ctx.bed.gz")"


# CollectLargeCNVSupportForCPX
# ---------------------------------------------------------------------------------------------------------------------

# GenerateCnvSegmentFromCpx
# ---------------------------------------------------------------------------------------------------------------------
GenerateCnvSegmentFromCpx_prefix=$(basename "${SplitCpxCtx_cpx_ctx_bed}" .bed.gz)

python /opt/sv-pipeline/scripts/generate_cnv_segment_from_cpx.py \
  --bed "${SplitCpxCtx_cpx_ctx_bed}" \
  --sample-batch-pe-map "${GetSampleBatchPEMap_sample_batch_pe_map}" \
  --prefix "${GenerateCnvSegmentFromCpx_prefix}"

bgzip "${GenerateCnvSegmentFromCpx_prefix}.lg_CNV.bed"

GenerateCnvSegmentFromCpx_cpx_lg_cnv="$(realpath "${GenerateCnvSegmentFromCpx_prefix}.lg_CNV.bed.gz")"

# ExtractCpxLgCnvByBatch
# ---------------------------------------------------------------------------------------------------------------------
zcat "${GenerateCnvSegmentFromCpx_cpx_lg_cnv}" \
  | awk -v batch="${batch_name}" '$4=="DEL" && $7==batch {print}' \
  > "${batch_name}.lg_cnv.DEL.bed"

zcat "${GenerateCnvSegmentFromCpx_cpx_lg_cnv}" \
  | awk -v batch="${batch_name}" '$4=="DUP" && $7==batch {print}' \
  > "${batch_name}.lg_cnv.DUP.bed"

bgzip "${batch_name}.lg_cnv.DEL.bed"
bgzip "${batch_name}.lg_cnv.DUP.bed"

ExtractCpxLgCnvByBatch_lg_cnv_del="$(realpath "${batch_name}.lg_cnv.DEL.bed.gz")"
ExtractCpxLgCnvByBatch_lg_cnv_dup="$(realpath "${batch_name}.lg_cnv.DUP.bed.gz")"

# SeekDepthSuppForCpx as seek_depth_supp_for_cpx_del
# ---------------------------------------------------------------------------------------------------------------------
cd "${working_dir}"
wd_seek_depth_supp_for_cpx_del=$(mktemp -d /wd_seek_depth_supp_for_cpx_del_XXXXXXXX)
wd_seek_depth_supp_for_cpx_del="$(realpath ${wd_seek_depth_supp_for_cpx_del})"
cd "${wd_seek_depth_supp_for_cpx_del}"
cpx_del_prefix="$(basename "${ExtractCpxLgCnvByBatch_lg_cnv_del}" .bed.gz)"

SeekDepthSuppForCpx "${ExtractCpxLgCnvByBatch_lg_cnv_del}" "${Depth_DEL_beds}" "${cpx_del_prefix}"

seek_depth_supp_for_cpx_del_cpx_cnv_depth_supp="$(realpath "${cpx_del_prefix}.depth_supp.bed.gz")"

cd "${working_dir}"

# SeekDepthSuppForCpx as seek_depth_supp_for_cpx_dup
# ---------------------------------------------------------------------------------------------------------------------
wd_seek_depth_supp_for_cpx_dup=$(mktemp -d /wd_seek_depth_supp_for_cpx_dup_XXXXXXXX)
wd_seek_depth_supp_for_cpx_dup="$(realpath ${wd_seek_depth_supp_for_cpx_dup})"
cd "${wd_seek_depth_supp_for_cpx_dup}"
cpx_dup_prefix="$(basename "${ExtractCpxLgCnvByBatch_lg_cnv_dup}" .bed.gz)"

SeekDepthSuppForCpx "${ExtractCpxLgCnvByBatch_lg_cnv_dup}" "${Depth_DUP_beds}" "${cpx_dup_prefix}"

seek_depth_supp_for_cpx_dup_cpx_cnv_depth_supp="$(realpath "${cpx_dup_prefix}.depth_supp.bed.gz")"

cd "${working_dir}"


# ConcatBeds as concat_beds_svtype
# ---------------------------------------------------------------------------------------------------------------------
head -n1 < <(zcat "${seek_depth_supp_for_cpx_del_cpx_cnv_depth_supp}") > header.txt

ConcatBeds_merged_bed_file="${batch_name}.lg_cnv.depth_supp.bed.gz"
{
  echo "${seek_depth_supp_for_cpx_del_cpx_cnv_depth_supp}"
  echo "${seek_depth_supp_for_cpx_dup_cpx_cnv_depth_supp}"
} > shard_bed_files.txt

while read SPLIT; do
  zcat $SPLIT
done < shard_bed_files.txt \
  | (grep -Ev "^#" || printf "") \
  | sort -Vk1,1 -k2,2n -k3,3n \
  | cat header.txt - \
  | bgzip -c \
  > "${ConcatBeds_merged_bed_file}"

tabix -p bed "${ConcatBeds_merged_bed_file}"

CollectLargeCNVSupportForCPX_lg_cnv_depth_supp="$(realpath "${ConcatBeds_merged_bed_file}")"


# GenerateCpxReviewScript
# ---------------------------------------------------------------------------------------------------------------------

GenerateCpxReviewScript_prefix="$(basename "${ExtractCpxLgCnvByBatch_lg_cnv_dup}" .bed.gz)"

cut -f1,3 "${GetSampleBatchPEMap_sample_batch_pe_map}" > sample_PE_metrics.tsv

bgzip -cd "${SplitCpxCtx_cpx_ctx_bed}" \
  | awk -F'\t' '/^#/{print; next} {$2++; print}' OFS='\t' - \
  | bgzip -c > shifted.bed.gz

python /opt/sv-pipeline/scripts/reformat_CPX_bed_and_generate_script.py \
-i shifted.bed.gz \
-s sample_PE_metrics.tsv \
-p CPX_CTX_disINS.PASS.PE_evidences \
-c collect_PE_evidences.CPX_CTX_disINS.PASS.sh \
-r "${GenerateCpxReviewScript_prefix}.svelter" \
-u "${GenerateCpxReviewScript_prefix}.unresolved_svids.txt"

GenerateCpxReviewScript_pe_evidence_collection_script="$(realpath "collect_PE_evidences.CPX_CTX_disINS.PASS.sh")"
GenerateCpxReviewScript_svelter="$(realpath "${GenerateCpxReviewScript_prefix}.svelter")"
GenerateCpxReviewScript_unresolved_svids="$(realpath "${GenerateCpxReviewScript_prefix}.unresolved_svids.txt")"


# CollectPEMetricsForCPX
# ---------------------------------------------------------------------------------------------------------------------

# CollectPEMetricsPerBatchCPX
# ---------------------------------------------------------------------------------------------------------------------

# SplitScripts
# ---------------------------------------------------------------------------------------------------------------------
# Note that in this implementation we're not splitting the files unlike the WDL implementation.
PE_collect_script="${GenerateCpxReviewScript_pe_evidence_collection_script}"
awk -v batch="${batch_name}" '$2 ~ batch {print}' "${PE_collect_script}" > collect_PE_evidences.sh

SplitScripts_script_splits="$(realpath collect_PE_evidences.sh)"



# CollectPEMetrics
# ---------------------------------------------------------------------------------------------------------------------
# Note that this step is a very long running step! (~1h on macbook)

mkdir PE_metrics/

# Note that the script generated above expects the PE file to be in the current working directory;
# hence, without modifying the above script generation method, I'm creating a symlink to the PE file in
# the current working directory. The WDL implementation re-localizes the file using gsutil, which is
# not ideal given the file size.
ln -s "${PE_metrics}" .
ln -s "${PE_metrics}.tbi" .

if [ $(wc -c < "${SplitScripts_script_splits}") -gt 0 ]; then
  bash "${SplitScripts_script_splits}"
fi

touch "${batch_name}.evidence"
touch "${batch_name}.0.PE_evidences"
for peEvFile in *.PE_evidences
do
  cat ${peEvFile} >> "${batch_name}.evidence"
done

bgzip "${batch_name}.evidence"

CollectPEMetrics_evidence="$(realpath "${batch_name}.evidence.gz")"

CollectPEMetricsPerBatchCPX_evidence="${CollectPEMetrics_evidence}"


# CalcuPEStat
# ---------------------------------------------------------------------------------------------------------------------

zcat "${CollectPEMetricsPerBatchCPX_evidence}" | cut -f3,6- | uniq -c > "${prefix}.evi_stat"
bgzip "${prefix}.evi_stat"

CalcuPEStat_evi_stat="$(realpath "${prefix}.evi_stat.gz")"

CollectPEMetricsForCPX_evidence="${CollectPEMetricsPerBatchCPX_evidence}"
CollectPEMetricsForCPX_evi_stat="${CalcuPEStat_evi_stat}"


# CalculateCpxEvidences
# ---------------------------------------------------------------------------------------------------------------------
cd "${working_dir}"
wd_CalculateCpxEvidences=$(mktemp -d /wd_CalculateCpxEvidences_XXXXXXXX)
wd_CalculateCpxEvidences="$(realpath ${wd_CalculateCpxEvidences})"
cd "${wd_CalculateCpxEvidences}"

awk '{print $6, $(NF-3)}' "${GenerateCpxReviewScript_pe_evidence_collection_script}" \
  | { grep -v "_CTX_" || true; } \
  | uniq > sample_SVID.tsv

python /opt/sv-pipeline/scripts/calculate_cpx_evidences.py \
  --sample-svid sample_SVID.tsv \
  --pe-supp "${CollectPEMetricsForCPX_evi_stat}" \
  --depth-supp "${CollectLargeCNVSupportForCPX_lg_cnv_depth_supp}" \
  --prefix "${prefix}" \
  --min-pe-cpx ${min_pe_cpx}

CalculateCpxEvidences_manual_revise_CPX_results="$(realpath "${prefix}.manual_revise.CPX_results")"


# CalculateCtxEvidences
# ---------------------------------------------------------------------------------------------------------------------
cd "${working_dir}"
wd_CalculateCtxEvidences=$(mktemp -d /wd_CalculateCtxEvidences_XXXXXXXX)
wd_CalculateCtxEvidences="$(realpath ${wd_CalculateCtxEvidences})"
cd "${wd_CalculateCtxEvidences}"

awk '{print $6, $(NF-3)}' "${GenerateCpxReviewScript_pe_evidence_collection_script}" \
  | { grep "_CTX_" || true; } \
  | uniq > sample_SVID.tsv

zcat "${CollectPEMetricsForCPX_evi_stat}" | { grep "_CTX_" || true; } | uniq > PE_supp.tsv

python /opt/sv-pipeline/scripts/calculate_ctx_evidences.py \
  --sample-svid "sample_SVID.tsv" \
  --pe-supp "PE_supp.tsv" \
  --prefix "${prefix}" \
  --min-pe-ctx ${min_pe_ctx}

CalculateCtxEvidences_manual_revise_CTX_results="$(realpath "${prefix}.manual_revise.CTX_results")"

cd "${working_dir}"
