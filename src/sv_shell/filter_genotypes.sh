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
  output_dir=$(mktemp -d /output_filter_genotypes_XXXXXXXX)
else
  mkdir -p "${output_dir}"
fi
output_dir="$(realpath ${output_dir})"

if [ -z "${output_json_filename}" ]; then
  output_json_filename="${output_dir}/output.json"
else
  output_json_filename="$(realpath ${output_json_filename})"
fi

working_dir=$(mktemp -d /wd_filter_genotypes_XXXXXXXX)
working_dir="$(realpath ${working_dir})"
cd "${working_dir}"
echo "Filter Genotypes working directory: ${working_dir}"

output_prefix=$(jq -r ".output_prefix" "$input_json")
vcf=$(jq -r ".vcf" "$input_json")
ploidy_table=$(jq -r ".ploidy_table" "$input_json")
no_call_rate_cutoff=$(jq -r ".no_call_rate_cutoff" "$input_json")
sl_cutoff_table=$(jq -r ".sl_cutoff_table" "$input_json")
sl_filter_args=$(jq -r ".sl_filter_args" "$input_json")
header_drop_fields=$(jq -r '.header_drop_fields // "FILTER/LOW_QUALITY,FORMAT/TRUTH_CN_EQUAL,FORMAT/GT_FILTER,FORMAT/CONC_ST,INFO/STATUS,INFO/TRUTH_AC,INFO/TRUTH_AN,INFO/TRUTH_AF,INFO/TRUTH_VID,INFO/CNV_CONCORDANCE,INFO/GENOTYPE_CONCORDANCE,INFO/HET_PPV,INFO/HET_SENSITIVITY,INFO/HOMVAR_PPV,INFO/HOMVAR_SENSITIVITY,INFO/MINSL,INFO/NON_REF_GENOTYPE_CONCORDANCE,INFO/SL_MAX,INFO/SL_MEAN,INFO/VAR_PPV,INFO/VAR_SENSITIVITY,INFO/VAR_SPECIFICITY"' "$input_json")


# -------------------------------------------------------
# ======================= Command =======================
# -------------------------------------------------------

# FilterVcf
# ---------------------------------------------------------------------------------------------------------------------
FilterVcf_output_prefix="${output_prefix}.filter_genotypes"
python /opt/sv-pipeline/scripts/apply_sl_filter.py \
  --vcf "${vcf}" \
  --out tmp.vcf.gz \
  --ploidy-table "${ploidy_table}" \
  --ncr-threshold "${no_call_rate_cutoff}" \
  --sl-cutoff-table "${sl_cutoff_table}" \
  ${sl_filter_args}

FilterVcf_out="$(realpath "${FilterVcf_output_prefix}.vcf.gz")"
bcftools +fill-tags tmp.vcf.gz -- -t AC,AN,AF \
  | bcftools view --no-update -i 'SVTYPE=="CNV" || AC>0' -Oz -o "${FilterVcf_out}"
tabix "${FilterVcf_out}"


# SanitizeHeader
# ---------------------------------------------------------------------------------------------------------------------

bcftools view --no-version -h "${FilterVcf_out}" > header.vcf

grep -v -e ^"##bcftools" header.vcf \
  -e ^"##source=depth" \
  -e ^"##source=cleanvcf" \
  -e ^"##ALT=<ID=UNR" \
  | sed 's/Split read genotype quality/Split-read genotype quality/g' \
  | sed 's/##ALT=<ID=BND,Description="Translocation">/##ALT=<ID=BND,Description="Breakend">/g' > newheader.vcf

SanitizeHeader_out="$(realpath "${output_prefix}.filter_genotypes.sanitized.vcf.gz")"

bcftools reheader -h newheader.vcf "${FilterVcf_out}" \
  | bcftools annotate -x "${header_drop_fields}" \
    --no-version \
    -O z \
    -o "${SanitizeHeader_out}"

tabix "${SanitizeHeader_out}"



# -------------------------------------------------------
# ======================= Output ========================
# -------------------------------------------------------

SanitizeHeader_output_dir="${output_dir}/$(basename "${SanitizeHeader_out}")"
mv "${SanitizeHeader_out}" "${SanitizeHeader_output_dir}"
mv "${SanitizeHeader_out}.tbi" "${SanitizeHeader_output_dir}.tbi"

jq -n \
  --arg filtered_vcf "${SanitizeHeader_output_dir}" \
  --arg filtered_vcf_index "${SanitizeHeader_output_dir}.tbi" \
  '{
      filtered_vcf: $filtered_vcf,
      filtered_vcf_index: $filtered_vcf_index
  }' > "${output_json_filename}"

echo "Finished FilterGenotypes successfully, output json filename: ${output_json_filename}"
