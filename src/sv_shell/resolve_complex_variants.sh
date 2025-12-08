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
  output_dir=$(mktemp -d /output_resolve_complex_variants_XXXXXXXX)
else
  mkdir -p "${output_dir}"
fi
output_dir="$(realpath ${output_dir})"

if [ -z "${output_json_filename}" ]; then
  output_json_filename="${output_dir}/output.json"
else
  output_json_filename="$(realpath ${output_json_filename})"
fi

working_dir="$(mktemp -d /wd_resolve_complex_variants_XXXXXXXX)"
working_dir="$(realpath ${working_dir})"
cd "${working_dir}"

cohort_name=$(jq -r ".cohort_name" "$input_json")
cluster_vcfs=$(jq -r ".cluster_vcfs" "$input_json")


# -------------------------------------------------------
# ======================= Command =======================
# -------------------------------------------------------

# SubsetInversions
# Subset inversions from PESR+RD VCF
# ---------------------------------------------------------------------------------------------------------------------

SubsetInversions_filtered_vcf="${cohort_name}.inversions_only.vcf.gz"
bcftools view --no-version --no-update -i 'INFO/SVTYPE="INV"' -O z -o "${SubsetInversions_filtered_vcf}" "${cluster_vcfs}"
tabix "${SubsetInversions_filtered_vcf}"