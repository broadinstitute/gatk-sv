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
  output_dir=$(mktemp -d /output_single_sample_XXXXXXXX)
else
  mkdir -p "${output_dir}"
fi
output_dir="$(realpath ${output_dir})"

if [ -z "${output_json_filename}" ]; then
  output_json_filename="${output_dir}/output.json"
else
  output_json_filename="$(realpath ${output_json_filename})"
fi

working_dir=$(mktemp -d /wd_single_sample_XXXXXXXX)
working_dir="$(realpath ${working_dir})"
cd "${working_dir}"
echo "Single-Sample Working directory: ${working_dir}"

batch=$(jq -r ".batch" "$input_json")
sample_id=$(jq -r ".sample_id" "$input_json")
run_vcf_qc=$(jq -r ".run_vcf_qc" "$input_json")
genome_file=$(jq -r ".genome_file" "$input_json")
wgd_scoring_mask=$(jq -r ".wgd_scoring_mask" "$input_json")
reference_dict=$(jq -r ".reference_dict" "$input_json")
ref_panel_bincov_matrix=$(jq -r ".ref_panel_bincov_matrix" "$input_json")





# TEMP --------- pipelining
# these inputs you get from gather sample evidence.
# WDL: # File case_counts_file_ = select_first([case_counts_file, GatherSampleEvidence.coverage_counts])
coverage_counts="/inputs/NA12878.counts.tsv.gz"



## TODO: the above file should be indexed.
## check if the output of the gather sample evidence pipeline already indexes it or not
## if not, then use the following solution
## The only solution I found that works is the following, but not sure if it is the best.
#SKIP_LINES=$(zcat "${coverage_counts}" | grep -c -E '^@|^CONTIG\s')
#tabix -S $SKIP_LINES -s 1 -b 2 -e 3 "${coverage_counts}"



# EvidenceQC
# ---------------------------------------------------------------------------------------------------------------------
evidence_qc_inputs_json_filename="${output_dir}/evidence_qc_inputs.json"
evidence_qc_outputs_json_filename="${output_dir}/evidence_qc_outputs.json"

jq -n \
  --arg batch "${batch}" \
  --arg samples "${sample_id}" \
  --arg run_vcf_qc "${run_vcf_qc}" \
  --arg genome_file "${genome_file}" \
  --arg count_files "${coverage_counts}" \
  --arg run_ploidy false \
  --arg wgd_scoring_mask "${wgd_scoring_mask}" \
  --arg reference_dict "${reference_dict}" \
  --arg bincov_matrix "${ref_panel_bincov_matrix}" \
  '{
      batch: $batch,
      samples: [$samples],
      run_vcf_qc: $run_vcf_qc,
      genome_file: $genome_file,
      count_files: [$count_files],
      run_ploidy: $run_ploidy,
      wgd_scoring_mask: $wgd_scoring_mask,
      reference_dict: $reference_dict,
      bincov_matrix: $bincov_matrix
  }' > "${evidence_qc_inputs_json_filename}"

bash /opt/sv_shell/evidence_qc.sh "${evidence_qc_inputs_json_filename}" "${evidence_qc_outputs_json_filename}"
