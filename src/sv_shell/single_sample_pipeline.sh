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
ref_panel_bincov_matrix=$(jq -r ".ref_panel_bincov_matrix" "$input_json")



# -------------------------------------------------------
# ======================= Command =======================
# -------------------------------------------------------

# GatherSampleEvidence
# ---------------------------------------------------------------------------------------------------------------------

gather_sample_evidence_output_dir=$(mktemp -d "/output_GatherSampleEvidence_XXXXXXXX")
gather_sample_evidence_output_dir="$(realpath ${gather_sample_evidence_output_dir})"
gather_sample_evidence_inputs_json="${gather_sample_evidence_output_dir}/inputs.json"
gather_sample_evidence_outputs_json="${gather_sample_evidence_output_dir}/outputs.json"

jq -n \
  --slurpfile inputs "${input_json}" \
  '{
    "sample_id": $inputs[0].sample_id,
    "bam_or_cram_file": $inputs[0].bam_or_cram_file,
    "bam_or_cram_index": $inputs[0].bam_or_cram_index,
    "reference_dict": $inputs[0].reference_dict,
    "reference_fasta": $inputs[0].reference_fasta,
    "reference_index": $inputs[0].reference_index,
    "primary_contigs_list": $inputs[0].primary_contigs_list,
    "primary_contigs_fai": $inputs[0].primary_contigs_fai,
    "preprocessed_intervals": $inputs[0].preprocessed_intervals,
    "manta_region_bed": $inputs[0].manta_region_bed,
    "manta_region_bed_index": $inputs[0].manta_region_bed_index,
    "sd_locs_vcf": $inputs[0].sd_locs_vcf,
    "mei_bed": $inputs[0].mei_bed,
    "wham_include_list_bed_file": $inputs[0].wham_include_list_bed_file,
    "reference_bwa_alt": $inputs[0].reference_bwa_alt,
    "reference_bwa_amb": $inputs[0].reference_bwa_amb,
    "reference_bwa_ann": $inputs[0].reference_bwa_ann,
    "reference_bwa_bwt": $inputs[0].reference_bwa_bwt,
    "reference_bwa_pac": $inputs[0].reference_bwa_pac,
    "reference_bwa_sa": $inputs[0].reference_bwa_sa,
    "collect_coverage": true,
    "run_scramble": true,
    "run_manta": true,
    "run_wham": true,
    "collect_pesr": true,
    "scramble_alignment_score_cutoff": 90,
    "run_module_metrics": $inputs[0].run_sampleevidence_metrics
  }' > "${gather_sample_evidence_inputs_json}"

bash /opt/sv_shell/gather_sample_evidence.sh \
  "${gather_sample_evidence_inputs_json}" \
  "${gather_sample_evidence_outputs_json}" \
  "${gather_sample_evidence_output_dir}"



# TODO: TEMP --------- pipelining
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
evidence_qc_output_dir=$(mktemp -d "/output_evidence_qc_XXXXXXXX")
evidence_qc_output_dir="$(realpath ${evidence_qc_output_dir})"
evidence_qc_inputs_json_filename="${evidence_qc_output_dir}/inputs.json"
evidence_qc_outputs_json_filename="${evidence_qc_output_dir}/outputs.json"

jq -n \
  --slurpfile inputs "${input_json}" \
  --arg samples "${sample_id}" \
  --arg count_files "${coverage_counts}" \
  '{
      batch: $inputs[0].batch,
      samples: [$samples],
      run_vcf_qc: $inputs[0].run_vcf_qc,
      genome_file: $inputs[0].genome_file,
      count_files: [$count_files],
      run_ploidy: false,
      wgd_scoring_mask: $inputs[0].wgd_scoring_mask,
      reference_dict: $inputs[0].reference_dict,
      bincov_matrix: $inputs[0].ref_panel_bincov_matrix
  }' > "${evidence_qc_inputs_json_filename}"

bash /opt/sv_shell/evidence_qc.sh \
  "${evidence_qc_inputs_json_filename}" \
  "${evidence_qc_outputs_json_filename}" \
  "${evidence_qc_output_dir}"
