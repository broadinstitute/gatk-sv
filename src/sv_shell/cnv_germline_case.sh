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
  output_dir=$(mktemp -d /output_cnv_germline_case_XXXXXXXX)
else
  mkdir -p "${output_dir}"
fi
output_dir="$(realpath ${output_dir})"

if [ -z "${output_json_filename}" ]; then
  output_json_filename="${output_dir}/output.json"
else
  output_json_filename="$(realpath ${output_json_filename})"
fi

working_dir=$(mktemp -d /wd_cnv_germline_case_XXXXXXXX)
working_dir="$(realpath ${working_dir})"
cd "${working_dir}"

contig_ploidy_model_tar=$(jq -r ".contig_ploidy_model_tar" "${input_json}")
counts=($(jq -r ".counts[]" "${input_json}"))
ploidy_mapping_error_rate=$(jq -r '.ploidy_mapping_error_rate // "0.01"' "${input_json}")
ploidy_sample_psi_scale=$(jq -r '.ploidy_sample_psi_scale // "0.0001"' "${input_json}")
gcnv_model_tars=($(jq -r ".gcnv_model_tars[]" "${input_json}"))
count_entity_ids=($(jq -r ".count_entity_ids[]" "${input_json}"))
allosomal_contigs=$(jq -c ".allosomal_contigs" "${input_json}")
ref_copy_number_autosomal_contigs=$(jq -r ".ref_copy_number_autosomal_contigs" "${input_json}")


cpu_core_count="$(nproc)"

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

# DetermineGermlineContigPloidyCaseMode
# ---------------------------------------------------------------------------------------------------------------------

source /opt/gatk_miniconda3/etc/profile.d/conda.sh
set +u
conda activate gatk
set -u

DetermineGermlineContigPloidyCaseMode_wd=$(mktemp -d /wd_DetermineGermlineContigPloidyCaseMode_XXXXXXXX)
DetermineGermlineContigPloidyCaseMode_wd="$(realpath ${DetermineGermlineContigPloidyCaseMode_wd})"
cd "${DetermineGermlineContigPloidyCaseMode_wd}"

export MKL_NUM_THREADS="${cpu_core_count}"
export OMP_NUM_THREADS="${cpu_core_count}"

mkdir input-contig-ploidy-model
tar xzf "${contig_ploidy_model_tar}" -C input-contig-ploidy-model

read_count_files_list="counts_list.tsv"
printf "%s\n" "${counts[@]}" > "${read_count_files_list}"


grep 'gz$' "${read_count_files_list}" | while IFS= read -r filename; do
    output_filename=$(basename "$filename" .gz)
    zcat "$filename" > "$output_filename"
done
awk -F'/' '{ sub(/\.gz$/, "", $NF); print "--input " $NF }' "${read_count_files_list}" > read_count_files.args


java "-Xmx${JVM_MAX_MEM}" -jar /opt/gatk.jar DetermineGermlineContigPloidy \
  --arguments_file read_count_files.args \
  --model input-contig-ploidy-model \
  --output out \
  --output-prefix case \
  --verbosity DEBUG \
  --mapping-error-rate "${ploidy_mapping_error_rate}" \
  --sample-psi-scale "${ploidy_sample_psi_scale}"

tar c -C out/case-calls . | gzip -1 > case-contig-ploidy-calls.tar.gz

DetermineGermlineContigPloidyCaseMode_contig_ploidy_calls_tar="$(realpath "case-contig-ploidy-calls.tar.gz")"

set +u
conda deactivate
set -u


# GermlineCNVCallerCaseMode
# ---------------------------------------------------------------------------------------------------------------------

temp_counter=0

gcnv_output_jsons=()
#calling_config_jsons=()
#denoising_config_jsons=()
#gcnvkernel_version_jsons=()
#sharded_interval_lists=()

for (( scatter_index=0; scatter_index<${#gcnv_model_tars[@]}; scatter_index++ )); do

  gcnv_case_scatter_wd=$(mktemp -d "/wd_germline_cnv_caller_case_mode_scatter_${scatter_index}_XXXXXXXX")
  gcnv_case_scatter_wd="$(realpath ${gcnv_case_scatter_wd})"
  gcnv_case_shard_inputs_json="${gcnv_case_scatter_wd}/inputs.json"
  gcnv_case_shard_outputs_json="${gcnv_case_scatter_wd}/outputs.json"

  jq -n \
    --slurpfile inputs "${input_json}" \
    --arg scatter_index "${scatter_index}" \
    --arg contig_ploidy_calls_tar "${DetermineGermlineContigPloidyCaseMode_contig_ploidy_calls_tar}" \
    --arg gcnv_model_tar "${gcnv_model_tars[scatter_index]}" \
    '{
      "scatter_index": $scatter_index,
      "read_count_files": $inputs[0].counts,
      "contig_ploidy_calls_tar": $contig_ploidy_calls_tar,
      "gcnv_model_tar": $gcnv_model_tar,
      "p_alt": $inputs[0].gcnv_p_alt,
      "cnv_coherence_length": $inputs[0].gcnv_cnv_coherence_length,
      "max_copy_number": $inputs[0].gcnv_max_copy_number,
      "mapping_error_rate": $inputs[0].gcnv_mapping_error_rate,
      "sample_psi_scale": $inputs[0].gcnv_sample_psi_scale,
      "depth_correction_tau": $inputs[0].gcnv_depth_correction_tau,
      "copy_number_posterior_expectation_mode": $inputs[0].gcnv_copy_number_posterior_expectation_mode,
      "active_class_padding_hybrid_mode": $inputs[0].gcnv_active_class_padding_hybrid_mode,
      "learning_rate": $inputs[0].gcnv_learning_rate,
      "adamax_beta_1": $inputs[0].gcnv_adamax_beta_1,
      "adamax_beta_2": $inputs[0].gcnv_adamax_beta_2,
      "log_emission_samples_per_round": $inputs[0].gcnv_log_emission_samples_per_round,
      "log_emission_sampling_median_rel_error": $inputs[0].gcnv_log_emission_sampling_median_rel_error,
      "log_emission_sampling_rounds": $inputs[0].gcnv_log_emission_sampling_rounds,
      "max_advi_iter_first_epoch": $inputs[0].gcnv_max_advi_iter_first_epoch,
      "max_advi_iter_subsequent_epochs": $inputs[0].gcnv_max_advi_iter_subsequent_epochs,
      "min_training_epochs": $inputs[0].gcnv_min_training_epochs,
      "max_training_epochs": $inputs[0].gcnv_max_training_epochs,
      "initial_temperature": $inputs[0].gcnv_initial_temperature,
      "num_thermal_advi_iters": $inputs[0].gcnv_num_thermal_advi_iters,
      "convergence_snr_averaging_window": $inputs[0].gcnv_convergence_snr_averaging_window,
      "convergence_snr_trigger_threshold": $inputs[0].gcnv_convergence_snr_trigger_threshold,
      "convergence_snr_countdown_window": $inputs[0].gcnv_convergence_snr_countdown_window,
      "max_calling_iters": $inputs[0].gcnv_max_calling_iters,
      "caller_update_convergence_threshold": $inputs[0].gcnv_caller_update_convergence_threshold,
      "caller_internal_admixing_rate": $inputs[0].gcnv_caller_internal_admixing_rate,
      "caller_external_admixing_rate": $inputs[0].gcnv_caller_external_admixing_rate,
      "disable_annealing": $inputs[0].gcnv_disable_annealing
    }' > "${gcnv_case_shard_inputs_json}"

    bash /opt/sv_shell/germline_cnv_caller_case_mode.sh \
      "${gcnv_case_shard_inputs_json}" \
      "${gcnv_case_shard_outputs_json}" \
      "${gcnv_case_scatter_wd}"

    gcnv_output_jsons+=("${gcnv_case_shard_outputs_json}")
#    calling_config_jsons+=($(jq -r ".calling_config_json" "${gcnv_case_shard_outputs_json}"))
#    denoising_config_jsons+=($(jq -r ".denoising_config_json" "${gcnv_case_shard_outputs_json}"))
#    gcnvkernel_version_jsons+=($(jq -r ".gcnvkernel_version_json" "${gcnv_case_shard_outputs_json}"))
#    sharded_interval_lists+=($(jq -r ".sharded_interval_list" "${gcnv_case_shard_outputs_json}"))

    # TODO: remove temp
    temp_counter=$((temp_counter + 1))
    if [ "$temp_counter" -ge 2 ]; then
        echo "DEBUG: Processed 2 samples. Breaking loop."
        break
    fi
done



# PostprocessGermlineCNVCalls_wd
# ---------------------------------------------------------------------------------------------------------------------

cd "${working_dir}"

# Implementation note:
# The wdl implementation uses wide scatters and uses the outputs of each iteration
# in a subsequent task; hence it creates many arrays, and relies on cromwell for managing
# these lists (each item in the list if an ABS file path). One way to implement a similar
# logic in sv-shell is to create an array for each of those outputs, and simply append to the
# array at each iteration, then to pass those arrays to the input JSON of `post_processing_germline_cnv_calls`,
# first combine items in each of those arrays as a JSON-list-formatted-string, then pass it to
# jq to create `post_processing_germline_cnv_calls` input.
# This approach becomes very convoluted and error-prune (ignoring its memory usage, and
# possibility of hitting max variable size when passing it between bash programs
# [you can check it using `getconf ARG_MAX` returns in bytes, using 2MB, which we most likely hit this limit]).
# Therefore, I'm implementing the following solution that relies on an intermediary json file
# and offloads possibly heavy array management to jq.
#
#
# note 1: `-s` (slurp) in the following that enables reading keys from an array of jsons
# note 2: regarding the transpose and map:
#         source:    [ [SampleA_Shard0, SampleB_Shard0], [SampleA_Shard1, SampleB_Shard1] ]
# transform into:    [ [SampleA_Shard0, SampleA_Shard1], [SampleB_Shard0, SampleB_Shard1] ]
#
# the transpose implementing the following wdl method:
# Array[Array[File]] call_tars_sample_by_shard = transpose(GermlineCNVCallerCaseMode.gcnv_call_tars)
aggregated_shards_json="${working_dir}/aggregated_shards.json"
jq -s '{
  calling_configs:        map(.calling_config_json),
  denoising_configs:      map(.denoising_config_json),
  gcnvkernel_versions:    map(.gcnvkernel_version_json),
  sharded_interval_lists: map(.sharded_interval_list),
  call_tars_matrix:       (map(.gcnv_call_tars) | transpose)
}' "${gcnv_output_jsons[@]}" > "${aggregated_shards_json}"


#call_tars_sample_by_shard_json=$(jq -s 'map(.gcnv_call_tars) | transpose' "${gcnv_output_jsons[@]}")
num_samples=${#counts[@]}
for (( sample_index=0; sample_index<num_samples; sample_index++ )); do

#  current_sample_tars_json=$(echo "${call_tars_sample_by_shard_json}" | jq -c ".[${sample_index}]")

  PostprocessGermlineCNVCalls_wd=$(mktemp -d "/wd_PostprocessGermlineCNVCalls_${sample_index}_XXXXXXXX")
  PostprocessGermlineCNVCalls_wd="$(realpath ${PostprocessGermlineCNVCalls_wd})"
  PostprocessGermlineCNVCalls_inputs_json="${PostprocessGermlineCNVCalls_wd}/inputs.json"
  PostprocessGermlineCNVCalls_outputs_json="${PostprocessGermlineCNVCalls_wd}/outputs.json"

  jq -n \
    --slurpfile inputs "${input_json}" \
    --slurpfile shards "${aggregated_shards_json}" \
    --arg entity_id "${count_entity_ids[sample_index]}" \
    --arg contig_ploidy_calls_tar "${DetermineGermlineContigPloidyCaseMode_contig_ploidy_calls_tar}" \
    --argjson allosomal_contigs "${allosomal_contigs}" \
    --arg ref_copy_number_autosomal_contigs "${ref_copy_number_autosomal_contigs}" \
    --argjson sample_index "${sample_index}" \
    '{
      "entity_id": $entity_id,
      "gcnv_calls_tars": $shards[0].call_tars_matrix[$sample_index],
      "gcnv_model_tars": $inputs[0].gcnv_model_tars,
      "calling_configs": $shards[0].calling_configs,
      "denoising_configs": $shards[0].denoising_configs,
      "gcnvkernel_version": $shards[0].gcnvkernel_versions,
      "sharded_interval_lists": $shards[0].sharded_interval_lists,
      "contig_ploidy_calls_tar": $contig_ploidy_calls_tar,
      "allosomal_contigs": $allosomal_contigs,
      "ref_copy_number_autosomal_contigs": $ref_copy_number_autosomal_contigs,
      "sample_index": $sample_index
    }' > "${PostprocessGermlineCNVCalls_inputs_json}"

  bash /opt/sv_shell/post_processing_germline_cnv_calls.sh \
    "${PostprocessGermlineCNVCalls_inputs_json}" \
    "${PostprocessGermlineCNVCalls_outputs_json}" \
    "${PostprocessGermlineCNVCalls_wd}"
done


# ExplodePloidyCalls
# ---------------------------------------------------------------------------------------------------------------------

cd "${working_dir}"

# Extract ploidy calls
mkdir calls
tar xzf "${DetermineGermlineContigPloidyCaseMode_contig_ploidy_calls_tar}" -C calls/

num_count_entity_ids=${#count_entity_ids[@]}
# Archive call files by sample, renaming so they will be glob'd in order
for (( i=0; i<num_count_entity_ids; i++ ))
do
  sample_id="${count_entity_ids[i]}"

  sample_no=$(printf "%06d" "$i")
  # note that the above part is implemented as the following in wdl (note the backticks), this causes
  # "unexpected EOF while looking for matching `"'" error, not sure how it works in cromwell (maybe old linux?!),
  # but the above syntax seems a working alternative.
  #  sample_no=`printf %06d $i`

  tar -czf "sample_${sample_no}.${sample_id}.contig_ploidy_calls.tar.gz" -C "calls/SAMPLE_${i}" .
done


# -------------------------------------------------------
# ======================= Output ========================
# -------------------------------------------------------

contig_ploidy_calls_tar_output_dir="${output_dir}/$(basename "${DetermineGermlineContigPloidyCaseMode_contig_ploidy_calls_tar}")"
mv "${DetermineGermlineContigPloidyCaseMode_contig_ploidy_calls_tar}" "${contig_ploidy_calls_tar_output_dir}"


sample_contig_ploidy_calls_tars