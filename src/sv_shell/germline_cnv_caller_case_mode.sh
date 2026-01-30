#!/bin/bash

set -Exeuo pipefail

function get_seeded_random() {
  SEED="$1"
  openssl enc -aes-256-ctr -pass pass:"$SEED" -nosalt </dev/zero 2>/dev/null
}

function run_gcnv_case() {
  java "-Xmx${JVM_MAX_MEM}" -jar /opt/gatk.jar GermlineCNVCaller \
    --run-mode CASE \
    --arguments_file read_count_files.args \
    --contig-ploidy-calls contig-ploidy-calls-dir \
    --model gcnv-model \
    --output "${working_dir}" \
    --output-prefix case \
    --verbosity DEBUG \
    --p-alt "${p_alt}" \
    --cnv-coherence-length "${cnv_coherence_length}" \
    --max-copy-number "${max_copy_number}" \
    --mapping-error-rate "${mapping_error_rate}" \
    --sample-psi-scale "${sample_psi_scale}" \
    --depth-correction-tau "${depth_correction_tau}" \
    --copy-number-posterior-expectation-mode "${copy_number_posterior_expectation_mode}" \
    --active-class-padding-hybrid-mode "${active_class_padding_hybrid_mode}" \
    --learning-rate "${learning_rate}" \
    --adamax-beta-1 "${adamax_beta_1}" \
    --adamax-beta-2 "${adamax_beta_2}" \
    --log-emission-samples-per-round "${log_emission_samples_per_round}" \
    --log-emission-sampling-median-rel-error "${log_emission_sampling_median_rel_error}" \
    --log-emission-sampling-rounds "${log_emission_sampling_rounds}" \
    --max-advi-iter-first-epoch "${max_advi_iter_first_epoch}" \
    --max-advi-iter-subsequent-epochs "${max_advi_iter_subsequent_epochs}" \
    --min-training-epochs "${min_training_epochs}" \
    --max-training-epochs "${max_training_epochs}" \
    --initial-temperature "${initial_temperature}" \
    --num-thermal-advi-iters "${num_thermal_advi_iters}" \
    --convergence-snr-averaging-window "${convergence_snr_averaging_window}" \
    --convergence-snr-trigger-threshold "${convergence_snr_trigger_threshold}" \
    --convergence-snr-countdown-window "${convergence_snr_countdown_window}" \
    --max-calling-iters "${max_calling_iters}" \
    --caller-update-convergence-threshold "${caller_update_convergence_threshold}" \
    --caller-internal-admixing-rate "${caller_internal_admixing_rate}" \
    --caller-external-admixing-rate "${caller_external_admixing_rate}" \
    --disable-annealing "${disable_annealing}"
}


# -------------------------------------------------------
# ==================== Input & Setup ====================
# -------------------------------------------------------

input_json=${1}
output_json_filename=${2:-""}
output_dir=${3:-""}

input_json="$(realpath ${input_json})"

if [ -z "${output_dir}" ]; then
  output_dir=$(mktemp -d ${SV_SHELL_BASE_DIR}/output_germline_cnv_caller_case_mode_XXXXXXXX)
else
  mkdir -p "${output_dir}"
fi
output_dir="$(realpath ${output_dir})"

if [ -z "${output_json_filename}" ]; then
  output_json_filename="${output_dir}/output.json"
else
  output_json_filename="$(realpath ${output_json_filename})"
fi

working_dir=$(mktemp -d ${SV_SHELL_BASE_DIR}/wd_germline_cnv_caller_case_mode_XXXXXXXX)
working_dir="$(realpath ${working_dir})"
cd "${working_dir}"


scatter_index=$(jq -r ".scatter_index" "${input_json}")
contig_ploidy_calls_tar=$(jq -r ".contig_ploidy_calls_tar" "${input_json}")
gcnv_model_tar=$(jq -r ".gcnv_model_tar" "${input_json}")
read_count_files=($(jq -r ".read_count_files[]" "${input_json}"))
p_alt=$(jq -r '.p_alt // "1e-6"' "${input_json}")
cnv_coherence_length=$(jq -r '.cnv_coherence_length // "10000.0"' "${input_json}")
max_copy_number=$(jq -r '.max_copy_number // "5"' "${input_json}")
mapping_error_rate=$(jq -r '.mapping_error_rate // "0.01"' "${input_json}")
sample_psi_scale=$(jq -r '.sample_psi_scale // "0.0001"' "${input_json}")
depth_correction_tau=$(jq -r '.depth_correction_tau // "10000.0"' "${input_json}")
copy_number_posterior_expectation_mode=$(jq -r '.copy_number_posterior_expectation_mode // "HYBRID"' "${input_json}")
active_class_padding_hybrid_mode=$(jq -r '.active_class_padding_hybrid_mode // "50000"' "${input_json}")
learning_rate=$(jq -r '.learning_rate // "0.05"' "${input_json}")
adamax_beta_1=$(jq -r '.adamax_beta_1 // "0.9"' "${input_json}")
adamax_beta_2=$(jq -r '.adamax_beta_2 // "0.99"' "${input_json}")
log_emission_samples_per_round=$(jq -r '.log_emission_samples_per_round // "50"' "${input_json}")
log_emission_sampling_median_rel_error=$(jq -r '.log_emission_sampling_median_rel_error // "0.005"' "${input_json}")
log_emission_sampling_rounds=$(jq -r '.log_emission_sampling_rounds // "10"' "${input_json}")
max_advi_iter_first_epoch=$(jq -r '.max_advi_iter_first_epoch // "5000"' "${input_json}")
max_advi_iter_subsequent_epochs=$(jq -r '.max_advi_iter_subsequent_epochs // "100"' "${input_json}")
min_training_epochs=$(jq -r '.min_training_epochs // "10"' "${input_json}")
max_training_epochs=$(jq -r '.max_training_epochs // "100"' "${input_json}")
initial_temperature=$(jq -r '.initial_temperature // "2.0"' "${input_json}")
num_thermal_advi_iters=$(jq -r '.num_thermal_advi_iters // "2500"' "${input_json}")
convergence_snr_averaging_window=$(jq -r '.convergence_snr_averaging_window // "500"' "${input_json}")
convergence_snr_trigger_threshold=$(jq -r '.convergence_snr_trigger_threshold // "0.1"' "${input_json}")
convergence_snr_countdown_window=$(jq -r '.convergence_snr_countdown_window // "10"' "${input_json}")
max_calling_iters=$(jq -r '.max_calling_iters // "10"' "${input_json}")
caller_update_convergence_threshold=$(jq -r '.caller_update_convergence_threshold // "0.001"' "${input_json}")
caller_internal_admixing_rate=$(jq -r '.caller_internal_admixing_rate // "0.75"' "${input_json}")
caller_external_admixing_rate=$(jq -r '.caller_external_admixing_rate // "1.00"' "${input_json}")
disable_annealing=$(jq -r '.disable_annealing // "false"' "${input_json}")


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


default_cpu_core_count="$(nproc)"
cpu_core_count=$(jq -r ".cpu // ${default_cpu_core_count}" "${input_json}")


# -------------------------------------------------------
# ======================= Command =======================
# -------------------------------------------------------


source /opt/gatk_miniconda3/etc/profile.d/conda.sh
set +u
conda activate gatk
set -u


export MKL_NUM_THREADS="${cpu_core_count}"
export OMP_NUM_THREADS="${cpu_core_count}"

mkdir contig-ploidy-calls-dir
tar xzf "${contig_ploidy_calls_tar}" -C contig-ploidy-calls-dir

mkdir gcnv-model
tar xzf "${gcnv_model_tar}" -C gcnv-model

read_count_files_list=$(mktemp)
printf "%s\n" "${read_count_files[@]}" > "${read_count_files_list}"

> read_count_files.args
for file_path in "${read_count_files[@]}"; do
  local_name=$(basename "$file_path" .gz)
  gunzip -c "$file_path" > "$local_name"
  echo "--input $local_name" >> read_count_files.args
done

{
  # Try to run gcnv case mode. Rarely, a bad initial seed results in NaN errors...
  run_gcnv_case
} || {
  # shuffle input arguments in a deterministic manner, resulting in a new seed
  shuf --random-source=<(get_seeded_random 42) --output=read_count_files.args read_count_files.args
  # run gcnv case mode one more time
  run_gcnv_case
}

tar c -C "${working_dir}/case-tracking" . | gzip -1 > "case-gcnv-tracking-${scatter_index}.tar.gz"

# tar output calls, ensuring output files are numbered in the same order as original sample list
num_samples=${#read_count_files[@]}
NUM_SAMPLES="${num_samples}"
NUM_DIGITS=${#NUM_SAMPLES}
CURRENT_SAMPLE=0
sed 's/\.gz$//' "$read_count_files_list" \
  | while read READ_COUNT_FILE; do
    SAMPLE_NAME=$(zgrep "^@RG" "$READ_COUNT_FILE" | cut -d: -f3)
    SAMPLE_PATH=$(dirname $(grep -lR -m1 "^$SAMPLE_NAME$" "${working_dir}/case-calls"))
    CURRENT_SAMPLE_WITH_LEADING_ZEROS=$(printf "%0${NUM_DIGITS}d" $CURRENT_SAMPLE)
    tar czf "case-gcnv-calls-shard-${scatter_index}-sample-$CURRENT_SAMPLE_WITH_LEADING_ZEROS.tar.gz" \
        -C "$SAMPLE_PATH" .
    ((++CURRENT_SAMPLE))
  done


set +u
conda deactivate
conda activate gatk-sv
set -u


# -------------------------------------------------------
# ======================= Output ========================
# -------------------------------------------------------


# the following part is (1) getting files matching the pattern, (2) move them all to the output dir,
# and (3) prepare them in a format that jq can write it as an array.
gcnv_call_tars=( ${working_dir}/case-gcnv-calls-shard-${scatter_index}-sample-*.tar.gz )
gcnv_call_tars_output_dir=()
for file_path_in_wd in "${gcnv_call_tars[@]}"; do
  if [ -f "${file_path_in_wd}" ]; then
    file_path_in_output_dir="${output_dir}/$(basename "${file_path_in_wd}")"
    mv "${file_path_in_wd}" "${file_path_in_output_dir}"
    gcnv_call_tars_output_dir+=("${file_path_in_output_dir}")
  fi
done
gcnv_call_tars_json=$(printf '%s\n' "${gcnv_call_tars_output_dir[@]}" | jq -R . | jq -s . -c)


gcnv_tracking_tar="$(realpath "case-gcnv-tracking-${scatter_index}.tar.gz")"
gcnv_tracking_tar_output_dir="${output_dir}/$(basename "${gcnv_tracking_tar}")"
mv "${gcnv_tracking_tar}" "${gcnv_tracking_tar_output_dir}"

calling_config_json="$(realpath "${working_dir}/case-calls/calling_config.json")"
calling_config_json_output_dir="${output_dir}/$(basename "${calling_config_json}")"
mv "${calling_config_json}" "${calling_config_json_output_dir}"

denoising_config_json="$(realpath "${working_dir}/case-calls/denoising_config.json")"
denoising_config_json_output_dir="${output_dir}/$(basename "${denoising_config_json}")"
mv "${denoising_config_json}" "${denoising_config_json_output_dir}"

gcnvkernel_version_json="$(realpath "${working_dir}/case-calls/gcnvkernel_version.json")"
gcnvkernel_version_json_output_dir="${output_dir}/$(basename "${gcnvkernel_version_json}")"
mv "${gcnvkernel_version_json}" "${gcnvkernel_version_json_output_dir}"

sharded_interval_list="$(realpath "${working_dir}/case-calls/interval_list.tsv")"
sharded_interval_list_output_dir="${output_dir}/$(basename "${sharded_interval_list}")"
mv "${sharded_interval_list}" "${sharded_interval_list_output_dir}"


jq -n \
  --argjson gcnv_call_tars "${gcnv_call_tars_json}" \
  --arg gcnv_tracking_tar "${gcnv_tracking_tar_output_dir}" \
  --arg calling_config_json "${calling_config_json_output_dir}" \
  --arg denoising_config_json "${denoising_config_json_output_dir}" \
  --arg gcnvkernel_version_json "${gcnvkernel_version_json_output_dir}" \
  --arg sharded_interval_list "${sharded_interval_list_output_dir}" \
  '{
      "gcnv_call_tars": $gcnv_call_tars,
      "gcnv_tracking_tar": $gcnv_tracking_tar,
      "calling_config_json": $calling_config_json,
      "denoising_config_json": $denoising_config_json,
      "gcnvkernel_version_json": $gcnvkernel_version_json,
      "sharded_interval_list": $sharded_interval_list
  }' > "${output_json_filename}"

echo "Successfully finished germline cnv caller case mode, output json filename: ${output_json_filename}"
