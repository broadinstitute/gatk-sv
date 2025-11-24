#!/bin/bash

set -Exeuo pipefail

# -------------------------------------------------------
# ==================== Input & Setup ====================
# -------------------------------------------------------

input_json=${1}
output_json_filename=${2:-""}
output_dir=${3:-""}

input_json="$(realpath "${input_json}")"

if [ -z "${output_dir}" ]; then
  output_dir=$(mktemp -d /output_post_process_germline_cnv_calls_XXXXXXXX)
else
  mkdir -p "${output_dir}"
fi
output_dir="$(realpath "${output_dir}")"

if [ -z "${output_json_filename}" ]; then
  output_json_filename="${output_dir}/output.json"
else
  output_json_filename="$(realpath "${output_json_filename}")"
fi

working_dir=$(mktemp -d /wd_post_process_germline_cnv_calls_XXXXXXXX)
working_dir="$(realpath "${working_dir}")"
cd "${working_dir}"

gcnv_calls_tars=($(jq -r ".gcnv_calls_tars[]" "${input_json}"))
gcnv_model_tars=($(jq -r ".gcnv_model_tars[]" "${input_json}"))
calling_configs=($(jq -r ".calling_configs[]" "${input_json}"))
denoising_configs=($(jq -r ".denoising_configs[]" "${input_json}"))
gcnvkernel_version=($(jq -r ".gcnvkernel_version[]" "${input_json}"))
sharded_interval_lists=($(jq -r ".sharded_interval_lists[]" "${input_json}"))
allosomal_contigs=($(jq -r '(.allosomal_contigs // [])[]' "${input_json}"))

entity_id=$(jq -r ".entity_id" "${input_json}")
contig_ploidy_calls_tar=$(jq -r ".contig_ploidy_calls_tar" "${input_json}")
ref_copy_number_autosomal_contigs=$(jq -r ".ref_copy_number_autosomal_contigs" "${input_json}")
sample_index=$(jq -r ".sample_index" "${input_json}")


function getJavaMem() {
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


source /opt/gatk_miniconda3/etc/profile.d/conda.sh
set +u
conda activate gatk
set -u

touch calls_and_model_args.txt

for (( CALL_INDEX=0; CALL_INDEX<${#gcnv_calls_tars[@]}; CALL_INDEX++ )); do

  GCNV_CALLS_TAR="${gcnv_calls_tars[$CALL_INDEX]}"
  CALLING_CONFIG="${calling_configs[$CALL_INDEX]}"
  DENOISING_CONFIG="${denoising_configs[$CALL_INDEX]}"
  GCNVKERNEL_VERSION="${gcnvkernel_version[$CALL_INDEX]}"
  SHARDED_INTERVAL_LIST="${sharded_interval_lists[$CALL_INDEX]}"
  GCNV_MODEL_TAR="${gcnv_model_tars[$CALL_INDEX]}"

  CALL_DIR="CALLS_$CALL_INDEX"
  mkdir -p "${CALL_DIR}/SAMPLE_${sample_index}"

  tar xzf "$GCNV_CALLS_TAR" -C "$CALL_DIR/SAMPLE_${sample_index}"

  ln -rsf "$CALLING_CONFIG" "${CALL_DIR}"
  ln -rsf "$DENOISING_CONFIG" "${CALL_DIR}"
  ln -rsf "$GCNVKERNEL_VERSION" "${CALL_DIR}"
  ln -rsf "$SHARDED_INTERVAL_LIST" "${CALL_DIR}"

  echo "--calls-shard-path ${CALL_DIR}" >> calls_and_model_args.txt

  MODEL_DIR="MODEL_${CALL_INDEX}"
  mkdir "${MODEL_DIR}"

  tar xzf "${GCNV_MODEL_TAR}" -C "${MODEL_DIR}"

  echo "--model-shard-path ${MODEL_DIR}" >> calls_and_model_args.txt

done


mkdir -p contig-ploidy-calls
tar xzf "${contig_ploidy_calls_tar}" -C contig-ploidy-calls


genotyped_intervals_vcf_filename="$(realpath "genotyped-intervals-${entity_id}.vcf.gz")"
genotyped_segments_vcf_filename="$(realpath "genotyped-segments-${entity_id}.vcf.gz")"
denoised_copy_ratios_filename="$(realpath "denoised_copy_ratios-${entity_id}.tsv")"

allosomal_contigs_args=""
if [ ${#allosomal_contigs[@]} -gt 0 ]; then
    allosomal_contigs_args=$(printf -- "--allosomal-contig %s " "${allosomal_contigs[@]}")
fi

java "-Xmx${JVM_MAX_MEM}" -jar /opt/gatk.jar PostprocessGermlineCNVCalls \
    --arguments_file calls_and_model_args.txt \
    ${allosomal_contigs_args} \
    --autosomal-ref-copy-number "${ref_copy_number_autosomal_contigs}" \
    --contig-ploidy-calls contig-ploidy-calls \
    --sample-index "${sample_index}" \
    --output-genotyped-intervals "${genotyped_intervals_vcf_filename}" \
    --output-genotyped-segments "${genotyped_segments_vcf_filename}" \
    --output-denoised-copy-ratios "${denoised_copy_ratios_filename}"


set +u
conda deactivate
set -u

# -------------------------------------------------------
# ======================= Output ========================
# -------------------------------------------------------

genotyped_intervals_vcf_filename_output_dir="${output_dir}/$(basename "${genotyped_intervals_vcf_filename}")"
mv "${genotyped_intervals_vcf_filename}" "${genotyped_intervals_vcf_filename_output_dir}"

genotyped_segments_vcf_filename_output_dir="${output_dir}/$(basename "${genotyped_segments_vcf_filename}")"
mv "${genotyped_segments_vcf_filename}" "${genotyped_segments_vcf_filename_output_dir}"

denoised_copy_ratios_filename_output_dir="${output_dir}/$(basename "${denoised_copy_ratios_filename}")"
mv "${denoised_copy_ratios_filename}" "${denoised_copy_ratios_filename_output_dir}"

jq -n \
  --arg genotyped_intervals_vcf "${genotyped_intervals_vcf_filename_output_dir}" \
  --arg genotyped_segments_vcf "${genotyped_segments_vcf_filename_output_dir}" \
  --arg denoised_copy_ratios "${denoised_copy_ratios_filename_output_dir}" \
  '{
      "genotyped_intervals_vcf": $genotyped_intervals_vcf,
      "genotyped_segments_vcf": $genotyped_segments_vcf,
      "denoised_copy_ratios": $denoised_copy_ratios
  }' > "${output_json_filename}"

echo "Successfully finished post processing germline CNV calls, output json filename: ${output_json_filename}"
