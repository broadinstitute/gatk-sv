#!/bin/bash

# -------------------------------------------------------
# ==================== Input & Setup ====================
# -------------------------------------------------------

set -Exeuo pipefail

input_json=${1}
output_json_filename=${2-""}
output_dir=${3:-""}

input_json="$(realpath "${input_json}")"

if [ -z "${output_dir}" ]; then
  output_dir=$(mktemp -d "${SV_SHELL_BASE_DIR}/output_ploidy_estimation_XXXXXXXX")
else
  mkdir -p "${output_dir}"
fi
output_dir="$(realpath "${output_dir}")"

if [ -z "${output_json_filename}" ]; then
  output_json_filename="${output_dir}/output.json"
else
  output_json_filename="$(realpath "${output_json_filename}")"
fi

working_dir=$(mktemp -d "${SV_SHELL_BASE_DIR}/wd_ploidy_estimation_XXXXXXXX")
working_dir="$(realpath "${working_dir}")"
cd "${working_dir}"

batch=$(jq -r ".batch" "${input_json}")
bincov_matrix=$(jq -r ".bincov_matrix" "${input_json}")
reference_dict=$(jq -r ".reference_dict" "${input_json}")

poor_regions=$(jq -r '.poor_regions // empty' "${input_json}")
truth_json=$(jq -r '.truth_json // empty' "${input_json}")
preprocess_args=$(jq -r '.preprocess_args // empty' "${input_json}")
polyploidy_args=$(jq -r '.polyploidy_args // empty' "${input_json}")
infer_args=$(jq -r '.infer_args // empty' "${input_json}")
ppd_args=$(jq -r '.ppd_args // empty' "${input_json}")
call_args=$(jq -r '.call_args // empty' "${input_json}")
plot_args=$(jq -r '.plot_args // empty' "${input_json}")
enable_ppd=$(jq -r '.enable_ppd // false' "${input_json}")
use_callq20=$(jq -r '.use_callq20 // false' "${input_json}")
max_interval_size=$(jq -r '.max_interval_size // 1000000' "${input_json}")
min_interval_size=$(jq -r '.min_interval_size // 1000000' "${input_json}")
java_mem_fraction=$(jq -r '.java_mem_fraction // ""' "${input_json}")

mapfile -t subset_sd_files < <(jq -r '.subset_sd_files // [] | .[]' "${input_json}")


function getJavaMem() {
  # get JVM memory in MiB by getting total memory from /proc/meminfo
  # and multiplying by java_mem_fraction
  local mem_fraction=${java_mem_fraction:=0.85}
  awk -v MEM_FIELD="$1" -v frac="${mem_fraction}" '{
    f[substr($1, 1, length($1)-1)] = $2
  } END {
    printf "%dM", f[MEM_FIELD] * frac / 1024
  }' /proc/meminfo
}
JVM_MAX_MEM=$(getJavaMem MemTotal)
echo "JVM memory: ${JVM_MAX_MEM}"


# -------------------------------------------------------
# ======================= Command =======================
# -------------------------------------------------------

# --- Condense depth matrix (parity with PloidyEstimation.wdl CondenseDepthMatrix)
condensed_prefix="${batch}_condensed_depth"
ploidy_matrix="$(realpath "${condensed_prefix}.rd.txt.gz")"
ploidy_matrix_index="${ploidy_matrix}.tbi"

java "-Xmx${JVM_MAX_MEM}" -jar /opt/gatk.jar CondenseDepthEvidence \
  -F "${bincov_matrix}" \
  -O "${ploidy_matrix}" \
  --sequence-dictionary "${reference_dict}" \
  --max-interval-size "${max_interval_size}" \
  --min-interval-size "${min_interval_size}"


# --- Build site-depth list if subset SD files were provided
sd_args=()
if (( ${#subset_sd_files[@]} > 0 )); then
  sd_list="$(realpath subset_sd_files.list)"
  : > "${sd_list}"
  for f in "${subset_sd_files[@]}"; do
    echo "${f}" >> "${sd_list}"
  done
  sd_args=(--site-depth-list "${sd_list}")
fi


# --- Run packaged gatk-sv-ploidy pipeline
OUTDIR="${working_dir}/${batch}_ploidy"
mkdir -p "${OUTDIR}"

extra_args=()
[[ -n "${poor_regions}" ]]     && extra_args+=(--poor-regions   "${poor_regions}")
[[ -n "${truth_json}" ]]       && extra_args+=(--truth-json     "${truth_json}")
[[ -n "${preprocess_args}" ]]  && extra_args+=(--preprocess-args "${preprocess_args}")
[[ -n "${polyploidy_args}" ]]  && extra_args+=(--polyploidy-args "${polyploidy_args}")
[[ -n "${infer_args}" ]]       && extra_args+=(--infer-args     "${infer_args}")
[[ -n "${ppd_args}" ]]         && extra_args+=(--ppd-args       "${ppd_args}")
[[ -n "${call_args}" ]]        && extra_args+=(--call-args      "${call_args}")
[[ -n "${plot_args}" ]]        && extra_args+=(--plot-args      "${plot_args}")
[[ "${enable_ppd}"  == "true" ]] && extra_args+=(--ppd)
[[ "${use_callq20}" == "true" ]] && extra_args+=(--use-callq20)

/opt/gatk-sv-ploidy/run_ploidy.sh \
  --input-depth "${ploidy_matrix}" \
  --work-dir "${OUTDIR}" \
  "${sd_args[@]}" \
  "${extra_args[@]}"


# --- Package all artifacts into a single tarball (parity with PloidyScore.ploidy_plots)
ploidy_plots="$(realpath "${batch}_ploidy_plots.tar.gz")"
tar -zcf "${ploidy_plots}" -C "${working_dir}" "${batch}_ploidy"


# -------------------------------------------------------
# ======================= Output ========================
# -------------------------------------------------------

ploidy_matrix_out="${output_dir}/$(basename "${ploidy_matrix}")"
ploidy_matrix_index_out="${output_dir}/$(basename "${ploidy_matrix_index}")"
mv "${ploidy_matrix}"       "${ploidy_matrix_out}"
mv "${ploidy_matrix_index}" "${ploidy_matrix_index_out}"

ploidy_plots_out="${output_dir}/$(basename "${ploidy_plots}")"
mv "${ploidy_plots}" "${ploidy_plots_out}"

chromosome_stats_src="${OUTDIR}/infer/chromosome_stats.tsv"
bin_stats_src="${OUTDIR}/infer/bin_stats.tsv.gz"
sex_assignments_src="${OUTDIR}/call/sex_assignments.txt.gz"
aneuploidy_src="${OUTDIR}/call/aneuploidy_type_predictions.tsv"

chromosome_stats_out="${output_dir}/$(basename "${chromosome_stats_src}")"
bin_stats_out="${output_dir}/$(basename "${bin_stats_src}")"
sample_sex_assignments_out="${output_dir}/$(basename "${sex_assignments_src}")"
aneuploidy_out="${output_dir}/$(basename "${aneuploidy_src}")"

cp "${chromosome_stats_src}" "${chromosome_stats_out}"
cp "${bin_stats_src}"        "${bin_stats_out}"
cp "${sex_assignments_src}"  "${sample_sex_assignments_out}"
cp "${aneuploidy_src}"       "${aneuploidy_out}"

jq -n \
  --arg ploidy_matrix "${ploidy_matrix_out}" \
  --arg ploidy_matrix_index "${ploidy_matrix_index_out}" \
  --arg ploidy_plots "${ploidy_plots_out}" \
  --arg chromosome_stats "${chromosome_stats_out}" \
  --arg bin_stats "${bin_stats_out}" \
  --arg sample_sex_assignments "${sample_sex_assignments_out}" \
  --arg aneuploidy_type_predictions "${aneuploidy_out}" \
  '{
     "ploidy_matrix": $ploidy_matrix,
     "ploidy_matrix_index": $ploidy_matrix_index,
     "ploidy_plots": $ploidy_plots,
     "chromosome_stats": $chromosome_stats,
     "bin_stats": $bin_stats,
     "sample_sex_assignments": $sample_sex_assignments,
     "aneuploidy_type_predictions": $aneuploidy_type_predictions
   }' > "${output_json_filename}"
