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
  output_dir=$(mktemp -d /output_resolve_complex_sv_XXXXXXXX)
else
  mkdir -p "${output_dir}"
fi
output_dir="$(realpath ${output_dir})"

if [ -z "${output_json_filename}" ]; then
  output_json_filename="${output_dir}/output.json"
else
  output_json_filename="$(realpath ${output_json_filename})"
fi

working_dir="$(mktemp -d /wd_resolve_complex_sv_XXXXXXXX)"
working_dir="$(realpath ${working_dir})"
cd "${working_dir}"

vcf=$(jq -r ".vcf" "$input_json")
prefix=$(jq -r ".prefix" "$input_json")
precluster_distance=$(jq -r ".precluster_distance" "$input_json")
precluster_overlap_frac=$(jq -r ".precluster_overlap_frac" "$input_json")
max_shard_size=$(jq -r ".max_shard_size" "$input_json")
rf_cutoff_files=($(jq -r ".rf_cutoff_files[]" "$input_json"))



# -------------------------------------------------------
# ======================= Command =======================
# -------------------------------------------------------


# ShardVcfCpx
# Shard vcf for complex resolution
# ---------------------------------------------------------------------------------------------------------------------
bcftools view -G "${vcf}" -Oz -o sites_only.vcf.gz
svtk vcfcluster <(echo "sites_only.vcf.gz") "${prefix}.vcf" \
  -d "${precluster_distance}" \
  -f "${precluster_overlap_frac}" \
  --single-end \
  -p candidate_complex_clusters \
  --svtypes DEL,DUP,INS,INV,BND \
  --ignore-svtypes \
  -o 0 \
  --preserve-header \
  --preserve-ids \
  --skip-merge
bgzip "${prefix}.vcf"

ShardVcfCpx_out="${prefix}.vcf.gz"


# ShardVidsForClustering
# ---------------------------------------------------------------------------------------------------------------------
python /shard_vids_for_clustering.py \
  --clustered_vcf "${ShardVcfCpx_out}" \
  --prefix "${prefix}" \
  --records_per_shard "${max_shard_size}"

shopt -s nullglob
ShardVidsForClustering_out=( "${prefix}.vids.shard_"*.list )
shopt -u nullglob


# ---------------------------------------------------------------------------------------------------------------------
if [ "${#ShardVidsForClustering_out[@]}" > 0 ]; then

  # GetSeCutoff
  # Get SR count cutoff from RF metrics to use in single-ender rescan procedure
  # Get SE cutoff: first quartile of PE cutoff from SR random forest across all batches
  # Defaults to 4 if first quartile < 4
  # -------------------------------------------------------------------------------------------------------------------
  rf_cutoffs_files="rf_cutoffs_files"
  printf "%s\n" "${rf_cutoff_files[@]}" > "${rf_cutoffs_files}"

  while read FILE; do
    /opt/sv-pipeline/04_variant_resolution/scripts/convert_poisson_p.py \
      $( awk -F '\t' '{if ( $5=="PE_log_pval") print $2 }' "${FILE}" | head -n1 )
  done < "${rf_cutoffs_files}" \
    | Rscript -e "cat(max(c(4,floor(quantile(as.numeric(scan('stdin',quiet=T)),probs=0.25)))),sep='\n')" \
    > median_cutoff.txt

  median_PE_cutoff=$(< median_cutoff.txt)


  # PullVcfShard
  # -------------------------------------------------------------------------------------------------------------------
  vids="${ShardVidsForClustering_out[0]}"
  bcftools view --no-version --include ID=@"${vids}" "${vcf}" -O z -o "${prefix}.shard_0.vcf.gz"
  tabix "${prefix}.shard_0.vcf.gz"
  wc -l < "${vids}" > count.txt

  PullVcfShard_out="${prefix}.shard_0.vcf.gz"
  PullVcfShard_out_index="${prefix}.shard_0.vcf.gz.tbi"
  PullVcfShard_count=$(< count.txt)


  # ResolvePrep
  # Prep files for svtk resolve using bucket streaming
  # -------------------------------------------------------------------------------------------------------------------

fi
