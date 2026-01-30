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
  output_dir=$(mktemp -d ${SV_SHELL_BASE_DIR}/output_resolve_complex_sv_XXXXXXXX)
else
  mkdir -p "${output_dir}"
fi
output_dir="$(realpath ${output_dir})"

if [ -z "${output_json_filename}" ]; then
  output_json_filename="${output_dir}/output.json"
else
  output_json_filename="$(realpath ${output_json_filename})"
fi

working_dir="$(mktemp -d ${SV_SHELL_BASE_DIR}/wd_resolve_complex_sv_XXXXXXXX)"
working_dir="$(realpath ${working_dir})"
cd "${working_dir}"

vcf=$(jq -r ".vcf" "$input_json")
vcf_idx="${vcf}.tbi"
prefix=$(jq -r ".prefix" "$input_json")
variant_prefix=$(jq -r ".variant_prefix" "$input_json")
precluster_distance=$(jq -r ".precluster_distance" "$input_json")
precluster_overlap_frac=$(jq -r ".precluster_overlap_frac" "$input_json")
max_shard_size=$(jq -r ".max_shard_size" "$input_json")
rf_cutoff_files=($(jq -r ".rf_cutoff_files[]" "$input_json"))
disc_files=($(jq -r ".disc_files[]" "$input_json"))
ref_dict=$(jq -r ".ref_dict" "$input_json")
mei_bed=$(jq -r ".mei_bed" "$input_json")
cytobands=$(jq -r ".cytobands" "$input_json")
pe_exclude_list=$(jq -r ".pe_exclude_list" "$input_json")


function getJavaMem() {
  # get JVM memory in MiB by getting total memory from /proc/meminfo
  # and multiplying by java_mem_fraction

  local mem_fraction=${java_mem_fraction:=0.85}
  cat /proc/meminfo | \
    awk -v MEM_FIELD="$1" -v frac="${mem_fraction}" '{
      f[substr($1, 1, length($1)-1)] = $2
    } END {
      printf "%dM", f[MEM_FIELD] * frac / 1024
    }'
}
JVM_MAX_MEM=$(getJavaMem MemTotal)
echo "JVM memory: $JVM_MAX_MEM"

if [[ ! -f "${cytobands}.tbi" ]]; then
    echo "Missing required file: '${cytobands}.tbi'"
    exit 1
fi

if [[ ! -f "${pe_exclude_list}.tbi" ]]; then
    echo "Missing required file: '${pe_exclude_list}.tbi'"
    exit 1
fi



# -------------------------------------------------------
# ======================= Command =======================
# -------------------------------------------------------

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

GetSeCutoff_median_PE_cutoff=$(< median_cutoff.txt)


# ResolvePrep
# Prep files for svtk resolve using bucket streaming
# -------------------------------------------------------------------------------------------------------------------
# Use GATK to pull down the discfile chunks within ±2kb of all
# INVERSION breakpoints, and bgzip / tabix
echo "Forming regions.bed"
svtk vcf2bed "${vcf}" input.bed --no-samples --no-header
cat input.bed \
  | (fgrep INV || printf "") \
  | awk -v OFS="\t" -v buffer=2000 \
    '{ print $1, $2-buffer, $2+buffer"\n"$1, $3-buffer, $3+buffer }' \
  | awk -v OFS="\t" '{ if ($2<1) $2=1; print $1, $2, $3 }' \
  | sort -Vk1,1 -k2,2n -k3,3n \
  | bedtools merge -i - \
  > regions.bed

disc_files_list_filename="disc_files_list"
printf "%s\n" "${disc_files[@]}" > "${disc_files_list_filename}"

echo "Localizing all discfile shards..."
DISC_FILE_NUM=0
while read GS_PATH_TO_DISC_FILE; do
  ((++DISC_FILE_NUM))
  SLICE="disc"$DISC_FILE_NUM"shard"

  if [ -s regions.bed ]; then
    java "-Xmx${JVM_MAX_MEM}" -jar /opt/gatk.jar PrintSVEvidence \
          --sequence-dictionary "${ref_dict}" \
          --evidence-file $GS_PATH_TO_DISC_FILE \
          -L regions.bed \
          -O ${SLICE}.PE.txt
  else
    touch ${SLICE}.PE.txt
  fi

  cat ${SLICE}.PE.txt \
    | awk '{ if ($1==$4 && $3==$6) print }' \
    | bgzip -c \
    > ${SLICE}.PE.txt.gz
  rm ${SLICE}.PE.txt
done < "${disc_files_list_filename}"

# Merge PE files and add one artificial pair corresponding to the chromosome of interest
# This makes it so that svtk doesn't break downstream
echo "Merging PE files"
{
  zcat disc*shard.PE.txt.gz;
} \
  | sort -Vk1,1 -k2,2n -k5,5n -k7,7 \
  | bgzip -c \
  > discfile.PE.txt.gz

rm disc*shard.PE.txt.gz

tabix -0 -s 1 -b 2 -e 2 -f discfile.PE.txt.gz

ResolvePrep_merged_discfile="$(realpath discfile.PE.txt.gz)"
ResolvePrep_merged_discfile_idx="$(realpath discfile.PE.txt.gz.tbi)"


# SvtkResolve
# Run svtk resolve on variants after all-ref exclusion
# -------------------------------------------------------------------------------------------------------------------
SvtkResolve_resolved_vcf="${prefix}.svtk_resolve.shard_0.resolved.unmerged.vcf"
SvtkResolve_unresolved_vcf="${prefix}.svtk_resolve.shard_0.unresolved.unmerged.vcf"
mkdir tmp
svtk resolve \
  "${vcf}" \
  "${SvtkResolve_resolved_vcf}" \
  -p "${variant_prefix}_shard_0_" \
  -u "${SvtkResolve_unresolved_vcf}" \
  -t ./tmp/ \
  --discfile "${ResolvePrep_merged_discfile}" \
  --mei-bed "${mei_bed}" \
  --cytobands "${cytobands}" \
  --min-rescan-pe-support ${GetSeCutoff_median_PE_cutoff} \
  -x "${pe_exclude_list}"

bgzip "${SvtkResolve_resolved_vcf}"
bgzip "${SvtkResolve_unresolved_vcf}"

SvtkResolve_rs_vcf="$(realpath "${SvtkResolve_resolved_vcf}.gz")"
SvtkResolve_un_vcf="$(realpath "${SvtkResolve_unresolved_vcf}.gz")"


# RestoreUnresolvedCnv
# -------------------------------------------------------------------------------------------------------------------
resolved_plus_cnv="${prefix}.restore_unresolved.shard_0.vcf.gz"

# get unresolved records
bcftools view --no-header "${SvtkResolve_un_vcf}" -Oz -o unresolved_records.vcf.gz
# rm "${SvtkResolve_un_vcf}"

# avoid possible obliteration of input file during later processing by writing
# to temporary file (and postCPX_cleanup.py writing final result to output name)
cp "${SvtkResolve_rs_vcf}" "${resolved_plus_cnv}.tmp.gz"

# Add unresolved CNVs to resolved VCF and wipe unresolved status
zcat unresolved_records.vcf.gz \
  | (fgrep -e "<DEL>" -e "<DUP>" -e "SVTYPE=DEL" -e "SVTYPE=DUP" -e "SVTYPE=CNV" -e "SVTYPE=MCNV" || printf "") \
  | sed -r -e 's/;EVENT=[^;]*;/;/' -e 's/;UNRESOLVED[^;]*;/;/g' \
  | sed -r -e 's/;UNRESOLVED_TYPE[^;]*;/;/g' -e 's/;UNRESOLVED_TYPE[^\t]*\t/\t/g' \
  | bgzip \
  >> "${resolved_plus_cnv}.tmp.gz"

# Add other unresolved variants & retain unresolved status (except for inversion single enders)
zcat unresolved_records.vcf.gz \
  | (fgrep -v -e "<DEL>" -e "<DUP>" -e "SVTYPE=DEL" -e "SVTYPE=DUP" -e "SVTYPE=CNV" -e "SVTYPE=MCNV" \
              -e "INVERSION_SINGLE_ENDER" || printf "") \
  | bgzip \
  >> "${resolved_plus_cnv}.tmp.gz"

#Add inversion single enders as SVTYPE=BND
zcat unresolved_records.vcf.gz \
  | (fgrep -v -e "<DEL>" -e "<DUP>" -e "SVTYPE=DEL" -e "SVTYPE=DUP" -e "SVTYPE=CNV" -e "SVTYPE=MCNV" || printf "") \
  | (fgrep -e "INVERSION_SINGLE_ENDER" || printf "") \
  | sed -e 's/SVTYPE=INV/SVTYPE=BND/g' \
  | sed -e 's/END=\([0-9]*\)/END=\1;END2=\1/' \
  | bgzip \
  >> "${resolved_plus_cnv}.tmp.gz"
rm unresolved_records.vcf.gz

#Sort, clean, and compress
zcat "${resolved_plus_cnv}.tmp.gz" \
  | vcf-sort -c \
  | /opt/sv-pipeline/04_variant_resolution/scripts/postCPX_cleanup.py \
    /dev/stdin /dev/stdout \
  | bgzip \
  > "${resolved_plus_cnv}"
tabix "${resolved_plus_cnv}"

RestoreUnresolvedCnv_res="$(realpath "${resolved_plus_cnv}")"
RestoreUnresolvedCnv_res_idx="$(realpath "${resolved_plus_cnv}.tbi")"

# ConcatResolvedPerShard
# -------------------------------------------------------------------------------------------------------------------
cp "${RestoreUnresolvedCnv_res}" "${prefix}.resolved.vcf.gz"
cp "${RestoreUnresolvedCnv_res_idx}" "${prefix}.resolved.vcf.gz.tbi"

resolved_vcf_merged="$(realpath "${prefix}.resolved.vcf.gz")"
resolved_vcf_merged_idx="$(realpath "${prefix}.resolved.vcf.gz.tbi")"


# -------------------------------------------------------
# ======================= Output ========================
# -------------------------------------------------------


resolved_vcf_merged_output_dir="${output_dir}/$(basename "${resolved_vcf_merged}")"
cp "${resolved_vcf_merged}" "${resolved_vcf_merged_output_dir}"

resolved_vcf_merged_idx_output_dir="${output_dir}/$(basename "${resolved_vcf_merged_idx}")"
cp "${resolved_vcf_merged_idx}" "${resolved_vcf_merged_idx_output_dir}"

jq -n \
  --arg resolved_vcf_merged "${resolved_vcf_merged_output_dir}" \
  --arg resolved_vcf_merged_idx "${resolved_vcf_merged_idx_output_dir}" \
  '{
      "resolved_vcf_merged": $resolved_vcf_merged,
      "resolved_vcf_merged_idx": $resolved_vcf_merged_idx
  }' > "${output_json_filename}"

echo "Successfully finished resolve complex sv, output json filename: ${output_json_filename}"
