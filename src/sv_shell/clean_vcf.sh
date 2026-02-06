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
  output_dir=$(mktemp -d ${SV_SHELL_BASE_DIR}/output_clean_vcf_XXXXXXXX)
else
  mkdir -p "${output_dir}"
fi
output_dir="$(realpath ${output_dir})"

if [ -z "${output_json_filename}" ]; then
  output_json_filename="${output_dir}/output.json"
else
  output_json_filename="$(realpath ${output_json_filename})"
fi

working_dir=$(mktemp -d ${SV_SHELL_BASE_DIR}/wd_clean_vcf_XXXXXXXX)
working_dir="$(realpath ${working_dir})"
cd "${working_dir}"
echo "Clean VCF working directory: ${working_dir}"


ped_file=$(jq -r ".ped_file" "$input_json")
cohort_name=$(jq -r ".cohort_name" "$input_json")
contig_list=$(jq -r ".contig_list" "$input_json")
complex_genotype_vcf=$(jq -r ".complex_genotype_vcf" "$input_json")
complex_resolve_bothside_pass_list=$(jq -r ".complex_resolve_bothside_pass_list" "$input_json")
complex_resolve_background_fail_list=$(jq -r ".complex_resolve_background_fail_list" "$input_json")
chr_x=$(jq -r ".chr_x" "$input_json")
chr_y=$(jq -r ".chr_y" "$input_json")
HERVK_reference=$(jq -r ".HERVK_reference" "$input_json")
LINE1_reference=$(jq -r ".LINE1_reference" "$input_json")
intron_reference=$(jq -r ".intron_reference" "$input_json")


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



# CreatePloidyTableFromPed
# ---------------------------------------------------------------------------------------------------------------------
# note that since the CleanVCF workflow calls this task by setting 'retain_female_chr_y=false',
# then the following implements only the related execution path
CreatePloidyTableFromPed_out="$(realpath "${cohort_name}.ploidy.tsv")"
python /opt/sv-pipeline/scripts/ploidy_table_from_ped.py \
  --ped "${ped_file}" \
  --out "${CreatePloidyTableFromPed_out}" \
  --contigs "${contig_list}" \
  --chr-x "${chr_x}" \
  --chr-y "${chr_y}"

# CleanVcfChromosome
# =====================================================================================================================

# FormatVcfToClean
# ---------------------------------------------------------------------------------------------------------------------
# Note: the WDL version re-runs CreatePloidyTableFromPed here,
# though seems unnecessary duplicated run, so skipping it here.

# FormatVcf
# ---------
# Convert format
FormatVcf_out="${cohort_name}.formatted.vcf.gz"
python /opt/sv-pipeline/scripts/format_svtk_vcf_for_gatk.py \
  --vcf "${complex_genotype_vcf}" \
  --out "${FormatVcf_out}" \
  --ploidy-table "${CreatePloidyTableFromPed_out}" \
  --bothside-pass-list "${complex_resolve_bothside_pass_list}" \
  --background-fail-list "${complex_resolve_background_fail_list}"

tabix "${FormatVcf_out}"
FormatVcfToClean_gatk_formatted_vcf="${FormatVcf_out}"


# CleanVcfPreprocess
# ---------------------------------------------------------------------------------------------------------------------
CleanVcfPreprocess_out="${cohort_name}.preprocess.vcf.gz"
_vcf="${FormatVcfToClean_gatk_formatted_vcf}"

bcftools view --header-only "${_vcf}" | grep '^##' > header.txt

cat <<EOF >> header.txt
##FILTER=<ID=UNRESOLVED,Description="Variant is unresolved">
##INFO=<ID=HIGH_SR_BACKGROUND,Number=0,Type=Flag,Description="Variant has high number of SR splits in background samples">
##INFO=<ID=BOTHSIDES_SUPPORT,Number=0,Type=Flag,Description="Variant has read-level support for both sides of breakpoint">
##INFO=<ID=REVISED_EVENT,Number=0,Type=Flag,Description="Variant has been revised due to a copy number mismatch">
EOF

bcftools view --header-only "${_vcf}" | grep '^#CHROM' >> header.txt

bcftools view "${_vcf}" | bcftools reheader -h header.txt | bgzip -c > processed.reheader.vcf.gz

rm header.txt

python /opt/sv-pipeline/04_variant_resolution/scripts/cleanvcf_preprocess.py \
  -V processed.reheader.vcf.gz \
  -O "${CleanVcfPreprocess_out}" \
  --chrX "${chr_x}" \
  --chrY "${chr_y}" \
  --fail-list "${complex_resolve_background_fail_list}" \
  --pass-list "${complex_resolve_bothside_pass_list}" \
  --ped-file "${ped_file}"

tabix -p vcf "${CleanVcfPreprocess_out}"


# CleanVcfReviseOverlappingCnvs
# ---------------------------------------------------------------------------------------------------------------------
CleanVcfReviseOverlappingCnvs_out="${cohort_name}.revise_overlapping_cnvs.vcf.gz"
java "-Xmx${JVM_MAX_MEM}" -jar /opt/gatk.jar SVReviseOverlappingCnvs \
  -V "${CleanVcfPreprocess_out}" \
  -O "${CleanVcfReviseOverlappingCnvs_out}"


# CleanVcfReviseMultiallelicCnvs
# ---------------------------------------------------------------------------------------------------------------------
CleanVcfReviseMultiallelicCnvs_out="${cohort_name}.revise_multiallelic_cnvs.vcf.gz"
java "-Xmx${JVM_MAX_MEM}" -jar /opt/gatk.jar SVReviseMultiallelicCnvs \
  -V "${CleanVcfReviseOverlappingCnvs_out}" \
  -O "${CleanVcfReviseMultiallelicCnvs_out}"


# CleanVcfReviseOverlappingMultiallelics
# ---------------------------------------------------------------------------------------------------------------------
CleanVcfReviseOverlappingMultiallelics_out="${cohort_name}.revise_overlapping_multiallelics.vcf.gz"
java "-Xmx${JVM_MAX_MEM}" -jar /opt/gatk.jar SVReviseOverlappingMultiallelics \
  -V "${CleanVcfReviseMultiallelicCnvs_out}" \
  -O "${CleanVcfReviseOverlappingMultiallelics_out}"


# CleanVcfPostprocess
# ---------------------------------------------------------------------------------------------------------------------
CleanVcfPostprocess_out="${cohort_name}.postprocess.vcf.gz"
if [ ! -f "${CleanVcfReviseOverlappingMultiallelics_out}.tbi" ]; then
  tabix -p vcf "${CleanVcfReviseOverlappingMultiallelics_out}"
fi

python /opt/sv-pipeline/04_variant_resolution/scripts/cleanvcf_postprocess.py \
  -V "${CleanVcfReviseOverlappingMultiallelics_out}" \
  -O processed.vcf.gz \
  --ped-file "${ped_file}"

bcftools annotate -x INFO/MULTIALLELIC,INFO/UNRESOLVED,INFO/EVENT,INFO/REVISED_EVENT,INFO/MULTI_CNV,INFO/varGQ processed.vcf.gz -o processed.annotated.vcf.gz -O z

bcftools view -h processed.annotated.vcf.gz | grep "^##" | \
  grep -v -E "CIPOS|CIEND|RMSSTD|source|bcftools|GATKCommandLine|##FORMAT=<ID=EV>|##ALT=<ID=UNR>|##INFO=<ID=(MULTIALLELIC|UNRESOLVED|EVENT|REVISED_EVENT|MULTI_CNV|varGQ)" > temp_header.txt
echo '##INFO=<ID=UNRESOLVED_TYPE,Number=1,Type=String,Description="Class of unresolved variant.">' >> temp_header.txt
echo '##ALT=<ID=CNV,Description="Copy Number Polymorphism">' >> temp_header.txt

bcftools view -h processed.annotated.vcf.gz | grep "^#CHROM" > chrom_header.txt

cat temp_header.txt chrom_header.txt > header.txt

bcftools reheader -h header.txt processed.annotated.vcf.gz -o "${CleanVcfPostprocess_out}"

tabix -p vcf "${CleanVcfPostprocess_out}"


# DropRedundantCnvs
# ---------------------------------------------------------------------------------------------------------------------
DropRedundantCnvs_out="${cohort_name}.drop_redundant_cnvs.vcf.gz"
/opt/sv-pipeline/04_variant_resolution/scripts/resolve_cpx_cnv_redundancies.py \
  "${CleanVcfPostprocess_out}" "${DropRedundantCnvs_out}" --temp-dir ./tmp


# SortDropRedundantCnvs
# ---------------------------------------------------------------------------------------------------------------------
SortDropRedundantCnvs_out="${cohort_name}.drop_redundant_cnvs.sorted.vcf.gz"
mkdir temp
bcftools sort \
    --temp-dir temp \
    --output-type z \
    --output-file "${SortDropRedundantCnvs_out}" \
    "${DropRedundantCnvs_out}"
tabix "${SortDropRedundantCnvs_out}"


# StitchFragmentedCnvs
# ---------------------------------------------------------------------------------------------------------------------
StitchFragmentedCnvs_out="${cohort_name}.stitch_fragmented_cnvs.vcf.gz"

echo "First pass..."
_vcf="stich_fragmented_cnv_tmp.tsv.gz"
cp "${SortDropRedundantCnvs_out}" "${_vcf}"
java "-Xmx${JVM_MAX_MEM}" -jar ${STITCH_JAR} 0.2 200000 0.2 "${_vcf}" \
  | bgzip \
  > tmp.vcf.gz
rm "${_vcf}"  # not sure why delete, but keeping it same as the wdl version
echo "Second pass..."
java "-Xmx${JVM_MAX_MEM}" -jar ${STITCH_JAR} 0.2 200000 0.2 tmp.vcf.gz \
  | bgzip \
  > "${StitchFragmentedCnvs_out}"


# RescueMobileElementDeletions
# ---------------------------------------------------------------------------------------------------------------------
_prefix="${cohort_name}.rescue_me_dels"
RescueMobileElementDeletions_out="${_prefix}.vcf.gz"

python <<CODE
import os
import pysam
fin=pysam.VariantFile("${StitchFragmentedCnvs_out}")
fo=pysam.VariantFile("${_prefix}.bnd_del.vcf.gz", 'w', header = fin.header)
for record in fin:
    if record.info['SVTYPE'] in ['BND'] and record.info['STRANDS']=="+-" and record.chrom == record.info['CHR2'] and record.info['END2'] - record.start < 10000:
        record.info['SVLEN'] = record.info['END2'] - record.start
        fo.write(record)
fin.close()
fo.close()
CODE

tabix -p vcf "${_prefix}.bnd_del.vcf.gz"

svtk vcf2bed "${_prefix}.bnd_del.vcf.gz" -i ALL --include-filters "${_prefix}.bnd_del.bed"
bgzip "${_prefix}.bnd_del.bed"

bedtools coverage -wo -a "${_prefix}.bnd_del.bed.gz" -b "${LINE1_reference}" | awk '{if ($NF>.5) print}' | cut -f4 | sed -e 's/$/\tDEL\tPASS\toverlap_LINE1/' > manual_revise.MEI_DEL_from_BND.SVID_SVTYPE_FILTER_INFO.tsv
bedtools coverage -wo -a "${_prefix}.bnd_del.bed.gz" -b "${HERVK_reference}" | awk '{if ($NF>.5) print}' | cut -f4 | sed -e 's/$/\tDEL\tPASS\toverlap_HERVK/' >> manual_revise.MEI_DEL_from_BND.SVID_SVTYPE_FILTER_INFO.tsv

python <<CODE
import pysam
def SVID_MEI_DEL_readin(MEI_DEL_reset):
    out={}
    fin=open(MEI_DEL_reset)
    for line in fin:
        pin=line.strip().split()
        if not pin[0] in out.keys():
            out[pin[0]] = pin[3]
    fin.close()
    return out

hash_MEI_DEL_reset = SVID_MEI_DEL_readin("manual_revise.MEI_DEL_from_BND.SVID_SVTYPE_FILTER_INFO.tsv")
fin=pysam.VariantFile("${StitchFragmentedCnvs_out}")
fo=pysam.VariantFile("${_prefix}.vcf.gz", 'w', header = fin.header)
for record in fin:
    if record.id in hash_MEI_DEL_reset.keys():
        del record.filter['UNRESOLVED']
        record.info['SVTYPE'] = 'DEL'
        record.info['SVLEN'] = record.info['END2'] - record.start
        record.stop = record.info['END2']
        record.info.pop("CHR2")
        record.info.pop("END2")
        record.info.pop("UNRESOLVED_TYPE")
        if hash_MEI_DEL_reset[record.id] == 'overlap_LINE1':
            record.alts = ('<DEL:ME:LINE1>',)
        if hash_MEI_DEL_reset[record.id] == 'overlap_HERVK':
            record.alts = ('<DEL:ME:HERVK>',)
    fo.write(record)
fin.close()
fo.close()
CODE


# AddHighFDRFilters
# ---------------------------------------------------------------------------------------------------------------------
_prefix="${cohort_name}.high_fdr_filtered"
_vcf="${RescueMobileElementDeletions_out}"
AddHighFDRFilters_out="${_prefix}.vcf.gz"

python <<CODE
import pysam
with pysam.VariantFile("${_vcf}", 'r') as fin:
  header = fin.header
  header.add_line("##FILTER=<ID=HIGH_ALGORITHM_FDR,Description=\"Categories of variants with low precision including Wham-only deletions and certain Scramble SVAs\">")
  with pysam.VariantFile("${_prefix}.vcf.gz", 'w', header=header) as fo:
    for record in fin:
        if (record.info['ALGORITHMS'] == ('wham',) and record.info['SVTYPE'] == 'DEL') or \
          (record.info['ALGORITHMS'] == ('scramble',) and record.info['HIGH_SR_BACKGROUND'] and record.alts == ('<INS:ME:SVA>',)):
            record.filter.add('HIGH_ALGORITHM_FDR')
        fo.write(record)
CODE


# AddRetroDelFilters
# ---------------------------------------------------------------------------------------------------------------------
AddRetroDelFilters_out="${cohort_name}.retro_del_filtered.vcf.gz"

python /opt/sv-pipeline/04_variant_resolution/scripts/add_retro_del_filters.py \
  "${AddHighFDRFilters_out}" \
  "${intron_reference}" \
  "${AddRetroDelFilters_out}"


# FinalCleanup
# ---------------------------------------------------------------------------------------------------------------------
_prefix="${cohort_name}.final_cleanup"
FinalCleanup_final_cleaned_shard="${_prefix}.vcf.gz"

/opt/sv-pipeline/04_variant_resolution/scripts/rename_after_vcfcluster.py \
  --prefix "${_prefix}" \
  "${AddRetroDelFilters_out}" stdout \
  | bcftools annotate --no-version -e 'SVTYPE=="CNV" && SVLEN<5000' -x INFO/MEMBERS -Oz -o "${FinalCleanup_final_cleaned_shard}"
tabix "${FinalCleanup_final_cleaned_shard}"


# FormatVcfToOutput
# ---------------------------------------------------------------------------------------------------------------------
# FormatVcf
# ---------
FormatVcfToOutput_out="${cohort_name}.final.vcf.gz"
python /opt/sv-pipeline/scripts/format_svtk_vcf_for_gatk.py \
  --vcf "${FinalCleanup_final_cleaned_shard}" \
  --out "${FormatVcfToOutput_out}" \
  --ploidy-table "${CreatePloidyTableFromPed_out}" \
  --bothside-pass-list "${complex_resolve_bothside_pass_list}" \
  --background-fail-list "${complex_resolve_background_fail_list}"

tabix "${FormatVcfToOutput_out}"


# -------------------------------------------------------
# ======================= Output ========================
# -------------------------------------------------------

cleaned_vcf="${output_dir}/$(basename "${FormatVcfToOutput_out}")"
cleaned_vcf_index="${output_dir}/$(basename "${FormatVcfToOutput_out}.tbi")"

mv "${FormatVcfToOutput_out}" "${cleaned_vcf}"
mv "${FormatVcfToOutput_out}.tbi" "${cleaned_vcf}.tbi"


jq -n \
  --arg cleaned_vcf "${cleaned_vcf}" \
  --arg cleaned_vcf_index "${cleaned_vcf}.tbi" \
  '{
     "cleaned_vcf": $cleaned_vcf,
     "cleaned_vcf_index": $cleaned_vcf_index,
   }' > "${output_json_filename}"

echo "Finished Clean VCF, output json filename: ${output_json_filename}"
