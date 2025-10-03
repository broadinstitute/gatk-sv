#!/bin/bash

set -Eeuo pipefail


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


function CNSampleNormal() {
  local _chr=$1
  local _mode=$2
  local _r=$3

  echo "----------- Starting CN Sample Normal -------------"
  echo "chr: ${_chr}"
  echo "mode: ${_mode}"
  echo "r: ${_r}"
  echo "---------------------------------------------------"

  java "-Xmx${JVM_MAX_MEM}" -jar /opt/gatk.jar PrintSVEvidence \
    --sequence-dictionary "${ref_dict}" \
    --evidence-file "${bincov_matrix}" \
    -L "${_chr}" \
    -O "${_chr}.RD.txt"

  if [ "${_mode}" == "normal" ]; then
    mv "${_chr}.RD.txt" "${_chr}.${_mode}.RD.txt"
  else
    awk -v sex="${_mode}" '$5==sex' "${ped_file}" | cut -f2 > ids.to.include
    col=$(head -n 1 "${_chr}.RD.txt" | tr '\t' '\n'|cat -n| grep -wf ids.to.include | awk -v ORS="," '{print $1}' | sed 's/,$//g' | sed 's:\([0-9]\+\):$&:g')
    col_a="{print \$1,\$2,\$3,$col}"
    awk -f <(echo "$col_a") "${_chr}.RD.txt" | tr ' ' '\t' > "${_chr}.${_mode}.RD.txt"
  fi

  # redirect stdout and stderr to cnmops.out so that EMPTY_OUTPUT_ERROR can be detected, but use tee to also output them to
  # terminal so that errors can be debugged
  EMPTY_OUTPUT_ERROR="No CNV regions in result object. Rerun cn.mops with different parameters!"
  set +e
  echo "Starting to run cnMOPS_workflow"
  bash /opt/WGD/bin/cnMOPS_workflow.sh -S "${exclude_list}" -x "${exclude_list}" -r "${_r}" -o . -M "${_chr}.${_mode}.RD.txt" </dev/null wor2>&1 | tee cnmops.out
  echo "Finished running cnMOPS_workflow"
  RC=$?
  set -e
  if [ ! $RC -eq 0 ]; then
    if grep -q "$EMPTY_OUTPUT_ERROR" "cnmops.out"; then
      touch calls/cnMOPS.cnMOPS.gff
    else
      echo "cnMOPS_workflow.sh returned a non-zero code that was not due to an empty call file."
      exit $RC
    fi
  fi

  echo "----------- Finished CN Sample Normal -------------"
}


function CleanCNMops() {
  local _sample_list=$1
  local _cnmops_gff=$2

  echo "----------- Starting Clearn CN Mops -------------"
  echo "sample_list: ${_sample_list}"
  echo "cnmops_gff: ${_cnmops_gff}"
  echo "---------------------------------------------------"

  cut -f2 "${_sample_list}" > sample.list

  mkdir calls
  grep -v "#" "${_cnmops_gff}" > cnmops.gff1
  echo "./cnmops.gff1" > GFF.list
  /opt/WGD/bin/cleancnMOPS.sh -z -o calls/ -S "${exclude_list}" sample.list GFF.list

  zcat calls/*/*.cnMOPS.DEL.bed.gz > DELS.bed
  awk -v batch="${batch}_DEL_" 'BEGIN{OFS="\t"} {print $1,$2,$3,batch,$4,"cnmops"}' DELS.bed | cat -n |\
  awk 'BEGIN{OFS="\t"} {print $2,$3,$4,$5$1,$6,"DEL",$7}' | sort -k1,1V -k2,2n > "${batch}.DEL.bed"

  cat <(echo -e "#chr\tstart\tend\tname\tsample\tsvtype\tsources") "${batch}.DEL.bed"  > "${batch}.DEL.${prefix}.bed"

  zcat calls/*/*.cnMOPS.DUP.bed.gz > DUPS.bed
  awk -v batch="${batch}_DUP_" 'BEGIN{OFS="\t"} {print $1,$2,$3,batch,$4,"cnmops"}' DUPS.bed | cat -n |\
  awk 'BEGIN{OFS="\t"} {print $2,$3,$4,$5$1,$6,"DUP",$7}' | sort -k1,1V -k2,2n > "${batch}.DUP.bed"

  cat <(echo -e "#chr\tstart\tend\tname\tsample\tsvtype\tsources") "${batch}.DUP.bed"  > "${batch}.DUP.${prefix}.bed"

  if [[ "${stitch_and_clean_large_events}" == "true" ]]; then
    echo "stitching and cleaning large events."
    mv "${batch}.DUP.${prefix}.bed" "${batch}.DUP.${prefix}.prestitch.bed"
    mv "${batch}.DEL.${prefix}.bed" "${batch}.DEL.${prefix}.prestitch.bed"
    cat "${chrom_file}" "${allo_file}" > contig.fai
    svtk rdtest2vcf --contigs contig.fai "${batch}.DUP.${prefix}.prestitch.bed" sample.list dup.vcf.gz
    svtk rdtest2vcf --contigs contig.fai "${batch}.DEL.${prefix}.prestitch.bed" sample.list del.vcf.gz
    bash /opt/sv-pipeline/04_variant_resolution/scripts/stitch_fragmented_CNVs.sh -d dup.vcf.gz dup1.vcf.gz
    bash /opt/sv-pipeline/04_variant_resolution/scripts/stitch_fragmented_CNVs.sh -d del.vcf.gz del1.vcf.gz
    svtk vcf2bed dup1.vcf.gz dup1.bed
    cat <(echo -e "#chr\tstart\tend\tname\tsample\tsvtype\tsources") \
        <(awk -v OFS="\t" -v minsize="${min_size}" '{if($3-$2>minsize)print $1,$2,$3,$4,$6,$5,"cnmops_large"}' dup1.bed) > "${batch}.DUP.${prefix}.bed"
    svtk vcf2bed del1.vcf.gz del1.bed
    cat <(echo -e "#chr\tstart\tend\tname\tsample\tsvtype\tsources") \
        <(awk -v OFS="\t" -v minsize="${min_size}" '{if($3-$2>minsize)print $1,$2,$3,$4,$6,$5,"cnmops_large"}' del1.bed) > "${batch}.DEL.${prefix}.bed"
    echo "Finished stitching."
  fi
  bgzip -f "${batch}.DEL.${prefix}.bed"
  tabix -f "${batch}.DEL.${prefix}.bed.gz"

  bgzip -f "${batch}.DUP.${prefix}.bed"
  tabix -f "${batch}.DUP.${prefix}.bed.gz"
}


# -------------------------------------------------------
# ==================== Input & Setup ====================
# -------------------------------------------------------


input_json=${1}
output_json_filename=${2-""}
output_dir=${3:-""}

input_json="$(realpath ${input_json})"

if [ -z "${output_dir}" ]; then
  output_dir=$(mktemp -d /output_batch_evidence_merging_XXXXXXXX)
else
  mkdir -p "${output_dir}"
fi
output_dir="$(realpath ${output_dir})"

if [ -z "${output_json_filename}" ]; then
  output_json_filename="${output_dir}/output.json"
else
  output_json_filename="$(realpath ${output_json_filename})"
fi

working_dir=$(mktemp -d /wd_batch_evidence_merging_XXXXXXXX)
working_dir="$(realpath ${working_dir})"
cd "${working_dir}"
echo "cnMOPS Working directory: ${working_dir}"

batch=($(jq -r '.batch' "$input_json"))
allo_file=($(jq -r '.allo_file' "$input_json"))
chrom_file=($(jq -r '.chrom_file' "$input_json"))
exclude_list=($(jq -r '.exclude_list' "$input_json"))
ped_file=($(jq -r '.ped_file' "$input_json"))
r1=($(jq -r '.r1' "$input_json"))
r2=($(jq -r '.r2' "$input_json"))
ref_dict=($(jq -r '.ref_dict' "$input_json"))
bincov_matrix=($(jq -r '.bincov_matrix' "$input_json"))
stitch_and_clean_large_events=($(jq -r '.stitch_and_clean_large_events' "$input_json"))
min_size=($(jq -r '.min_size' "$input_json"))
prefix=($(jq -r '.prefix' "$input_json"))


# -------------------------------------------------------
# ======================= Command =======================
# -------------------------------------------------------


male_r1_gff=()
male_r2_gff=()
allos=($(awk '{print $1}' "${allo_file}"))
for allo in "${allos[@]}"; do
  # Male R2
  working_dir=$(mktemp -d /wd_cn_sample_normal_${allo}_${r2}_XXXXXXXX)
  working_dir="$(realpath ${working_dir})"
  cd "${working_dir}"
  CNSampleNormal "${allo}" "1" "${r2}"
  male_r2_gff+=("${working_dir}/calls/cnMOPS.cnMOPS.gff")

  # Male R1
  working_dir=$(mktemp -d /wd_cn_sample_normal_${allo}_${r1}_XXXXXXXX)
  working_dir="$(realpath ${working_dir})"
  cd "${working_dir}"
  CNSampleNormal "${allo}" "1" "${r1}"
  male_r1_gff+=("${working_dir}/calls/cnMOPS.cnMOPS.gff")
done

normal_r1_gff=()
normal_r2_gff=()
chroms=($(awk '{print $1}' "${chrom_file}"))
for chrom in "${chroms[@]}"; do
  # Normal R2
  working_dir=$(mktemp -d /wd_cn_sample_normal_${chrom}_${r2}_XXXXXXXX)
  working_dir="$(realpath ${working_dir})"
  cd "${working_dir}"
  CNSampleNormal "${chrom}" "normal" "${r2}"
  normal_r2_gff+=("${working_dir}/calls/cnMOPS.cnMOPS.gff")

  # Normal R1
  working_dir=$(mktemp -d /wd_cn_sample_normal_${chrom}_${r1}_XXXXXXXX)
  working_dir="$(realpath ${working_dir})"
  cd "${working_dir}"
  CNSampleNormal "${chrom}" "normal" "${r1}"
  normal_r1_gff+=("${working_dir}/calls/cnMOPS.cnMOPS.gff")
done


# Female R2
working_dir=$(mktemp -d /wd_cn_sample_normal_chrX_${r2}_XXXXXXXX)
working_dir="$(realpath ${working_dir})"
cd "${working_dir}"
CNSampleNormal "chrX" "2" "${r2}"
female_r2_gff=("${working_dir}/calls/cnMOPS.cnMOPS.gff")

# Female R1
working_dir=$(mktemp -d /wd_cn_sample_normal_chrX_${r1}_XXXXXXXX)
working_dir="$(realpath ${working_dir})"
cd "${working_dir}"
CNSampleNormal "chrX" "2" "${r1}"
female_r1_gff=("${working_dir}/calls/cnMOPS.cnMOPS.gff")


# TODO: TEMP TEST
#mapfile -t normal_r1_gff < <(find /inputs_sv_shell/ -name "NR1_*_cnMOPS.cnMOPS.gff" | sort -V)
#mapfile -t normal_r2_gff < <(find /inputs_sv_shell/ -name "NR2_*_cnMOPS.cnMOPS.gff" | sort -V)
#mapfile -t male_r1_gff < <(find /inputs_sv_shell/ -name "MR1_*_cnMOPS.cnMOPS.gff" | sort -V)
#mapfile -t male_r2_gff < <(find /inputs_sv_shell/ -name "MR2_*_cnMOPS.cnMOPS.gff" | sort -V)
#mapfile -t female_r1_gff < <(find /inputs_sv_shell/ -name "FR1_*_cnMOPS.cnMOPS.gff" | sort -V)
#mapfile -t female_r2_gff < <(find /inputs_sv_shell/ -name "FR2_*_cnMOPS.cnMOPS.gff" | sort -V)
#sample_list="/inputs_sv_shell/combined_ped_file.ped"
#exclude_list="/inputs_sv_shell/GRCh38_Nmask.bed"



# CleanCNMops
working_dir=$(mktemp -d /wd_clean_cnmops_XXXXXXXX)
working_dir="$(realpath ${working_dir})"
cd "${working_dir}"

#echo "--------------- ${working_dir}"
#echo "Bash script is using these files:"
#echo "${normal_r1_gff[@]}" "${normal_r2_gff[@]}" "${male_r1_gff[@]}" "${male_r2_gff[@]}" "${female_r1_gff[@]}" "${female_r2_gff[@]}"
#echo "********************************************"
cnmops_gff_filename="$(realpath ${working_dir}/cnmops.gff)"
cat \
  "${normal_r1_gff[@]}" \
  "${normal_r2_gff[@]}" \
  "${male_r1_gff[@]}" \
  "${male_r2_gff[@]}" \
  "${female_r1_gff[@]}" \
  "${female_r2_gff[@]}" > "${cnmops_gff_filename}"

CleanCNMops "${ped_file}" "${cnmops_gff_filename}"

echo "Finished running CleanCNMops."
#first_row_string="${Allos[0]}"

#echo "${first_row_string}"

# Split that string into a new array called 'fields'
# read -ra splits on whitespace (tabs and spaces)
#read -ra fields <<< "$first_row_string"

# Now access an element from the 'fields' array (e.g., the 2nd field is at index 1)
#echo "${fields[1]}"


# -------------------------------------------------------
# ======================= Output ========================
# -------------------------------------------------------

clean_cnmops_del="${working_dir}/${batch}.DEL.${prefix}.bed.gz"
clean_cnmops_del_index="${working_dir}/${batch}.DEL.${prefix}.bed.gz.tbi"
clean_cnmops_dup="${working_dir}/${batch}.DUP.${prefix}.bed.gz"
clean_cnmops_dup_index="${working_dir}/${batch}.DUP.${prefix}.bed.gz.tbi"

del_filename="${output_dir}/$(basename "${clean_cnmops_del}")"
del_index_filename="${output_dir}/$(basename "${clean_cnmops_del_index}")"
dup_filename="${output_dir}/$(basename "${clean_cnmops_dup}")"
dup_index_filename="${output_dir}/$(basename "${clean_cnmops_dup_index}")"

mv "${clean_cnmops_del}" "${del_filename}"
mv "${clean_cnmops_del_index}" "${del_index_filename}"
mv "${clean_cnmops_dup}" "${dup_filename}"
mv "${clean_cnmops_dup_index}" "${dup_index_filename}"

outputs_json=$(jq -n \
  --arg del "${del_filename}" \
  --arg del_index "${del_index_filename}" \
  --arg dup "${dup_filename}" \
  --arg dup_index "${dup_index_filename}" \
  '{Del: $del, Del_idx: $del_index, Dup: $dup, Dup_idx: $dup_index}' )
echo "${outputs_json}" > "${output_json_filename}"

echo "Finished cnmops successfully, output json filename: ${output_json_filename}"
