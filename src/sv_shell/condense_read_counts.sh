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
  output_dir=$(mktemp -d /output_condense_read_counts_XXXXXXXX)
else
  mkdir -p "${output_dir}"
fi
output_dir="$(realpath ${output_dir})"

if [ -z "${output_json_filename}" ]; then
  output_json_filename="${output_dir}/output.json"
else
  output_json_filename="$(realpath ${output_json_filename})"
fi

working_dir=$(mktemp -d /wd_condense_read_counts_XXXXXXXX)
working_dir="$(realpath ${working_dir})"
cd "${working_dir}"
echo "Condense Read Counts Working directory: ${working_dir}"

sample=($(jq -r '.sample' "$input_json"))
counts=($(jq -r '.counts' "$input_json"))


gatk4_jar_override=$(jq -r --arg default_value "/root/gatk.jar" '.gatk4_jar_override // $default_value' "${input_json}")

min_interval_size=$(jq -r --arg default_value "101" '.min_interval_size // $default_value' "${input_json}")
max_interval_size=$(jq -r --arg default_value "2000" '.max_interval_size // $default_value' "${input_json}")


# -------------------------------------------------------
# ======================= Command =======================
# -------------------------------------------------------

export GATK_LOCAL_JAR="${gatk4_jar_override}"

zcat "${counts}" | grep '^@' | grep -v '@RG' > ref.dict
zcat "${counts}" | grep -v '^@' | sed -e 1d | \
    awk 'BEGIN{FS=OFS="\t";print "#Chr\tStart\tEnd\tNA21133"}{print $1,$2-1,$3,$4}' | bgzip > in.rd.txt.gz
tabix -0 -s1 -b2 -e3 in.rd.txt.gz
java -Xmx2g -jar /opt/gatk.jar CondenseDepthEvidence \
    -F in.rd.txt.gz \
    -O out.rd.txt.gz \
    --sequence-dictionary ref.dict \
    --max-interval-size "${max_interval_size}" \
    --min-interval-size "${min_interval_size}"

cat ref.dict <(zcat out.rd.txt.gz | \
    awk 'BEGIN{FS=OFS="\t";print "@RG\tID:GATKCopyNumber\tSM:${sample}\nCONTIG\tSTART\tEND\tCOUNT"}{if(NR>1)print $1,$2+1,$3,$4}') | \
    bgzip > condensed_counts."${sample}".tsv.gz


# -------------------------------------------------------
# ======================= Output ========================
# -------------------------------------------------------

out_file_wd="${working_dir}/condensed_counts.${sample}.tsv.gz"
out_file_output="${output_dir}/$(basename "${out_file_wd}")"
mv "${out_file_wd}" "${out_file_output}"

outputs_json=$(jq -n \
  --arg out "${out_file_output}" \
  '{out: $out}' )
echo "${outputs_json}" > "${output_json_filename}"

echo "Finished condense read counts successfully, output json filename: ${output_json_filename}"
