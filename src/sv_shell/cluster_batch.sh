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
  output_dir=$(mktemp -d /output_cluster_batch_XXXXXXXX)
else
  mkdir -p "${output_dir}"
fi
output_dir="$(realpath ${output_dir})"

if [ -z "${output_json_filename}" ]; then
  output_json_filename="${output_dir}/output.json"
else
  output_json_filename="$(realpath ${output_json_filename})"
fi

working_dir=$(mktemp -d /wd_cluster_batch_XXXXXXXX)
working_dir="$(realpath ${working_dir})"
cd "${working_dir}"

batch=$(jq -r ".batch" "${input_json}")
contig_list=$(jq -r ".contig_list" "${input_json}")
ped_file=$(jq -r ".ped_file" "${input_json}")
dragen_vcf_tar=$(jq -r '.dragen_vcf_tar // ""' "${input_json}")
manta_vcf_tar=$(jq -r '.manta_vcf_tar // ""' "${input_json}")
melt_vcf_tar=$(jq -r '.melt_vcf_tar // ""' "${input_json}")
wham_vcf_tar=$(jq -r '.wham_vcf_tar // ""' "${input_json}")
chr_x=$(jq -r '.chr_x // ""' "${input_json}")
chr_y=$(jq -r '.chr_y // ""' "${input_json}")
retain_female_chr_y=$(jq -r '.retain_female_chr_y' "${input_json}")


# -------------------------------------------------------
# ======================= Command =======================
# -------------------------------------------------------


# GetSampleIdsFromVcfTar
# ---------------------------------------------------------------------------------------------------------------------
if [ -n "${dragen_vcf_tar}" ]; then
  vcf_tar="${dragen_vcf_tar}"
elif [ -n "${manta_vcf_tar}" ]; then
  vcf_tar="${manta_vcf_tar}"
elif [ -n "${wham_vcf_tar}" ]; then
  vcf_tar="${wham_vcf_tar}"
else
  echo "Error: No VCF tar file provided" >&2
  exit 1
fi
mkdir vcfs
tar xzf "${vcf_tar}" -C vcfs/
ls vcfs/*.vcf.gz | xargs -n1 bcftools query -l | sort -u > "${batch}.samples.txt"



# CreatePloidyTableFromPed
# ---------------------------------------------------------------------------------------------------------------------

if [[ "${retain_female_chr_y}" == "true" ]]; then
  ploidy_table_output_file="${batch}.ploidy.FEMALE_chrY_1.tsv"
else
  ploidy_table_output_file="${batch}.ploidy.tsv"
fi
ploidy_table_output_file="$(realpath ${ploidy_table_output_file})"

# parameter expansion: use value after :+ if chr_x provided, empty string otherwise
python /opt/sv-pipeline/scripts/ploidy_table_from_ped.py \
  --ped "${ped_file}" \
  --out tmp.tsv \
  --contigs "${contig_list}" \
  ${chr_x:+--chr-x $chr_x} \
  ${chr_y:+--chr-y $chr_y}

# TODO : For now we retain female Y genotypes for metric generation
if [ "${retain_female_chr_y}" = true ]; then
    sed -e 's/\t0/\t1/g' tmp.tsv > "${ploidy_table_output_file}"
else
    mv tmp.tsv "${ploidy_table_output_file}"
fi
