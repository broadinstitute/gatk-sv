#!/bin/bash

# -------------------------------------------------------
# ==================== Input & Setup ====================
# -------------------------------------------------------

set -Eeuo pipefail

input_json=${1}
output_json_filename=${2-""}
output_dir=${3:-""}

input_json="$(realpath ${input_json})"

if [ -z "${output_dir}" ]; then
  output_dir=$(mktemp -d output_ploidy_estimation_XXXXXXXX)
else
  mkdir -p "${output_dir}"
fi
output_dir="$(realpath ${output_dir})"

if [ -z "${output_json_filename}" ]; then
  output_json_filename="${output_dir}/output.json"
else
  output_json_filename="$(realpath ${output_json_filename})"
fi

working_dir=$(mktemp -d wd_ploidy_estimation_XXXXXXXX)
working_dir="$(realpath ${working_dir})"
cd "${working_dir}"

batch=$(jq -r ".batch" "${input_json}")
bincov_matrix=$(jq -r ".bincov_matrix" "${input_json}")
bin_size=1000000


# -------------------------------------------------------
# ======================= Command =======================
# -------------------------------------------------------

# --- Build ploidy matrix
zcat "${bincov_matrix}" \
 | awk -v bin_size="${bin_size}" ' \
   function printRow() \
     {printf "%s\t%d\t%d",chr,start,stop; \
      for(i=4;i<=nf;++i) {printf "\t%d",vals[i]; vals[i]=0}; \
      print ""} \
   BEGIN {binSize=bin_size} \
   NR==1 {print substr($0,2)} \
   NR==2 {chr=$1; start=$2; stop=start+binSize; nf=NF; for(i=4;i<=nf;++i) {vals[i]=$i}} \
   NR>2  {if($1!=chr){printRow(); chr=$1; start=$2; stop=start+binSize} \
          else if($2>=stop) {printRow(); while($2>=stop) {start=stop; stop=start+binSize}} \
          for(i=4;i<=nf;++i) {vals[i]+=$i}} \
   END   {if(nf!=0)printRow()}' \
 | bgzip > "${batch}_ploidy_matrix.bed.gz"

ploidy_matrix="$(realpath "${batch}_ploidy_matrix.bed.gz")"


# --- Ploidy Score
mkdir ploidy_est
Rscript /opt/WGD/bin/estimatePloidy.R -z -O ./ploidy_est "${ploidy_matrix}"


python /opt/sv-pipeline/02_evidence_assessment/estimated_CN_denoising.py \
  --binwise-copy-number ./ploidy_est/binwise_estimated_copy_numbers.bed.gz \
  --estimated-copy-number ./ploidy_est/estimated_copy_numbers.txt.gz \
  --output-stats cn_denoising_stats.tsv \
  --output-pdf cn_denoising_plots.pdf

cp cn_denoising_stats.tsv ./ploidy_est/
cp cn_denoising_plots.pdf ./ploidy_est/

tar -zcf ./ploidy_est.tar.gz ./ploidy_est
mv ploidy_est.tar.gz "${batch}_ploidy_plots.tar.gz"
