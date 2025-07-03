#!/bin/bash

set -Exeuo pipefail

sample_id=${1}
vcf=${2}
caller=${3}
contig_index=${4}
min_size=${5}

svtk standardize --min-size "${min_size}" --contigs "${contig_index}" "${vcf}" "${sample_id}.${caller}.std.vcf" "${caller}"
bgzip "${sample_id}.${caller}.std.vcf"
