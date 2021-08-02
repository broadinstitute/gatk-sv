#!/usr/bin/env bash
# Author: asanchis@broadinstitute.org
export SV_DIR=/data/talkowski/an436/software/svtoolkit
mx="-Xmx64g"
classpath="${SV_DIR}/lib/SVToolkit.jar:${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar:${SV_DIR}/lib/gatk/Queue.jar"
while getopts s:o:r:d:g:a: flag
do
    case "${flag}" in
        s) sites=${OPTARG};;
        o) output=${OPTARG};;
        r) report=${OPTARG};;
        d) discovery=${OPTARG};;
        g) genome=${OPTARG};;
        a) array=${OPTARG};;
    esac
done
echo "sites: $sites";
echo "output: $output";
echo "report: $report";
echo "discovery: $discovery";
echo "genome: $genome";
echo "array: $array";
if [[ $genome = "37" ]]; then
  genome_path=/data/talkowski/an436/resources/genomes/human_g1k_v37_chrom/human_g1k_v37.chr.canonic.fasta
elif [[ $genome = "38" ]]; then
  genome_path=/data/talkowski/xuefang/data/reference/GRCh38.1KGP/GRCh38_full_analysis_set_plus_decoy_hla.fa
fi
echo "genome_path: $genome_path"
# Run discovery
java ${mx} -cp ${classpath} \
     org.broadinstitute.sv.main.SVAnnotator \
     -A IntensityRankSum \
     -R $genome_path \
     -vcf $sites \
     -O $output \
     -arrayIntensityFile $array \
     -sample $discovery \
     -irsSampleTag SAMPLES \
     -writeReport true \
     -reportFile $report
