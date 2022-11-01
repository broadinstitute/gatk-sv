##########################################################################################

## Github commit: talkowski-lab/gatk-sv-v1:<ENTER HASH HERE IN FIRECLOUD>

##########################################################################################

## Copyright Broad Institute, 2020
## 
## This WDL pipeline implements Duphold 
##
##
## LICENSING : 
## This script is released under the WDL source code license (BSD-3) (see LICENSE in 
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may 
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker 
## page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
## licensing information pertaining to the included programs.

version 1.0

import "Structs.wdl"
import "TasksBenchmark.wdl" as mini_tasks
import "AnnotateTrainingFeaturesPerSample.wdl" as annotate_training_feature_per_sample

workflow AnnotateTrainingFeatures{
    input{
        Array[String] samples
        Array[File] IL_anno_files

        File ref_fasta
        File ref_fai
        File ref_dict
        Int min_shard_size

        Array[File?] vapor_result_files

        Boolean requester_pays_crams = false
        Boolean run_genomic_context_anno = false
        Boolean run_extract_algo_evi = false
        Boolean run_duphold = false
        Boolean run_extract_gt_gq = true
        Boolean run_versus_raw_vcf = true
        Boolean run_rdpesr_anno = true
        Boolean run_vapor = true

        String vapor_docker
        String sv_base_mini_docker
        String sv_benchmark_docker
        String sv_pipeline_docker

        RuntimeAttr? runtime_attr_vapor 
        RuntimeAttr? runtime_attr_bcf2vcf
        RuntimeAttr? runtime_attr_vcf2bed
        RuntimeAttr? runtime_attr_SplitVcf
        RuntimeAttr? runtime_attr_ConcatBeds
        RuntimeAttr? runtime_attr_bed_vs_hgsv
        RuntimeAttr? runtime_attr_bed_vs_array
        RuntimeAttr? runtime_attr_bed_vs_pacbio
        RuntimeAttr? runtime_attr_bed_vs_bionano
        RuntimeAttr? runtime_attr_bed_vs_ont
        RuntimeAttr? runtime_attr_override_inte_anno
        RuntimeAttr? runtime_attr_override_format_ref 
    }

    scatter (i in range(length(samples))){
        call annotate_training_feature_per_sample.AnnotateTrainingFeaturesPerSample as AnnotateTrainingFeaturesPerSample{
            input:
                sample = samples[i],
                IL_anno = IL_anno_files[i],

                vapor_result = vapor_result_files[i],

                run_vapor = run_vapor,
                run_duphold = run_duphold,
                run_rdpesr_anno = run_rdpesr_anno ,
                run_extract_gt_gq = run_extract_gt_gq,
                run_versus_raw_vcf = run_versus_raw_vcf,
                requester_pays_crams = requester_pays_crams,
                run_extract_algo_evi = run_extract_algo_evi,
                run_genomic_context_anno = run_genomic_context_anno,

                ref_fasta = ref_fasta,
                ref_fai = ref_fai,
                ref_dict = ref_dict,
                min_shard_size = min_shard_size,

                vapor_docker = vapor_docker,
                sv_base_mini_docker = sv_base_mini_docker,
                sv_pipeline_docker = sv_pipeline_docker,
                sv_benchmark_docker = sv_benchmark_docker,

                runtime_attr_vapor = runtime_attr_vapor,
                runtime_attr_bcf2vcf = runtime_attr_bcf2vcf,
                runtime_attr_vcf2bed = runtime_attr_vcf2bed,
                runtime_attr_SplitVcf = runtime_attr_SplitVcf,
                runtime_attr_ConcatBeds = runtime_attr_ConcatBeds,
                runtime_attr_bed_vs_ont = runtime_attr_bed_vs_ont,
                runtime_attr_bed_vs_hgsv = runtime_attr_bed_vs_hgsv,
                runtime_attr_bed_vs_array = runtime_attr_bed_vs_array,
                runtime_attr_bed_vs_pacbio = runtime_attr_bed_vs_pacbio,
                runtime_attr_bed_vs_bionano = runtime_attr_bed_vs_bionano,
                runtime_attr_override_inte_anno = runtime_attr_override_inte_anno,
                runtime_attr_override_format_ref = runtime_attr_override_format_ref

        }
    }
    output{
        Array[File] annotated_files = AnnotateTrainingFeaturesPerSample.annotated_file
    }
}

