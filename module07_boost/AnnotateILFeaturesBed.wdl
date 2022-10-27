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
import "AnnotateILFeaturesPerSampleBed.wdl" as annotate_il_features_per_sample


workflow AnnotateILFeatures{
    input{
        String prefix
        Array[String] samples

        Array[File] raw_mantas
        Array[File] raw_whams
        Array[File] raw_melts
        Array[File] raw_depths
        Array[File] gtgqs
        Array[File] beds
        #Array[File?] array_queries

        File? ref_fasta
        File? ref_fai
        File? ref_dict

        Boolean requester_pays_crams = false
        Boolean run_genomic_context_anno = false
        Boolean run_extract_algo_evi = false
        Boolean run_duphold = false
        Boolean run_rdpesr_anno = false
        Boolean run_extract_gt_gq = true
        Boolean run_versus_raw_vcf = true

        String sv_benchmark_docker
        String duphold_docker
        String sv_base_mini_docker
        String sv_pipeline_docker

        RuntimeAttr? runtime_attr_duphold
        RuntimeAttr? runtime_attr_rdpesr
        RuntimeAttr? runtime_attr_bcf2vcf
        RuntimeAttr? runtime_attr_LocalizeCram
        RuntimeAttr? runtime_attr_vcf2bed
        RuntimeAttr? runtime_attr_SplitVcf
        RuntimeAttr? runtime_attr_ConcatBeds
        RuntimeAttr? runtime_attr_ConcatVcfs
        RuntimeAttr? runtime_inte_anno
        RuntimeAttr? runtime_attr_split_vcf
        RuntimeAttr? runtime_attr_override_generate_tarball
    }

    scatter (i in range(length(samples))){
        call annotate_il_features_per_sample.AnnotateILFeaturesPerSample as AnnotateILFeaturesPerSample{
            input:
                sample = samples[i],
                raw_manta = raw_mantas[i],
                raw_wham = raw_whams[i],
                raw_melt = raw_melts[i],
                raw_depth = raw_depths[i],
                gtgq  = gtgqs[i],
                bed = beds[i],
                #array_query = array_queries[i],

                ref_fasta = ref_fasta,
                ref_fai = ref_fai, 
                ref_dict = ref_dict,

                requester_pays_crams = requester_pays_crams,
                run_genomic_context_anno = run_genomic_context_anno,
                run_extract_algo_evi = run_extract_algo_evi,
                run_duphold = run_duphold,
                run_extract_gt_gq = run_extract_gt_gq,
                run_versus_raw_vcf = run_versus_raw_vcf,
                run_rdpesr_anno = run_rdpesr_anno,

                sv_benchmark_docker = sv_benchmark_docker,
                duphold_docker = duphold_docker,
                sv_base_mini_docker = sv_base_mini_docker,
                sv_pipeline_docker = sv_pipeline_docker,

                runtime_attr_duphold = runtime_attr_duphold,
                runtime_attr_rdpesr = runtime_attr_rdpesr,
                runtime_attr_bcf2vcf = runtime_attr_bcf2vcf,
                runtime_attr_LocalizeCram = runtime_attr_LocalizeCram,
                runtime_attr_vcf2bed = runtime_attr_vcf2bed,
                runtime_attr_SplitVcf = runtime_attr_SplitVcf,
                runtime_attr_ConcatBeds = runtime_attr_ConcatBeds,
                runtime_attr_ConcatVcfs = runtime_attr_ConcatVcfs,
                runtime_inte_anno = runtime_inte_anno,
                runtime_attr_split_vcf = runtime_attr_split_vcf
        }
    }

    call mini_tasks.FilesToTarredFolder as FilesToTarredFolder{
        input:
            in_files = AnnotateILFeaturesPerSample.annotated_file,
            folder_name = "merged",
            tarball_prefix = prefix,
            sv_base_mini_docker = sv_base_mini_docker,
            runtime_attr_override = runtime_attr_override_generate_tarball
    }

    output{
        Array[File] annotation = AnnotateILFeaturesPerSample.annotated_file
        File anno_tar = FilesToTarredFolder.tarball
    }
}



