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
import "VaPoR.wdl" as vapor
import "AnnotateILPBFeaturesPerSamplePerBed.wdl" as annotate_il_pb_featuers_per_sample_per_bed

workflow AnnotateILPBFeaturesPerSampleBed{
    input{
        String sample
        File raw_manta
        File raw_wham
        File raw_melt
        File raw_depth
        File gtgq
        File bed

        File pacbio_seq
        File pacbio_index

        File hgsv_query
        File pacbio_query
        File bionano_query

        File ref_fasta
        File ref_fai
        File ref_dict
        File contig_list
        Int min_shard_size

        Boolean requester_pays_crams = false
        Boolean run_genomic_context_anno = false
        Boolean run_extract_algo_evi = false
        Boolean run_duphold = false
        Boolean run_extract_gt_gq = true
        Boolean run_versus_raw_vcf = true
        Boolean run_rdpesr_anno = true
        Boolean run_vapor = true

        String rdpesr_benchmark_docker
        String vapor_docker
        String duphold_docker
        String sv_base_mini_docker
        String sv_pipeline_docker

        RuntimeAttr? runtime_attr_vapor 
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
        RuntimeAttr? runtime_attr_bed_vs_hgsv
        RuntimeAttr? runtime_attr_bed_vs_pacbio
        RuntimeAttr? runtime_attr_bed_vs_bionano
        RuntimeAttr? runtime_attr_bed_vs_array
    }


        call annotate_il_pb_featuers_per_sample_per_bed.AnnotateILPBFeaturesPerSamplePerBed as AnnotateILPBFeaturesPerSamplePerBed{
            input:
                sample = sample,
                raw_manta = raw_manta,
                raw_wham = raw_wham,
                raw_melt = raw_melt,
                raw_depth = raw_depth,
                gtgq = gtgq,
                bed = bed,

                pacbio_seq = pacbio_seq,
                pacbio_index = pacbio_index,

                hgsv_query = hgsv_query,
                pacbio_query = pacbio_query,
                bionano_query = bionano_query,

                ref_fasta = ref_fasta,
                ref_fai = ref_fai,
                ref_dict = ref_dict,
                contig_list = contig_list,
                min_shard_size = min_shard_size,

                requester_pays_crams = requester_pays_crams,
                run_genomic_context_anno = run_genomic_context_anno,
                run_extract_algo_evi = run_extract_algo_evi,
                run_duphold = run_duphold,
                run_extract_gt_gq = run_extract_gt_gq,
                run_versus_raw_vcf = run_versus_raw_vcf,
                run_rdpesr_anno = run_rdpesr_anno,
                run_vapor = run_vapor,

                rdpesr_benchmark_docker = rdpesr_benchmark_docker,
                vapor_docker = vapor_docker, 
                duphold_docker = duphold_docker,
                sv_base_mini_docker = sv_base_mini_docker,
                sv_pipeline_docker = sv_pipeline_docker,

                runtime_attr_vapor = runtime_attr_vapor,
                runtime_attr_duphold = runtime_attr_duphold, 
                runtime_attr_rdpesr = runtime_attr_rdpesr,
                runtime_attr_bcf2vcf = runtime_attr_bcf2vcf,
                runtime_attr_LocalizeCram = runtime_attr_LocalizeCram,
                runtime_attr_vcf2bed = runtime_attr_vcf2bed,
                runtime_attr_SplitVcf = runtime_attr_SplitVcf,
                runtime_attr_ConcatBeds = runtime_attr_ConcatBeds,
                runtime_attr_ConcatVcfs = runtime_attr_ConcatVcfs,
                runtime_inte_anno = runtime_inte_anno,
                runtime_attr_split_vcf = runtime_attr_split_vcf,
                runtime_attr_bed_vs_hgsv = runtime_attr_bed_vs_hgsv,
                runtime_attr_bed_vs_pacbio = runtime_attr_bed_vs_pacbio,
                runtime_attr_bed_vs_bionano = runtime_attr_bed_vs_bionano,
                runtime_attr_bed_vs_array = runtime_attr_bed_vs_array
        }


    output{
        File annotated_file = AnnotateILPBFeaturesPerSamplePerBed.integrated_file
        File vapor_plot = AnnotateILPBFeaturesPerSamplePerBed.vapor_plots
    }
}

