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
import "AnnotateILFeaturesPerSample.wdl" as anno_il
import "TasksBenchmark.wdl" as mini_tasks
import "VaPoR.wdl" as vapor
import "AnnotateILPBFeaturesPerSamplePerBed.wdl" as annotate_il_pb_featuers_per_sample_per_bed

workflow AnnotateILPBFeatures{
    input{
        Array[String] samples
        #Array[String?] il_bams
        #Array[String?] il_bam_bais
        
        #Array[File?] pe_metrics
        #Array[File?] pe_indexes
        #Array[File?] sr_metrics
        #Array[File?] sr_indexes
        #Array[File?] rd_metrics
        #Array[File?] rd_indexes
        
        #Array[File?] ref_SegDups
        #Array[File?] ref_SimpReps
        #Array[File?] ref_RepMasks


        Array[File] raw_mantas
        Array[File] raw_whams
        Array[File] raw_melts
        Array[File] raw_depths
        Array[File] gtgqs
        Array[File] beds


        Array[File] pacbio_seqs
        Array[File] pacbio_indexes

        Array[File] hgsv_queries
        Array[File] pacbio_queries
        Array[File] bionano_queries

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
        RuntimeAttr? runtime_attr_ConcatBeds
        RuntimeAttr? runtime_attr_ConcatVcfs
        RuntimeAttr? runtime_inte_anno
        RuntimeAttr? runtime_attr_SplitVcf
        RuntimeAttr? runtime_attr_bed_vs_hgsv
        RuntimeAttr? runtime_attr_bed_vs_pacbio
        RuntimeAttr? runtime_attr_bed_vs_bionano
        RuntimeAttr? runtime_attr_bed_vs_array
    }

    scatter (i in range(length(samples))){
        call annotate_il_pb_featuers_per_sample_per_bed.AnnotateILPBFeaturesPerSamplePerBed as AnnotateILPBFeaturesPerSamplePerBed{
            input:
                sample = samples[i],
                raw_manta = raw_mantas[i],
                raw_wham = raw_whams[i],
                raw_melt = raw_melts[i],
                raw_depth = raw_depths[i],
                gtgq = gtgqs[i],
                bed = beds[i],

                pacbio_seq = pacbio_seqs[i],
                pacbio_index = pacbio_indexes[i],

                hgsv_query = hgsv_queries[i],
                pacbio_query = pacbio_queries[i],
                bionano_query = bionano_queries[i],

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
                runtime_attr_bed_vs_hgsv = runtime_attr_bed_vs_hgsv,
                runtime_attr_bed_vs_pacbio = runtime_attr_bed_vs_pacbio,
                runtime_attr_bed_vs_bionano = runtime_attr_bed_vs_bionano,
                runtime_attr_bed_vs_array = runtime_attr_bed_vs_array
        }
    }
    output{
        Array[File] annotated_files = AnnotateILPBFeaturesPerSamplePerBed.integrated_file
        Array[File] vapor_plots = AnnotateILPBFeaturesPerSamplePerBed.vapor_plots
    }
}
