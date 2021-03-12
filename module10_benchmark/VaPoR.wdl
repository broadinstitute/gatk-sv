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
import "VaPoRBedPerContig.wdl" as vapor_bed
import "VaPoRVcf.wdl" as vapor_vcf
workflow VaPoR{
    input{
        String prefix
        String bam_or_cram_file
        String bam_or_cram_index
        File? vcf_file
        File? bed_file
        File ref_fasta
        File ref_fai
        File ref_dict
        Array[String] contigs
        String vapor_docker
        String sv_base_mini_docker
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_Vapor 
        RuntimeAttr? runtime_attr_bcf2vcf
        RuntimeAttr? runtime_attr_vcf2bed
        RuntimeAttr? runtime_attr_SplitVcf
        RuntimeAttr? runtime_attr_ConcatBeds
    }

    if (defined(vcf_file)) {
        call vapor_vcf.VaPoRVcf as VaPoR_vcf{
            input:
                prefix = prefix,
                bam_or_cram_file = bam_or_cram_file,
                bam_or_cram_index = bam_or_cram_index,
                vcf_file = vcf_file,
                ref_fasta = ref_fasta,
                ref_fai = ref_fai,
                ref_dict = ref_dict,
                contigs = contigs,
                vapor_docker = vapor_docker,
                sv_base_mini_docker = sv_base_mini_docker,
                sv_pipeline_docker = sv_pipeline_docker,
                runtime_attr_Vapor    = runtime_attr_Vapor,
                runtime_attr_bcf2vcf = runtime_attr_bcf2vcf,
                runtime_attr_vcf2bed = runtime_attr_vcf2bed,
                runtime_attr_SplitVcf = runtime_attr_SplitVcf,
                runtime_attr_ConcatBeds = runtime_attr_ConcatBeds
        }
    }

    if (defined(bed_file)) {
        call vapor_bed.VaPoRBed as VaPoR_bed{
            input:
                prefix = prefix,
                bam_or_cram_file = bam_or_cram_file,
                bam_or_cram_index = bam_or_cram_index,
                bed_file = bed_file,
                ref_fasta = ref_fasta,
                ref_fai = ref_fai,
                ref_dict = ref_dict,
                contigs = contigs,
                vapor_docker = vapor_docker,
                sv_base_mini_docker = sv_base_mini_docker,
                sv_pipeline_docker = sv_pipeline_docker,
                runtime_attr_Vapor    = runtime_attr_Vapor,
                runtime_attr_bcf2vcf = runtime_attr_bcf2vcf,
                runtime_attr_vcf2bed = runtime_attr_vcf2bed,
                runtime_attr_SplitVcf = runtime_attr_SplitVcf,
                runtime_attr_ConcatBeds = runtime_attr_ConcatBeds
        }
    }
    output{
            File? vcf_out = VaPoR_vcf.bed
            File? bed_out = VaPoR_bed.bed
            Array[File]? vcf_plots = VaPoR_vcf.plots
            Array[File]? bed_plots = VaPoR_bed.plots
    }
}


