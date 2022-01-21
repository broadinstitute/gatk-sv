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
import "SplitPerSampleGTGQPerVcf.wdl" as split_per_sample_gtgq_per_vcf
import "SplitPerSampleGTGQPerSampleList.wdl" as split_per_sample_gtgq_per_samplelist

workflow SplitPerSampleGTGQ{
    input{
        Array[File] cleanVcfs
        Array[File] cleanVcfIdxes
        Array[File] cleanBeds

        Array[File] SampleLists
        Array[String] prefixes

        String rdpesr_benchmark_docker
        String sv_base_mini_docker

        RuntimeAttr? runtime_split_per_sample_gtgq
    }

    scatter(SampleList in SampleLists){
        call split_per_sample_gtgq_per_samplelist.SplitPerSampleGTGQPerSampleList as SplitPerSampleGTGQPerSampleList{
            input:
                cleanVcfs = cleanVcfs,
                cleanVcfIdxes = cleanVcfIdxes,
                cleanBeds = cleanBeds,
                SampleList = SampleList,
                prefixes = prefixes,

                rdpesr_benchmark_docker = rdpesr_benchmark_docker,
                sv_base_mini_docker = sv_base_mini_docker,

                runtime_split_per_sample_gtgq = runtime_split_per_sample_gtgq
        }
    }

    output{
        Array[Array[File]] gtgq = SplitPerSampleGTGQPerSampleList.gtgq_out
        Array[Array[File]] bed = SplitPerSampleGTGQPerSampleList.bed_out
    }
}

