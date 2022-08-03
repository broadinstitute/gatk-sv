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
import "SplitPerSampleGTGQandBEDPerSampleList.wdl" as split_per_sample_gtgq_and_bed_per_samplelist

workflow SplitPerSampleGTGQandBED{
    input{
        Array[File] cleanVcfs
        Array[File] cleanVcfIdxes
        Array[File]? cleanBeds

        Array[File] SampleLists
        Array[String] prefixes

        String rdpesr_benchmark_docker
        String sv_base_mini_docker
        String sv_pipeline_docker

        RuntimeAttr? runtime_split_per_sample_gtgq
        RuntimeAttr? runtime_attr_concat_gtgqs
        RuntimeAttr? runtime_attr_concat_refs

    }

    scatter(SampleList in SampleLists){
        call split_per_sample_gtgq_and_bed_per_samplelist.SplitPerSampleGTGQandBEDPerSampleList as SplitPerSampleGTGQandBEDPerSampleList{
            input:
                cleanVcfs = cleanVcfs,
                cleanVcfIdxes = cleanVcfIdxes,
                cleanBeds = cleanBeds,
                SampleList = SampleList,
                prefixes = prefixes,

                rdpesr_benchmark_docker = rdpesr_benchmark_docker,
                sv_base_mini_docker = sv_base_mini_docker,
                sv_pipeline_docker = sv_pipeline_docker,

                runtime_split_per_sample_gtgq = runtime_split_per_sample_gtgq,
                runtime_attr_concat_gtgqs = runtime_attr_concat_gtgqs,
                runtime_attr_concat_refs = runtime_attr_concat_refs
        }
    }

    output{
        Array[Array[File]] gtgqs = SplitPerSampleGTGQPerSampleList.gtgq_out
        Array[Array[File]] beds = SplitPerSampleGTGQPerSampleList.bed_out
    }
}

