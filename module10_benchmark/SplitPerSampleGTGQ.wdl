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


workflow SplitPerSampleGTGQ{
    input{
        Array[File] cleanVcfs
        Array[File] cleanVcfIdxes
        Array[String] prefixes

        Array[String] samples

        String sv_pipeline_docker

        RuntimeAttr? runtime_split_per_sample_gtgq
    }

    scatter(i in range(length(cleanVcfs))){
        call mini_tasks.split_per_sample_gtgq{
            input:
                clean_vcf = cleanVcfs[i],
                clean_vcf_idx = cleanVcfIdxes[i],
                samples = samples,
                prefix = prefixes[i],
                sv_pipeline_docker = sv_pipeline_docker,
                runtime_attr_override = runtime_split_per_sample_gtgq
        }
    }

    output{
        Array[Array[File]] gtgq = split_per_sample_gtgq.gtgq_file
    }
}

