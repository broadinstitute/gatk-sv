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
import "AnnotateInheritancePerSample.wdl" as annotate_inheritance_per_sample

workflow AnnotateInheritance{
    input{
        Array[File] cleanVcfs
        Array[String] contigs

        Array[String] probands_ids
        Array[String] father_ids
        Array[String] mother_ids

        String sv_pipeline_docker
        String sv_base_mini_docker
        String rdpesr_benchmark_docker

        RuntimeAttr? runtime_attr_split_vcf
        RuntimeAttr? runtime_attr_vcf2bed
        RuntimeAttr? runtime_attr_bed_comp
    }

    scatter (i in range(length(probands_ids))) {
        call annotate_inheritance_per_sample.AnnotateInheritancePerSample as AnnotateInheritancePerSample{
            input:
                cleanVcfs = cleanVcfs,
                contigs = contigs,

                probands_id = probands_ids[i],
                father_id = father_ids[i],
                mother_id = mother_ids[i],

                sv_pipeline_docker = sv_pipeline_docker,
                sv_base_mini_docker = sv_base_mini_docker,
                rdpesr_benchmark_docker = rdpesr_benchmark_docker,

                runtime_attr_split_vcf = runtime_attr_split_vcf,
                runtime_attr_vcf2bed = runtime_attr_vcf2bed,
                runtime_attr_bed_comp = runtime_attr_bed_comp
        }
    }

    output{
        Array[File] inheritance = AnnotateInheritancePerSample.inheritance
    }
}

