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
import "AnnotateInheritancePerSamplePerVcf.wdl" as annotate_inheritance_per_sample_per_vcf

workflow AnnotateInheritancePerSample{
    input{
        Array[File] cleanVcfs
        Array[String] contigs

        String probands_id
        String father_id
        String mother_id


        String sv_pipeline_docker
        String sv_base_mini_docker
        String rdpesr_benchmark_docker

        RuntimeAttr? runtime_attr_split_vcf
        RuntimeAttr? runtime_attr_vcf2bed
        RuntimeAttr? runtime_attr_bed_comp
    }

    scatter (i in range(length(cleanVcfs))) {
        call annotate_inheritance_per_sample_per_vcf.AnnotateInheritancePerSamplePerVcf as AnnotateInheritancePerSamplePerVcf {
            input:
                cleanVcf = cleanVcfs[i],
                contig = contigs[i],

                probands_id = probands_id,
                father_id = father_id,
                mother_id = mother_id,

                sv_pipeline_docker = sv_pipeline_docker,
                sv_base_mini_docker = sv_base_mini_docker,
                rdpesr_benchmark_docker = rdpesr_benchmark_docker,

                runtime_attr_split_vcf = runtime_attr_split_vcf,
                runtime_attr_vcf2bed = runtime_attr_vcf2bed,
                runtime_attr_bed_comp = runtime_attr_bed_comp
        }
    }

    call mini_tasks.ConcatBeds as concat_inheri{
        input:
            shard_bed_files = AnnotateInheritancePerSamplePerVcf.inheritance,
            prefix = probands_id,
            sv_base_mini_docker = sv_base_mini_docker
    }

    output{
        File inheritance = concat_inheri.merged_bed_file
    }
}

