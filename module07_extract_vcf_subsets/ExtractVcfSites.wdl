## Copyright Broad Institute, 2022
## 
##
## Consolidate boost scores per sample across all batches and write those scores
## directly into an input VCF
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

workflow ExtractVcfSites {
    input {
        Array[File] vcf_list
        Array[File] vcf_idx_list
        Array[String] contig_list
        Boolean shard_vcf
        String midfix
        String sv_pipeline_docker
        String sv_base_mini_docker
        String sv_pipeline_hail_docker
        RuntimeAttr? runtime_attr_override_extract_subset_vcf
    }

    scatter(i in range(length(vcf_list))){
        if (!shard_vcf){
            call ExtractSubsetSamples{
                input:
                    vcf = vcf_list[i],
                    vcf_idx = vcf_idx_list[i],
                    sample_list = sample_list
                    midfix = contig_list[i],
                    sv_pipeline_docker = sv_pipeline_docker,
                    runtime_attr_override = runtime_attr_override_extract_subset_vcf
           }
        }

     }

    output{
        Array[File?] subset_vcf_list = ExtractSubsetSamples.out_vcf
        Array[File?] subset_vcf_idx_list = ExtractSubsetSamples.out_vcf_idx
   }
}



task ExtractSubsetSamples {
    input {
        File vcf
        File vcf_idx
        File sample_list
        String midfix
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }


    Float input_size = size(vcf, "GB")
    Float base_disk_gb = 10.0
    RuntimeAttr runtime_default = object {
            mem_gb: 3,
            disk_gb: ceil(base_disk_gb + (input_size * 2.0)),
            cpu_cores: 1,
            preemptible_tries: 3,
            max_retries: 1,
            boot_disk_gb: 10
    }
    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
    runtime {
            memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
            disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
            cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
            preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
            maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
            docker: sv_pipeline_docker
            bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    String prefix = basename(vcf, '.vcf.gz')
    command <<<
        set -eu -o pipefail

        bcftools view -S ~{sample_list} ~{vcf} \
        | bgzip > ~{prefix}.~{midfix}.vcf.gz

        tabix -p vcf ~{prefix}.~{midfix}.vcf.gz

    >>>

    output {
        File out_vcf = "~{prefix}.~{midfix}.vcf.gz"
        File out_vcf_idx = "~{prefix}.~{midfix}.vcf.gz.tbi"
    }
}

