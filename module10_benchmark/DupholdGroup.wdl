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
import "Duphold.wdl" as duphold
import "DupholdRequestPay.wdl" as DupholdRequestPay
workflow DupholdGroup{
  input{
    Array[String] prefixes
    Array[String] samples
    Array[String] bam_or_cram_files
    Array[String] bam_or_cram_indexes
    File vcf_file
    File ref_fasta
    File ref_fai
    File ref_dict
    File contig_list
    String duphold_docker
    String sv_base_mini_docker
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_split_vcf_per_sample
    RuntimeAttr? runtime_attr_duphold
    RuntimeAttr? runtime_attr_bcf2vcf
    RuntimeAttr? runtime_attr_LocalizeCram
    RuntimeAttr? runtime_attr_SplitVcf
    RuntimeAttr? runtime_attr_ConcatVcfs
  }

  scatter (i in range(length(samples))){

  }

    scatter (i in range(length(samples))){
        call split_per_sample_vcf{
            input:
                vcf = vcf_file,
                sample = samples[i],
                sv_pipeline_docker = sv_pipeline_docker,
                runtime_attr_override = runtime_attr_split_vcf_per_sample
        }

        call duphold.Duphold as duphold{
            input:
                prefix = prefixes[i],
                bam_or_cram_file = bam_or_cram_files[i],
                bam_or_cram_index = bam_or_cram_indexes[i],
                vcf_file = split_per_sample_vcf.vcf_file,
                ref_fasta = ref_fasta,
                ref_fai = ref_fai,
                ref_dict = ref_dict,
                contig_list = contig_list,
                duphold_docker = duphold_docker,
                sv_base_mini_docker = sv_base_mini_docker,
                sv_pipeline_docker = sv_pipeline_docker,
                runtime_attr_duphold = runtime_attr_duphold,
                runtime_attr_bcf2vcf = runtime_attr_bcf2vcf,
                runtime_attr_LocalizeCram = runtime_attr_LocalizeCram,
                runtime_attr_SplitVcf = runtime_attr_SplitVcf,
                runtime_attr_ConcatVcfs = runtime_attr_ConcatVcfs

        }
    }

    output{
      Array[File] vcf = duphold.vcf
      Array[File] vcf_idx = duphold.vcf_idx
    }
  }



task split_per_sample_vcf{
    input{
        File vcf
        String sample
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
        }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 3.75, 
        disk_gb: 10,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])


    output {
        File vcf_file = "~{sample}.vcf.gz"
        File vcf_idx ="~{sample}.vcf.gz.tbi"
    }
    command <<<

        set -Eeuo pipefail
        
        bcftools view -s ~{sample} ~{vcf} | grep -v "0/0" | bgzip > ~{sample}.vcf.gz
        tabix -p vcf ~{sample}.vcf.gz

    >>>
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: sv_pipeline_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}




