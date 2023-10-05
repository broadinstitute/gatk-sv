##########################################################################################

## Copyright Broad Institute, 2022
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
import "Utils.wdl" as Utils
import "HailMerge.wdl" as HailMerge
import "TasksMakeCohortVcf.wdl" as MiniTasks
import "RemoveDuplicateEvents.wdl" as remove_duplicate_events
import "EnforceMinNoCallRate.wdl" as enforce_min_no_call_rate
import "ApplyNoCallRateCutoffs.wdl" as apply_no_call_rate_cutoffs

workflow CombineVcfWithCtx{
    input{ 
        File vcf
        File vcf_idx
        File ctx_vcf
        File ctx_vcf_idx
        String contig


        String sv_pipeline_docker

        RuntimeAttr? runtime_attr_combine_ctx
        }

    call CombineWithCtx{
        input:
            vcf = vcf,
            vcf_idx = vcf_idx,
            contig = contig,
            ctx_vcf = ctx_vcf,
            ctx_vcf_idx = ctx_vcf_idx,
            sv_pipeline_docker = sv_pipeline_docker,
            runtime_attr_override = runtime_attr_combine_ctx

    }
    output{
        File output_vcf = CombineWithCtx.combined_vcf
        File output_vcf_idx = CombineWithCtx.combined_vcf_idx
    }
}


task CombineWithCtx{
    input{
        File vcf
        File vcf_idx
        File ctx_vcf
        File ctx_vcf_idx
        String contig
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 5, 
        disk_gb: 200,
        boot_disk_gb: 30,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    String prefix = basename(vcf,".vcf.gz")

    command<<<
        set -euo pipefail

        bcftools view ~{ctx_vcf} ~{contig} -O z -o ~{contig}.ctx.vcf.gz
        tabix -p vcf ~{contig}.ctx.vcf.gz
        
        echo "~{vcf}" > vcf_list
        echo "~{contig}.ctx.vcf.gz" >> vcf_list

        bcftools concat -a --allow-overlaps --output-type z --file-list vcf_list --output  "~{prefix}.with_ctx.vcf.gz"
        tabix -p vcf  "~{prefix}.with_ctx.vcf.gz"

    >>>

    output{
        File out_vcf_list = "vcf_list"
        File combined_vcf = "~{prefix}.with_ctx.vcf.gz"
        File combined_vcf_idx = "~{prefix}.with_ctx.vcf.gz.tbi"
    }

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


