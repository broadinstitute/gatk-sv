#########################################################################################

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
import "HailMerge.wdl" as HailMerge
import "TasksMakeCohortVcf.wdl" as MiniTasks
import "RemoveDuplicateEvents.wdl" as remove_duplicate_events
import "CollectPEMetricsForCPX.wdl" as collect_pe_metrics_for_cpx
import "CollectLgCnvSupportForCPX.wdl" as collect_lg_cnv_supp_for_cpx
import "ReviseVcfWithManualResults.wdl" as revise_vcf_with_manual_results

### this script is designed to split CTX event from vcf, so they can be manually reformatted and inserted back
workflow SplitCtxForReformat{
    input{ 
        File vcf
        File vcf_idx
        String contig #use chromosome name if the input vcf is per-contig; else, put "whole_genome"


        String prefix
        Int n_per_split

        Boolean use_hail = false
        String? gcs_project

        String sv_benchmark_docker
        String sv_base_mini_docker
        String sv_pipeline_docker
        String sv_pipeline_hail_docker

        RuntimeAttr? runtime_attr_remove_duplicate_events_task 
        RuntimeAttr? runtime_attr_concat_sharded_cluster
        RuntimeAttr? runtime_attr_scatter_vcf
    }


    call MiniTasks.ScatterVcf{
        input:
          vcf = vcf,
          prefix = prefix,
          records_per_shard = n_per_split,
          sv_pipeline_docker = sv_pipeline_docker,
          runtime_attr_override = runtime_attr_scatter_vcf
      }

    scatter (i in range(length(ScatterVcf.shards))) {

        call ExtractCtxCalls{
            input:
                rerun_vcf = ScatterVcf.shards[i],
                rerun_vcf_idx = ScatterVcf.shards_idx[i],
                sv_pipeline_docker = sv_pipeline_docker,
                runtime_attr_override = runtime_attr_remove_duplicate_events_task
        }
    }

    call MiniTasks.ConcatVcfs {
        input:
            vcfs=ExtractCtxCalls.ctx_vcf,
            vcfs_idx=ExtractCtxCalls.ctx_vcf_idx,
            allow_overlaps=true,
            outfile_prefix="~{prefix}.ctx",
            sv_base_mini_docker=sv_base_mini_docker,
            runtime_attr_override=runtime_attr_concat_sharded_cluster
    }

    output{
        File ctx_vcf = ConcatVcfs.concat_vcf
        File ctx_vcf_idx = ConcatVcfs.concat_vcf_idx
    }
}



task ExtractCtxCallsV2{

    input{
        File rerun_vcf
        File rerun_vcf_idx
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 5, 
        disk_gb: 20,
        boot_disk_gb: 30,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])


    String prefix  = basename(rerun_vcf, ".vcf.gz")

    command<<<
        set -euo pipefail

        bcftools view -Oz --include 'INFO/SVTYPE="CTX"' -o ~{prefix}.ctx.vcf.gz ~{rerun_vcf}

        tabix -p vcf "~{prefix}.ctx.vcf.gz"

    >>>

    output{
        File ctx_vcf = "~{prefix}.ctx.vcf.gz"
        File ctx_vcf_idx ="~{prefix}.ctx.vcf.gz.tbi"
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

task ExtractCtxCalls{
    input{
        File rerun_vcf
        File rerun_vcf_idx
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 5, 
        disk_gb: 20,
        boot_disk_gb: 30,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])


    String prefix  = basename(rerun_vcf, ".vcf.gz")

    command<<<
        set -euo pipefail

        python <<CODE

        import os
        import pysam
        fin2=pysam.VariantFile("~{rerun_vcf}")
        fo2=pysam.VariantFile("~{prefix}.ctx.vcf.gz", 'w', header = fin2.header)
        for record in fin2:
            if record.info['SVTYPE']=="CTX":
                fo2.write(record)
        fin2.close()
        fo2.close()
        CODE

        tabix -p vcf "~{prefix}.ctx.vcf.gz"

    >>>

    output{
        File ctx_vcf = "~{prefix}.ctx.vcf.gz"
        File ctx_vcf_idx ="~{prefix}.ctx.vcf.gz.tbi"
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


