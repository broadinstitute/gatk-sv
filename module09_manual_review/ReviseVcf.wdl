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
import "TasksBenchmark.wdl" as mini_tasks
import "ReviseVcfWithManualResults.wdl" as revise_vcf_with_manual_results
workflow ReviseVcf{
    input{ 
        Array[File] vcf_files
        Array[File] vcf_indexes
        Array[String] chr_names
        File raw_SVs
        File SVID_to_Remove
        File MEI_DEL_Rescue
        File CPX_manual
        File CTX_manual
        File duplicated_SVID_manual

        String prefix

        String sv_benchmark_docker
        String sv_base_mini_docker
        String sv_pipeline_docker
        String sv_pipeline_hail_docker

        RuntimeAttr? runtime_attr_sort_merged_vcf
        RuntimeAttr? runtime_attr_override_revise_vcf
        RuntimeAttr? runtime_attr_override_add_raw_SVs
        RuntimeAttr? runtime_attr_override_concat_vcfs
        RuntimeAttr? runtime_attr_override_split_raw_SVs_per_chr
        RuntimeAttr? runtime_attr_concat_sharded_cluster
    }

    scatter(i in range(length(vcf_files))){
        call SplitRawSVsPerChr{
            input:
                raw_SVs = raw_SVs,
                chr_name = chr_names[i],
                sv_base_mini_docker = sv_base_mini_docker,
                runtime_attr_override = runtime_attr_override_split_raw_SVs_per_chr
        }
        call  revise_vcf_with_manual_results.ReviseVcfWithManualResults as ReviseVcfWithManualResults{
            input:
                vcf_file = vcf_files[i],
                vcf_index = vcf_indexes[i],
                chr_name = chr_names[i],
                raw_SVs = SplitRawSVsPerChr.raw_SV_per_chr,
                SVID_to_Remove = SVID_to_Remove,
                MEI_DEL_Rescue = MEI_DEL_Rescue,
                CPX_manual = CPX_manual,
                CTX_manual = CTX_manual,
                duplicated_SVID_manual = duplicated_SVID_manual,
                prefix = prefix,
                sv_benchmark_docker = sv_benchmark_docker,
                sv_base_mini_docker = sv_base_mini_docker,
                sv_pipeline_docker = sv_pipeline_docker,
                sv_pipeline_hail_docker = sv_pipeline_hail_docker,
                runtime_attr_override_revise_vcf = runtime_attr_override_revise_vcf,
                runtime_attr_override_add_raw_SVs = runtime_attr_override_add_raw_SVs,
                runtime_attr_override_concat_vcfs = runtime_attr_override_concat_vcfs,
                runtime_attr_sort_merged_vcf = runtime_attr_sort_merged_vcf,
                runtime_attr_concat_sharded_cluster = runtime_attr_concat_sharded_cluster
        }
    }

    output{
        Array[File] revised_vcf_list = ReviseVcfWithManualResults.revised_vcf
        Array[File] revised_idx_list = ReviseVcfWithManualResults.revised_vcf_idx
    }
 }
   

task SplitRawSVsPerChr{
    input{
        File raw_SVs
        String chr_name
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 5, 
        disk_gb: 10,
        boot_disk_gb: 30,
        preemptible_tries: 1,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    command<<<
        set -euo pipefail
        touch "raw_SVs_to_add.~{chr_name}.tsv"
        awk '{if ($1=="~{chr_name}") print}' ~{raw_SVs} >> "raw_SVs_to_add.~{chr_name}.tsv"
    >>>

    output{
        File raw_SV_per_chr = "raw_SVs_to_add.~{chr_name}.tsv"
    }

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: sv_base_mini_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}




