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
import "ShardedManualRevise.wdl" as sharded_manual_revise
import "ManualRevise.wdl" as manual_revise
workflow ReviseVcfWithManualResults{
    input{ 
        File vcf_file
        File vcf_index
        File? raw_SVs
        File SVID_to_Remove
        File MEI_DEL_Rescue
        File CPX_manual
        File CTX_manual
        #File duplicated_SVID_manual

        String prefix
        String chr_name

        String sv_benchmark_docker
        String sv_base_mini_docker
        String sv_pipeline_docker
        String sv_pipeline_hail_docker

        RuntimeAttr? runtime_attr_override_revise_vcf
        RuntimeAttr? runtime_attr_override_add_raw_SVs
        RuntimeAttr? runtime_attr_override_concat_vcfs
        RuntimeAttr? runtime_attr_sort_merged_vcf
        RuntimeAttr? runtime_attr_concat_sharded_cluster
        RuntimeAttr? runtime_attr_override_split_cpx_per_contig
    }

    if (chr_name != "whole_genome"){
        call SplitFilePerContig{
            input:
                input_file = CPX_manual,
                contig = chr_name,
                sv_base_mini_docker = sv_base_mini_docker,
                runtime_attr_override = runtime_attr_override_split_cpx_per_contig
            }
    }

    File CPX_manual_output = select_first([SplitFilePerContig.output_file, CPX_manual])
    call manual_revise.ManualRevise as ManualRevise{
        input:
            vcf = vcf_file,
            prefix = prefix,
            SVID_to_Remove = SVID_to_Remove,
            MEI_DEL_Rescue = MEI_DEL_Rescue,
            CPX_manual = CPX_manual_output,
            CTX_manual = CTX_manual,
            use_hail = false,
            sv_benchmark_docker = sv_benchmark_docker,
            sv_pipeline_docker = sv_pipeline_docker,
            sv_base_mini_docker = sv_base_mini_docker,
            sv_pipeline_hail_docker = sv_pipeline_hail_docker,
            runtime_attr_sort_merged_vcf = runtime_attr_sort_merged_vcf,
            runtime_attr_concat_sharded_cluster = runtime_attr_concat_sharded_cluster
        }


    if(defined(raw_SVs)){
        call AddRawSVs{
            input:
                prefix = prefix,
                chr_name = chr_name,
                vcf_file = ManualRevise.cpx_ctx_vcf,
                raw_SVs = raw_SVs,
                sv_benchmark_docker = sv_benchmark_docker,
                runtime_attr_override = runtime_attr_override_add_raw_SVs
        }
    }

    File Cpx_Ctx_vcf = select_first([ManualRevise.cpx_ctx_vcf, AddRawSVs.vcf_with_raw_SVs])

    call mini_tasks.ConcatVcfs as ConcatVcfs{
        input:
            vcfs = [ManualRevise.manual_revised_vcf, Cpx_Ctx_vcf],
            outfile_prefix = "~{prefix}.~{chr_name}.manually_revised",
            sv_base_mini_docker = sv_base_mini_docker,
            runtime_attr_override = runtime_attr_override_concat_vcfs
    }

    output{
        File revised_vcf = ConcatVcfs.concat_vcf
        File revised_vcf_idx = ConcatVcfs.concat_vcf_idx

    }
 }
   
task SplitFilePerContig{
    input{
        File input_file
        String contig
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_override
        }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 7.5, 
        disk_gb: 15,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    String prefix = basename(input_file,'.tsv')
    command<<<
        grep ".~{contig}." ~{input_file} > ~{prefix}.~{contig}
    >>>
    output{
        File output_file = "~{prefix}.~{contig}"
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

task ReviseVcf{
    input{
        File vcf_file
        File vcf_index
        File SVID_to_Remove
        File MEI_DEL_Rescue
        File CPX_manual
        File CTX_manual
        File duplicated_SVID_manual
        String sv_benchmark_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 7.5, 
        disk_gb: ceil(5.0 +  size(vcf_file, "GB")*3),
        boot_disk_gb: 30,
        preemptible_tries: 1,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])


    String prefix  = basename(vcf_file, ".vcf.gz")
    command<<<
        set -euo pipefail

        python /src/revise_vcf_with_manual_review_results.py \
         ~{vcf_file} ~{prefix}.Manual_Revised.vcf.gz \
          --cpx_vcf ~{prefix}.CPX_CTX.vcf.gz \
          --SVID_to_remove  ~{SVID_to_Remove} \
          --MEI_DEL_rescue  ~{MEI_DEL_Rescue} \
          --CPX_manual  ~{CPX_manual} \
          --CTX_manual  ~{CTX_manual} \
          --duplicated_SVID_manual ~{duplicated_SVID_manual}
    >>>

    output{
        File manual_revised_vcf = "~{prefix}.Manual_Revised.vcf.gz"
        File cpx_ctx_vcf = "~{prefix}.CPX_CTX.vcf.gz"
    }

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: sv_benchmark_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}


task AddRawSVs{
    input{
        String prefix
        String chr_name
        File vcf_file
        File? raw_SVs
        String sv_benchmark_docker
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
        tabix -p vcf ~{vcf_file}
        python /src/add_manual_review_ctx_to_vcf.py \
            --vcf ~{vcf_file} \
            --reviewed-events-file ~{raw_SVs} \
            --cohort-name ~{prefix} \
            --batch-name ~{chr_name} \
            --out ~{prefix}.with_raw.vcf.gz 
    >>>

    output{
        File vcf_with_raw_SVs = "~{prefix}.with_raw.vcf.gz"
    }

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: sv_benchmark_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}


