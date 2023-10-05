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
import "HailMerge.wdl" as HailMerge
import "TasksMakeCohortVcf.wdl" as MiniTasks
import "RemoveDuplicateEvents.wdl" as remove_duplicate_events
import "CollectPEMetricsForCPX.wdl" as collect_pe_metrics_for_cpx
import "CollectLgCnvSupportForCPX.wdl" as collect_lg_cnv_supp_for_cpx
import "ReviseVcfWithManualResults.wdl" as revise_vcf_with_manual_results

workflow ShardeManualReviewConcatVcfRerun{
    input{ 
        File vcf
        File vcf_idx
        String contig #use chromosome name if the input vcf is per-contig; else, put "whole_genome"
        String prefix
        File? raw_SVs

        Array[File] sharded_annotated_vcf
        Array[File] sharded_annotated_vcf_idx

        Boolean use_hail = false
        String? gcs_project

        String sv_benchmark_docker
        String sv_base_mini_docker
        String sv_pipeline_docker
        String sv_pipeline_hail_docker

        RuntimeAttr? runtime_attr_vcf2bed 
        RuntimeAttr? runtime_attr_vcf2bed_sm
        RuntimeAttr? runtime_attr_bnd_vs_mei
        RuntimeAttr? runtime_attr_collect_pe
        RuntimeAttr? runtime_attr_scatter_vcf
        RuntimeAttr? runtime_attr_split_script
        RuntimeAttr? runtime_attr_calcu_pe_stat
        RuntimeAttr? runtime_attr_extract_cpx_ctx
        RuntimeAttr? runtime_attr_extract_bnd_del
        RuntimeAttr? runtime_attr_concat_evidence
        RuntimeAttr? runtime_attr_split_raw_SVs_per_chr
        RuntimeAttr? runtime_attr_concat_sharded_cluster
        RuntimeAttr? runtime_attr_preconcat_sharded_cluster
        RuntimeAttr? runtime_attr_hail_merge_sharded_cluster
        RuntimeAttr? runtime_attr_fix_header_sharded_cluster
        RuntimeAttr? runtime_attr_generate_cpx_review_script
        RuntimeAttr? runtime_attr_remove_duplicate_events_task
        RuntimeAttr? runtime_attr_generate_cnv_segments_from_cpx
        RuntimeAttr? runtime_attr_get_vcf_header_with_members_info_line
        RuntimeAttr? runtime_attr_generate_cnv_segments_from_cpx
        RuntimeAttr? runtime_attr_extract_cpx_lg_cnv_by_batch
        RuntimeAttr? runtime_attr_seek_depth_supp_for_cpx
        RuntimeAttr? runtime_attr_concat_bed_Step1
        RuntimeAttr? runtime_attr_concat_bed_Step2
        RuntimeAttr? runtime_attr_add_raw_SVs
    }

    if (length(sharded_annotated_vcf) == 0) {
        call MiniTasks.GetVcfHeaderWithMembersInfoLine as GetVcfHeader_annotated {
            input:
                vcf_gz=vcf,
                prefix="~{prefix}.manual_review",
                sv_base_mini_docker=sv_base_mini_docker,
                runtime_attr_override=runtime_attr_get_vcf_header_with_members_info_line
        }
    }

    if (length(sharded_annotated_vcf) > 0) {
        if (use_hail) {
            call HailMerge.HailMerge as ConcatVcfsHail_annotated {
                input:
                    vcfs=sharded_annotated_vcf,
                    prefix="~{prefix}.manual_review",
                    gcs_project=gcs_project,
                    sv_base_mini_docker=sv_base_mini_docker,
                    sv_pipeline_docker=sv_pipeline_docker,
                    sv_pipeline_hail_docker=sv_pipeline_hail_docker,
                    runtime_attr_preconcat=runtime_attr_preconcat_sharded_cluster,
                    runtime_attr_hail_merge=runtime_attr_hail_merge_sharded_cluster,
                    runtime_attr_fix_header=runtime_attr_fix_header_sharded_cluster
            }
        }

        if (!use_hail) {
            call MiniTasks.ConcatVcfs as ConcatVcfs_annotated {
                input:
                    vcfs=sharded_annotated_vcf,
                    vcfs_idx=sharded_annotated_vcf_idx,
                    allow_overlaps=true,
                    outfile_prefix="~{prefix}.manual_review",
                    sv_base_mini_docker=sv_base_mini_docker,
                    runtime_attr_override=runtime_attr_concat_sharded_cluster
            }
        }
    }

    File reviewed_vcf = select_first([GetVcfHeader_annotated.out, ConcatVcfs_annotated.concat_vcf, ConcatVcfsHail_annotated.merged_vcf])
    File reviewed_vcf_idx = select_first([GetVcfHeader_annotated.out_idx, ConcatVcfs_annotated.concat_vcf_idx, ConcatVcfsHail_annotated.merged_vcf_index])



    if (defined(raw_SVs)){

        if (contig!="whole_genome"){
            call SplitRawSVsPerChr{
                input:
                    raw_SVs = raw_SVs,
                    contig = contig,
                    sv_base_mini_docker = sv_base_mini_docker,
                    runtime_attr_override = runtime_attr_split_raw_SVs_per_chr
            }
        }

        File split_raw_SVs = select_first([SplitRawSVsPerChr.raw_SV_per_chr, raw_SVs])

        call revise_vcf_with_manual_results.AddRawSVs{
            input:
                prefix = prefix,
                batch_name = contig,
                vcf_file = reviewed_vcf,
                raw_SVs = split_raw_SVs,
                sv_benchmark_docker = sv_benchmark_docker,
                runtime_attr_override = runtime_attr_add_raw_SVs
        }
    }

    output{
        File revised_output_vcf = select_first([AddRawSVs.vcf_with_raw_SVs, reviewed_vcf])
        File revised_output_vcf_idx = select_first([AddRawSVs.vcf_idx_with_raw_SVs, reviewed_vcf_idx])
    }
}

task SplitRawSVsPerChr{
    input{
        File? raw_SVs
        String contig
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
        touch "raw_SVs_to_add.~{contig}.tsv"
        awk '{if ($1=="~{contig}") print}' ~{raw_SVs} >> "raw_SVs_to_add.~{contig}.tsv"
    >>>

    output{
        File raw_SV_per_chr = "raw_SVs_to_add.~{contig}.tsv"
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


