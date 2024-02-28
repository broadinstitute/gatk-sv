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
import "ShardedManualReview.wdl" as sharded_manual_review

workflow ManualReview{
    input{ 
        File SVID_to_Remove #large CNVs, mCNVs that failed manual review, or redundant large CNVs that overlaps with Seg Dup and mCNVs
        File SVID_to_Remove_DEL_bump
        File CTX_manual
        #File CPX_manual
        #File duplicated_SVID_manual
        
        File LINE1_Ref
        File HERVK_Ref

        Array[File] vcf_list
        Array[File] vcf_idx_list
        Array[String] contig_list     
        Array[String] prefix_list

        Array[String] batch_name_list
        Array[File] PE_metrics
        Array[File] PE_metrics_idxes
        Array[File] Depth_DEL_beds
        Array[File] Depth_DUP_beds

        Int n_per_split
        File sample_PE_metrics
        File sample_depth_calls
        File? raw_SVs

        Boolean use_hail = false
        String? gcs_project
        Boolean run_fix_ends
        Boolean clean_del_bump

        String sv_base_mini_docker
        String sv_pipeline_docker
        String sv_pipeline_hail_docker

        RuntimeAttr? runtime_attr_vcf2bed 
        RuntimeAttr? runtime_attr_vcf2bed_sm
        RuntimeAttr? runtime_attr_extract_cpx_ctx
        RuntimeAttr? runtime_attr_extract_bnd_del
        RuntimeAttr? runtime_attr_bnd_vs_mei
        RuntimeAttr? runtime_attr_collect_pe
        RuntimeAttr? runtime_attr_split_script
        RuntimeAttr? runtime_attr_calcu_pe_stat
        RuntimeAttr? runtime_attr_concat_evidence
        RuntimeAttr? runtime_attr_generate_cpx_review_script
        RuntimeAttr? runtime_attr_split_raw_SVs_per_chr
        RuntimeAttr? runtime_attr_scatter_vcf
        RuntimeAttr? runtime_attr_get_vcf_header_with_members_info_line
        RuntimeAttr? runtime_attr_concat_sharded_cluster
        RuntimeAttr? runtime_attr_preconcat_sharded_cluster
        RuntimeAttr? runtime_attr_hail_merge_sharded_cluster
        RuntimeAttr? runtime_attr_fix_header_sharded_cluster
    }

    scatter (i in range(length(vcf_list))){
        call sharded_manual_review.ShardedManualReview as ShardedManualReview{
            input:
                SVID_to_Remove = SVID_to_Remove,
                SVID_to_Remove_DEL_bump = SVID_to_Remove_DEL_bump,
                CTX_manual = CTX_manual,
                #CPX_manual = CPX_manual,
                #duplicated_SVID_manual = duplicated_SVID_manual,

                LINE1_Ref = LINE1_Ref,
                HERVK_Ref = HERVK_Ref,

                vcf = vcf_list[i],
                vcf_idx = vcf_idx_list[i],
                contig = contig_list[i],

                batch_name_list = batch_name_list,
                PE_metrics = PE_metrics,
                PE_metrics_idxes  = PE_metrics_idxes,
                Depth_DEL_beds = Depth_DEL_beds,
                Depth_DUP_beds = Depth_DUP_beds,

                prefix = prefix_list[i],
                n_per_split = n_per_split,
                sample_PE_metrics = sample_PE_metrics,
                sample_depth_calls = sample_depth_calls,
                raw_SVs = raw_SVs,

                use_hail  = use_hail,
                gcs_project = gcs_project,
                run_fix_ends = run_fix_ends,
                clean_del_bump = clean_del_bump,

                sv_base_mini_docker = sv_base_mini_docker,
                sv_pipeline_docker = sv_pipeline_docker,

                runtime_attr_vcf2bed = runtime_attr_vcf2bed,
                runtime_attr_vcf2bed_sm  = runtime_attr_vcf2bed_sm,
                runtime_attr_extract_cpx_ctx = runtime_attr_extract_cpx_ctx,
                runtime_attr_extract_bnd_del = runtime_attr_extract_bnd_del,
                runtime_attr_bnd_vs_mei  = runtime_attr_bnd_vs_mei,
                runtime_attr_collect_pe  = runtime_attr_collect_pe,
                runtime_attr_split_script = runtime_attr_split_script,
                runtime_attr_calcu_pe_stat = runtime_attr_calcu_pe_stat,
                runtime_attr_concat_evidence = runtime_attr_concat_evidence,
                runtime_attr_generate_cpx_review_script  = runtime_attr_generate_cpx_review_script,
                runtime_attr_split_raw_SVs_per_chr = runtime_attr_split_raw_SVs_per_chr,
                runtime_attr_scatter_vcf = runtime_attr_scatter_vcf,
                runtime_attr_get_vcf_header_with_members_info_line = runtime_attr_get_vcf_header_with_members_info_line,
                runtime_attr_concat_sharded_cluster = runtime_attr_concat_sharded_cluster,
                sv_pipeline_hail_docker = sv_pipeline_hail_docker,
                runtime_attr_preconcat_sharded_cluster = runtime_attr_preconcat_sharded_cluster,
                runtime_attr_hail_merge_sharded_cluster = runtime_attr_hail_merge_sharded_cluster,
                runtime_attr_fix_header_sharded_cluster = runtime_attr_fix_header_sharded_cluster

        }
    }

    output{
        Array[File] reviewed_vcfs = ShardedManualReview.revised_output_vcf
        Array[File] reviewed_vcf_idxes = ShardedManualReview.revised_output_vcf_idx
        Array[File] cpx_evidences = ShardedManualReview.cpx_evidences
        Array[File] ctx_evidences = ShardedManualReview.ctx_evidences
    }


}
