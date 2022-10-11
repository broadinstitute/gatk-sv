version 1.0

import "CombineSRBothsidePass.wdl" as CombineSRBothsidePass
import "VcfClusterSingleChromsome.wdl" as VcfClusterContig
import "TasksMakeCohortVcf.wdl" as MiniTasks
import "HailMerge.wdl" as HailMerge
import "HarmonizeHeaders.wdl" as HarmonizeHeaders
import "MergePesrDepth.wdl" as MergePesrDepth
import "Utils.wdl" as Utils

workflow CombineBatches {
  input {
    String cohort_name
    Array[String] batches

    Boolean merge_vcfs = false

    Array[File] pesr_vcfs
    Array[File] depth_vcfs

    Array[File] raw_sr_bothside_pass_files
    Array[File] raw_sr_background_fail_files

    File contig_list
    Int localize_shard_size = 100000
    File pe_exclude_list
    File depth_exclude_list
    Float min_sr_background_fail_batches

    File empty_file

    Boolean use_hail = false
    String? gcs_project

    String sv_base_mini_docker
    String sv_pipeline_docker
    String sv_pipeline_hail_docker

    # overrides for local tasks
    RuntimeAttr? runtime_override_update_sr_list
    RuntimeAttr? runtime_override_merge_pesr_depth
    RuntimeAttr? runtime_override_reheader
    RuntimeAttr? runtime_override_pull_header

    # overrides for mini tasks
    RuntimeAttr? runtime_attr_get_non_ref_vids
    RuntimeAttr? runtime_attr_calculate_support_frac
    RuntimeAttr? runtime_override_clean_background_fail
    RuntimeAttr? runtime_override_concat
    RuntimeAttr? runtime_override_concat_pesr_depth
    RuntimeAttr? runtime_override_update_fix_pesr_header
    RuntimeAttr? runtime_override_count_samples
    RuntimeAttr? runtime_override_preconcat_pesr_depth
    RuntimeAttr? runtime_override_hail_merge_pesr_depth
    RuntimeAttr? runtime_override_fix_header_pesr_depth
    RuntimeAttr? runtime_override_concat_large_pesr_depth

    # overrides for VcfClusterContig
    RuntimeAttr? runtime_override_localize_vcfs
    RuntimeAttr? runtime_override_join_vcfs
    RuntimeAttr? runtime_override_fix_multiallelic
    RuntimeAttr? runtime_override_fix_ev_tags
    RuntimeAttr? runtime_override_subset_bothside_pass
    RuntimeAttr? runtime_override_subset_background_fail
    RuntimeAttr? runtime_override_subset_sv_type
    RuntimeAttr? runtime_override_shard_clusters
    RuntimeAttr? runtime_override_shard_vids
    RuntimeAttr? runtime_override_pull_vcf_shard
    RuntimeAttr? runtime_override_svtk_vcf_cluster
    RuntimeAttr? runtime_override_get_vcf_header_with_members_info_line
    RuntimeAttr? runtime_override_concat_vcf_cluster
    RuntimeAttr? runtime_override_concat_svtypes
    RuntimeAttr? runtime_override_concat_sharded_cluster
    RuntimeAttr? runtime_override_make_sites_only
    RuntimeAttr? runtime_override_sort_merged_vcf_cluster
    RuntimeAttr? runtime_override_preconcat_sharded_cluster
    RuntimeAttr? runtime_override_hail_merge_sharded_cluster
    RuntimeAttr? runtime_override_fix_header_sharded_cluster

    # overerides for merge pesr depth
    RuntimeAttr? runtime_override_shard_clusters_mpd
    RuntimeAttr? runtime_override_shard_vids_mpd
    RuntimeAttr? runtime_override_pull_vcf_shard_mpd
    RuntimeAttr? runtime_override_merge_pesr_depth_mpd

    RuntimeAttr? runtime_override_sort_merged_vcf_mpd
    RuntimeAttr? runtime_override_subset_small_mpd
    RuntimeAttr? runtime_override_subset_large_mpd
    RuntimeAttr? runtime_override_make_sites_only_mpd
    RuntimeAttr? runtime_override_concat_large_pesr_depth_mpd
    RuntimeAttr? runtime_override_concat_shards_mpd

    RuntimeAttr? runtime_override_preconcat_large_pesr_depth_mpd
    RuntimeAttr? runtime_override_hail_merge_large_pesr_depth_mpd
    RuntimeAttr? runtime_override_fix_header_large_pesr_depth_mpd

    RuntimeAttr? runtime_override_preconcat_pesr_depth_shards_mpd
    RuntimeAttr? runtime_override_hail_merge_pesr_depth_shards_mpd
    RuntimeAttr? runtime_override_fix_header_pesr_depth_shards_mpd

  }

  # Preprocess some inputs
  call CombineSRBothsidePass.CombineSRBothsidePass {
    input:
      pesr_vcfs=pesr_vcfs,
      raw_sr_bothside_pass_files=raw_sr_bothside_pass_files,
      prefix="~{cohort_name}.sr_bothside_pass",
      sv_base_mini_docker=sv_base_mini_docker,
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_get_non_ref_vids=runtime_attr_get_non_ref_vids,
      runtime_attr_calculate_support_frac=runtime_attr_calculate_support_frac
  }

  Float min_background_fail_first_col = min_sr_background_fail_batches * length(raw_sr_background_fail_files)
  call MiniTasks.CatUncompressedFiles as CleanBackgroundFail {
    input:
      shards=raw_sr_background_fail_files,
      filter_command="sort | uniq -c | awk -v OFS='\\t' '{if($1 >= ~{min_background_fail_first_col}) print $2}'",
      outfile_name="~{cohort_name}.background_fail.txt",
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_override_clean_background_fail
  }

  call Utils.CountSamples {
    input:
      vcf=depth_vcfs[0],
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_override_count_samples
  }

  #Scatter per chromosome
  Array[String] contigs = transpose(read_tsv(contig_list))[0]
  scatter ( contig in contigs ) {

    #Subset PESR VCFs to single chromosome & cluster
    #Note: also subsets bothside_pass and background_fail files to variants
    #present on chromosome of interest
    call VcfClusterContig.VcfClusterSingleChrom as ClusterPesr {
      input:
        vcfs=pesr_vcfs,
        num_samples=CountSamples.num_samples,
        batches=batches,
        prefix="~{cohort_name}.~{contig}.pesr",
        dist=300,
        frac=0.1,
        sample_overlap=0.5,
        exclude_list=pe_exclude_list,
        sv_size=50,
        sv_types=["DEL","DUP","INV","BND","INS"],
        contig=contig,
        evidence_type="pesr",
        cohort_name=cohort_name,
        localize_shard_size=localize_shard_size,
        subset_sr_lists=true,
        bothside_pass=CombineSRBothsidePass.out,
        background_fail=CleanBackgroundFail.outfile,
        empty_file=empty_file,
        use_hail=use_hail,
        gcs_project=gcs_project,
        sv_pipeline_docker=sv_pipeline_docker,
        sv_pipeline_hail_docker=sv_pipeline_hail_docker,
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_override_localize_vcfs = runtime_override_localize_vcfs,
        runtime_override_join_vcfs = runtime_override_join_vcfs,
        runtime_override_fix_multiallelic = runtime_override_fix_multiallelic,
        runtime_override_fix_ev_tags = runtime_override_fix_ev_tags,
        runtime_override_subset_bothside_pass=runtime_override_subset_bothside_pass,
        runtime_override_subset_background_fail=runtime_override_subset_background_fail,
        runtime_override_subset_sv_type=runtime_override_subset_sv_type,
        runtime_override_shard_clusters=runtime_override_shard_clusters,
        runtime_override_shard_vids=runtime_override_shard_vids,
        runtime_override_pull_vcf_shard=runtime_override_pull_vcf_shard,
        runtime_override_svtk_vcf_cluster=runtime_override_svtk_vcf_cluster,
        runtime_override_get_vcf_header_with_members_info_line=runtime_override_get_vcf_header_with_members_info_line,
        runtime_override_concat_vcf_cluster=runtime_override_concat_vcf_cluster,
        runtime_override_concat_svtypes=runtime_override_concat_svtypes,
        runtime_override_concat_sharded_cluster=runtime_override_concat_sharded_cluster,
        runtime_override_make_sites_only=runtime_override_make_sites_only,
        runtime_override_sort_merged_vcf=runtime_override_sort_merged_vcf_cluster,
        runtime_override_preconcat_sharded_cluster=runtime_override_preconcat_sharded_cluster,
        runtime_override_hail_merge_sharded_cluster=runtime_override_hail_merge_sharded_cluster,
        runtime_override_fix_header_sharded_cluster=runtime_override_fix_header_sharded_cluster
    }

    #Subset RD VCFs to single chromosome & cluster
    call VcfClusterContig.VcfClusterSingleChrom as ClusterDepth {
      input:
        vcfs=depth_vcfs,
        num_samples=CountSamples.num_samples,
        batches=batches,
        prefix="~{cohort_name}.~{contig}.depth",
        dist=500000,
        frac=0.5,
        sample_overlap=0.5,
        exclude_list=depth_exclude_list,
        sv_size=5000,
        sv_types=["DEL","DUP"],
        contig=contig,
        evidence_type="depth",
        cohort_name=cohort_name,
        localize_shard_size=localize_shard_size,
        subset_sr_lists=false,
        bothside_pass=CombineSRBothsidePass.out,
        background_fail=CleanBackgroundFail.outfile,
        empty_file=empty_file,
        use_hail=use_hail,
        gcs_project=gcs_project,
        sv_pipeline_docker=sv_pipeline_docker,
        sv_pipeline_hail_docker=sv_pipeline_hail_docker,
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_override_localize_vcfs = runtime_override_localize_vcfs,
        runtime_override_join_vcfs = runtime_override_join_vcfs,
        runtime_override_fix_multiallelic = runtime_override_fix_multiallelic,
        runtime_override_fix_ev_tags = runtime_override_fix_ev_tags,
        runtime_override_subset_bothside_pass=runtime_override_subset_bothside_pass,
        runtime_override_subset_background_fail=runtime_override_subset_background_fail,
        runtime_override_subset_sv_type=runtime_override_subset_sv_type,
        runtime_override_shard_clusters=runtime_override_shard_clusters,
        runtime_override_shard_vids=runtime_override_shard_vids,
        runtime_override_svtk_vcf_cluster=runtime_override_svtk_vcf_cluster,
        runtime_override_get_vcf_header_with_members_info_line=runtime_override_get_vcf_header_with_members_info_line,
        runtime_override_concat_vcf_cluster=runtime_override_concat_vcf_cluster,
        runtime_override_concat_svtypes=runtime_override_concat_svtypes,
        runtime_override_concat_sharded_cluster=runtime_override_concat_sharded_cluster,
        runtime_override_make_sites_only=runtime_override_make_sites_only,
        runtime_override_sort_merged_vcf=runtime_override_sort_merged_vcf_cluster,
        runtime_override_preconcat_sharded_cluster=runtime_override_preconcat_sharded_cluster,
        runtime_override_hail_merge_sharded_cluster=runtime_override_hail_merge_sharded_cluster,
        runtime_override_fix_header_sharded_cluster=runtime_override_fix_header_sharded_cluster
    }

    call MiniTasks.ConcatVcfs as ConcatPesrSitesOnly {
      input:
        vcfs=ClusterPesr.clustered_vcfs,
        vcfs_idx=ClusterPesr.clustered_vcf_indexes,
        naive=true,
        generate_index=false,
        sites_only=true,
        outfile_prefix="~{cohort_name}.clustered_pesr.sites_only",
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_attr_override=runtime_override_concat
    }

    #Update SR background fail & bothside pass files (1)
    call MiniTasks.UpdateSrList as UpdateBackgroundFailFirst {
      input:
        vcf=ConcatPesrSitesOnly.concat_vcf,
        original_list=ClusterPesr.filtered_background_fail,
        outfile="~{cohort_name}.~{contig}.sr_background_fail.updated.txt",
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_override_update_sr_list
    }
    call MiniTasks.UpdateSrList as UpdateBothsidePassFirst {
      input:
        vcf=ConcatPesrSitesOnly.concat_vcf,
        original_list=ClusterPesr.filtered_bothside_pass,
        outfile="~{cohort_name}.~{contig}.sr_bothside_pass.updated.txt",
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_override_update_sr_list
    }

    call HarmonizeHeaders.HarmonizeHeaders {
      input:
        header_vcf=ClusterDepth.clustered_vcfs[0],
        vcfs=ClusterPesr.clustered_vcfs,
        prefix="~{cohort_name}.~{contig}.harmonize_headers",
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_override_reheader=runtime_override_reheader,
        runtime_override_pull_header=runtime_override_pull_header
    }

    call MergePesrDepth.MergePesrDepth as MergeDeletions {
      input:
        subtyped_pesr_vcf=HarmonizeHeaders.out[0],
        subtyped_depth_vcf=ClusterDepth.clustered_vcfs[0],
        svtype="DEL",
        num_samples=CountSamples.num_samples,
        prefix="~{cohort_name}.~{contig}.merge_del",
        cohort_name=cohort_name,
        contig=contig,
        use_hail=use_hail,
        gcs_project=gcs_project,
        sv_base_mini_docker=sv_base_mini_docker,
        sv_pipeline_docker=sv_pipeline_docker,
        sv_pipeline_hail_docker=sv_pipeline_hail_docker,
        runtime_override_shard_clusters=runtime_override_shard_clusters_mpd,
        runtime_override_shard_vids=runtime_override_shard_vids_mpd,
        runtime_override_pull_vcf_shard=runtime_override_pull_vcf_shard_mpd,
        runtime_override_merge_pesr_depth=runtime_override_merge_pesr_depth_mpd,
        runtime_override_sort_merged_vcf=runtime_override_sort_merged_vcf_mpd,
        runtime_override_subset_small=runtime_override_subset_small_mpd,
        runtime_override_subset_large=runtime_override_subset_large_mpd,
        runtime_override_make_sites_only=runtime_override_make_sites_only_mpd,
        runtime_override_concat_large_pesr_depth=runtime_override_concat_large_pesr_depth_mpd,
        runtime_override_concat_shards=runtime_override_concat_shards_mpd,
        runtime_override_preconcat_large_pesr_depth=runtime_override_preconcat_large_pesr_depth_mpd,
        runtime_override_hail_merge_large_pesr_depth=runtime_override_hail_merge_large_pesr_depth_mpd,
        runtime_override_fix_header_large_pesr_depth=runtime_override_fix_header_large_pesr_depth_mpd,
        runtime_override_preconcat_pesr_depth_shards=runtime_override_preconcat_pesr_depth_shards_mpd,
        runtime_override_hail_merge_pesr_depth_shards=runtime_override_hail_merge_pesr_depth_shards_mpd,
        runtime_override_fix_header_pesr_depth_shards=runtime_override_fix_header_pesr_depth_shards_mpd
    }

    call MergePesrDepth.MergePesrDepth as MergeDuplications {
      input:
        subtyped_pesr_vcf=HarmonizeHeaders.out[1],
        subtyped_depth_vcf=ClusterDepth.clustered_vcfs[1],
        svtype="DUP",
        num_samples=CountSamples.num_samples,
        prefix="~{cohort_name}.~{contig}.merge_dup",
        cohort_name=cohort_name,
        contig=contig,
        use_hail=use_hail,
        gcs_project=gcs_project,
        sv_base_mini_docker=sv_base_mini_docker,
        sv_pipeline_docker=sv_pipeline_docker,
        sv_pipeline_hail_docker=sv_pipeline_hail_docker,
        runtime_override_shard_clusters=runtime_override_shard_clusters_mpd,
        runtime_override_shard_vids=runtime_override_shard_vids_mpd,
        runtime_override_pull_vcf_shard=runtime_override_pull_vcf_shard_mpd,
        runtime_override_merge_pesr_depth=runtime_override_merge_pesr_depth_mpd,
        runtime_override_sort_merged_vcf=runtime_override_sort_merged_vcf_mpd,
        runtime_override_subset_small=runtime_override_subset_small_mpd,
        runtime_override_subset_large=runtime_override_subset_large_mpd,
        runtime_override_make_sites_only=runtime_override_make_sites_only_mpd,
        runtime_override_concat_large_pesr_depth=runtime_override_concat_large_pesr_depth_mpd,
        runtime_override_concat_shards=runtime_override_concat_shards_mpd,
        runtime_override_preconcat_large_pesr_depth=runtime_override_preconcat_large_pesr_depth_mpd,
        runtime_override_hail_merge_large_pesr_depth=runtime_override_hail_merge_large_pesr_depth_mpd,
        runtime_override_fix_header_large_pesr_depth=runtime_override_fix_header_large_pesr_depth_mpd,
        runtime_override_preconcat_pesr_depth_shards=runtime_override_preconcat_pesr_depth_shards_mpd,
        runtime_override_hail_merge_pesr_depth_shards=runtime_override_hail_merge_pesr_depth_shards_mpd,
        runtime_override_fix_header_pesr_depth_shards=runtime_override_fix_header_pesr_depth_shards_mpd
    }

    #Merge PESR & RD VCFs
    if (use_hail) {
      call HailMerge.HailMerge as ConcatPesrDepthHail {
        input:
          vcfs=[MergeDeletions.out, MergeDuplications.out, HarmonizeHeaders.out[2], HarmonizeHeaders.out[3], HarmonizeHeaders.out[4]],
          prefix="~{cohort_name}.~{contig}.concat_pesr_depth",
          gcs_project=gcs_project,
          sv_base_mini_docker=sv_base_mini_docker,
          sv_pipeline_docker=sv_pipeline_docker,
          sv_pipeline_hail_docker=sv_pipeline_hail_docker,
          runtime_override_preconcat=runtime_override_preconcat_pesr_depth,
          runtime_override_hail_merge=runtime_override_hail_merge_pesr_depth,
          runtime_override_fix_header=runtime_override_fix_header_pesr_depth
      }
    }
    if (!use_hail) {
      call MiniTasks.ConcatVcfs as ConcatPesrDepth {
        input:
          vcfs=[MergeDeletions.out, MergeDuplications.out, HarmonizeHeaders.out[2], HarmonizeHeaders.out[3], HarmonizeHeaders.out[4]],
          vcfs_idx=[MergeDeletions.out+".tbi", MergeDuplications.out+".tbi", HarmonizeHeaders.out[2]+".tbi", HarmonizeHeaders.out[3]+".tbi", HarmonizeHeaders.out[4]+".tbi"],
          allow_overlaps=true,
          outfile_prefix="~{cohort_name}.~{contig}.concat_pesr_depth",
          sv_base_mini_docker=sv_base_mini_docker,
          runtime_attr_override=runtime_override_concat_large_pesr_depth
      }
    }

    #Update SR background fail & bothside pass files (2)
    call MiniTasks.UpdateSrList as UpdateBackgroundFailSecond {
      input:
        vcf=select_first([ConcatPesrDepth.concat_vcf, ConcatPesrDepthHail.merged_vcf]),
        original_list=UpdateBackgroundFailFirst.updated_list,
        outfile="~{cohort_name}.~{contig}.sr_background_fail.updated2.txt",
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_override_update_sr_list
    }
    call MiniTasks.UpdateSrList as UpdateBothsidePassSecond {
      input:
        vcf=select_first([ConcatPesrDepth.concat_vcf, ConcatPesrDepthHail.merged_vcf]),
        original_list=UpdateBothsidePassFirst.updated_list,
        outfile="~{cohort_name}.~{contig}.sr_bothside_pass.updated2.txt",
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_override_update_sr_list
    }

    File vcfs_out_ = select_first([ConcatPesrDepth.concat_vcf, ConcatPesrDepthHail.merged_vcf])
    File vcf_indexes_out_ = select_first([ConcatPesrDepth.concat_vcf_idx, ConcatPesrDepthHail.merged_vcf_index])
  }

  #Merge resolved vcfs for QC
  if (merge_vcfs) {
    call MiniTasks.ConcatVcfs {
      input:
        vcfs=vcfs_out_,
        vcfs_idx=vcf_indexes_out_,
        naive=true,
        outfile_prefix="~{cohort_name}.combine_batches",
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_attr_override=runtime_override_concat
    }
  }

  #Final outputs
  output {
    Array[File] combined_vcfs = vcfs_out_
    Array[File] combined_vcf_indexes = vcf_indexes_out_
    Array[File] cluster_bothside_pass_lists = UpdateBothsidePassSecond.updated_list
    Array[File] cluster_background_fail_lists = UpdateBackgroundFailSecond.updated_list
    File? combine_batches_merged_vcf = ConcatVcfs.concat_vcf
    File? combine_batches_merged_vcf_index = ConcatVcfs.concat_vcf_idx
  }
}
