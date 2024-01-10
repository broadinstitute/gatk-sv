version 1.0

# Author: Ryan Collins <rlcollins@g.harvard.edu>

import "TasksMakeCohortVcf.wdl" as MiniTasks
import "ShardedCluster.wdl" as ShardedCluster
import "HailMerge.wdl" as HailMerge

# Workflow to perform sharding & clustering of a vcf for a single chromosome
workflow ClusterSingleChrom {
  input {
    File vcf
    File vcf_index
    Int num_samples
    String contig
    String cohort_name
    String evidence_type
    String prefix
    Int dist
    Float frac
    Float sample_overlap
    File? exclude_list
    Int sv_size
    Array[String] sv_types
    File empty_file

    Boolean use_hail
    String? gcs_project

    String sv_pipeline_docker
    String sv_pipeline_hail_docker
    String sv_base_mini_docker

    # overrides for MiniTasks
    RuntimeAttr? runtime_override_subset_sv_type
    RuntimeAttr? runtime_override_cat_vid_lists_chrom

    # overrides for ShardedCluster
    RuntimeAttr? runtime_override_shard_clusters
    RuntimeAttr? runtime_override_shard_vids
    RuntimeAttr? runtime_override_pull_vcf_shard
    RuntimeAttr? runtime_override_svtk_vcf_cluster
    RuntimeAttr? runtime_override_get_vcf_header_with_members_info_line
    RuntimeAttr? runtime_override_concat_svtypes
    RuntimeAttr? runtime_override_concat_sharded_cluster
    RuntimeAttr? runtime_override_cat_vid_lists_sharded
    RuntimeAttr? runtime_override_make_sites_only
    RuntimeAttr? runtime_override_sort_merged_vcf

    RuntimeAttr? runtime_override_preconcat_sharded_cluster
    RuntimeAttr? runtime_override_hail_merge_sharded_cluster
    RuntimeAttr? runtime_override_fix_header_sharded_cluster
  }

  #Scatter over svtypes
  scatter ( sv_type in sv_types ) {
    #Subset vcf to only contain records for that svtype

    call MiniTasks.FilterVcf as SubsetSvType {
      input:
        vcf=vcf,
        vcf_index=vcf_index,
        records_filter='INFO/SVTYPE="~{sv_type}"',
        outfile_prefix="~{prefix}.~{sv_type}",
        use_ssd=true,
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_attr_override=runtime_override_subset_sv_type
    }

    #For each svtype, intelligently shard VCF for clustering
    call ShardedCluster.ShardedCluster {
      input:
        vcf=SubsetSvType.filtered_vcf,
        num_samples=num_samples,
        dist=dist,
        frac=frac,
        prefix="~{prefix}.~{sv_type}",
        cohort_name=cohort_name,
        contig=contig,
        evidence_type=evidence_type,
        sv_type=sv_type,
        sample_overlap=sample_overlap,
        exclude_list=exclude_list,
        sv_size=sv_size,
        sv_types=sv_types,
        empty_file=empty_file,
        use_hail=use_hail,
        gcs_project=gcs_project,
        sv_pipeline_docker=sv_pipeline_docker,
        sv_pipeline_hail_docker=sv_pipeline_hail_docker,
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_override_shard_clusters=runtime_override_shard_clusters,
        runtime_override_shard_vids=runtime_override_shard_vids,
        runtime_override_pull_vcf_shard=runtime_override_pull_vcf_shard,
        runtime_override_svtk_vcf_cluster=runtime_override_svtk_vcf_cluster,
        runtime_override_get_vcf_header_with_members_info_line=runtime_override_get_vcf_header_with_members_info_line,
        runtime_override_concat_sharded_cluster=runtime_override_concat_sharded_cluster,
        runtime_override_cat_vid_lists_sharded=runtime_override_cat_vid_lists_sharded,
        runtime_override_make_sites_only=runtime_override_make_sites_only,
        runtime_override_sort_merged_vcf=runtime_override_sort_merged_vcf,
        runtime_override_preconcat_sharded_cluster=runtime_override_preconcat_sharded_cluster,
        runtime_override_hail_merge_sharded_cluster=runtime_override_hail_merge_sharded_cluster,
        runtime_override_fix_header_sharded_cluster=runtime_override_fix_header_sharded_cluster
    }
  }

  #Output clustered vcf
  output {
    Array[File] clustered_vcfs = ShardedCluster.clustered_vcf
    Array[File] clustered_vcf_indexes = ShardedCluster.clustered_vcf_idx
  }
}
