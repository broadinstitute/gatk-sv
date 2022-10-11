version 1.0

import "Structs.wdl"
import "TasksMakeCohortVcf.wdl" as MiniTasks
import "HailMerge.wdl" as HailMerge
import "ShardedCluster.wdl" as ShardedCluster
import "Utils.wdl" as utils

workflow MergePesrDepth {
    input {
        File subtyped_pesr_vcf
        File subtyped_depth_vcf
        Int num_samples

        String prefix
        String cohort_name
        String svtype
        String contig
        Float merging_shard_scale_factor = 30000000

        Boolean use_hail = false
        String? gcs_project

        String sv_pipeline_docker
        String sv_pipeline_hail_docker
        String sv_base_mini_docker

        # overrides for local tasks
        RuntimeAttr? runtime_override_shard_clusters
        RuntimeAttr? runtime_override_shard_vids
        RuntimeAttr? runtime_override_pull_vcf_shard
        RuntimeAttr? runtime_override_merge_pesr_depth

        # overrides for MiniTasks
        RuntimeAttr? runtime_override_sort_merged_vcf
        RuntimeAttr? runtime_override_subset_small
        RuntimeAttr? runtime_override_subset_large
        RuntimeAttr? runtime_override_make_sites_only
        RuntimeAttr? runtime_override_concat_large_pesr_depth
        RuntimeAttr? runtime_override_concat_shards

        RuntimeAttr? runtime_override_preconcat_large_pesr_depth
        RuntimeAttr? runtime_override_hail_merge_large_pesr_depth
        RuntimeAttr? runtime_override_fix_header_large_pesr_depth

        RuntimeAttr? runtime_override_preconcat_pesr_depth_shards
        RuntimeAttr? runtime_override_hail_merge_pesr_depth_shards
        RuntimeAttr? runtime_override_fix_header_pesr_depth_shards
    }

    # Pull out CNVs too small to cluster (less than reciprocal_overlap_fraction * min_depth_only_length)
    call MiniTasks.FilterVcf as SubsetSmall {
        input:
            vcf=subtyped_pesr_vcf,
            vcf_index=subtyped_pesr_vcf + ".tbi",
            outfile_prefix="~{prefix}.subset_small",
            records_filter='INFO/SVLEN<2500',
            use_ssd=true,
            sv_base_mini_docker=sv_base_mini_docker,
            runtime_attr_override=runtime_override_subset_small
    }

    call MiniTasks.FilterVcf as SubsetLarge {
        input:
            vcf=subtyped_pesr_vcf,
            vcf_index=subtyped_pesr_vcf + ".tbi",
            outfile_prefix="~{prefix}.subset_large",
            records_filter='INFO/SVLEN>=2500',
            use_ssd=true,
            sv_base_mini_docker=sv_base_mini_docker,
            runtime_attr_override=runtime_override_subset_large
    }

    if (use_hail) {
        call HailMerge.HailMerge as ConcatLargePesrDepthHail {
            input:
                vcfs=[SubsetLarge.filtered_vcf, subtyped_depth_vcf],
                prefix="~{prefix}.large_pesr_depth",
                gcs_project=gcs_project,
                sv_base_mini_docker=sv_base_mini_docker,
                sv_pipeline_docker=sv_pipeline_docker,
                sv_pipeline_hail_docker=sv_pipeline_hail_docker,
                runtime_override_preconcat=runtime_override_preconcat_large_pesr_depth,
                runtime_override_hail_merge=runtime_override_hail_merge_large_pesr_depth,
                runtime_override_fix_header=runtime_override_fix_header_large_pesr_depth
        }
    }
    if (!use_hail) {
        call MiniTasks.ConcatVcfs as ConcatLargePesrDepth {
            input:
                vcfs=[SubsetLarge.filtered_vcf, subtyped_depth_vcf],
                vcfs_idx=[SubsetLarge.filtered_vcf + ".tbi", subtyped_depth_vcf + ".tbi"],
                allow_overlaps=true,
                outfile_prefix="~{prefix}.large_pesr_depth",
                sv_base_mini_docker=sv_base_mini_docker,
                runtime_attr_override=runtime_override_concat_large_pesr_depth
        }
    }

    call MiniTasks.MakeSitesOnlyVcf {
        input:
            vcf=select_first([ConcatLargePesrDepth.concat_vcf, ConcatLargePesrDepthHail.merged_vcf]),
            vcf_index=select_first([ConcatLargePesrDepth.concat_vcf_idx, ConcatLargePesrDepthHail.merged_vcf_index]),
            prefix="~{prefix}.large_pesr_depth.sites_only",
            sv_base_mini_docker=sv_base_mini_docker,
            runtime_attr_override=runtime_override_make_sites_only
    }

    # Fast cluster without sample overlap linkage for sharding
    Int merge_shard_size = ceil(merging_shard_scale_factor / num_samples)
    call ShardedCluster.ShardClusters {
        input:
            vcf=MakeSitesOnlyVcf.out,
            prefix="~{prefix}.shard_clusters",
            dist=1000000000,
            frac=0.5,
            svsize=0,
            sv_types=[svtype],
            sv_pipeline_docker=sv_pipeline_docker,
            runtime_attr_override=runtime_override_shard_clusters
    }

    call MiniTasks.ShardVidsForClustering {
        input:
            clustered_vcf=ShardClusters.out,
            prefix=prefix,
            records_per_shard=merge_shard_size,
            sv_pipeline_docker=sv_pipeline_docker,
            runtime_attr_override=runtime_override_shard_vids
    }

    scatter (i in range(length(ShardVidsForClustering.out))) {
        call MiniTasks.PullVcfShard {
            input:
                vcf=select_first([ConcatLargePesrDepth.concat_vcf, ConcatLargePesrDepthHail.merged_vcf]),
                vids=ShardVidsForClustering.out[i],
                prefix="~{prefix}.unclustered.shard_${i}",
                sv_base_mini_docker=sv_base_mini_docker,
                runtime_attr_override=runtime_override_pull_vcf_shard
        }
        call MergePesrDepthShard {
            input:
                vcf=PullVcfShard.out,
                vcf_index=PullVcfShard.out_index,
                prefix="~{prefix}.merge_pesr_depth.shard_~{i}",
                vid_prefix="~{cohort_name}_~{contig}_mpd~{i}",
                sv_pipeline_docker=sv_pipeline_docker,
                runtime_attr_override=runtime_override_merge_pesr_depth
        }
        call MiniTasks.SortVcf {
            input:
                vcf = MergePesrDepthShard.out,
                outfile_prefix = "~{prefix}.sorted.shard_${i}",
                sv_base_mini_docker = sv_base_mini_docker,
                runtime_attr_override = runtime_override_sort_merged_vcf
        }
    }

    if (use_hail) {
        call HailMerge.HailMerge as ConcatShardsHail {
            input:
                vcfs=flatten([[SubsetSmall.filtered_vcf], SortVcf.out]),
                prefix="~{prefix}.concat_shards",
                gcs_project=gcs_project,
                sv_base_mini_docker=sv_base_mini_docker,
                sv_pipeline_docker=sv_pipeline_docker,
                sv_pipeline_hail_docker=sv_pipeline_hail_docker,
                runtime_override_preconcat=runtime_override_preconcat_pesr_depth_shards,
                runtime_override_hail_merge=runtime_override_hail_merge_pesr_depth_shards,
                runtime_override_fix_header=runtime_override_fix_header_pesr_depth_shards
        }
    }
    if (!use_hail) {
        call MiniTasks.ConcatVcfs as ConcatShards {
            input:
                vcfs=flatten([[SubsetSmall.filtered_vcf], SortVcf.out]),
                vcfs_idx=flatten([[SubsetSmall.filtered_vcf_idx], SortVcf.out_index]),
                allow_overlaps=true,
                outfile_prefix="~{prefix}.concat_shards",
                sv_base_mini_docker=sv_base_mini_docker,
                runtime_attr_override=runtime_override_concat_shards
        }
    }

    output {
        File out = select_first([ConcatShards.concat_vcf, ConcatShardsHail.merged_vcf])
        File out_index = select_first([ConcatShards.concat_vcf_idx, ConcatShardsHail.merged_vcf_index])
    }
}


task MergePesrDepthShard {
    input {
        File vcf
        File vcf_index
        String prefix
        String vid_prefix
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    String output_file = prefix + ".vcf.gz"

    # when filtering/sorting/etc, memory usage will likely go up (much of the data will have to
    # be held in memory or disk while working, potentially in a form that takes up more space)
    Float input_size = size(vcf, "GiB")
    RuntimeAttr runtime_default = object {
                                      mem_gb: 2.0 + 0.6 * input_size,
                                      disk_gb: ceil(10.0 + 6 * input_size),
                                      cpu_cores: 1,
                                      preemptible_tries: 3,
                                      max_retries: 1,
                                      boot_disk_gb: 10
                                  }
    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GiB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} SSD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: sv_pipeline_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command <<<
        set -euo pipefail
        /opt/sv-pipeline/04_variant_resolution/scripts/merge_pesr_depth.py \
            --prefix ~{vid_prefix} \
            ~{vcf} \
            ~{output_file}
    >>>

    output {
        File out = output_file
    }
}
