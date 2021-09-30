version 1.0

# Author: Ryan Collins <rlcollins@g.harvard.edu>

import "TasksMakeCohortVcf.wdl" as MiniTasks
import "ShardedCluster.wdl" as ShardedCluster

# Workflow to perform sharding & clustering of a vcf for a single chromosome
workflow ClusterSingleChrom {
  input {
    File vcf
    File vcf_index
    String contig
    String prefix
    Int dist
    Float frac
    Float sample_overlap
    File? exclude_list
    Int sv_size
    Array[String] sv_types

    String sv_pipeline_docker
    String sv_base_mini_docker

    # overrides for local tasks
    RuntimeAttr? runtime_override_concat_svtypes

    # overrides for MiniTasks
    RuntimeAttr? runtime_override_subset_sv_type

    # overrides for ShardedCluster
    RuntimeAttr? runtime_override_shard_vcf_precluster
    RuntimeAttr? runtime_override_pull_vcf_shard
    RuntimeAttr? runtime_override_svtk_vcf_cluster
    RuntimeAttr? runtime_override_get_vcf_header_with_members_info_line
    RuntimeAttr? runtime_override_concat_svtypes
    RuntimeAttr? runtime_override_concat_sharded_cluster
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
    call ShardedCluster.ShardedCluster as ShardedCluster {
      input:
        vcf=SubsetSvType.filtered_vcf,
        dist=dist,
        frac=frac,
        prefix="~{prefix}.~{sv_type}",
        contig=contig,
        sv_type=sv_type,
        sample_overlap=sample_overlap,
        exclude_list=exclude_list,
        sv_size=sv_size,
        sv_types=sv_types,
        sv_pipeline_docker=sv_pipeline_docker,
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_override_shard_vcf_precluster=runtime_override_shard_vcf_precluster,
        runtime_override_pull_vcf_shard=runtime_override_pull_vcf_shard,
        runtime_override_svtk_vcf_cluster=runtime_override_svtk_vcf_cluster,
        runtime_override_get_vcf_header_with_members_info_line=runtime_override_get_vcf_header_with_members_info_line,
        runtime_override_concat_sharded_cluster=runtime_override_concat_sharded_cluster
    }
    call RenameVariants {
      input:
        vcf=ShardedCluster.clustered_vcf,
        vcf_index=ShardedCluster.clustered_vcf_idx,
        prefix=prefix,
        contig=contig,
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_override_concat_svtypes
    }
  }

  #Merge svtypes
  call MiniTasks.ConcatVcfs as ConcatSvTypes {
    input:
      vcfs=RenameVariants.out,
      vcfs_idx=RenameVariants.out_index,
      allow_overlaps=true,
      outfile_prefix="~{prefix}.~{contig}.precluster_concat",
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_override_concat_svtypes
  }

  #Output clustered vcf
  output {
    File clustered_vcf = ConcatSvTypes.concat_vcf
    File clustered_vcf_idx = ConcatSvTypes.concat_vcf_idx
  }
}

task RenameVariants {
  input {
    File vcf
    File vcf_index
    String prefix
    String contig

    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  String vcf_name = prefix + "." + contig + ".renamed.vcf.gz"

  # when filtering/sorting/etc, memory usage will likely go up (much of the data will have to
  # be held in memory or disk while working, potentially in a form that takes up more space)
  Float input_size = size(vcf, "GiB")
  RuntimeAttr runtime_default = object {
    mem_gb: 3.75,
    disk_gb: ceil(10.0 + 2.0 * input_size),
    cpu_cores: 1,
    preemptible_tries: 3,
    max_retries: 1,
    boot_disk_gb: 10
  }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GiB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euo pipefail

    /opt/sv-pipeline/04_variant_resolution/scripts/rename_after_vcfcluster.py \
      --chrom ~{contig} \
      --prefix ~{prefix} \
      ~{vcf} - \
      | bgzip -c \
      > ~{vcf_name}

    tabix -p vcf -f ~{vcf_name}
  >>>

  output {
    File out = vcf_name
    File out_index = vcf_name + ".tbi"
  }
}
