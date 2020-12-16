version 1.0

# Author: Ryan Collins <rlcollins@g.harvard.edu>

import "Tasks0506.wdl" as MiniTasks
import "ShardedCluster.wdl" as ShardedCluster

# Workflow to perform sharding & clustering of a vcf for a single chromosome
workflow ClusterSingleChrom {
  input {
    File vcf
    File vcf_index
    String contig
    String prefix
    Int max_shards
    Int min_per_shard
    Int dist
    Float frac
    Float sample_overlap
    File? exclude_list
    Int sv_size
    Array[String] sv_types

    String sv_pipeline_docker
    String sv_base_mini_docker

    # overrides for local tasks
    RuntimeAttr? runtime_override_concat_sv_types

    # overrides for MiniTasks
    RuntimeAttr? runtime_override_subset_sv_type

    # overrides for ShardedCluster
    RuntimeAttr? runtime_override_shard_vcf_precluster
    RuntimeAttr? runtime_override_svtk_vcf_cluster
    RuntimeAttr? runtime_override_get_vcf_header_with_members_info_line
    RuntimeAttr? runtime_override_concat_shards
  }

  String contig_prefix = prefix + "." + contig

  #Scatter over svtypes
  scatter ( sv_type in sv_types ) {
    #Subset vcf to only contain records for that svtype

    call MiniTasks.FilterVcf as SubsetSvType {
      input:
        vcf=vcf,
        vcf_index=vcf_index,
        records_filter='INFO/SVTYPE="~{sv_type}"',
        outfile_prefix=contig_prefix + ".~{sv_type}",
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_attr_override=runtime_override_subset_sv_type
    }

    #For each svtype, intelligently shard VCF for clustering
    call ShardedCluster.ShardedCluster as ShardedCluster {
      input:
        vcf=SubsetSvType.filtered_vcf,
        dist=dist,
        frac=frac,
        max_shards=max_shards,
        min_per_shard=min_per_shard,
        prefix=prefix,
        contig=contig,
        sv_type=sv_type,
        sample_overlap=sample_overlap,
        exclude_list=exclude_list,
        sv_size=sv_size,
        sv_types=sv_types,
        sv_pipeline_docker=sv_pipeline_docker,
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_override_shard_vcf_precluster=runtime_override_shard_vcf_precluster,
        runtime_override_svtk_vcf_cluster=runtime_override_svtk_vcf_cluster,
        runtime_override_get_vcf_header_with_members_info_line=runtime_override_get_vcf_header_with_members_info_line,
        runtime_override_concat_shards=runtime_override_concat_shards
    }
  }

  #Merge svtypes
  call ConcatAndRenameVcfs as ConcatSvTypes {
    input:
      vcfs=ShardedCluster.clustered_vcf,
      prefix=prefix,
      contig=contig,
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_override_concat_sv_types
  }

  #Output clustered vcf
  output {
    File clustered_vcf = ConcatSvTypes.concat_vcf
    File clustered_vcf_idx = ConcatSvTypes.concat_vcf_idx
  }
}


#Merge multiple vcfs
task ConcatAndRenameVcfs {
  input {
    Array[File] vcfs
    String prefix
    String contig

    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  String raw_vcf_name = prefix + "." + contig + ".raw.vcf.gz"
  String vcf_name = prefix + "." + contig + ".vcf.gz"

  # when filtering/sorting/etc, memory usage will likely go up (much of the data will have to
  # be held in memory or disk while working, potentially in a form that takes up more space)
  Float input_size = size(vcfs, "GiB")
  RuntimeAttr runtime_default = object {
    mem_gb: 2.0 + 5.0 * input_size,
    disk_gb: ceil(10.0 + 40.0 * input_size),
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
    set -eu -o pipefail
    
    vcf-concat -f ~{write_lines(vcfs)} \
      | vcf-sort -c \
      | bgzip -c \
      > ~{raw_vcf_name}

    /opt/sv-pipeline/04_variant_resolution/scripts/rename_after_vcfcluster.py \
      --chrom ~{contig} \
      --prefix ~{prefix} \
      ~{raw_vcf_name} - \
      | bgzip -c \
      > ~{vcf_name}

    tabix -p vcf -f ~{vcf_name}
  >>>

  output {
    File concat_vcf = vcf_name
    File concat_vcf_idx = vcf_name + ".tbi"
  }
}
