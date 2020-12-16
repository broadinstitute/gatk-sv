version 1.0

# Author: Ryan Collins <rlcollins@g.harvard.edu>

import "Structs.wdl"
import "TasksMakeCohortVcf.wdl" as MiniTasks
import "HailMerge.wdl" as HailMerge

# Workflow to shard a filtered vcf & run vcfcluster (sub-sub-sub workflow)
workflow ShardedCluster {
  input {
    File vcf
    Int num_samples
    Int dist
    Float frac
    String prefix
    String contig
    String cohort_name
    String evidence_type
    String sv_type
    Float sample_overlap
    File? exclude_list
    File empty_file
    Int sv_size
    Array[String] sv_types
    Float merging_shard_scale_factor = 30000000

    File hail_script
    String project

    String sv_pipeline_docker
    String sv_base_mini_docker

    # overrides for local tasks
    RuntimeAttr? runtime_override_shard_clusters
    RuntimeAttr? runtime_override_shard_vids
    RuntimeAttr? runtime_override_pull_vcf_shard
    RuntimeAttr? runtime_override_svtk_vcf_cluster
    RuntimeAttr? runtime_override_get_vcf_header_with_members_info_line

    RuntimeAttr? runtime_override_preconcat_sharded_cluster
    RuntimeAttr? runtime_override_hail_merge_sharded_cluster
    RuntimeAttr? runtime_override_fix_header_sharded_cluster

    # overrides for merge subworkflow
    RuntimeAttr? runtime_override_merge_clusters
    RuntimeAttr? runtime_override_concat_inner_shards

    # overrides for MiniTasks
    RuntimeAttr? runtime_override_concat_sharded_cluster
    RuntimeAttr? runtime_override_sort_merged_vcf
    RuntimeAttr? runtime_override_count_samples
    RuntimeAttr? runtime_override_get_vids
    RuntimeAttr? runtime_override_cat_vid_lists_sharded
    RuntimeAttr? runtime_override_make_sites_only
  }


  File vcf_idx = vcf + ".tbi"
  if (defined(exclude_list)) {
    File exclude_list_idx = exclude_list + ".tbi"
  }

  call MiniTasks.MakeSitesOnlyVcf {
    input:
      vcf=vcf,
      vcf_index=vcf + ".tbi",
      prefix="~{prefix}.sites_only",
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_override_make_sites_only
  }

  Int merge_shard_size = ceil(merging_shard_scale_factor / num_samples)

  call ShardClusters {
    input:
      vcf=MakeSitesOnlyVcf.out,
      prefix="~{prefix}.sites_only.shard_clusters",
      dist=dist,
      frac=frac,
      exclude_list=exclude_list,
      exclude_list_idx=exclude_list_idx,
      svsize=sv_size,
      sv_types=sv_types,
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_override_shard_clusters
  }

  call MiniTasks.ShardVidsForClustering {
    input:
      clustered_vcf=ShardClusters.out,
      prefix="~{prefix}.sites_only.clustered",
      records_per_shard=merge_shard_size,
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_override_shard_vids
  }

  #Run vcfcluster per shard
  scatter (i in range(length(ShardVidsForClustering.out))) {
    call MiniTasks.PullVcfShard {
      input:
        vcf=vcf,
        vids=ShardVidsForClustering.out[i],
        prefix="~{prefix}.unclustered.shard_${i}",
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_attr_override=runtime_override_pull_vcf_shard
    }
    call SvtkVcfCluster {
      input:
        vcf=PullVcfShard.out,
        num_samples=num_samples,
        num_vids=PullVcfShard.count,
        prefix="~{prefix}.clustered.shard_${i}",
        vid_prefix="~{cohort_name}_~{contig}_~{evidence_type}_~{sv_type}_~{i}",
        dist=dist,
        frac=frac,
        exclude_list=exclude_list,
        exclude_list_idx=exclude_list_idx,
        svsize=sv_size,
        sample_overlap=sample_overlap,
        sv_types=sv_types,
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_override_svtk_vcf_cluster
    }
    call MiniTasks.SortVcf {
      input:
        vcf = SvtkVcfCluster.out,
        outfile_prefix = "~{prefix}.sorted.shard_${i}",
        sv_base_mini_docker = sv_base_mini_docker,
        runtime_attr_override = runtime_override_sort_merged_vcf
    }
  }

  if (length(SvtkVcfCluster.out) == 0) {
    call GetVcfHeaderWithMembersInfoLine {
      input:
        vcf_gz=vcf,
        prefix="~{prefix}.clustered",
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_attr_override=runtime_override_get_vcf_header_with_members_info_line
    }
  }
  if (length(SvtkVcfCluster.out) > 0) {
    call HailMerge.HailMerge as ConcatVcfs {
      input:
        vcfs=SortVcf.out,
        prefix="~{prefix}.clustered",
        hail_script=hail_script,
        project=project,
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_override_preconcat=runtime_override_preconcat_sharded_cluster,
        runtime_override_hail_merge=runtime_override_hail_merge_sharded_cluster,
        runtime_override_fix_header=runtime_override_fix_header_sharded_cluster
    }
  }

  #Output
  output {
    File clustered_vcf = select_first([GetVcfHeaderWithMembersInfoLine.out, ConcatVcfs.merged_vcf])
    File clustered_vcf_idx = select_first([GetVcfHeaderWithMembersInfoLine.out_idx, ConcatVcfs.merged_vcf_index])
  }
}

# Adds MEMBERS definition to header (workaround for when VIDs_list is empty)
task GetVcfHeaderWithMembersInfoLine {
  input {
    File vcf_gz
    String prefix
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr runtime_default = object {
    mem_gb: 1,
    disk_gb: 10,
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
    docker: sv_base_mini_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euxo pipefail
    zgrep "^##" ~{vcf_gz} > ~{prefix}.vcf
    echo "##INFO=<ID=MEMBERS,Number=.,Type=String,Description=\"IDs of cluster's constituent records.\">" >> ~{prefix}.vcf
    zgrep "^#CHROM" ~{vcf_gz} >> ~{prefix}.vcf
    bgzip ~{prefix}.vcf
    tabix ~{prefix}.vcf.gz
  >>>

  output {
    File out = "~{prefix}.vcf.gz"
    File out_idx = "~{prefix}.vcf.gz.tbi"
  }
}

#Do fast cluster on vcf (sample_overlap = 0) to generate shards
task ShardClusters {
  input {
    File vcf
    String prefix
    Int dist
    Float frac
    File? exclude_list
    File? exclude_list_idx
    Int svsize
    Array[String] sv_types
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size(vcf, "GiB")
  Float base_disk_gb = 10.0
  Float input_disk_scale = 1.0
  RuntimeAttr runtime_default = object {
    mem_gb: 2.0,
    disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
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
    ~{if defined(exclude_list) && !defined(exclude_list_idx) then "tabix -p bed ~{exclude_list}" else ""}

    #Run clustering
    svtk vcfcluster <(echo "~{vcf}") ~{prefix}.vcf.gz \
      -d ~{dist} \
      -f ~{frac} \
      ~{if defined(exclude_list) then "-x ~{exclude_list}" else ""} \
      -z ~{svsize} \
      -p ~{prefix} \
      -t ~{sep=',' sv_types} \
      -o 0 \
      --preserve-header \
      --preserve-ids \
      --skip-merge
  >>>

  output {
    File out = "~{prefix}.vcf.gz"
  }
}

task SvtkVcfCluster {
  input {
    File vcf
    String prefix
    String vid_prefix
    Int num_vids
    Int num_samples
    Int dist
    Float frac
    Float sample_overlap
    File? exclude_list
    File? exclude_list_idx
    Int svsize
    Array[String] sv_types
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float default_mem_gb = 3.75 + (120.0 * (num_vids / 19000.0) * (num_samples / 140000.0))
  RuntimeAttr runtime_default = object {
                                  mem_gb: default_mem_gb,
                                  disk_gb: ceil(10.0 + size(vcf, "GiB") * 2.0),
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
    ~{if defined(exclude_list) && !defined(exclude_list_idx) then "tabix -p bed ~{exclude_list}" else ""}

    #Run clustering
    svtk vcfcluster <(echo "~{vcf}") - \
        -d ~{dist} \
        -f ~{frac} \
        ~{if defined(exclude_list) then "-x ~{exclude_list}" else ""} \
        -z ~{svsize} \
        -p ~{vid_prefix} \
        -t ~{sep=',' sv_types} \
        -o ~{sample_overlap} \
        --preserve-ids \
        --preserve-genotypes \
        --preserve-header \
      | gzip > ~{prefix}.vcf.gz
  >>>

  output {
    File out = "~{prefix}.vcf.gz"
  }
}