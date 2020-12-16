version 1.0

# Author: Ryan Collins <rlcollins@g.harvard.edu>

import "Structs.wdl"
import "Tasks0506.wdl" as MiniTasks
import "Utils.wdl" as utils

# Workflow to shard a filtered vcf & run vcfcluster (sub-sub-sub workflow)
workflow ShardedCluster {
  input {
    File vcf
    Int dist
    Float frac
    Int max_shards
    Int min_per_shard
    String prefix
    String contig
    String sv_type
    Float sample_overlap
    File? exclude_list
    Int sv_size
    Array[String] sv_types
    Float merging_shard_scale_factor = 30000000

    String sv_pipeline_docker
    String sv_base_mini_docker

    # Do not use
    File? NONE_FILE_

    # overrides for local tasks
    RuntimeAttr? runtime_override_shard_vcf_precluster
    RuntimeAttr? runtime_override_svtk_vcf_cluster
    RuntimeAttr? runtime_override_get_vcf_header_with_members_info_line

    # overrides for merge subworkflow
    RuntimeAttr? runtime_override_merge_clusters
    RuntimeAttr? runtime_override_concat_inner_shards

    # overrides for MiniTasks
    RuntimeAttr? runtime_override_concat_shards
    RuntimeAttr? runtime_override_sort_merged_vcf
    RuntimeAttr? runtime_override_count_samples
  }

  File vcf_idx = vcf + ".tbi"
  if (defined(exclude_list)) {
    File exclude_list_idx = exclude_list + ".tbi"
  }

  call utils.CountSamples {
    input:
    vcf=vcf,
    sv_base_mini_docker=sv_base_mini_docker,
    runtime_attr_override=runtime_override_count_samples
  }
  Int merge_shard_size = ceil(merging_shard_scale_factor / CountSamples.num_samples)

  call ShardClusters {
    input:
      vcf=vcf,
      prefix=prefix,
      dist=dist,
      frac=frac,
      exclude_list=exclude_list,
      exclude_list_idx=exclude_list_idx,
      svsize=sv_size,
      records_per_shard=merge_shard_size,
      sv_types=sv_types,
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_override_svtk_vcf_cluster
  }

  #Run vcfcluster per shard
  scatter (i in range(length(ShardClusters.out))) {
    call CountLines {
      input:
        file=ShardClusters.out[i],
        sv_base_mini_docker=sv_base_mini_docker
    }
    call SvtkVcfCluster {
      input:
        vcf=vcf,
        vids=ShardClusters.out[i],
        num_samples=CountSamples.num_samples,
        num_vids=CountLines.out,
        prefix="~{prefix}.~{contig}.~{sv_type}.shard_${i}.clustered",
        vid_prefix="~{prefix}_~{contig}_~{sv_type}_~{i}",
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
        outfile_prefix = "~{prefix}.~{contig}.~{sv_type}.shard_${i}.sorted",
        sv_base_mini_docker = sv_base_mini_docker,
        runtime_attr_override = runtime_override_sort_merged_vcf
    }
  }

  if (length(SvtkVcfCluster.out) == 0) {
    call GetVcfHeaderWithMembersInfoLine {
      input:
        vcf_gz=vcf,
        prefix="~{prefix}.~{contig}.~{sv_type}",
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_attr_override=runtime_override_get_vcf_header_with_members_info_line
    }
  }
  if (length(SvtkVcfCluster.out) > 0) {
    call MiniTasks.ConcatVcfs {
      input:
        vcfs=SortVcf.out,
        vcfs_idx=SortVcf.out_index,
        merge_sort=true,
        outfile_prefix="~{prefix}.~{contig}.~{sv_type}.clustered",
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_attr_override=runtime_override_concat_shards
    }
  }

  #Output
  output {
    File clustered_vcf = select_first([GetVcfHeaderWithMembersInfoLine.out, ConcatVcfs.concat_vcf])
    File clustered_vcf_idx = select_first([GetVcfHeaderWithMembersInfoLine.out_idx, ConcatVcfs.concat_vcf_idx])
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
    set -euo pipefail
    gunzip -c ~{vcf_gz} | grep "^##" > header
    echo "##INFO=<ID=MEMBERS,Number=.,Type=String,Description=\"IDs of cluster's constituent records.\">" >> header
    gunzip -c ~{vcf_gz} | grep "^#" | grep -v "^##" >> header
    bgzip -c header > ~{prefix}.members.vcf.gz
    tabix ~{prefix}.members.vcf.gz
  >>>

  output {
    File out = "~{prefix}.members.vcf.gz"
    File out_idx = "~{prefix}.members.vcf.gz.tbi"
  }
}

#Do fast cluster on sites-only vcf (sample_overlap = 0) to generate shards
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
    Int records_per_shard
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
    bcftools view -G ~{vcf} -o sites_only.vcf
    ~{if defined(exclude_list) && !defined(exclude_list_idx) then "tabix -p bed ~{exclude_list}" else ""}
    #Run clustering
    svtk vcfcluster <(echo "sites_only.vcf") unmerged_clusters.vcf \
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

    # Shard output
    python3 <<CODE
    import sys
    import os
    import pysam

    vcf = pysam.VariantFile("unmerged_clusters.vcf")

    current_cluster = None
    current_cluster_vids = []
    current_shard = 0
    current_shard_size = 0
    shard_path_format = "~{prefix}.vids.shard_{}.list"
    shard_path = shard_path_format.format(current_shard)
    fout = open(shard_path, 'w')
    if fout is None:
      raise IOError("Could not open '{}'".format(shard_path))
      sys.exit(1)

    for record in vcf.fetch():
      cluster_id = record.info['CLUSTER']
      if cluster_id == current_cluster:
        current_cluster_vids.append(record.id)
      else:
        for vid in current_cluster_vids:
          fout.write(vid + '\n')
        current_shard_size += len(current_cluster_vids)
        if current_shard_size >= ~{records_per_shard}:
          current_shard += 1
          current_shard_size = 0
          fout.close()
          shard_path = shard_path_format.format(current_shard)
          fout = open(shard_path, 'w')
          if fout is None:
            raise IOError("Could not open '{}'".format(shard_path))
            sys.exit(1)
        current_cluster_vids = [record.id]
        current_cluster = cluster_id

    # Write last cluster
    for vid in current_cluster_vids:
      fout.write(vid + '\n')
    current_shard_size += len(current_cluster_vids)
    fout.close()

    # Delete trailing empty shard
    if current_shard > 0 and current_shard_size == 0:
      os.remove(shard_path)
    CODE
  >>>

  output {
    Array[File] out = glob("~{prefix}.vids.shard_*.list")
  }
}

task SvtkVcfCluster {
  input {
    File vcf
    File vids
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
  String output_prefix = "~{prefix}"
  RuntimeAttr runtime_default = object {
                                  mem_gb: default_mem_gb,
                                  disk_gb: ceil(10.0 + size(vcf, "GiB") * 10.0),
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
    bcftools view --no-version --include ID=@~{vids} ~{vcf} -O z -o unclustered.vcf.gz
    ~{if defined(exclude_list) && !defined(exclude_list_idx) then "tabix -p bed ~{exclude_list}" else ""}
    #Run clustering
    svtk vcfcluster <(echo "unclustered.vcf.gz") ~{output_prefix}.vcf \
      -d ~{dist} \
      -f ~{frac} \
      ~{if defined(exclude_list) then "-x ~{exclude_list}" else ""} \
      -z ~{svsize} \
      -p ~{vid_prefix} \
      -t ~{sep=',' sv_types} \
      -o ~{sample_overlap} \
      --preserve-ids \
      --preserve-genotypes \
      --preserve-header
    gzip ~{output_prefix}.vcf
  >>>

  output {
    File out = "~{output_prefix}.vcf.gz"
  }
}

task CountLines {
  input {
    File file
    String sv_base_mini_docker
  }
  runtime {
    memory: "0.9 GiB"
    disks: "local-disk 10 HDD"
    cpu: 1
    preemptible: 3
    maxRetries: 1
    docker: sv_base_mini_docker
    bootDiskSizeGb: 10
  }

  command <<<
    set -euo pipefail
    wc -l < ~{file} > count.txt
  >>>

  output {
    Int out = read_int("count.txt")
  }
}