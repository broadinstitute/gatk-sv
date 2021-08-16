version 1.0

# Author: Ryan Collins <rlcollins@g.harvard.edu>

# Workflow to gather SV VCF summary stats for an input VCF

import "TasksMakeCohortVcf.wdl" as MiniTasks

workflow ShardedQcCollection {
  input {
    File vcf
    String contig
    Int sv_per_shard
    String prefix

    File? vcf_idx

    String sv_base_mini_docker
    String sv_pipeline_docker

    # overrides for local tasks
    RuntimeAttr? runtime_override_collect_sharded_vcf_stats
    RuntimeAttr? runtime_override_svtk_vcf_2_bed

    # overrides for MiniTasks
    RuntimeAttr? runtime_override_split_vcf_to_qc
    RuntimeAttr? runtime_override_merge_subvcf_stat_shards
    RuntimeAttr? runtime_override_merge_svtk_vcf_2_bed
  }

  # Tabix to chromosome of interest, and shard input VCF for stats collection
  call MiniTasks.SplitVcf as SplitVcfToQc {
    input:
      vcf=vcf,
      vcf_idx=vcf_idx,
      contig=contig,
      min_vars_per_shard=sv_per_shard,
      prefix="vcf.shard.",
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_override_split_vcf_to_qc
  }

  # Scatter over VCF shards
  scatter (shard in SplitVcfToQc.vcf_shards) {
    # Collect VCF-wide summary stats
    call CollectShardedVcfStats {
      input:
        vcf=shard,
        prefix=prefix + ".shard_",
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_override_collect_sharded_vcf_stats
      }
    # Run vcf2bed_subworkflow for record purposes
    call SvtkVcf2bed {
      input:
        vcf=shard,
        prefix="~{prefix}.shard",
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_override_svtk_vcf_2_bed
      }
    }

  # Merge shards into single VCF stats file
  call MiniTasks.ConcatBeds as MergeSubvcfStatShards {
    input:
      shard_bed_files=CollectShardedVcfStats.vcf_stats,
      prefix=prefix + ".VCF_sites.stats",
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_override_merge_subvcf_stat_shards
  }

  # Merge vcf2bed_subworkflow output
  call MiniTasks.ConcatBeds as MergeSvtkVcf2bed {
    input:
      shard_bed_files=SvtkVcf2bed.vcf2bed_subworkflow_out,
      prefix=prefix + ".vcf2bed_subworkflow",
      index_output=false,
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_override_merge_svtk_vcf_2_bed
  }

  # Final output
  output {
    File vcf_stats=MergeSubvcfStatShards.merged_bed_file
    File vcf_stats_idx=MergeSubvcfStatShards.merged_bed_idx
    File samples_list=CollectShardedVcfStats.samples_list[0]
    File vcf2bed_out=MergeSvtkVcf2bed.merged_bed_file
  }
}

# Task to collect VCF-wide QC stats
task CollectShardedVcfStats {
  input {
    File vcf
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  # when filtering/sorting/etc, memory usage will likely go up (much of the data will have to
  # be held in memory or disk while working, potentially in a form that takes up more space)
  Float input_size = size(vcf, "GiB")
  Float compression_factor = 5.0
  Float base_disk_gb = 5.0
  Float base_mem_gb = 2.0
  RuntimeAttr runtime_default = object {
    mem_gb: base_mem_gb + compression_factor * input_size,
    disk_gb: ceil(base_disk_gb + input_size * (2.0 + 2.0 * compression_factor)),
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
  
    # Run QC script
    /opt/sv-pipeline/scripts/vcf_qc/collectQC.vcf_wide.sh \
      ~{vcf} \
      /opt/sv-pipeline/scripts/vcf_qc/SV_colors.txt \
      collectQC_vcfwide_output/
    
    # Prep outputs
    cp collectQC_vcfwide_output/data/VCF_sites.stats.bed.gz \
      ~{prefix}.VCF_sites.stats.bed.gz
    cp collectQC_vcfwide_output/data/VCF_sites.stats.bed.gz.tbi \
      ~{prefix}.VCF_sites.stats.bed.gz.tbi
    cp collectQC_vcfwide_output/analysis_samples.list \
      ~{prefix}.analysis_samples.list
    tar -czvf ~{prefix}.collectQC_vcfwide_output.tar.gz \
      collectQC_vcfwide_output
  >>>

  output {
    File vcf_stats = "~{prefix}.VCF_sites.stats.bed.gz"
    File vcf_stats_idx = "~{prefix}.VCF_sites.stats.bed.gz.tbi"
    File samples_list = "~{prefix}.analysis_samples.list"
    File vcfwide_tarball = "~{prefix}.collectQC_vcfwide_output.tar.gz"
  }
}


# Run vcf2bed_subworkflow on an input vcf
task SvtkVcf2bed {
  input {
    File vcf
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  String output_file = prefix + ".vcf2bed_subworkflow.bed.gz"
  
  # simple record-by-record processing, overhead should be O(1), with disk space usage increased because the operation
  # is copying input into new format
  Float input_size = size(vcf, "GiB")
  Float base_disk_gb = 5.0
  Float base_mem_gb = 2.0
  RuntimeAttr runtime_default = object {
    mem_gb: base_mem_gb,
    disk_gb: ceil(base_disk_gb + input_size * 2.0),
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
    
    svtk vcf2bed --info ALL ~{vcf} stdout \
      | bgzip -c \
      > "~{output_file}"
  >>>

  output {
    File vcf2bed_subworkflow_out = output_file
  }
}
