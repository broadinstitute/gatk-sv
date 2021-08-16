version 1.0

# Author: Ryan Collins <rlcollins@g.harvard.edu>

import "TasksMakeCohortVcf.wdl" as MiniTasks

# Workflow to perform per-sample benchmarking from an SV VCF vs an external dataset
workflow PerSampleExternalBenchmark {
  input {
    File vcf_stats
    File samples_list
    File per_sample_tarball
    File comparison_tarball
    String prefix
    String comparison_set_name
    Int samples_per_shard
    Int? random_seed

    String sv_base_mini_docker
    String sv_pipeline_docker

    # overrides for local tasks
    RuntimeAttr? runtime_override_benchmark_samples

    # overrides for mini tasks
    RuntimeAttr? runtime_override_split_shuffled_list
    RuntimeAttr? runtime_override_merge_and_tar_shard_benchmarks
  }

  call MiniTasks.SplitUncompressed as SplitShuffledList {
      input:
        whole_file=samples_list,
        lines_per_shard=samples_per_shard,
        shard_prefix=prefix + ".list_shard.",
        shuffle_file=true,
        random_seed=random_seed,
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_override_split_shuffled_list
    }

  # Collect benchmarking results per sample list shard
  scatter (sublist in SplitShuffledList.shards) {
    call BenchmarkSamples {
      input:
        vcf_stats=vcf_stats,
        samples_list=sublist,
        per_sample_tarball=per_sample_tarball,
        comparison_tarball=comparison_tarball,
        prefix=prefix,
        comparison_set_name=comparison_set_name,
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_override_benchmark_samples
    }
  }

  call MiniTasks.FilesToTarredFolder as MergeAndTarShardBenchmarks {
    input:
      in_files=flatten(BenchmarkSamples.benchmarking_results),
      folder_name="~{prefix}_~{comparison_set_name}_results_merged",
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_override_merge_and_tar_shard_benchmarks
  }

  # Return tarball of results
  output {
    File benchmarking_results_tarball = MergeAndTarShardBenchmarks.tarball
  }
}


# Task to collect per-sample benchmarking stats
task BenchmarkSamples {
  input {
    File vcf_stats
    File samples_list
    File per_sample_tarball
    File comparison_tarball
    String prefix
    String comparison_set_name
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  String output_folder = "~{prefix}_~{comparison_set_name}_perSample_results"

  # when filtering/sorting/etc, memory usage will likely go up (much of the data will have to
  # be held in memory or disk while working, potentially in a form that takes up more space)
  Float input_size = size([vcf_stats, samples_list, per_sample_tarball, comparison_tarball], "GiB")
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
    
    # Run benchmarking script
    mkdir ~{output_folder}
    /opt/sv-pipeline/scripts/vcf_qc/collectQC.perSample_benchmarking.sh \
      -p ~{comparison_set_name} \
      ~{vcf_stats} \
      ~{samples_list} \
      ~{per_sample_tarball} \
      ~{comparison_tarball} \
      ~{output_folder}/
  >>>

  output {
    Array[File] benchmarking_results = flatten([
      glob("~{output_folder}/*.sensitivity.bed.gz"),
      glob("~{output_folder}/*.specificity.bed.gz")
    ])
  }
}

