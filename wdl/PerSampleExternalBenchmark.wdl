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
    Array[String] contigs
    String comparison_set_name
    Int samples_per_shard
    Int? random_seed

    String sv_base_mini_docker
    String sv_pipeline_docker
    String sv_pipeline_qc_docker

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
        contigs=contigs,
        comparison_set_name=comparison_set_name,
        sv_pipeline_qc_docker=sv_pipeline_qc_docker,
        runtime_attr_override=runtime_override_benchmark_samples
    }
  }

  call MergeTarballs as MergeTarredResults {
    input:
      in_tarballs=BenchmarkSamples.benchmarking_results,
      folder_name=prefix + "_vs_" + comparison_set_name,
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_override_merge_and_tar_shard_benchmarks
  }

  # Return tarball of results
  output {
    File benchmarking_results_tarball = MergeTarredResults.tarball
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
    Array[String] contigs
    String comparison_set_name
    String sv_pipeline_qc_docker
    RuntimeAttr? runtime_attr_override
  }

  String output_folder = "~{prefix}_~{comparison_set_name}_perSample_results"

  # Scale disk dynamically w/r/t input size
  Float input_size = size([vcf_stats, samples_list, per_sample_tarball, comparison_tarball], "GiB")
  Float compression_factor = 3.5
  Float base_disk_gb = 5.0
  RuntimeAttr runtime_default = object {
    mem_gb: 7.5,
    disk_gb: ceil(base_disk_gb + (input_size * compression_factor)),
    cpu_cores: 1,
    preemptible_tries: 1,
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
    docker: sv_pipeline_qc_docker
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
      ~{write_lines(contigs)} \
      ~{per_sample_tarball} \
      ~{comparison_tarball} \
      ~{output_folder}/

    # Tar benchmarking results for easier caching of downstream steps
    tar -czvf ~{output_folder}.tar.gz ~{output_folder}
  >>>

  output {
    File benchmarking_results = "~{output_folder}.tar.gz"
  }
}


# Task to merge benchmarking results across shards
task MergeTarballs {
  input {
    Array[File] in_tarballs
    String? folder_name
    String? tarball_prefix
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  String tar_folder_name = select_first([folder_name, "merged"])
  String outfile_name = select_first([tarball_prefix, tar_folder_name]) + ".tar.gz"

  # Since the input files are often/always compressed themselves, assume compression factor for tarring is 1.0
  Float input_size = size(in_tarballs, "GB")
  Float base_disk_gb = 10.0
  Float base_mem_gb = 2.0
  RuntimeAttr runtime_default = object {
    mem_gb: base_mem_gb,
    disk_gb: ceil(base_disk_gb + input_size * 2.0),
    cpu_cores: 1,
    preemptible_tries: 1,
    max_retries: 1,
    boot_disk_gb: 10
  }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_base_mini_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    # Create final output directory
    mkdir "~{tar_folder_name}"

    while read tarball_path; do
      tar -xzvf "$tarball_path" --directory ~{tar_folder_name}/
    done < ~{write_lines(in_tarballs)}

    # Compress final output directory
    tar -czvf "~{outfile_name}" "~{tar_folder_name}"
  >>>

  output {
    File tarball = outfile_name
  }
}
