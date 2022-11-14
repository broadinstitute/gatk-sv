version 1.0

import "TasksMakeCohortVcf.wdl" as MiniTasks

# Workflow to perform per-sample benchmarking from an SV VCF vs an external dataset
workflow CollectPerSampleBenchmarking {
  input {
    File vcf_stats
    File samples_list
    File per_sample_tarball
    File comparison_tarball
    File? sample_renaming_tsv
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

  String output_prefix = "~{prefix}.sample_benchmark"

  if(defined(sample_renaming_tsv)) {
    # rename samples in external benchmarking data
    call RenameBenchmarkTarfileSamples {
      input:
        benchmark_tarfile=comparison_tarball,
        sample_renaming_tsv=select_first([sample_renaming_tsv]),
        sv_base_mini_docker=sv_base_mini_docker
    }
  }
  File comparison_tarball_=select_first([RenameBenchmarkTarfileSamples.benchmark_renamed_samples, comparison_tarball])

  call MiniTasks.SplitUncompressed as SplitShuffledList {
      input:
        whole_file=samples_list,
        lines_per_shard=samples_per_shard,
        shard_prefix="~{output_prefix}.list_shard.",
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
        comparison_tarball=comparison_tarball_,
        prefix=output_prefix,
        contigs=contigs,
        comparison_set_name=comparison_set_name,
        sv_pipeline_qc_docker=sv_pipeline_qc_docker,
        runtime_attr_override=runtime_override_benchmark_samples
    }
  }

  call MergeTarballs as MergeTarredResults {
    input:
      in_tarballs=BenchmarkSamples.benchmarking_results,
      folder_name="~{output_prefix}_vs_~{comparison_set_name}",
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_override_merge_and_tar_shard_benchmarks
  }

  # Return tarball of results
  output {
    File benchmarking_results_tarball = MergeTarredResults.tarball
  }
}


# rename samples in comparison benchmark to match the cohort VCF
task RenameBenchmarkTarfileSamples {
  input {
    File benchmark_tarfile
    File sample_renaming_tsv  # TSV file with original sample IDs in 1st column, and desired sample IDs in 2nd column
    String sv_base_mini_docker
    String renamed_suffix = "renamed_samples"
    RuntimeAttr? runtime_attr_override
  }

  String benchmark_renamed_samples_filename = basename(benchmark_tarfile, ".tar.gz") + "_~{renamed_suffix}.tar.gz"

  # Disk must be scaled proportionally to the size of the archive. Since the archive contains compressed files, the
  # extracted files should be sized comparably to the archive
  Float input_size = size(benchmark_tarfile, "GiB")
  RuntimeAttr default_attr = object {
    mem_gb: 2.0,
    disk_gb: ceil(10.0 + 3 * input_size),
    cpu_cores: 1,
    preemptible_tries: 3,
    max_retries: 1,
    boot_disk_gb: 10
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail
    # find the folder the benchmark tarfile will extract to
    ARCHIVE_DIR=$(basename $(tar --list -f "~{benchmark_tarfile}" | head -n1))
    # extrct archive
    tar -xf "~{benchmark_tarfile}"

    # define function to rename the sample IDs
    function rename_file() {
      filename=$1
      original_id=$2
      renamed_id=$3
      if [ "$original_id" != "$renamed_id" ]; then
        new_filename=$(echo "$filename" | sed "s/$original_id/$renamed_id/")
        mv "$filename" "$new_filename"
      fi
    }
    export -f rename_file

    # loop over sample IDs to rename and do the renaming
    cd "$ARCHIVE_DIR"
    while read ORIGINAL_SAMPLE_ID RENAMED_SAMPLE_ID; do
      find . -name "$ORIGINAL_SAMPLE_ID.*.bed.gz*" \
        | xargs -I{} bash -c "rename_file {} $ORIGINAL_SAMPLE_ID $RENAMED_SAMPLE_ID"
    done < "~{sample_renaming_tsv}"
    cd ..

    # archive the new files
    tar cz -f "~{benchmark_renamed_samples_filename}" "$ARCHIVE_DIR"

    # remove the unarchived files (in case this is being run in local mode)
    rm -r "$ARCHIVE_DIR"
  >>>

  output {
    File benchmark_renamed_samples = benchmark_renamed_samples_filename
  }

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_base_mini_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
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
  RuntimeAttr runtime_default = object {
    mem_gb: 3.75,
    disk_gb: ceil(10.0 + input_size * 15),
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
  RuntimeAttr runtime_default = object {
    mem_gb: 2.0,
    disk_gb: ceil(10.0 + input_size * 2.0),
    cpu_cores: 1,
    preemptible_tries: 3,
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
