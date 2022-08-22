version 1.0

# Workflow to scatter site-level benchmarking vs. an external dataset by chromosome

import "TasksMakeCohortVcf.wdl" as MiniTasks

workflow CollectSiteLevelBenchmarking {
  input {
    File vcf_stats
    String prefix
    Array[String] contigs
    String benchmark_url
    String benchmark_name
    String sv_pipeline_qc_docker
    String sv_base_mini_docker
    RuntimeAttr? runtime_override_site_level_benchmark
    RuntimeAttr? runtime_override_merge_site_level_benchmark
  }

  String output_prefix="~{prefix}.site_benchmark"
	
  # Collect site-level external benchmarking data per chromosome
  scatter ( contig in contigs ) {
    call VcfExternalBenchmarkSingleChrom {
      input:
        vcf_stats=vcf_stats,
        prefix=output_prefix,
        contig=contig,
        benchmark_url=benchmark_url,
        benchmark_name=benchmark_name,
        sv_pipeline_qc_docker=sv_pipeline_qc_docker,
        runtime_attr_override=runtime_override_site_level_benchmark
    }
  }

  # Merge results across chromosomes
  call MergeContigBenchmarks {
    input:
      in_tarballs=VcfExternalBenchmarkSingleChrom.benchmarking_results_tarball,
      prefix=output_prefix,
      benchmark_name=benchmark_name,
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_override_merge_site_level_benchmark
  }

  output {
    File benchmarking_results_tarball = MergeContigBenchmarks.merged_results_tarball
  }
}


# Task to collect external benchmarking data for a single chromosome
task VcfExternalBenchmarkSingleChrom {
  input {
    File vcf_stats
    String benchmark_url
    String prefix
    String contig
    String benchmark_name
    String sv_pipeline_qc_docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr runtime_default = object {
    mem_gb: 3.75,
    disk_gb: 40,
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

    # Copy benchmarking BED files to local directory
    mkdir benchmarks
    gsutil -m cp ~{benchmark_url}/*.bed.gz benchmarks/
    
    # Run benchmarking script
    echo ~{contig} > contigs.list
    /opt/sv-pipeline/scripts/vcf_qc/collectQC.external_benchmarking.sh \
      ~{vcf_stats} \
      /opt/sv-pipeline/scripts/vcf_qc/SV_colors.txt \
      contigs.list \
      benchmarks \
      collectQC_benchmarking_~{benchmark_name}_~{contig}_output/
    
    # Prep outputs
    tar -czvf ~{prefix}.collectQC_benchmarking_~{benchmark_name}_~{contig}_output.tar.gz \
      collectQC_benchmarking_~{benchmark_name}_~{contig}_output
  >>>

  output {
    File benchmarking_results_tarball = "~{prefix}.collectQC_benchmarking_~{benchmark_name}_~{contig}_output.tar.gz"
  }
}


# Task to merge external benchmarking data across chromosomes
task MergeContigBenchmarks {
  input {
    Array[File] in_tarballs
    String benchmark_name
    String prefix
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr runtime_default = object {
    mem_gb: 3.75,
    disk_gb: 40,
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
    # Untar all shards
    mkdir sharded_results
    while read tarball_path; do
      tar -xzvf "$tarball_path" --directory sharded_results/
    done < ~{write_lines(in_tarballs)}

    # Create final output directory
    mkdir ~{prefix}_collectQC_benchmarking_~{benchmark_name}_output
    mkdir ~{prefix}_collectQC_benchmarking_~{benchmark_name}_output/data

    # Merge each unique BED
    find sharded_results/ -name "*.bed.gz" | xargs -I {} basename {} | sort -V | uniq > bed_filenames.list
    while read fname; do
      find sharded_results/ -name $fname > matching_beds.list
      sed -n '1p' matching_beds.list | xargs -I {} zcat {} | sed -n '1p' > header.bed
      cat matching_beds.list | xargs -I {} zcat {} | fgrep -v "#" \
      | sort -Vk1,1 -k2,2n -k3,3n | cat header.bed - | bgzip -c \
      > ~{prefix}_collectQC_benchmarking_~{benchmark_name}_output/data/$fname
      tabix -f ~{prefix}_collectQC_benchmarking_~{benchmark_name}_output/data/$fname
    done < bed_filenames.list

    # Compress final output directory
    tar -czvf \
      ~{prefix}_collectQC_benchmarking_~{benchmark_name}_output.tar.gz \
      ~{prefix}_collectQC_benchmarking_~{benchmark_name}_output
  >>>

  output {
    File merged_results_tarball = "~{prefix}_collectQC_benchmarking_~{benchmark_name}_output.tar.gz"
  }
}
