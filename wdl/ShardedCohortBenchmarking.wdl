version 1.0

# Author: Ryan Collins <rlcollins@g.harvard.edu>

# Workflow to scatter site-level benchmarking vs. an external dataset by chromosome

import "Tasks0506.wdl" as MiniTasks

workflow ShardedCohortBenchmarking {
  input {
    File vcf_stats
    String prefix
    Array[String] contigs
    String benchmarking_bucket
    String comparator
    String sv_pipeline_qc_docker
    String sv_base_mini_docker
    RuntimeAttr? runtime_override_site_level_benchmark
    RuntimeAttr? runtime_override_merge_site_level_benchmark
  }
	
  # Collect site-level external benchmarking data per chromosome
  scatter ( contig in contigs ) {
    call VcfExternalBenchmarkSingleChrom as CollectSiteLevelBenchmarking {
      input:
        vcf_stats=vcf_stats,
        prefix=prefix,
        contig=contig,
        benchmarking_bucket=benchmarking_bucket,
        comparator=comparator,
        sv_pipeline_qc_docker=sv_pipeline_qc_docker,
        runtime_attr_override=runtime_override_site_level_benchmark
    }
  }

  # Merge results across chromosomes
  call MergeContigBenchmarks as MergeBenchmarking {
    input:
      in_tarballs=CollectSiteLevelBenchmarking.benchmarking_results_tarball,
      comparator=comparator,
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_override_merge_site_level_benchmark
  }

  output {
    File benchmarking_results_tarball = MergeBenchmarking.merged_results_tarball
  }
}


# Task to collect external benchmarking data for a single chromosome
task VcfExternalBenchmarkSingleChrom {
  input {
    File vcf_stats
    String benchmarking_bucket
    String prefix
    String contig
    String comparator
    String sv_pipeline_qc_docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr runtime_default = object {
    mem_gb: 3.75,
    disk_gb: 40,
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

    # Copy benchmarking BED files to local directory
    mkdir benchmarks
    gsutil -m cp ~{benchmarking_bucket}/*.bed.gz benchmarks/
    
    # Run benchmarking script
    echo ~{contig} > contigs.list
    /opt/sv-pipeline/scripts/vcf_qc/collectQC.external_benchmarking.sh \
      ~{vcf_stats} \
      /opt/sv-pipeline/scripts/vcf_qc/SV_colors.txt \
      contigs.list \
      benchmarks \
      collectQC_benchmarking_~{comparator}_~{contig}_output/
    
    # Prep outputs
    tar -czvf ~{prefix}.collectQC_benchmarking_~{comparator}_~{contig}_output.tar.gz \
      collectQC_benchmarking_~{comparator}_~{contig}_output
  >>>

  output {
    File benchmarking_results_tarball = "~{prefix}.collectQC_benchmarking_~{comparator}_~{contig}_output.tar.gz"
  }
}


# Task to merge external benchmarking data across chromosomes
task MergeContigBenchmarks {
  input {
    Array[File] in_tarballs
    String comparator
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr runtime_default = object {
    mem_gb: 3.75,
    disk_gb: 40,
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
    mkdir collectQC_benchmarking_~{comparator}_output
    mkdir collectQC_benchmarking_~{comparator}_output/data

    # Merge each unique BED
    find sharded_results/ -name "*.bed.gz" | xargs -I {} basename {} | sort -V | uniq > bed_filenames.list
    while read fname; do
      find sharded_results/ -name $fname > matching_beds.list
      sed -n '1p' matching_beds.list | xargs -I {} zcat {} | sed -n '1p' > header.bed
      cat matching_beds.list | xargs -I {} zcat {} | fgrep -v "#" \
      | sort -Vk1,1 -k2,2n -k3,3n | cat header.bed - | bgzip -c \
      > collectQC_benchmarking_~{comparator}_output/data/$fname
      tabix -f collectQC_benchmarking_~{comparator}_output/data/$fname
    done < bed_filenames.list

    # Compress final output directory
    tar -czvf \
      collectQC_benchmarking_~{comparator}_output.tar.gz \
      collectQC_benchmarking_~{comparator}_output
  >>>

  output {
    File merged_results_tarball = "collectQC_benchmarking_~{comparator}_output.tar.gz"
  }
}
