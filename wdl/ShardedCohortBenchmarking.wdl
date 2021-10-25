version 1.0

# Author: Ryan Collins <rlcollins@g.harvard.edu>

# Workflow to scatter site-level benchmarking vs. an external dataset by chromosome

workflow ShardedCohortBenchmarking {
  input {
    File vcf_stats
    String prefix
    Array[String] contigs
    String benchmarking_bucket
    String comparator
    String sv_pipeline_qc_docker
    RuntimeAttr? runtime_override_site_level_benchmark

  }
	
  # Collect site-level external benchmarking data per chromosome
  scatter ( contig in contigs ) {
    call VcfExternalBenchmark as CollectSiteLevelBenchmarking {
      input:
        vcf_stats=vcf_stats,
        prefix=prefix,
        contigs=[contig],
        benchmarking_bucket=benchmarking_bucket,
        comparator=comparator,
        sv_pipeline_qc_docker=sv_pipeline_qc_docker,
        runtime_attr_override=runtime_override_site_level_benchmark
    }
  }
  

}


# Task to collect external benchmarking data
task VcfExternalBenchmark {
  input {
    File vcf_stats
    String benchmarking_bucket
    String prefix
    Array[String] contigs
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
    /opt/sv-pipeline/scripts/vcf_qc/collectQC.external_benchmarking.sh \
      ~{vcf_stats} \
      /opt/sv-pipeline/scripts/vcf_qc/SV_colors.txt \
      ~{write_lines(contigs)} \
      benchmarks \
      collectQC_benchmarking_~{comparator}_output/
    
    # Prep outputs
    tar -czvf ~{prefix}.collectQC_benchmarking_~{comparator}_output.tar.gz \
      collectQC_benchmarking_~{comparator}_output
  >>>

  output {
    File benchmarking_results_tarball = "~{prefix}.collectQC_benchmarking_~{comparator}_output.tar.gz"
  }
}