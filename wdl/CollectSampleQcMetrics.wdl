version 1.0

import "Structs.wdl"

workflow CollectSampleQcMetrics {
  input {
    String sample_id
    File bam_or_cram_file
    File? bam_or_cram_index

    File ref_fasta
    File? ref_index
    File? ref_dict

    File? wgs_metrics_intervals
    File? multiple_metrics_intervals

    Int? read_length_override

    Boolean run_wgs_metrics = true
    Boolean run_raw_wgs_metrics = true
    Boolean run_alignment_summary_metrics = true
    Boolean run_insert_size_metrics = true
    Boolean run_quality_score_distribution = true
    Boolean run_mean_quality_by_cycle = true
    Boolean run_base_distribution_by_cycle = true
    Boolean run_gc_bias_metrics = true
    Boolean run_sequencing_artifact_metrics = true
    Boolean run_quality_yield_metrics = true
    Boolean run_multiple_metrics = false

    String genomes_in_the_cloud_docker
    String gatk_docker

    RuntimeAttr? runtime_attr_wgs_metrics
    RuntimeAttr? runtime_attr_raw_wgs_metrics
    RuntimeAttr? runtime_attr_alignment_summary_metrics
    RuntimeAttr? runtime_attr_insert_size_metrics
    RuntimeAttr? runtime_attr_quality_score_distribution
    RuntimeAttr? runtime_attr_mean_quality_by_cycle
    RuntimeAttr? runtime_attr_base_distribution_by_cycle
    RuntimeAttr? runtime_attr_gc_bias_metrics
    RuntimeAttr? runtime_attr_sequencing_artifact_metrics
    RuntimeAttr? runtime_attr_quality_yield_metrics
    RuntimeAttr? runtime_attr_multiple_metrics

  }

  Boolean is_bam_ = basename(bam_or_cram_file, ".bam") + ".bam" == basename(bam_or_cram_file)
  String index_ext_ = if is_bam_ then ".bai" else ".crai"
  File bam_or_cram_index_ = if defined(bam_or_cram_index) then select_first([bam_or_cram_index]) else bam_or_cram_file + index_ext_

  File ref_index_ = select_first([ref_index, ref_fasta + ".fai"])
  File ref_dict_ = select_first([ref_dict, sub(ref_fasta, "\\.fasta$", ".dict")])

  if (run_wgs_metrics) {
    call CollectWgsMetrics {
      input:
        sample_id = sample_id,
        bam_or_cram_file = bam_or_cram_file,
        bam_or_cram_index = bam_or_cram_index_,
        read_length = read_length_override,
        reference_fasta = ref_fasta,
        reference_index = ref_index_,
        intervals = wgs_metrics_intervals,
        genomes_in_the_cloud_docker = genomes_in_the_cloud_docker,
        runtime_attr_override = runtime_attr_wgs_metrics
    }
  }

  if (run_raw_wgs_metrics) {
    call CollectRawWgsMetrics {
      input:
        sample_id = sample_id,
        bam_or_cram_file = bam_or_cram_file,
        bam_or_cram_index = bam_or_cram_index_,
        read_length = read_length_override,
        reference_fasta = ref_fasta,
        reference_index = ref_index_,
        intervals = wgs_metrics_intervals,
        genomes_in_the_cloud_docker = genomes_in_the_cloud_docker,
        runtime_attr_override = runtime_attr_raw_wgs_metrics
    }
  }

  if (run_alignment_summary_metrics) {
    call CollectAlignmentSummaryMetrics {
      input:
        sample_id = sample_id,
        bam_or_cram_file = bam_or_cram_file,
        bam_or_cram_index = bam_or_cram_index_,
        reference_fasta = ref_fasta,
        reference_index = ref_index_,
        reference_dict = ref_dict_,
        genomes_in_the_cloud_docker = genomes_in_the_cloud_docker,
        runtime_attr_override = runtime_attr_raw_wgs_metrics
    }
  }

  if (run_insert_size_metrics) {
    call CollectInsertSizeMetrics {
      input:
        sample_id = sample_id,
        bam_or_cram_file = bam_or_cram_file,
        bam_or_cram_index = bam_or_cram_index_,
        reference_fasta = ref_fasta,
        reference_index = ref_index_,
        reference_dict = ref_dict_,
        genomes_in_the_cloud_docker = genomes_in_the_cloud_docker,
        runtime_attr_override = runtime_attr_insert_size_metrics
    }
  }

  if (run_quality_score_distribution) {
    call QualityScoreDistribution {
      input:
        sample_id = sample_id,
        bam_or_cram_file = bam_or_cram_file,
        bam_or_cram_index = bam_or_cram_index_,
        reference_fasta = ref_fasta,
        reference_index = ref_index_,
        reference_dict = ref_dict_,
        genomes_in_the_cloud_docker = genomes_in_the_cloud_docker,
        runtime_attr_override = runtime_attr_quality_score_distribution
    }
  }

  if (run_mean_quality_by_cycle) {
    call MeanQualityByCycle {
      input:
        sample_id = sample_id,
        bam_or_cram_file = bam_or_cram_file,
        bam_or_cram_index = bam_or_cram_index_,
        reference_fasta = ref_fasta,
        reference_index = ref_index_,
        reference_dict = ref_dict_,
        genomes_in_the_cloud_docker = genomes_in_the_cloud_docker,
        runtime_attr_override = runtime_attr_mean_quality_by_cycle
    }
  }

  if (run_base_distribution_by_cycle) {
    call CollectBaseDistributionByCycle {
      input:
        sample_id = sample_id,
        bam_or_cram_file = bam_or_cram_file,
        bam_or_cram_index = bam_or_cram_index_,
        reference_fasta = ref_fasta,
        reference_index = ref_index_,
        reference_dict = ref_dict_,
        genomes_in_the_cloud_docker = genomes_in_the_cloud_docker,
        runtime_attr_override = runtime_attr_base_distribution_by_cycle
    }
  }

  if (run_gc_bias_metrics) {
    call CollectGcBiasMetrics {
      input:
        sample_id = sample_id,
        bam_or_cram_file = bam_or_cram_file,
        bam_or_cram_index = bam_or_cram_index_,
        reference_fasta = ref_fasta,
        reference_index = ref_index_,
        reference_dict = ref_dict_,
        genomes_in_the_cloud_docker = genomes_in_the_cloud_docker,
        runtime_attr_override = runtime_attr_gc_bias_metrics
    }
  }

  if (run_sequencing_artifact_metrics) {
    call CollectSequencingArtifactMetrics {
      input:
        sample_id = sample_id,
        bam_or_cram_file = bam_or_cram_file,
        bam_or_cram_index = bam_or_cram_index_,
        reference_fasta = ref_fasta,
        reference_index = ref_index_,
        reference_dict = ref_dict_,
        intervals = multiple_metrics_intervals,
        genomes_in_the_cloud_docker = genomes_in_the_cloud_docker,
        runtime_attr_override = runtime_attr_sequencing_artifact_metrics
    }
  }

  if (run_quality_yield_metrics) {
    call CollectQualityYieldMetrics {
      input:
        sample_id = sample_id,
        bam_or_cram_file = bam_or_cram_file,
        bam_or_cram_index = bam_or_cram_index_,
        reference_fasta = ref_fasta,
        reference_index = ref_index_,
        reference_dict = ref_dict_,
        genomes_in_the_cloud_docker = genomes_in_the_cloud_docker,
        runtime_attr_override = runtime_attr_quality_yield_metrics
    }
  }

  if (run_multiple_metrics) {
    call CollectMultipleMetrics {
      input:
        sample_id = sample_id,
        bam_or_cram_file = bam_or_cram_file,
        bam_or_cram_index = bam_or_cram_index_,
        reference_fasta = ref_fasta,
        reference_index = ref_index_,
        reference_dict = ref_dict_,
        intervals = multiple_metrics_intervals,
        gatk_docker = gatk_docker,
        runtime_attr_override = runtime_attr_multiple_metrics
    }
  }

  output {
    File? wgs_metrics = CollectWgsMetrics.metrics_file
    File? raw_wgs_metrics = CollectRawWgsMetrics.metrics_file
    File? alignment_summary_metrics = CollectAlignmentSummaryMetrics.metrics_file
    File? insert_size_metrics = CollectInsertSizeMetrics.metrics_file
    File? insert_size_histogram = CollectInsertSizeMetrics.insert_size_histogram
    File? quality_score_table = QualityScoreDistribution.table
    File? quality_score_chart = QualityScoreDistribution.chart
    File? mean_quality_by_cycle_table = MeanQualityByCycle.table
    File? mean_quality_by_cycle_chart = MeanQualityByCycle.chart
    File? base_distribution_by_cycle_table = CollectBaseDistributionByCycle.table
    File? base_distribution_by_cycle_chart = CollectBaseDistributionByCycle.chart
    File? gc_bias_metrics = CollectGcBiasMetrics.metrics_file
    File? gc_bias_summary_metrics = CollectGcBiasMetrics.summary_metrics_file
    File? gc_bias_chart = CollectGcBiasMetrics.chart
    File? sequencing_artifact_detail_metrics = CollectSequencingArtifactMetrics.metrics_file
    File? sequencing_artifact_summary_metrics = CollectSequencingArtifactMetrics.summary_metrics_file
    File? quality_yield_metrics = CollectQualityYieldMetrics.metrics_file
    Array[File]? multiple_metrics = CollectMultipleMetrics.metrics_files
  }
}


task CollectWgsMetrics {
  input {
    String sample_id
    File bam_or_cram_file
    File bam_or_cram_index
    Int? read_length
    File reference_fasta
    File reference_index
    File? intervals
    Boolean? use_fast_algorithm
    String genomes_in_the_cloud_docker
    RuntimeAttr? runtime_attr_override
  }

  Boolean fast_algorithm = select_first([use_fast_algorithm, true])

  String metrics_file_name = "~{sample_id}.wgs_metrics.txt"

  Int num_cpu = if defined(runtime_attr_override) then select_first([select_first([runtime_attr_override]).cpu_cores, 1]) else 1
  Float mem_size_gb = num_cpu * 8.0
  Int java_mem_mb = round(mem_size_gb * 0.8 * 1000)
  Float bam_or_cram_size = size(bam_or_cram_file, "GiB")
  Float disk_overhead = 20.0
  Float ref_size = size(reference_fasta, "GiB")
  Int vm_disk_size = ceil(bam_or_cram_size + ref_size + disk_overhead)

  RuntimeAttr default_attr = object {
    cpu_cores: num_cpu,
    mem_gb: mem_size_gb,
    disk_gb: vm_disk_size,
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: genomes_in_the_cloud_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

  command <<<

    set -Eeuo pipefail

    # calculate coverage
    java -Xms~{java_mem_mb}m -jar /usr/gitc/picard.jar \
      CollectWgsMetrics \
      INPUT=~{bam_or_cram_file} \
      VALIDATION_STRINGENCY=SILENT \
      REFERENCE_SEQUENCE=~{reference_fasta} \
      ~{"READ_LENGTH=" + read_length} \
      INCLUDE_BQ_HISTOGRAM=true \
      ~{if defined(intervals) then "INTERVALS=~{intervals}" else ""} \
      OUTPUT="~{metrics_file_name}" \
      USE_FAST_ALGORITHM=~{fast_algorithm}

  >>>

  output {
    File metrics_file = metrics_file_name
  }
}

task CollectRawWgsMetrics {
  input {
    String sample_id
    File bam_or_cram_file
    File bam_or_cram_index
    Int? read_length
    File reference_fasta
    File reference_index
    File? intervals
    Boolean? use_fast_algorithm
    String genomes_in_the_cloud_docker
    RuntimeAttr? runtime_attr_override
  }

  Boolean fast_algorithm = select_first([use_fast_algorithm, true])

  String metrics_file_name = "~{sample_id}.raw_wgs_metrics.txt"

  Int num_cpu = if defined(runtime_attr_override) then select_first([select_first([runtime_attr_override]).cpu_cores, 1]) else 1
  Float mem_size_gb = num_cpu * 8.0
  Int java_mem_mb = round(mem_size_gb * 0.8 * 1000)
  Float bam_or_cram_size = size(bam_or_cram_file, "GiB")
  Float disk_overhead = 20.0
  Float ref_size = size(reference_fasta, "GiB")
  Int vm_disk_size = ceil(bam_or_cram_size + ref_size + disk_overhead)

  RuntimeAttr default_attr = object {
    cpu_cores: num_cpu,
    mem_gb: mem_size_gb,
    disk_gb: vm_disk_size,
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: genomes_in_the_cloud_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

  command <<<

    set -Eeuo pipefail

    # calculate coverage
    java -Xms~{java_mem_mb}m -jar /usr/gitc/picard.jar \
      CollectRawWgsMetrics \
      INPUT=~{bam_or_cram_file} \
      VALIDATION_STRINGENCY=SILENT \
      REFERENCE_SEQUENCE=~{reference_fasta} \
      ~{"READ_LENGTH=" + read_length} \
      INCLUDE_BQ_HISTOGRAM=true \
      ~{if defined(intervals) then "INTERVALS=~{intervals}" else ""} \
      OUTPUT="~{metrics_file_name}" \
      USE_FAST_ALGORITHM=~{fast_algorithm}

  >>>

  output {
    File metrics_file = metrics_file_name
  }
}

task CollectAlignmentSummaryMetrics {
  input {
    String sample_id
    File bam_or_cram_file
    File bam_or_cram_index
    File reference_fasta
    File reference_index
    File reference_dict
    File? intervals
    String genomes_in_the_cloud_docker
    RuntimeAttr? runtime_attr_override
  }

  String metrics_file_name = "~{sample_id}.alignment_summary_metrics.txt"

  Int num_cpu = if defined(runtime_attr_override) then select_first([select_first([runtime_attr_override]).cpu_cores, 1]) else 1
  Float mem_size_gb = num_cpu * 8.0
  Int java_mem_mb = round(mem_size_gb * 0.8 * 1000)
  Float bam_or_cram_size = size(bam_or_cram_file, "GiB")
  Float disk_overhead = 20.0
  Float ref_size = size(reference_fasta, "GiB")
  Int vm_disk_size = ceil(bam_or_cram_size + ref_size + disk_overhead)

  RuntimeAttr default_attr = object {
    cpu_cores: num_cpu,
    mem_gb: mem_size_gb,
    disk_gb: vm_disk_size,
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: genomes_in_the_cloud_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

  command <<<

    set -Eeuo pipefail

    # calculate coverage
    java -Xms~{java_mem_mb}m -jar /usr/gitc/picard.jar \
      CollectAlignmentSummaryMetrics \
      INPUT=~{bam_or_cram_file} \
      REFERENCE_SEQUENCE=~{reference_fasta} \
      OUTPUT="~{metrics_file_name}" 
  >>>

  output {
    File metrics_file = metrics_file_name
  }
}

task CollectInsertSizeMetrics {
  input {
    String sample_id
    File bam_or_cram_file
    File? bam_or_cram_index
    File reference_fasta
    File reference_index
    File reference_dict
    String genomes_in_the_cloud_docker
    RuntimeAttr? runtime_attr_override
  }

  String metrics_file_name = "~{sample_id}.insert_size_metrics.txt"
  String metrics_histogram_name = "~{sample_id}.insert_size_histogram.pdf"

  Int num_cpu = if defined(runtime_attr_override) then select_first([select_first([runtime_attr_override]).cpu_cores, 1]) else 1
  Float mem_size_gb = num_cpu * 8.0
  Int java_mem_mb = round(mem_size_gb * 0.8 * 1000)
  Float bam_or_cram_size = size(bam_or_cram_file, "GiB")
  Float disk_overhead = 20.0
  Float ref_size = size(reference_fasta, "GiB")
  Int vm_disk_size = ceil(bam_or_cram_size + ref_size + disk_overhead)

  RuntimeAttr default_attr = object {
    cpu_cores: num_cpu,
    mem_gb: mem_size_gb,
    disk_gb: vm_disk_size,
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: genomes_in_the_cloud_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

  command <<<

    set -Eeuo pipefail

    # calculate coverage
    java -Xms~{java_mem_mb}m -jar /usr/gitc/picard.jar \
      CollectInsertSizeMetrics \
      INPUT=~{bam_or_cram_file} \
      REFERENCE_SEQUENCE=~{reference_fasta} \
      HISTOGRAM_FILE="~{metrics_histogram_name}" \
      OUTPUT="~{metrics_file_name}" 
  >>>

  output {
    File metrics_file = metrics_file_name
    File insert_size_histogram = metrics_histogram_name
  }
}

task QualityScoreDistribution {
  input {
    String sample_id
    File bam_or_cram_file
    File? bam_or_cram_index
    File reference_fasta
    File reference_index
    File reference_dict
    String genomes_in_the_cloud_docker
    RuntimeAttr? runtime_attr_override
  }

  String table_name = "~{sample_id}.quality_score_distribution.txt"
  String chart_name = "~{sample_id}.quality_score_chart.pdf"

  Int num_cpu = if defined(runtime_attr_override) then select_first([select_first([runtime_attr_override]).cpu_cores, 1]) else 1
  Float mem_size_gb = num_cpu * 8.0
  Int java_mem_mb = round(mem_size_gb * 0.8 * 1000)
  Float bam_or_cram_size = size(bam_or_cram_file, "GiB")
  Float disk_overhead = 20.0
  Float ref_size = size(reference_fasta, "GiB")
  Int vm_disk_size = ceil(bam_or_cram_size + ref_size + disk_overhead)

  RuntimeAttr default_attr = object {
    cpu_cores: num_cpu,
    mem_gb: mem_size_gb,
    disk_gb: vm_disk_size,
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: genomes_in_the_cloud_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

  command <<<

    set -Eeuo pipefail

    # calculate coverage
    java -Xms~{java_mem_mb}m -jar /usr/gitc/picard.jar \
      QualityScoreDistribution \
      INPUT=~{bam_or_cram_file} \
      REFERENCE_SEQUENCE=~{reference_fasta} \
      CHART_OUTPUT="~{chart_name}" \
      OUTPUT="~{table_name}" 
  >>>

  output {
    File table = table_name
    File chart = chart_name
  }
}

task MeanQualityByCycle {
  input {
    String sample_id
    File bam_or_cram_file
    File? bam_or_cram_index
    File reference_fasta
    File reference_index
    File reference_dict
    String genomes_in_the_cloud_docker
    RuntimeAttr? runtime_attr_override
  }

  String table_name = "~{sample_id}.mean_quality_by_cycle_table.txt"
  String chart_name = "~{sample_id}.mean_quality_by_cycle_chart.pdf"

  Int num_cpu = if defined(runtime_attr_override) then select_first([select_first([runtime_attr_override]).cpu_cores, 1]) else 1
  Float mem_size_gb = num_cpu * 8.0
  Int java_mem_mb = round(mem_size_gb * 0.8 * 1000)
  Float bam_or_cram_size = size(bam_or_cram_file, "GiB")
  Float disk_overhead = 20.0
  Float ref_size = size(reference_fasta, "GiB")
  Int vm_disk_size = ceil(bam_or_cram_size + ref_size + disk_overhead)

  RuntimeAttr default_attr = object {
    cpu_cores: num_cpu,
    mem_gb: mem_size_gb,
    disk_gb: vm_disk_size,
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: genomes_in_the_cloud_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

  command <<<

    set -Eeuo pipefail

    # calculate coverage
    java -Xms~{java_mem_mb}m -jar /usr/gitc/picard.jar \
      MeanQualityByCycle \
      INPUT=~{bam_or_cram_file} \
      REFERENCE_SEQUENCE=~{reference_fasta} \
      CHART_OUTPUT="~{chart_name}" \
      OUTPUT="~{table_name}" 
  >>>

  output {
    File table = table_name
    File chart = chart_name
  }
}


task CollectBaseDistributionByCycle {
  input {
    String sample_id
    File bam_or_cram_file
    File? bam_or_cram_index
    File reference_fasta
    File reference_index
    File reference_dict
    String genomes_in_the_cloud_docker
    RuntimeAttr? runtime_attr_override
  }

  String table_name = "~{sample_id}.base_distribution_by_cycle_table.txt"
  String chart_name = "~{sample_id}.base_distribution_by_cycle_chart.pdf"

  Int num_cpu = if defined(runtime_attr_override) then select_first([select_first([runtime_attr_override]).cpu_cores, 1]) else 1
  Float mem_size_gb = num_cpu * 8.0
  Int java_mem_mb = round(mem_size_gb * 0.8 * 1000)
  Float bam_or_cram_size = size(bam_or_cram_file, "GiB")
  Float disk_overhead = 20.0
  Float ref_size = size(reference_fasta, "GiB")
  Int vm_disk_size = ceil(bam_or_cram_size + ref_size + disk_overhead)

  RuntimeAttr default_attr = object {
    cpu_cores: num_cpu,
    mem_gb: mem_size_gb,
    disk_gb: vm_disk_size,
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: genomes_in_the_cloud_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

  command <<<

    set -Eeuo pipefail

    # calculate coverage
    java -Xms~{java_mem_mb}m -jar /usr/gitc/picard.jar \
      CollectBaseDistributionByCycle \
      INPUT=~{bam_or_cram_file} \
      REFERENCE_SEQUENCE=~{reference_fasta} \
      CHART_OUTPUT="~{chart_name}" \
      OUTPUT="~{table_name}" 
  >>>

  output {
    File table = table_name
    File chart = chart_name
  }
}

task CollectGcBiasMetrics {
  input {
    String sample_id
    File bam_or_cram_file
    File? bam_or_cram_index
    File reference_fasta
    File reference_index
    File reference_dict
    String genomes_in_the_cloud_docker
    RuntimeAttr? runtime_attr_override
  }

  String metrics_file_name = "~{sample_id}.gc_bias_metrics.txt"
  String summary_metrics_file_name = "~{sample_id}.gc_bias_summary_metrics.txt"
  String chart_name = "~{sample_id}.gc_bias_chart.pdf"

  Int num_cpu = if defined(runtime_attr_override) then select_first([select_first([runtime_attr_override]).cpu_cores, 1]) else 1
  Float mem_size_gb = num_cpu * 8.0
  Int java_mem_mb = round(mem_size_gb * 0.8 * 1000)
  Float bam_or_cram_size = size(bam_or_cram_file, "GiB")
  Float disk_overhead = 20.0
  Float ref_size = size(reference_fasta, "GiB")
  Int vm_disk_size = ceil(bam_or_cram_size + ref_size + disk_overhead)

  RuntimeAttr default_attr = object {
    cpu_cores: num_cpu,
    mem_gb: mem_size_gb,
    disk_gb: vm_disk_size,
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: genomes_in_the_cloud_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

  command <<<

    set -Eeuo pipefail

    # calculate coverage
    java -Xms~{java_mem_mb}m -jar /usr/gitc/picard.jar \
      CollectGcBiasMetrics \
      INPUT=~{bam_or_cram_file} \
      REFERENCE_SEQUENCE=~{reference_fasta} \
      CHART_OUTPUT="~{chart_name}" \
      SUMMARY_OUTPUT="~{summary_metrics_file_name}" \
      OUTPUT="~{metrics_file_name}" 
  >>>

  output {
    File metrics_file = metrics_file_name
    File summary_metrics_file = summary_metrics_file_name
    File chart = chart_name
  }
}

task CollectSequencingArtifactMetrics {
  input {
    String sample_id
    File bam_or_cram_file
    File bam_or_cram_index
    File reference_fasta
    File reference_index
    File reference_dict
    File? intervals
    String genomes_in_the_cloud_docker
    RuntimeAttr? runtime_attr_override
  }

  String metrics_file_base = "~{sample_id}"

  Int num_cpu = if defined(runtime_attr_override) then select_first([select_first([runtime_attr_override]).cpu_cores, 1]) else 1
  Float mem_size_gb = num_cpu * 8.0
  Int java_mem_mb = round(mem_size_gb * 0.8 * 1000)
  Float bam_or_cram_size = size(bam_or_cram_file, "GiB")
  Float disk_overhead = 20.0
  Float ref_size = size(reference_fasta, "GiB")
  Int vm_disk_size = ceil(bam_or_cram_size + ref_size + disk_overhead)

  RuntimeAttr default_attr = object {
    cpu_cores: num_cpu,
    mem_gb: mem_size_gb,
    disk_gb: vm_disk_size,
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: genomes_in_the_cloud_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

  command <<<

    set -Eeuo pipefail

    # calculate coverage
    java -Xms~{java_mem_mb}m -jar /usr/gitc/picard.jar \
      CollectSequencingArtifactMetrics \
      INPUT=~{bam_or_cram_file} \
      REFERENCE_SEQUENCE=~{reference_fasta} \
      FILE_EXTENSION=".txt" \
      ~{if defined(intervals) then "--INTERVALS ~{intervals}" else ""} \
      OUTPUT="~{metrics_file_base}" 
  >>>

  output {
    File metrics_file = "~{metrics_file_base}.pre_adapter_detail_metrics.txt"
    File summary_metrics_file = "~{metrics_file_base}.pre_adapter_summary_metrics.txt"
  }
}

task CollectQualityYieldMetrics {
  input {
    String sample_id
    File bam_or_cram_file
    File? bam_or_cram_index
    File reference_fasta
    File reference_index
    File reference_dict
    String genomes_in_the_cloud_docker
    RuntimeAttr? runtime_attr_override
  }

  String metrics_file_name = "~{sample_id}.quality_yield_metrics.txt"

  Int num_cpu = if defined(runtime_attr_override) then select_first([select_first([runtime_attr_override]).cpu_cores, 1]) else 1
  Float mem_size_gb = num_cpu * 8.0
  Int java_mem_mb = round(mem_size_gb * 0.8 * 1000)
  Float bam_or_cram_size = size(bam_or_cram_file, "GiB")
  Float disk_overhead = 20.0
  Float ref_size = size(reference_fasta, "GiB")
  Int vm_disk_size = ceil(bam_or_cram_size + ref_size + disk_overhead)

  RuntimeAttr default_attr = object {
    cpu_cores: num_cpu,
    mem_gb: mem_size_gb,
    disk_gb: vm_disk_size,
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: genomes_in_the_cloud_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

  command <<<

    set -Eeuo pipefail

    # calculate coverage
    java -Xms~{java_mem_mb}m -jar /usr/gitc/picard.jar \
      CollectQualityYieldMetrics \
      INPUT=~{bam_or_cram_file} \
      REFERENCE_SEQUENCE=~{reference_fasta} \
      OUTPUT="~{metrics_file_name}" 
  >>>

  output {
    File metrics_file = metrics_file_name
  }
}


task CollectMultipleMetrics {
  input {
    String sample_id
    File bam_or_cram_file
    File bam_or_cram_index
    File reference_fasta
    File? reference_index
    File? reference_dict
    File? intervals
    String gatk_docker
    RuntimeAttr? runtime_attr_override
  }

  Int num_cpu = if defined(runtime_attr_override)
    then select_first([select_first([runtime_attr_override]).cpu_cores, 1])
    else 1

  Float mem_use_gb = num_cpu * 6.0
  Float java_mem_pad_gb = num_cpu * 1.0
  Float mem_size_gb = mem_use_gb + java_mem_pad_gb
  Float bam_or_cram_size = size(bam_or_cram_file, "GiB")
  Float disk_overhead = 20.0
  Float ref_size = size(reference_fasta, "GiB")
  Int vm_disk_size = ceil(bam_or_cram_size + ref_size + disk_overhead)

  String metrics_base = "~{sample_id}.multiple_metrics"
  String alignment_summary_filename = metrics_base + ".alignment_summary_metrics"
  String insert_size_filename = metrics_base + ".insert_size_metrics"
  String gc_bias_filename = metrics_base + ".gc_bias.summary_metrics"
  String metrics_table_filename=metrics_base + "_table.tsv"

  RuntimeAttr default_attr = object {
    cpu_cores: num_cpu,
    mem_gb: mem_size_gb,
    disk_gb: vm_disk_size,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: gatk_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
  Int java_mem_mb = round(1024 * (select_first([runtime_attr.mem_gb, default_attr.mem_gb]) - java_mem_pad_gb))

  command <<<

    set -Eeuo pipefail

    gatk --java-options -Xmx~{java_mem_mb}m CollectMultipleMetrics \
      -I "~{bam_or_cram_file}" \
      -O "~{metrics_base}" \
      -R "~{reference_fasta}" \
      ~{if defined(intervals) then "--INTERVALS ~{intervals}" else ""} \
      --ASSUME_SORTED true \
      --PROGRAM null \
      --PROGRAM CollectAlignmentSummaryMetrics \
      --PROGRAM CollectInsertSizeMetrics \
      --PROGRAM CollectSequencingArtifactMetrics \
      --PROGRAM CollectGcBiasMetrics \
      --PROGRAM QualityScoreDistribution \
      --METRIC_ACCUMULATION_LEVEL null \
      --METRIC_ACCUMULATION_LEVEL SAMPLE

    
  >>>

  output {
    Array[File] metrics_files = glob("~{metrics_base}.*")
  }
}
