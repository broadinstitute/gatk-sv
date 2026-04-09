version 1.0

import "CollectCoverage.wdl" as cov
import "Structs.wdl"

workflow CallGenomicDisorderCNVs {
  input {
    String batch

    File baf_matrix
    File baf_matrix_index
    File high_res_rd_matrix
    File high_res_rd_matrix_index
    File reference_dict

    File gd_table
    File segdup_bed
    File centromere_bed
    File acrocentric_arm_bed
    Array[File]? flank_exclusion_intervals
    File gaps_bed
    File gtf
    File transition_matrix
    File breakpoint_transition_matrix
    File? truth_table

    Int rebinned_interval_size = 100000
    String sv_pipeline_docker
    String gatk_docker
    String? preprocess_args
    String? infer_args
    String? call_args
    String? eval_args
    RuntimeAttr? runtime_attr_condense
    RuntimeAttr? runtime_attr_gd
  }

  call cov.CondenseReadCounts {
    input:
      counts = high_res_rd_matrix,
      output_prefix = batch + ".low_res",
      sequence_dictionary = reference_dict,
      min_interval_size = rebinned_interval_size,
      max_interval_size = rebinned_interval_size,
      gatk_docker = gatk_docker,
      runtime_attr_override = runtime_attr_condense
  }

  call RunGenomicDisorderCNVs {
    input:
      batch = batch,
      low_res_counts = CondenseReadCounts.out,
      baf_matrix = baf_matrix,
      baf_matrix_index = baf_matrix_index,
      high_res_counts = high_res_rd_matrix,
      high_res_counts_index = high_res_rd_matrix_index,
      gd_table = gd_table,
      segdup_bed = segdup_bed,
      centromere_bed = centromere_bed,
      acrocentric_arm_bed = acrocentric_arm_bed,
      flank_exclusion_intervals = flank_exclusion_intervals,
      gaps_bed = gaps_bed,
      gtf = gtf,
      transition_matrix = transition_matrix,
      breakpoint_transition_matrix = breakpoint_transition_matrix,
      truth_table = truth_table,
      sv_pipeline_docker = sv_pipeline_docker,
      preprocess_args = preprocess_args,
      infer_args = infer_args,
      call_args = call_args,
      eval_args = eval_args,
      runtime_attr_override = runtime_attr_gd
  }

  output {
    File gd_output_tarball = RunGenomicDisorderCNVs.output_tarball
  }
}

task RunGenomicDisorderCNVs {
  input {
    String batch
    File low_res_counts
    File baf_matrix
    File baf_matrix_index
    File high_res_counts
    File high_res_counts_index
    File gd_table
    File segdup_bed
    File centromere_bed
    File acrocentric_arm_bed
    Array[File]? flank_exclusion_intervals
    File gaps_bed
    File gtf
    File transition_matrix
    File breakpoint_transition_matrix
    File? truth_table
    String sv_pipeline_docker
    String? preprocess_args
    String? infer_args
    String? call_args
    String? eval_args
    RuntimeAttr? runtime_attr_override
  }

  String output_tarball_name = batch + ".gatk_sv_gd.tar.gz"
  Int disk_gb = ceil(50 + size([
    low_res_counts,
    baf_matrix,
    high_res_counts,
  ], "GiB"))

  RuntimeAttr default_attr = object {
    cpu_cores: 2,
    mem_gb: 16.0,
    disk_gb: disk_gb,
    boot_disk_gb: 20,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail

    WORK_DIR="gd_outputs"
    bash /opt/gatk-sv-gd/run_gd.sh \
      --work-dir "${WORK_DIR}" \
      --input-depth ~{low_res_counts} \
      --high-res-depth ~{high_res_counts} \
      --baf-table ~{baf_matrix} \
      --gd-table ~{gd_table} \
      --segdup-bed ~{segdup_bed} \
      --centromere-bed ~{centromere_bed} \
      --acrocentric-arm-bed ~{acrocentric_arm_bed} \
      ~{if length(select_first([flank_exclusion_intervals, []])) > 0 then "--flank-exclusion-intervals" else ""} ~{sep=" " select_first([flank_exclusion_intervals, []])} \
      --gaps-bed ~{gaps_bed} \
      --gtf ~{gtf} \
      --transition-matrix ~{transition_matrix} \
      --breakpoint-transition-matrix ~{breakpoint_transition_matrix} \
      --preprocess-args '~{default="" preprocess_args}' \
      --infer-args '~{default="" infer_args}' \
      --call-args '~{default="" call_args}' \
      --eval-args '~{default="" eval_args}' \
      ~{if defined(truth_table) then "--truth-table '" + select_first([truth_table]) + "'" else ""}

    tar -czf ~{output_tarball_name} "${WORK_DIR}"
  >>>

  runtime {
    docker: sv_pipeline_docker
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    noAddress: true
  }

  output {
    File output_tarball = output_tarball_name
  }
}
