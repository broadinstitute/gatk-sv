version 1.0

import "Structs.wdl"

workflow SVConcordance {
  input {
    # Vcfs must be formatted using FormatVcfForGatk (if unsure, check for ECN FORMAT field)
    # Per-contig vcfs, ordered to match contig_list
    Array[File] eval_vcfs
    Array[File] truth_vcfs
    String output_prefix

    File contig_list
    File reference_dict

    # Stratification parameters
    File? clustering_config
    File? stratification_config
    Array[String]? track_names
    Array[File]? track_intervals

    String gatk_docker

    Float? java_mem_fraction

    RuntimeAttr? runtime_attr_sv_concordance
  }

  Array[String] contigs = read_lines(contig_list)

  scatter (i in range(length(contigs))) {
    call SVConcordanceTask {
      input:
        eval_vcf=eval_vcfs[i],
        truth_vcf=truth_vcfs[i],
        output_prefix="~{output_prefix}.concordance.~{contigs[i]}",
        contig=contigs[i],
        clustering_config=clustering_config,
        stratification_config=stratification_config,
        track_names=track_names,
        track_intervals=track_intervals,
        reference_dict=reference_dict,
        java_mem_fraction=java_mem_fraction,
        gatk_docker=gatk_docker,
        runtime_attr_override=runtime_attr_sv_concordance
    }
  }

  output {
    Array[File] concordance_vcfs = SVConcordanceTask.out
    Array[File] concordance_vcf_indexes = SVConcordanceTask.out_index
  }
}

task SVConcordanceTask {
  input {
    File truth_vcf
    File eval_vcf
    String output_prefix
    File reference_dict
    String? contig
    File? clustering_config
    File? stratification_config
    Array[String]? track_names
    Array[File]? track_intervals
    String? additional_args

    Float? java_mem_fraction
    String gatk_docker
    RuntimeAttr? runtime_attr_override
  }

  parameter_meta {
    truth_vcf: {
                 localization_optional: true
               }
    eval_vcf:  {
                 localization_optional: true
               }
  }

  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 16,
                               disk_gb: ceil(10 + size(eval_vcf, "GB") * 4 + size(truth_vcf, "GB")),
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File out = "~{output_prefix}.vcf.gz"
    File out_index = "~{output_prefix}.vcf.gz.tbi"
  }
  command <<<
    set -euo pipefail

    function getJavaMem() {
    # get JVM memory in MiB by getting total memory from /proc/meminfo
    # and multiplying by java_mem_fraction
      cat /proc/meminfo \
        | awk -v MEM_FIELD="$1" '{
          f[substr($1, 1, length($1)-1)] = $2
        } END {
          printf "%dM", f[MEM_FIELD] * ~{default="0.85" java_mem_fraction} / 1024
        }'
    }
    JVM_MAX_MEM=$(getJavaMem MemTotal)
    echo "JVM memory: $JVM_MAX_MEM"

    gatk --java-options "-Xmx${JVM_MAX_MEM}" SVConcordance \
      ~{"-L " + contig} \
      --sequence-dictionary ~{reference_dict} \
      --eval ~{eval_vcf} \
      --truth ~{truth_vcf} \
      -O ~{output_prefix}.vcf.gz \
      ~{if defined(clustering_config) then "--clustering-config " + clustering_config else ""} \
      ~{if defined(stratification_config) then "--stratify-config "  + stratification_config else ""} \
      ~{if length(select_first([track_names, []])) > 0 then "--track-name" else ""} ~{sep=" --track-name " select_first([track_names, []])} \
      ~{if length(select_first([track_intervals, []])) > 0 then "--track-intervals" else ""} ~{sep=" --track-intervals " select_first([track_intervals, []])} \
      ~{additional_args}
  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: gatk_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}
