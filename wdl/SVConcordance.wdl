version 1.0

import "Structs.wdl"
import "TasksMakeCohortVcf.wdl" as tasks_cohort

workflow SVConcordance {
  input {
    # Vcfs must be formatted using FormatVcfForGatk (if unsure, check for ECN FORMAT field)
    File eval_vcf
    File truth_vcf
    String output_prefix

    File reference_dict

    # Records per eval_vcf shard; set to -1 to skip sharding
    Int n_per_split = -1

    # Stratification parameters
    File? clustering_config
    File? stratification_config
    Array[String]? track_names
    Array[File]? track_intervals

    String? additional_args

    String gatk_docker
    String sv_base_mini_docker

    Float? java_mem_fraction

    RuntimeAttr? runtime_attr_scatter_vcf
    RuntimeAttr? runtime_attr_sv_concordance
    RuntimeAttr? runtime_attr_concat
  }

  File truth_vcf_index = truth_vcf + ".tbi"

  if (n_per_split != -1) {
    call tasks_cohort.ScatterVcf {
      input:
        vcf=eval_vcf,
        prefix=output_prefix,
        records_per_shard=n_per_split,
        index_shards=true,
        sv_pipeline_docker=sv_base_mini_docker,
        runtime_attr_override=runtime_attr_scatter_vcf
    }
  }

  Array[File] eval_vcf_shards = select_first([ScatterVcf.shards, [eval_vcf]])
  Array[File] eval_vcf_shard_indexes = select_first([ScatterVcf.shard_indexes, [eval_vcf + ".tbi"]])

  scatter (i in range(length(eval_vcf_shards))) {
    call SVConcordanceTask {
      input:
        eval_vcf=eval_vcf_shards[i],
        eval_vcf_index=eval_vcf_shard_indexes[i],
        truth_vcf=truth_vcf,
        truth_vcf_index=truth_vcf_index,
        output_prefix="~{output_prefix}.concordance.shard_~{i}",
        clustering_config=clustering_config,
        stratification_config=stratification_config,
        track_names=track_names,
        track_intervals=track_intervals,
        additional_args=additional_args,
        reference_dict=reference_dict,
        java_mem_fraction=java_mem_fraction,
        gatk_docker=gatk_docker,
        runtime_attr_override=runtime_attr_sv_concordance
    }
  }

  call tasks_cohort.ConcatVcfs {
    input:
      vcfs=SVConcordanceTask.out,
      vcfs_idx=SVConcordanceTask.out_index,
      naive=true,
      outfile_prefix="~{output_prefix}.concordance",
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_attr_concat
  }

  output {
    File concordance_vcf = ConcatVcfs.concat_vcf
    File concordance_vcf_index = ConcatVcfs.concat_vcf_idx
  }
}

task SVConcordanceTask {
  input {
    File truth_vcf
    File truth_vcf_index
    File eval_vcf
    File? eval_vcf_index
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

    # The index may be localized to a different directory than the vcf (e.g. separate glob outputs),
    # so link both into the working directory where GATK can find the index
    EVAL_VCF="~{eval_vcf}"
    if ~{defined(eval_vcf_index)}; then
      if [ "$(basename "~{eval_vcf_index}")" != "$(basename "~{eval_vcf}").tbi" ]; then
        echo "eval_vcf_index ~{eval_vcf_index} does not match eval_vcf ~{eval_vcf}" >&2
        exit 1
      fi
      EVAL_VCF="eval.$(basename "~{eval_vcf}")"
      ln -s "~{eval_vcf}" "$EVAL_VCF"
      ln -s "~{eval_vcf_index}" "$EVAL_VCF.tbi"
    fi

    gatk --java-options "-Xmx${JVM_MAX_MEM}" SVConcordance \
      --sequence-dictionary ~{reference_dict} \
      --eval "$EVAL_VCF" \
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
    memory: "640 GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " SSD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: gatk_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}
