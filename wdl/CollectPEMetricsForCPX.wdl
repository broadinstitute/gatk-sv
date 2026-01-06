version 1.0

import "Structs.wdl"
import "CollectPEMetricsPerBatchCPX.wdl" as per_batch

workflow CollectPEMetricsForCPX {
    input {
        Array[String] batch_name_list
        Array[File] PE_metrics
        Array[File] PE_metrics_idxes
        File PE_collect_script
        String prefix
        Int n_per_split
        String sv_pipeline_docker
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_collect_pe
        RuntimeAttr? runtime_attr_split_script
        RuntimeAttr? runtime_attr_calcu_pe_stat
        RuntimeAttr? runtime_attr_concat_evidence
        }

    scatter (i in range(length(batch_name_list))) {
        call per_batch.CollectPEMetricsPerBatchCPX {
            input:
                n_per_split = n_per_split,
                prefix = "~{prefix}.~{batch_name_list[i]}",
                batch_name = batch_name_list[i],
                PE_metric = PE_metrics[i],
                PE_metrics_idx = PE_metrics_idxes[i],
                PE_collect_script = PE_collect_script,
                sv_pipeline_docker = sv_pipeline_docker,
                sv_base_mini_docker = sv_base_mini_docker,
                runtime_attr_collect_pe = runtime_attr_collect_pe,
                runtime_attr_split_script = runtime_attr_split_script,
                runtime_attr_calcu_pe_stat = runtime_attr_calcu_pe_stat,
                runtime_attr_concat_evidence = runtime_attr_concat_evidence
        }
     }

    call per_batch.ConcatEvidences {
        input:
            evidences = CollectPEMetricsPerBatchCPX.evidence,
            prefix = prefix,
            sv_base_mini_docker = sv_base_mini_docker,
            runtime_attr_override = runtime_attr_concat_evidence
    }

    call CalcuPEStat {
        input:
            evidence = ConcatEvidences.concat_evidence,
            prefix = prefix,
            sv_base_mini_docker = sv_base_mini_docker,
            runtime_attr_override = runtime_attr_calcu_pe_stat
    }

    output {
        File evidence = ConcatEvidences.concat_evidence
        File evi_stat = CalcuPEStat.evi_stat
    }
}


task CalcuPEStat {
  input {
    File evidence
    String prefix
    String sv_base_mini_docker
   RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 3.75,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail

    zcat ~{evidence} | cut -f3,6- | uniq -c > ~{prefix}.evi_stat
    bgzip ~{prefix}.evi_stat
  >>>

  output {
    File evi_stat = "~{prefix}.evi_stat.gz"
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
