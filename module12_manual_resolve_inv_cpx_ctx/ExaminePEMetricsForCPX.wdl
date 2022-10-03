##########################
## EXPERIMENTAL WORKFLOW
##########################


version 1.0

import "Structs.wdl"
workflow CollectPEMetricsForCPX {
    input{
        Array[String] batch_name_list
        Array[File] PE_metrics
        Array[File] PE_metrics_idxes
        File PE_collect_script
        String prefix
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_override_collect_pe
        RuntimeAttr? runtime_attr_override_concat_evidence
        RuntimeAttr? runtime_attr_override_calcu_pe_stat
        }

    scatter (i in range(length(batch_name_list))){
        call CollectPEMetrics{
            input:
                batch_name = batch_name_list[i],
                PE_metric = PE_metrics[i],
                PE_metrics_idx = PE_metrics_idxes[i],
                PE_collect_script = PE_collect_script,
                sv_base_mini_docker = sv_base_mini_docker,
                runtime_attr_override = runtime_attr_override_collect_pe
        }
    }

    call ConcatEvidences{
        input:
            evidences = CollectPEMetrics.evidence,
            prefix = prefix,
            sv_base_mini_docker = sv_base_mini_docker,
            runtime_attr_override = runtime_attr_override_concat_evidence
    }

    call CalcuPEStat{
        input:
            evidence = ConcatEvidences.concat_evidence,
            prefix = prefix,
            sv_base_mini_docker = sv_base_mini_docker,
            runtime_attr_override = runtime_attr_override_calcu_pe_stat
    }

    output{
        File evidence = ConcatEvidences.concat_evidence
        File evi_stat = CalcuPEStat.evi_stat
    }
}



# collect PE metrics
task CollectPEMetrics{
  input{
    String batch_name
    File PE_metric
    File PE_metrics_idx
    File PE_collect_script
    String sv_base_mini_docker
   RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 5,
    disk_gb: ceil(10.0 + size(PE_metric, "GiB") * 2),
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail
    mkdir PE_metrics/
    gsutil cp ~{PE_metric} ./PE_metrics/
    gsutil cp ~{PE_metrics_idx} ./PE_metrics/
    grep -w ~{batch_name} ~{PE_collect_script} > tmp_metrics.sh
    bash tmp_metrics.sh
    mv out.evidence ~{batch_name}.evidence
    bgzip ~{batch_name}.evidence

  >>>

  output {
    File evidence = "~{batch_name}.evidence.gz"
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

task ConcatEvidences{
  input{
    Array[File] evidences
    String prefix
    String sv_base_mini_docker
   RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 5,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail

    while read SPLIT; do
      zcat $SPLIT
    done < ~{write_lines(evidences)} \
      | bgzip -c \
      > ~{prefix}.evidence.gz

  >>>

  output {
    File concat_evidence = "~{prefix}.evidence.gz"
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

task CalcuPEStat{
  input{
    File evidence
    String prefix
    String sv_base_mini_docker
   RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 5,
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




 
