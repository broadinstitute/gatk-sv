version 1.0

import "Structs.wdl"
workflow CollectPEMetricsPerBatchCPX {
    input {
        Int n_per_split
        String batch_name
        File PE_metric
        File PE_metrics_idx
        File PE_collect_script
        String prefix
        String sv_base_mini_docker
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_collect_pe
        RuntimeAttr? runtime_attr_concat_evidence
        RuntimeAttr? runtime_attr_calcu_pe_stat
        RuntimeAttr? runtime_attr_split_script
        }

    call SplitScripts {
        input:
            script = PE_collect_script,
            n_per_split = n_per_split,
            batch_name = batch_name,
            sv_base_mini_docker = sv_base_mini_docker,
            runtime_attr_override = runtime_attr_split_script
    }

    scatter (script in SplitScripts.script_splits ) {
        call CollectPEMetrics {
            input:
                batch_name = batch_name,
                PE_metric = PE_metric,
                PE_metrics_idx = PE_metrics_idx,
                PE_collect_script = script,
                sv_pipeline_docker = sv_pipeline_docker,
                runtime_attr_override = runtime_attr_collect_pe
        }
    }

    call ConcatEvidences {
        input:
            evidences = CollectPEMetrics.evidence,
            prefix = prefix,
            sv_base_mini_docker = sv_base_mini_docker,
            runtime_attr_override = runtime_attr_concat_evidence
    }

    output {
        File evidence = ConcatEvidences.concat_evidence
    }
}



# collect PE metrics
task CollectPEMetrics {
  input {
    String batch_name
    File PE_metric
    File PE_metrics_idx
    File PE_collect_script
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 3.75,
    disk_gb: ceil(30.0 + size(PE_metric, "GiB") * 3),
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail
    mkdir PE_metrics/
    gsutil cp ~{PE_metric} ./
    gsutil cp ~{PE_metrics_idx} ./
    echo "Metrics succesfully downloaded."
    echo "Starting the evidence collection ..."

    if [ $(wc -c < "~{PE_collect_script}") -gt 0 ]; then
      bash ~{PE_collect_script}
    fi

    touch ~{batch_name}.evidence
    touch ~{batch_name}.0.PE_evidences
    for peEvFile in *.PE_evidences
    do
       cat ${peEvFile} >> ~{batch_name}.evidence
    done

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
    docker: sv_pipeline_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

task ConcatEvidences {
  input {
    Array[File] evidences
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

task SplitScripts {
  input {
    File script
    String batch_name
    Int n_per_split
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

  output {
    Array[File] script_splits = glob("collect_PE_evidences.*")
  }
  command <<<

    set -e pipefail

    awk '{if ($2~"~{batch_name}") print}' ~{script} > tmp_metrics.sh

    if [ $(wc -c < "tmp_metrics.sh") -gt 0 ]; then
      split --additional-suffix ".sh" -l ~{n_per_split} -a 6 tmp_metrics.sh collect_PE_evidences.
    else
      touch collect_PE_evidences.aaaaaa.sh
    fi

  >>>
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
