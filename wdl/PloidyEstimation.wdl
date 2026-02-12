version 1.0

import "Structs.wdl"

workflow Ploidy {
  input {
    File merged_depth_file
    File bincov_matrix
    File reference_dict
    String batch
    String? plot_highlight_sample
    String? model_args
    String? plot_args
    String gatk_docker
    String sv_pipeline_qc_docker
    RuntimeAttr? runtime_attr_score
    RuntimeAttr? runtime_attr_build
  }

  call CondenseDepthMatrix {
    input:
      merged_depth_file = merged_depth_file,
      merged_depth_file_index = merged_depth_file + ".tbi",
      prefix = "~{batch}_condensed_depth",
      reference_dict = reference_dict,
      gatk_docker = gatk_docker,
      runtime_attr_override = runtime_attr_build
  }

  call PloidyScore {
    input:
      ploidy_matrix = CondenseDepthMatrix.out,
      batch = batch,
      plot_highlight_sample = plot_highlight_sample,
      model_args = model_args,
      plot_args = plot_args,
      sv_pipeline_qc_docker = sv_pipeline_qc_docker,
      runtime_attr_override = runtime_attr_score
  }

  output {
    File ploidy_matrix = CondenseDepthMatrix.out
    File ploidy_matrix_index = CondenseDepthMatrix.out_index
    File ploidy_plots = PloidyScore.ploidy_plots
  }
}

task CondenseDepthMatrix {
  input {
    File merged_depth_file
    File merged_depth_file_index
    File reference_dict
    String prefix
    Int? max_interval_size
    Int? min_interval_size

    # Runtime parameters
    Float? java_mem_fraction
    String gatk_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 3.0,
                               disk_gb: 10,
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

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

    gatk --java-options "-Xmx${JVM_MAX_MEM}" CondenseDepthEvidence -F ~{merged_depth_file} -O ~{prefix}.rd.txt.gz --sequence-dictionary ~{reference_dict} \
      --max-interval-size ~{default=1000000 max_interval_size} --min-interval-size ~{default=1000000 min_interval_size}
  >>>

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: gatk_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    noAddress: true
  }

  output {
    File out = "~{prefix}.rd.txt.gz"
    File out_index = "~{prefix}.rd.txt.gz.tbi"
  }
}

task PloidyScore {
  input {
    File ploidy_matrix
    String batch
    String? plot_highlight_sample
    String? model_args
    String? plot_args
    String sv_pipeline_qc_docker
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
    File ploidy_plots = "~{batch}_ploidy.tar.gz"
  }
  
  command <<<
    set -euo pipefail

    mkdir
    mkdir "~{batch}_ploidy" "~{batch}_ploidy/model" "~{batch}_ploidy/results"
    
    # Run aneuploidy detection
    python /opt/sv-pipeline/scripts/aneuploidy_pyro.py \
      --input ~{ploidy_matrix} \
      --output-dir ./model \
       ~{model_args}
    
    # Aggregate results
    python /opt/sv-pipeline/scripts/aggregate_ploidy_output.py \
      --chrom-stats ./model/chromosome_stats.tsv \
      --bin-stats model/bin_stats.tsv.gz \
      --training-loss model/training_loss.tsv \
      --output-dir ./results \
      ~{"--highlight-sample " + plot_highlight_sample} \
      ~{plot_args}

    # Package all outputs
    tar -zcf ~{batch}_ploidy.tar.gz ./~{batch}_ploidy/model ./~{batch}_ploidy/results
  >>>
  
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_qc_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    noAddress: true
  }
}

