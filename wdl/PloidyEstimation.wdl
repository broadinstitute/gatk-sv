version 1.0

import "Structs.wdl"

workflow Ploidy {
  input {
    File merged_depth_file
    File bincov_matrix
    File reference_dict
    String batch
    String? plot_highlight_sample
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
      --max-interval-size ~{default=2000000 max_interval_size} --min-interval-size ~{default=1000000 min_interval_size}
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
    File ploidy_plots = "${batch}_ploidy_plots.tar.gz"
  }
  
  command <<<
    set -euo pipefail
    
    mkdir ploidy_est
    
    # Run aneuploidy detection
    python /opt/sv-pipeline/scripts/aneuploidy_pyro.py \
      --input ~{ploidy_matrix} \
      --output-dir ./ploidy_est \
      --max-iter 5000 \
      --early-stopping \
      --patience 50
    
    # Create file list for aggregation
    echo "./ploidy_est/chromosome_stats.tsv" > chromosome_stats_list.txt
    
    # Aggregate results
    python /opt/sv-pipeline/scripts/aggregate_ploidy_output.py \
      --input chromosome_stats_list.txt \
      --output-dir ./ploidy_est
    
    # Run CN denoising if bin_stats exists
    python /opt/sv-pipeline/02_evidence_assessment/estimated_CN_denoising.py \
      --binwise-copy-number ./ploidy_est/bin_stats.tsv.gz \
      --estimated-copy-number ./ploidy_est/chromosome_stats.tsv \
      --output-stats ./ploidy_est/cn_denoising_stats.tsv \
      --output-pdf ./ploidy_est/cn_denoising_plots.pdf
    
    # Package all outputs
    tar -zcf ./ploidy_est.tar.gz ./ploidy_est
    mv ploidy_est.tar.gz ~{batch}_ploidy_plots.tar.gz
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

