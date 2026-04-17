version 1.0

import "Structs.wdl"

workflow Ploidy {
  input {
    File merged_depth_file
    File reference_dict
    String batch
    Array[File] sparse_sd_files
    File? truth_json
    String? plot_highlight_sample
    String? preprocess_args
    String? model_args
    String? plot_args
    File? poor_regions
    Float min_poor_region_coverage = 0.5
    Boolean enable_ppd = false
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
      sparse_sd_files = sparse_sd_files,
      truth_json = truth_json,
      plot_highlight_sample = plot_highlight_sample,
      preprocess_args = preprocess_args,
      model_args = model_args,
      plot_args = plot_args,
      poor_regions = poor_regions,
      min_poor_region_coverage = min_poor_region_coverage,
      enable_ppd = enable_ppd,
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
    Array[File] sparse_sd_files
    File? truth_json
    String? plot_highlight_sample
    String? preprocess_args
    String? model_args
    String? plot_args
    File? poor_regions
    Float min_poor_region_coverage
    Boolean enable_ppd
    String sv_pipeline_qc_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 15,
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

    OUTDIR="./~{batch}_ploidy"
    mkdir -p "${OUTDIR}"

    # Build site-depth list if SD files were provided
    SD_ARGS=""
    if [[ ~{length(sparse_sd_files)} -gt 0 ]]; then
      SD_LIST="sparse_sd_files.list"
      for f in ~{sep=' ' sparse_sd_files}; do
        echo "$f" >> "${SD_LIST}"
      done
      SD_ARGS="--site-depth-list ${SD_LIST}"
    fi

    # Run packaged pipeline script which wraps the gatk-sv-ploidy subcommands
    /opt/gatk-sv-ploidy/run_ploidy.sh \
      --input-depth ~{ploidy_matrix} \
      --work-dir ${OUTDIR} \
      ${SD_ARGS} \
      ~{if defined(truth_json) then "--truth-json " + truth_json else ""} \
      ~{if defined(preprocess_args) then "--preprocess-args " + preprocess_args else ""} \
      ~{if defined(model_args) then "--model-args " + model_args else ""} \
      ~{if defined(plot_args) then "--plot-args " + plot_args else ""} \
      ~{if defined(poor_regions) then "--poor-regions " + poor_regions else ""} \
      --min-poor-region-coverage ~{min_poor_region_coverage} \
      ~{if enable_ppd then "--ppd" else ""} \
      ~{if defined(plot_highlight_sample) then "--highlight-sample " + plot_highlight_sample else ""}

    # Package all outputs
    tar -zcf ~{batch}_ploidy.tar.gz "${OUTDIR}"
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

