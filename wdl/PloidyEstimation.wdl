version 1.0

import "Structs.wdl"
import "SubsetSD.wdl" as subset

workflow Ploidy {
  input {
    File merged_depth_file
    File reference_dict
    String batch
    Array[File] sd_files = []
    File? ploidy_sd_locs_vcf
    File? poor_regions
    Int subset_sd_stride = 10
    File? truth_json
    String? preprocess_args
    String? polyploidy_args
    String? infer_args
    String? ppd_args
    String? call_args
    String? plot_args
    Boolean enable_ppd = false
    Boolean use_callq20 = false
    String gatk_docker
    String sv_pipeline_docker
    String sv_pipeline_qc_docker
    RuntimeAttr? runtime_attr_subset_sd
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

  if (length(sd_files) > 0) {
    scatter (i in range(length(sd_files))) {
      call subset.SubsetSDTask as SubsetSDForPloidy {
        input:
          sd_file = sd_files[i],
          sites_vcf = select_first([ploidy_sd_locs_vcf]),
          prefix = basename(sd_files[i], ".sd.txt.gz"),
          stride = subset_sd_stride,
          sv_pipeline_docker = sv_pipeline_docker,
          runtime_attr_override = runtime_attr_subset_sd
      }
    }
  }

  Array[File] subset_sd_files = select_first([SubsetSDForPloidy.subset_sd_file, []])

  call PloidyScore {
    input:
      ploidy_matrix = CondenseDepthMatrix.out,
      batch = batch,
      subset_sd_files = subset_sd_files,
      poor_regions = poor_regions,
      truth_json = truth_json,
      preprocess_args = preprocess_args,
      polyploidy_args = polyploidy_args,
      infer_args = infer_args,
      ppd_args = ppd_args,
      call_args = call_args,
      plot_args = plot_args,
      enable_ppd = enable_ppd,
      use_callq20 = use_callq20,
      sv_pipeline_qc_docker = sv_pipeline_qc_docker,
      runtime_attr_override = runtime_attr_score
  }

  output {
    File ploidy_matrix = CondenseDepthMatrix.out
    File ploidy_matrix_index = CondenseDepthMatrix.out_index
    File chromosome_stats = PloidyScore.chromosome_stats
    File bin_stats = PloidyScore.bin_stats
    File sample_sex_assignments = PloidyScore.sample_sex_assignments
    File aneuploidy_type_predictions = PloidyScore.aneuploidy_type_predictions
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
    Array[File] subset_sd_files
    File? poor_regions
    File? truth_json
    String? preprocess_args
    String? polyploidy_args
    String? infer_args
    String? ppd_args
    String? call_args
    String? plot_args
    Boolean enable_ppd
    Boolean use_callq20
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
    File chromosome_stats = "~{batch}_ploidy/infer/chromosome_stats.tsv"
    File bin_stats = "~{batch}_ploidy/infer/bin_stats.tsv.gz"
    File sample_sex_assignments = "~{batch}_ploidy/call/sex_assignments.txt.gz"
    File aneuploidy_type_predictions = "~{batch}_ploidy/call/aneuploidy_type_predictions.tsv"
    File ploidy_plots = "~{batch}_ploidy.tar.gz"
  }
  
  command <<<
    set -euo pipefail

    OUTDIR="./~{batch}_ploidy"
    mkdir -p "${OUTDIR}"

    # Build site-depth list if SD files were provided
    SD_ARGS=""
    if [[ ~{length(subset_sd_files)} -gt 0 ]]; then
      SD_LIST="subset_sd_files.list"
      for f in ~{sep=' ' subset_sd_files}; do
        echo "$f" >> "${SD_LIST}"
      done
      SD_ARGS="--site-depth-list ${SD_LIST}"
    fi

    # Run packaged pipeline script which wraps the gatk-sv-ploidy subcommands
    /opt/gatk-sv-ploidy/run_ploidy.sh \
      --input-depth ~{ploidy_matrix} \
      --work-dir ${OUTDIR} \
      ${SD_ARGS} \
      ~{if defined(poor_regions) then "--poor-regions " + poor_regions else ""} \
      ~{if defined(truth_json) then "--truth-json " + truth_json else ""} \
      ~{if defined(preprocess_args) then "--preprocess-args '" + preprocess_args + "'" else ""} \
      ~{if defined(polyploidy_args) then "--polyploidy-args '" + polyploidy_args + "'" else ""} \
      ~{if defined(infer_args) then "--infer-args '" + infer_args + "'" else ""} \
      ~{if defined(ppd_args) then "--ppd-args '" + ppd_args + "'" else ""} \
      ~{if defined(call_args) then "--call-args '" + call_args + "'" else ""} \
      ~{if defined(plot_args) then "--plot-args '" + plot_args + "'" else ""} \
      ~{if enable_ppd then "--ppd" else ""} \
      ~{if use_callq20 then "--use-callq20" else ""}

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

