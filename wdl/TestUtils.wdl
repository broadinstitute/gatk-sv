version 1.0
import "Structs.wdl"

task CatMetrics {
  input {
    String prefix
    Array[File] metric_files
    String? search_string
    String? replace_string
    String linux_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr runtime_default = object {
    mem_gb: 1.0,
    disk_gb: 10,
    cpu_cores: 1,
    preemptible_tries: 3,
    max_retries: 1,
    boot_disk_gb: 10
  }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

  output {
    File out = "~{prefix}.metrics.tsv"
  }
  command <<<

    set -eu
    cat ~{sep=" " metric_files} > ~{prefix}.metrics.tsv
    if ~{defined(search_string) && defined(replace_string)}; then
      sed -i 's/~{search_string}/~{replace_string}/g' ~{prefix}.metrics.tsv
    fi

  >>>
  runtime {
    memory: select_first([runtime_override.mem_gb, runtime_default.mem_gb]) + " GB"
    disks: "local-disk " + select_first([runtime_override.disk_gb, runtime_default.disk_gb]) + " HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: linux_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    noAddress: true
  }
}

task StandardizeVCF {
  input {
    File vcf
    String sample_id
    String caller
    File contig_index
    Int min_size
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr runtime_default = object {
    mem_gb: 1.0,
    disk_gb: 10,
    cpu_cores: 1,
    preemptible_tries: 3,
    max_retries: 1,
    boot_disk_gb: 10
  }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  output {
    File out = "~{sample_id}.~{caller}.std.vcf.gz"
  }
  command <<<

    set -eu
    svtk standardize --min-size ~{min_size} --contigs ~{contig_index} ~{vcf} ~{sample_id}.~{caller}.std.vcf ~{caller}
    bgzip ~{sample_id}.~{caller}.std.vcf

  >>>
  runtime {
    memory: select_first([runtime_override.mem_gb, runtime_default.mem_gb]) + " GB"
    disks: "local-disk " + select_first([runtime_override.disk_gb, runtime_default.disk_gb]) + " HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    noAddress: true
  }
}

task VCFMetrics {
  input {
    File vcf
    File? baseline_vcf
    Array[String] samples
    String types
    String prefix
    File contig_list
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  File samples_list = write_lines(samples)

  RuntimeAttr runtime_default = object {
    mem_gb: 3.75,
    disk_gb: 10,
    cpu_cores: 1,
    preemptible_tries: 3,
    max_retries: 1,
    boot_disk_gb: 10
  }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

  output {
    File out = "~{prefix}.vcf.tsv"
  }
  command <<<

    svtest vcf \
     ~{vcf} \
     ~{contig_list} \
     ~{samples_list} \
     ~{types} \
     ~{prefix} \
     ~{if defined(baseline_vcf) then "--baseline-vcf " + baseline_vcf else ""} \
     > ~{prefix}.vcf.tsv

  >>>
  runtime {
    memory: select_first([runtime_override.mem_gb, runtime_default.mem_gb]) + " GB"
    disks: "local-disk " + select_first([runtime_override.disk_gb, runtime_default.disk_gb]) + " HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    noAddress: true
  }
}

task BAFMetrics {
  input {
    File baf_file
    Array[String] samples
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  String prefix = if length(samples) > 1 then "merged" else samples[0]
  File samples_list = write_lines(samples)

  RuntimeAttr runtime_default = object {
    mem_gb: 3.75,
    disk_gb: 10,
    cpu_cores: 1,
    preemptible_tries: 3,
    max_retries: 1,
    boot_disk_gb: 10
  }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

  output {
    File out = "~{prefix}.baf-file.tsv"
  }
  command <<<

    svtest baf-file ~{baf_file} ~{samples_list} > ~{prefix}.baf-file.tsv

  >>>
  runtime {
    memory: select_first([runtime_override.mem_gb, runtime_default.mem_gb]) + " GB"
    disks: "local-disk " + select_first([runtime_override.disk_gb, runtime_default.disk_gb]) + " HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    noAddress: true
  }
}

task SRMetrics {
  input {
    File sr_file
    Array[String] samples
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  String prefix = if length(samples) > 1 then "merged" else samples[0]
  File samples_list = write_lines(samples)

  RuntimeAttr runtime_default = object {
    mem_gb: 3.75,
    disk_gb: 25,
    cpu_cores: 1,
    preemptible_tries: 3,
    max_retries: 1,
    boot_disk_gb: 10
  }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

  output {
    File out = "~{prefix}.sr-file.tsv"
  }
  command <<<

    svtest sr-file ~{sr_file} ~{samples_list} > ~{prefix}.sr-file.tsv

  >>>
  runtime {
    memory: select_first([runtime_override.mem_gb, runtime_default.mem_gb]) + " GB"
    disks: "local-disk " + select_first([runtime_override.disk_gb, runtime_default.disk_gb]) + " HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    noAddress: true
  }
}

task PEMetrics {
  input {
    File pe_file
    Array[String] samples
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  String prefix = if length(samples) > 1 then "merged" else samples[0]
  File samples_list = write_lines(samples)

  RuntimeAttr runtime_default = object {
    mem_gb: 3.75,
    disk_gb: 10,
    cpu_cores: 1,
    preemptible_tries: 3,
    max_retries: 1,
    boot_disk_gb: 10
  }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

  output {
    File out = "~{prefix}.pe-file.tsv"
  }
  command <<<

    svtest pe-file ~{pe_file} ~{samples_list} > ~{prefix}.pe-file.tsv

  >>>
  runtime {
    memory: select_first([runtime_override.mem_gb, runtime_default.mem_gb]) + " GB"
    disks: "local-disk " + select_first([runtime_override.disk_gb, runtime_default.disk_gb]) + " HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    noAddress: true
  }
}

task CountsMetrics {
  input {
    File counts_file
    String sample_id
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr runtime_default = object {
    mem_gb: 3.75,
    disk_gb: 10,
    cpu_cores: 1,
    preemptible_tries: 3,
    max_retries: 1,
    boot_disk_gb: 10
  }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

  output {
    File out = "~{sample_id}.raw-counts.tsv"
  }
  command <<<

    svtest raw-counts ~{counts_file} ~{sample_id} > ~{sample_id}.raw-counts.tsv

  >>>
  runtime {
    memory: select_first([runtime_override.mem_gb, runtime_default.mem_gb]) + " GB"
    disks: "local-disk " + select_first([runtime_override.disk_gb, runtime_default.disk_gb]) + " HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    noAddress: true
  }
}

task BincovMetrics {
  input {
    File bincov_matrix
    Array[String] samples
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  File samples_list = write_lines(samples)
  String low_mem_mode_arg = if length(samples) > 10 then "--low-mem-mode" else ""

  RuntimeAttr runtime_default = object {
    mem_gb: 15.0,
    disk_gb: 10,
    cpu_cores: 1,
    preemptible_tries: 3,
    max_retries: 1,
    boot_disk_gb: 10
  }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

  output {
    File out = "bincov-matrix.tsv"
  }
  command <<<

    svtest bincov-matrix ~{bincov_matrix} ~{samples_list} ~{low_mem_mode_arg} > bincov-matrix.tsv

  >>>
  runtime {
    memory: select_first([runtime_override.mem_gb, runtime_default.mem_gb]) + " GB"
    disks: "local-disk " + select_first([runtime_override.disk_gb, runtime_default.disk_gb]) + " HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    noAddress: true
  }
}

task MedcovMetrics {
  input {
    File medcov_file
    Array[String] samples
    File? baseline_medcov_file
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  File samples_list = write_lines(samples)

  RuntimeAttr runtime_default = object {
    mem_gb: 1.0,
    disk_gb: 10,
    cpu_cores: 1,
    preemptible_tries: 3,
    max_retries: 1,
    boot_disk_gb: 10
  }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

  output {
    File out = "medcov.tsv"
  }
  command <<<

    svtest medcov \
      ~{medcov_file} \
      ~{samples_list} \
      ~{if defined(baseline_medcov_file) then "--baseline-file " + baseline_medcov_file else ""} \
      > medcov.tsv

  >>>
  runtime {
    memory: select_first([runtime_override.mem_gb, runtime_default.mem_gb]) + " GB"
    disks: "local-disk " + select_first([runtime_override.disk_gb, runtime_default.disk_gb]) + " HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    noAddress: true
  }
}

task MergedDepthMetricsWithBaseline {
  input {
    File bed
    File baseline_bed
    String type
    File contig_list
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr runtime_default = object {
    mem_gb: 3.75,
    disk_gb: 10,
    cpu_cores: 1,
    preemptible_tries: 3,
    max_retries: 1,
    boot_disk_gb: 10
  }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

  output {
    File out = "~{type}.merged-depth.tsv"
  }
  command <<<

    set -euo pipefail
    bedtools intersect -wa -u -f 0.5 -r -a ~{bed} -b ~{baseline_bed} | cut -f4 > overlap.test.list
    bedtools intersect -wa -u -f 0.5 -r -b ~{bed} -a ~{baseline_bed} | cut -f4 > overlap.base.list
    svtest merged-depth \
      ~{bed} \
      ~{contig_list} \
      ~{type} \
      --baseline-bed ~{baseline_bed} \
      --test-hits overlap.test.list \
      --baseline-hits overlap.base.list \
      > ~{type}.merged-depth.tsv

  >>>
  runtime {
    memory: select_first([runtime_override.mem_gb, runtime_default.mem_gb]) + " GB"
    disks: "local-disk " + select_first([runtime_override.disk_gb, runtime_default.disk_gb]) + " HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    noAddress: true
  }
}

task MergedDepthMetricsWithoutBaseline {
  input {
    File bed
    String type
    File contig_list
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr runtime_default = object {
    mem_gb: 3.75,
    disk_gb: 10,
    cpu_cores: 1,
    preemptible_tries: 3,
    max_retries: 1,
    boot_disk_gb: 10
  }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

  output {
    File out = "~{type}.merged-depth.tsv"
  }
  command <<<

    svtest merged-depth ~{bed} ~{contig_list} ~{type} > ~{type}.merged-depth.tsv

  >>>
  runtime {
    memory: select_first([runtime_override.mem_gb, runtime_default.mem_gb]) + " GB"
    disks: "local-disk " + select_first([runtime_override.disk_gb, runtime_default.disk_gb]) + " HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    noAddress: true
  }
}

task MetricsFileMetrics {
  input {
    File metrics_file
    File contig_list
    Boolean common
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  String common_arg = if common then "--common" else ""

  RuntimeAttr runtime_default = object {
    mem_gb: 3.75,
    disk_gb: 10,
    cpu_cores: 1,
    preemptible_tries: 3,
    max_retries: 1,
    boot_disk_gb: 10
  }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

  output {
    File out = "~{prefix}.metrics.tsv"
  }
  command <<<

    svtest metrics-file ~{metrics_file} ~{contig_list} ~{common_arg} > ~{prefix}.metrics.tsv

  >>>
  runtime {
    memory: select_first([runtime_override.mem_gb, runtime_default.mem_gb]) + " GB"
    disks: "local-disk " + select_first([runtime_override.disk_gb, runtime_default.disk_gb]) + " HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    noAddress: true
  }
}

task CutoffAndOutlierMetrics {
  input {
    File cutoffs
    File outlier_list
    File filtered_ped_file
    Array[String] samples
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  File samples_list = write_lines(samples)

  RuntimeAttr runtime_default = object {
    mem_gb: 1.0,
    disk_gb: 10,
    cpu_cores: 1,
    preemptible_tries: 3,
    max_retries: 1,
    boot_disk_gb: 10
  }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

  output {
    File out = "rf_cutoff.metrics.tsv"
  }
  command <<<

    set -euo pipefail
    svtest rf-cutoffs ~{cutoffs} > rf_cutoff.metrics.tsv
    svtest sample-list --prefix rf_outliers --valid-sample-list ~{samples_list} ~{outlier_list} >> rf_cutoff.metrics.tsv
    svtest ped-file ~{filtered_ped_file} >> rf_cutoff.metrics.tsv

  >>>
  runtime {
    memory: select_first([runtime_override.mem_gb, runtime_default.mem_gb]) + " GB"
    disks: "local-disk " + select_first([runtime_override.disk_gb, runtime_default.disk_gb]) + " HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    noAddress: true
  }
}

task GenotypingCutoffMetrics {
  input {
    File cutoffs
    String name
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr runtime_default = object {
    mem_gb: 1.0,
    disk_gb: 10,
    cpu_cores: 1,
    preemptible_tries: 3,
    max_retries: 1,
    boot_disk_gb: 10
  }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  output {
    File out = "cutoffs.~{name}.metrics.tsv"
  }
  command <<<

    svtest gt-cutoffs ~{cutoffs} ~{name} > cutoffs.~{name}.metrics.tsv

  >>>
  runtime {
    memory: select_first([runtime_override.mem_gb, runtime_default.mem_gb]) + " GB"
    disks: "local-disk " + select_first([runtime_override.disk_gb, runtime_default.disk_gb]) + " HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    noAddress: true
  }
}

task IdListMetrics {
  input {
    File id_list
    String name
    String linux_docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr runtime_default = object {
    mem_gb: 1.0,
    disk_gb: 10,
    cpu_cores: 1,
    preemptible_tries: 3,
    max_retries: 1,
    boot_disk_gb: 10
  }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  output {
    File out = "~{name}.id_list_metrics.tsv"
  }
  command <<<

    set -euo pipefail
    COUNT=$(cat ~{id_list} | wc -l)
    echo "~{name}	${COUNT}" > ~{name}.id_list_metrics.tsv

  >>>
  runtime {
    memory: select_first([runtime_override.mem_gb, runtime_default.mem_gb]) + " GB"
    disks: "local-disk " + select_first([runtime_override.disk_gb, runtime_default.disk_gb]) + " HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: linux_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    noAddress: true
  }
}

task PlotMetrics {
  input {
    String name
    File test_metrics
    File base_metrics
    Array[String] samples
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  File sample_list = write_lines(samples)

  RuntimeAttr runtime_default = object {
    mem_gb: 1.0,
    disk_gb: 10,
    cpu_cores: 1,
    preemptible_tries: 3,
    max_retries: 1,
    boot_disk_gb: 10
  }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

  output {
    File metrics_plot_pdf = "~{name}.plot_metrics.pdf"
    File metrics_plot_tsv = "~{name}.plot_metrics.tsv"
  }
  command <<<

    svtest plot-metrics --sample-list ~{sample_list} --metrics-out ~{name}.plot_metrics.tsv ~{test_metrics} ~{base_metrics} ~{name}.plot_metrics.pdf

  >>>
  runtime {
    memory: select_first([runtime_override.mem_gb, runtime_default.mem_gb]) + " GB"
    disks: "local-disk " + select_first([runtime_override.disk_gb, runtime_default.disk_gb]) + " HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    noAddress: true
  }
}