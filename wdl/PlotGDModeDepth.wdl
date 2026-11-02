version 1.0

import "CollectCoverage.wdl" as cov
import "Structs.wdl"

# Most-common normalized read depth per sample at every autosomal NAHR-mediated
# GD locus (chrX and chrY loci are excluded).
#
# Inputs are per-batch coverage matrices, parallel to batches: give either
# low_res_rd_matrices (used as is) or high_res_rd_matrices (condensed to low
# resolution with CondenseReadCounts using reference_dict and
# rebinned_interval_size; the low-res matrices are workflow outputs). Counts
# are normalized as plot.py normalizes its raw-depth overlay: each sample's
# counts are scaled so its autosomal median is CN 2. For each sample and locus,
# the depth of each bin between the locus breakpoints is rounded to the nearest
# 0.1, the most commonly observed rounded value is found, and the median of the
# actual depths that rounded to it is reported. Expected values: non-carriers
# ~2, germline het DEL carriers ~1 (het DUP ~3), mosaic carriers in between,
# and carriers of an event covering only part of the locus ~2 or ~carrier
# value.
#
# Preprocessing always runs and outputs the table needed to plot
# (mode_depth_table, one row per batch/locus/sample), plus per-batch tables.
# The PDF is drawn from that table only when make_plot is true. The plotting
# script accepts several preprocessed tables, stacked with a shared x-axis:
# add them with additional_mode_depth_tables / additional_dataset_labels.
#
# Loci are ordered as in gd_table (NAHR == "yes" only), with DEL and DUP at one
# locus sharing an x position. Only loci with probable mosaic carriers (median
# depth between 1.2 and 1.8, or between 2.2 and 2.8, by default) get an axis
# label.
#
# plot_mode_depth_script is plot_mode_depth.py from gatk-sv-gd; the gatk_sv_gd
# package must be installed in the image.
workflow PlotGDModeDepth {
  input {
    Array[String] batches
    Array[File]? low_res_rd_matrices
    Array[File]? high_res_rd_matrices
    File? reference_dict
    Int rebinned_interval_size = 10000
    File gd_table

    Boolean make_plot = true
    String dataset_label = "All batches"
    Array[File] additional_mode_depth_tables = []
    Array[String] additional_dataset_labels = []
    Int min_mosaic_samples = 1
    Float mosaic_del_min = 1.2
    Float mosaic_del_max = 1.8
    Float mosaic_dup_min = 2.2
    Float mosaic_dup_max = 2.8
    Float fig_width_in = 3.5
    Float fig_height_in = 1.8
    Float depth_axis_min = 0.0
    Float depth_axis_max = 4.0

    File plot_mode_depth_script

    String output_prefix = "gd_mode_depth"
    String sv_pipeline_docker
    String? gatk_docker
    RuntimeAttr? runtime_attr_condense
    RuntimeAttr? runtime_attr_summarize
    RuntimeAttr? runtime_attr_merge
    RuntimeAttr? runtime_attr_plot
  }

  Boolean use_low_res = defined(low_res_rd_matrices)
  Array[File] low_res_inputs = select_first([low_res_rd_matrices, []])
  Array[File] high_res_inputs = select_first([high_res_rd_matrices, []])

  scatter (i in range(length(batches))) {
    if (!use_low_res) {
      call cov.CondenseReadCounts {
        input:
          counts = high_res_inputs[i],
          output_prefix = batches[i] + ".low_res",
          sequence_dictionary = reference_dict,
          min_interval_size = rebinned_interval_size,
          max_interval_size = rebinned_interval_size,
          gatk_docker = select_first([gatk_docker]),
          runtime_attr_override = runtime_attr_condense
      }
    }

    File low_res_counts = if use_low_res then low_res_inputs[i] else select_first([CondenseReadCounts.out])

    call SummarizeBatchModeDepth {
      input:
        batch = batches[i],
        low_res_counts = low_res_counts,
        gd_table = gd_table,
        plot_mode_depth_script = plot_mode_depth_script,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_summarize
    }
  }

  call MergeModeDepth {
    input:
      tables = SummarizeBatchModeDepth.mode_depth,
      output_prefix = output_prefix,
      plot_mode_depth_script = plot_mode_depth_script,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_merge
  }

  if (make_plot) {
    call PlotModeDepth {
      input:
        mode_depth_tables = flatten([[MergeModeDepth.merged], additional_mode_depth_tables]),
        dataset_labels = flatten([[dataset_label], additional_dataset_labels]),
        gd_table = gd_table,
        min_mosaic_samples = min_mosaic_samples,
        mosaic_del_min = mosaic_del_min,
        mosaic_del_max = mosaic_del_max,
        mosaic_dup_min = mosaic_dup_min,
        mosaic_dup_max = mosaic_dup_max,
        fig_width_in = fig_width_in,
        fig_height_in = fig_height_in,
        depth_axis_min = depth_axis_min,
        depth_axis_max = depth_axis_max,
        output_prefix = output_prefix,
        plot_mode_depth_script = plot_mode_depth_script,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_plot
    }
  }

  output {
    File mode_depth_table = MergeModeDepth.merged
    Array[File] batch_mode_depth_tables = SummarizeBatchModeDepth.mode_depth
    # Low-res matrices made from high_res_rd_matrices (empty if low-res was given).
    Array[File] condensed_low_res_matrices = select_all(CondenseReadCounts.out)
    File? mode_depth_plot = PlotModeDepth.pdf
  }
}

# Normalizes one batch's low-res counts as plot.py does and reports the
# most-common depth (median of the values in the modal 0.1 bin) per sample and
# autosomal NAHR locus.
task SummarizeBatchModeDepth {
  input {
    String batch
    File low_res_counts
    File gd_table
    File plot_mode_depth_script
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Int disk_gb = ceil(20 + 3 * size(low_res_counts, "GiB"))

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 32.0,
    disk_gb: disk_gb,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail

    python ~{plot_mode_depth_script} summarize \
      --batch '~{batch}' \
      --gd-table ~{gd_table} \
      --low-res-counts ~{low_res_counts} \
      --out '~{batch}.mode_depth.tsv.gz'
  >>>

  runtime {
    docker: sv_pipeline_docker
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

  output {
    File mode_depth = batch + ".mode_depth.tsv.gz"
  }
}

task MergeModeDepth {
  input {
    Array[File] tables
    String output_prefix
    File plot_mode_depth_script
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 4.0,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail

    python ~{plot_mode_depth_script} merge \
      --inputs ~{sep=" " tables} \
      --out '~{output_prefix}.mode_depth.tsv.gz'
  >>>

  runtime {
    docker: sv_pipeline_docker
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

  output {
    File merged = output_prefix + ".mode_depth.tsv.gz"
  }
}

task PlotModeDepth {
  input {
    Array[File] mode_depth_tables
    Array[String] dataset_labels
    File gd_table
    Int min_mosaic_samples
    Float mosaic_del_min
    Float mosaic_del_max
    Float mosaic_dup_min
    Float mosaic_dup_max
    Float fig_width_in
    Float fig_height_in
    Float depth_axis_min
    Float depth_axis_max
    String output_prefix
    File plot_mode_depth_script
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 4.0,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail

    python ~{plot_mode_depth_script} plot \
      --gd-table ~{gd_table} \
      --data ~{sep=" " mode_depth_tables} \
      --labels ~{sep=" " dataset_labels} \
      --min-mosaic-samples ~{min_mosaic_samples} \
      --mosaic-del-min ~{mosaic_del_min} \
      --mosaic-del-max ~{mosaic_del_max} \
      --mosaic-dup-min ~{mosaic_dup_min} \
      --mosaic-dup-max ~{mosaic_dup_max} \
      --width ~{fig_width_in} \
      --height ~{fig_height_in} \
      --ymin ~{depth_axis_min} \
      --ymax ~{depth_axis_max} \
      --out-pdf '~{output_prefix}.pdf'
  >>>

  runtime {
    docker: sv_pipeline_docker
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

  output {
    File pdf = output_prefix + ".pdf"
  }
}
