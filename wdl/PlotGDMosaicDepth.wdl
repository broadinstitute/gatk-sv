version 1.0

import "CollectCoverage.wdl" as cov
import "Structs.wdl"

# Publication panel of normalized read depth across a GD locus (or its whole
# cluster), comparing mosaic carriers with random germline carriers and
# non-carriers drawn from the same batches.
#
# Inputs are per-batch outputs of CallGenomicDisorderCNVs. Arrays are indexed
# together: batches[i] goes with gd_output_tarballs[i] and high_res_rd_matrices[i].
# Use the same gd_table, reference_dict, rebinned_interval_size and
# high_res_rd_matrix given to CallGenomicDisorderCNVs so depth is normalized
# the way plot.py normalizes it. Only the condensed low-res counts are plotted;
# the high-res matrix is only the source for condensing.
#
# plot_mosaic_script, if given, replaces `gatk-sv-gd plot-mosaic` in the image
# with `python3 <script>` (development override; the gatk_sv_gd package must
# still be installed in the image).
workflow PlotGDMosaicDepth {
  input {
    String gd_id
    Array[String] mosaic_samples

    Array[String] batches
    Array[File] gd_output_tarballs
    Array[File] high_res_rd_matrices

    File gd_table
    File reference_dict
    File? segdup_bed
    Int rebinned_interval_size = 10000

    Int n_germline_carriers = 10
    Int n_non_carriers = 50
    Int random_seed = 42

    Float fig_width_in = 2.0
    Float fig_height_in = 2.1
    Float flank_scale = 0.15
    Float depth_axis_min = 0.0
    Float depth_axis_max = 4.0
    String? plot_title
    File plot_mosaic_script

    String output_prefix = "gd_mosaic_depth"
    String sv_pipeline_docker
    String gatk_docker
    RuntimeAttr? runtime_attr_condense
    RuntimeAttr? runtime_attr_classify
    RuntimeAttr? runtime_attr_select
    RuntimeAttr? runtime_attr_extract
    RuntimeAttr? runtime_attr_plot
  }

  scatter (i in range(length(batches))) {
    call ClassifyBatch {
      input:
        batch = batches[i],
        gd_output_tarball = gd_output_tarballs[i],
        gd_table = gd_table,
        gd_id = gd_id,
        plot_mosaic_script = plot_mosaic_script,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_classify
    }
  }

  call SelectSamples {
    input:
      mosaic_samples = mosaic_samples,
      classes = ClassifyBatch.classes,
      regions = ClassifyBatch.region,
      n_germline_carriers = n_germline_carriers,
      n_non_carriers = n_non_carriers,
      random_seed = random_seed,
      output_prefix = output_prefix,
      plot_mosaic_script = plot_mosaic_script,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_select
  }

  Array[String] has_mosaic = read_lines(SelectSamples.has_mosaic_flags)

  scatter (i in range(length(batches))) {
    if (has_mosaic[i] == "true") {
      call cov.CondenseReadCounts {
        input:
          counts = high_res_rd_matrices[i],
          output_prefix = batches[i] + ".low_res",
          sequence_dictionary = reference_dict,
          min_interval_size = rebinned_interval_size,
          max_interval_size = rebinned_interval_size,
          gatk_docker = gatk_docker,
          runtime_attr_override = runtime_attr_condense
      }

      call ExtractBatchDepth {
        input:
          batch = batches[i],
          low_res_counts = CondenseReadCounts.out,
          selected = SelectSamples.selected,
          region = SelectSamples.region,
          gd_table = gd_table,
          gd_id = gd_id,
          plot_mosaic_script = plot_mosaic_script,
          sv_pipeline_docker = sv_pipeline_docker,
          runtime_attr_override = runtime_attr_extract
      }
    }
  }

  call PlotMosaicDepth {
    input:
      depth_tables = select_all(ExtractBatchDepth.depth),
      selected = SelectSamples.selected,
      region = SelectSamples.region,
      gd_table = gd_table,
      gd_id = gd_id,
      segdup_bed = segdup_bed,
      fig_width_in = fig_width_in,
      fig_height_in = fig_height_in,
      flank_scale = flank_scale,
      depth_axis_min = depth_axis_min,
      depth_axis_max = depth_axis_max,
      plot_title = plot_title,
      output_prefix = output_prefix,
      plot_mosaic_script = plot_mosaic_script,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_plot
  }

  output {
    File mosaic_depth_plot = PlotMosaicDepth.pdf
    File selected_samples = SelectSamples.selected
  }
}

# Labels every sample at the locus's cluster as germline carrier of gd_id,
# non-carrier of all GDs in the cluster, or other.
task ClassifyBatch {
  input {
    String batch
    File gd_output_tarball
    File gd_table
    String gd_id
    File plot_mosaic_script
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Int disk_gb = ceil(20 + 3 * size(gd_output_tarball, "GiB"))

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 8.0,
    disk_gb: disk_gb,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail

    tar -xzf ~{gd_output_tarball} \
      gd_outputs/call/gd_cnv_calls.tsv.gz \
      gd_outputs/infer/cn_posteriors.tsv.gz
    PLOIDY_ARG=""
    if tar -xzf ~{gd_output_tarball} gd_outputs/preprocess/ploidy_estimates.tsv 2>/dev/null; then
      PLOIDY_ARG="--ploidy-table gd_outputs/preprocess/ploidy_estimates.tsv"
    fi

    python ~{plot_mosaic_script} classify \
      --batch '~{batch}' \
      --gd-table ~{gd_table} \
      --gd-id '~{gd_id}' \
      --calls gd_outputs/call/gd_cnv_calls.tsv.gz \
      --cn-posteriors gd_outputs/infer/cn_posteriors.tsv.gz \
      ${PLOIDY_ARG} \
      --out-prefix '~{batch}'
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
    File classes = batch + ".classes.tsv"
    File region = batch + ".region.json"
  }
}

# Finds the batches holding the mosaic carriers, then draws random germline
# carriers and non-carriers from those batches only.
task SelectSamples {
  input {
    Array[String] mosaic_samples
    Array[File] classes
    Array[File] regions
    Int n_germline_carriers
    Int n_non_carriers
    Int random_seed
    String output_prefix
    File plot_mosaic_script
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

    python ~{plot_mosaic_script} select \
      --classes ~{sep=" " classes} \
      --regions ~{sep=" " regions} \
      --mosaic-samples ~{write_lines(mosaic_samples)} \
      --n-germline ~{n_germline_carriers} \
      --n-non-carrier ~{n_non_carriers} \
      --seed ~{random_seed} \
      --out-prefix '~{output_prefix}'
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
    File selected = output_prefix + ".selected.tsv"
    File region = output_prefix + ".region.json"
    # One "true"/"false" per entry of regions, in the same order.
    File has_mosaic_flags = output_prefix + ".has_mosaic.txt"
  }
}

# Normalized low-res depth for the selected samples of one batch, built the
# way plot.py builds its raw-normalized column.
task ExtractBatchDepth {
  input {
    String batch
    File low_res_counts
    File selected
    File region
    File gd_table
    String gd_id
    File plot_mosaic_script
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Int disk_gb = ceil(20 + 2 * size(low_res_counts, "GiB"))

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 16.0,
    disk_gb: disk_gb,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail

    python ~{plot_mosaic_script} extract \
      --batch '~{batch}' \
      --gd-table ~{gd_table} \
      --gd-id '~{gd_id}' \
      --selected ~{selected} \
      --region ~{region} \
      --low-res-counts ~{low_res_counts} \
      --out-prefix '~{batch}'
  >>>

  runtime {
    docker: sv_pipeline_docker
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " SSD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

  output {
    File depth = batch + ".depth.tsv.gz"
  }
}

task PlotMosaicDepth {
  input {
    Array[File] depth_tables
    File selected
    File region
    File gd_table
    String gd_id
    File plot_mosaic_script
    File? segdup_bed
    Float fig_width_in
    Float fig_height_in
    Float flank_scale
    Float depth_axis_min
    Float depth_axis_max
    String? plot_title
    String output_prefix
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

    python ~{plot_mosaic_script} plot \
      --gd-table ~{gd_table} \
      --gd-id '~{gd_id}' \
      --selected ~{selected} \
      --region ~{region} \
      --depth ~{sep=" " depth_tables} \
      ~{"--segdup-bed " + segdup_bed} \
      --width ~{fig_width_in} \
      --height ~{fig_height_in} \
      --flank-scale ~{flank_scale} \
      --ymin ~{depth_axis_min} \
      --ymax ~{depth_axis_max} \
      ~{if defined(plot_title) then "--title '" + select_first([plot_title, ""]) + "'" else ""} \
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
