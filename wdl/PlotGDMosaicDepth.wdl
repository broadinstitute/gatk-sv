version 1.0

import "CollectCoverage.wdl" as cov
import "Structs.wdl"

# Normalized read depth across a GD locus (or its whole cluster) for mosaic
# carriers, random germline carriers and random non-carriers. Non-carriers come
# from the mosaic carriers' batches; germline carriers do too, plus the fewest
# extra batches needed to reach n_germline_carriers (or all that exist).
#
# Germline carriers include those of the DEL and DUP at the target locus (the
# partner GD_ID differs only in "DEL"/"DUP"), plotted in one color on one
# legend line with separate DEL and DUP counts. n_germline_carriers is the total
# over both types: a type with 3 or fewer carriers is included whole, otherwise
# the types share the total in proportion to their cohort frequency, each with
# at least 3. Mosaic carriers are typed DEL/DUP by median depth over the locus
# (below/above ploidy) and their counts are listed by type when both occur.
# When fewer germline carriers (of a type) or non-carriers are plotted than are
# eligible in the cohort, the legend reads e.g. "Germline DEL (n=20, selected)".
#
# include_samples (optional, one sample ID per line) restricts plotting to
# those samples: comparators are drawn only from the list, and mosaic samples
# outside it are dropped with a warning.
#
# Preprocessing and plotting are separate. Preprocessing always runs and
# outputs everything needed to draw the plot: the selected samples with their
# group (selected_samples), the plotted region (locus_region), and normalized
# low-res depth for the selected samples at the locus/cluster, one table per
# batch contributing a selected sample (normalized_depth_tables; batches_used
# lists them). The publication PDF is drawn from
# those outputs only when make_plot is true.
#
# Samples with abnormal depth in the flanks or in cluster intervals outside the
# target GD (at least flank_min_deviant_bins of any flank_window_bins
# consecutive usable bins differ from the batch median of same-ploidy samples
# by more than max_flank_deviation, in the same direction) are left out. Bins
# overlapping a segdup (segdup_bed) and bins where more than
# flank_max_population_fraction of the batch deviates (common polymorphism)
# are ignored of the germline and non-carrier comparators;
# a mosaic sample with such depth is kept, with a warning. This needs the
# low-res counts of every batch: give low_res_rd_matrices (condensed counts,
# e.g. condensed_low_res_matrices from a previous run) or leave it unset and
# every batch's high-res matrix is condensed.
#
# Inputs are per-batch outputs of CallGenomicDisorderCNVs. Arrays are indexed
# together: batches[i] goes with gd_output_tarballs[i], low_res_rd_matrices[i]
# and high_res_rd_matrices[i]. Use the same gd_table, reference_dict,
# rebinned_interval_size and high_res_rd_matrix given to
# CallGenomicDisorderCNVs so depth is normalized the way plot.py normalizes it.
# Only the condensed low-res counts are plotted; the high-res matrix is only
# the source for condensing, and is not needed when low_res_rd_matrices is set.
#
# Which gatk_sv_gd code runs (first match wins): gatk_sv_gd_src, a tarball of the
# package made with `tar czf gatk_sv_gd_src.tar.gz -C gatk-sv-gd/src gatk_sv_gd`
# (put on PYTHONPATH ahead of the image's copy, so no image rebuild is needed to
# test changes); else plot_mosaic_script, a single script (it still imports the rest of the
# package from the image); else the gatk_sv_gd package installed in the image.
workflow PlotGDMosaicDepth {
  input {
    String gd_id
    Array[String] mosaic_samples
    File? include_samples

    Array[String] batches
    Array[File] gd_output_tarballs
    Array[File]? low_res_rd_matrices
    Array[File]? high_res_rd_matrices

    File gd_table
    File reference_dict
    File? segdup_bed
    Int rebinned_interval_size = 10000

    Int n_germline_carriers = 10
    Int n_non_carriers = 50
    Int random_seed = 42
    Float max_flank_deviation = 0.4
    Int flank_window_bins = 5
    Int flank_min_deviant_bins = 3
    Float flank_max_population_fraction = 0.05

    Boolean make_plot = true
    Float fig_width_in = 2.0
    Float fig_height_in = 2.1
    Float flank_scale = 0.15
    Float depth_axis_min = 0.0
    Float depth_axis_max = 4.0
    String? plot_title
    File? plot_mosaic_script
    File? gatk_sv_gd_src

    String output_prefix = "gd_mosaic_depth"
    String sv_pipeline_docker
    String gatk_docker
    RuntimeAttr? runtime_attr_condense
    RuntimeAttr? runtime_attr_classify
    RuntimeAttr? runtime_attr_select
    RuntimeAttr? runtime_attr_extract
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
          gatk_docker = gatk_docker,
          runtime_attr_override = runtime_attr_condense
      }
    }

    File low_res_counts = if use_low_res then low_res_inputs[i] else select_first([CondenseReadCounts.out])

    call ClassifyBatch {
      input:
        batch = batches[i],
        gd_output_tarball = gd_output_tarballs[i],
        low_res_counts = low_res_counts,
        max_flank_deviation = max_flank_deviation,
        flank_window_bins = flank_window_bins,
        flank_min_deviant_bins = flank_min_deviant_bins,
        flank_max_population_fraction = flank_max_population_fraction,
        segdup_bed = segdup_bed,
        gd_table = gd_table,
        gd_id = gd_id,
        plot_mosaic_script = plot_mosaic_script,
        gatk_sv_gd_src = gatk_sv_gd_src,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_classify
    }
  }

  call SelectSamples {
    input:
      mosaic_samples = mosaic_samples,
      include_samples = include_samples,
      classes = ClassifyBatch.classes,
      regions = ClassifyBatch.region,
      n_germline_carriers = n_germline_carriers,
      n_non_carriers = n_non_carriers,
      random_seed = random_seed,
      output_prefix = output_prefix,
      plot_mosaic_script = plot_mosaic_script,
      gatk_sv_gd_src = gatk_sv_gd_src,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_select
  }

  Array[String] use_batch = read_lines(SelectSamples.use_batch_flags)

  scatter (i in range(length(batches))) {
    if (use_batch[i] == "true") {
      call ExtractBatchDepth {
        input:
          batch = batches[i],
          low_res_counts = low_res_counts[i],
          selected = SelectSamples.selected,
          region = SelectSamples.region,
          gd_table = gd_table,
          gd_id = gd_id,
          plot_mosaic_script = plot_mosaic_script,
          gatk_sv_gd_src = gatk_sv_gd_src,
          sv_pipeline_docker = sv_pipeline_docker,
          runtime_attr_override = runtime_attr_extract
      }
    }
  }

  Array[File] depth_tables = select_all(ExtractBatchDepth.depth)

  if (make_plot) {
    call PlotMosaicDepth {
      input:
        depth_tables = depth_tables,
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
        gatk_sv_gd_src = gatk_sv_gd_src,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_plot
    }
  }

  output {
    File selected_samples = SelectSamples.selected
    File selected_batches = SelectSamples.batches_used
    Array[File] condensed_low_res_matrices = if use_low_res then low_res_inputs else select_all(CondenseReadCounts.out)
    File locus_region = SelectSamples.region
    Array[File] normalized_depth_tables = depth_tables
    File? mosaic_depth_plot = PlotMosaicDepth.pdf
  }
}

# Labels every sample at the locus's cluster as germline carrier of gd_id,
# non-carrier of all GDs in the cluster, or other; also flags samples with
# abnormal depth in the flanks / other cluster intervals.
task ClassifyBatch {
  input {
    String batch
    File gd_output_tarball
    File low_res_counts
    Float max_flank_deviation
    Int flank_window_bins
    Int flank_min_deviant_bins
    Float flank_max_population_fraction
    File? segdup_bed
    File gd_table
    String gd_id
    File? plot_mosaic_script
    File? gatk_sv_gd_src
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Int disk_gb = ceil(20 + 3 * size(gd_output_tarball, "GiB") + 2 * size(low_res_counts, "GiB"))

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

    GD_SRC='~{if defined(gatk_sv_gd_src) then gatk_sv_gd_src else ""}'
    GD_SCRIPT='~{if defined(plot_mosaic_script) then plot_mosaic_script else ""}'
    if [ -n "$GD_SRC" ]; then
      mkdir -p gd_pkg
      tar xzf "$GD_SRC" -C gd_pkg
      test -d gd_pkg/gatk_sv_gd || { echo "ERROR: gatk_sv_gd_src must contain a top-level gatk_sv_gd/ directory" >&2; exit 1; }
      export PYTHONPATH="$PWD/gd_pkg${PYTHONPATH:+:$PYTHONPATH}"
      GD_PY=(python -m gatk_sv_gd.plot_mosaic)
    elif [ -n "$GD_SCRIPT" ]; then
      GD_PY=(python "$GD_SCRIPT")
    else
      GD_PY=(python -m gatk_sv_gd.plot_mosaic)
    fi

    "${GD_PY[@]}" classify \
      --batch '~{batch}' \
      --gd-table ~{gd_table} \
      --gd-id '~{gd_id}' \
      --calls gd_outputs/call/gd_cnv_calls.tsv.gz \
      --cn-posteriors gd_outputs/infer/cn_posteriors.tsv.gz \
      ${PLOIDY_ARG} \
      --low-res-counts ~{low_res_counts} \
      --max-flank-deviation ~{max_flank_deviation} \
      --flank-window-bins ~{flank_window_bins} \
      --flank-min-deviant-bins ~{flank_min_deviant_bins} \
      --flank-max-population-fraction ~{flank_max_population_fraction} \
      ~{if defined(segdup_bed) then "--segdup-bed " + segdup_bed else ""} \
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

# Finds the batches holding the mosaic carriers and draws random non-carriers
# from them. Draws min(n_germline_carriers, all available) germline carriers,
# from the mosaic batches first; if those hold too few, the fewest additional
# batches (random among ties) supply the rest. use_batch_flags has one
# true/false per entry of regions, in order.
task SelectSamples {
  input {
    Array[String] mosaic_samples
    File? include_samples
    Array[File] classes
    Array[File] regions
    Int n_germline_carriers
    Int n_non_carriers
    Int random_seed
    String output_prefix
    File? plot_mosaic_script
    File? gatk_sv_gd_src
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

    GD_SRC='~{if defined(gatk_sv_gd_src) then gatk_sv_gd_src else ""}'
    GD_SCRIPT='~{if defined(plot_mosaic_script) then plot_mosaic_script else ""}'
    if [ -n "$GD_SRC" ]; then
      mkdir -p gd_pkg
      tar xzf "$GD_SRC" -C gd_pkg
      test -d gd_pkg/gatk_sv_gd || { echo "ERROR: gatk_sv_gd_src must contain a top-level gatk_sv_gd/ directory" >&2; exit 1; }
      export PYTHONPATH="$PWD/gd_pkg${PYTHONPATH:+:$PYTHONPATH}"
      GD_PY=(python -m gatk_sv_gd.plot_mosaic)
    elif [ -n "$GD_SCRIPT" ]; then
      GD_PY=(python "$GD_SCRIPT")
    else
      GD_PY=(python -m gatk_sv_gd.plot_mosaic)
    fi

    "${GD_PY[@]}" select \
      --classes ~{sep=" " classes} \
      --regions ~{sep=" " regions} \
      --mosaic-samples ~{write_lines(mosaic_samples)} \
      ~{if defined(include_samples) then "--include-samples " + include_samples else ""} \
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
        File batches_used = output_prefix + ".batches.txt"
    File use_batch_flags = output_prefix + ".use_batch.txt"
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
    File? plot_mosaic_script
    File? gatk_sv_gd_src
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

    GD_SRC='~{if defined(gatk_sv_gd_src) then gatk_sv_gd_src else ""}'
    GD_SCRIPT='~{if defined(plot_mosaic_script) then plot_mosaic_script else ""}'
    if [ -n "$GD_SRC" ]; then
      mkdir -p gd_pkg
      tar xzf "$GD_SRC" -C gd_pkg
      test -d gd_pkg/gatk_sv_gd || { echo "ERROR: gatk_sv_gd_src must contain a top-level gatk_sv_gd/ directory" >&2; exit 1; }
      export PYTHONPATH="$PWD/gd_pkg${PYTHONPATH:+:$PYTHONPATH}"
      GD_PY=(python -m gatk_sv_gd.plot_mosaic)
    elif [ -n "$GD_SCRIPT" ]; then
      GD_PY=(python "$GD_SCRIPT")
    else
      GD_PY=(python -m gatk_sv_gd.plot_mosaic)
    fi

    "${GD_PY[@]}" extract \
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
    File? plot_mosaic_script
    File? gatk_sv_gd_src
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

    GD_SRC='~{if defined(gatk_sv_gd_src) then gatk_sv_gd_src else ""}'
    GD_SCRIPT='~{if defined(plot_mosaic_script) then plot_mosaic_script else ""}'
    if [ -n "$GD_SRC" ]; then
      mkdir -p gd_pkg
      tar xzf "$GD_SRC" -C gd_pkg
      test -d gd_pkg/gatk_sv_gd || { echo "ERROR: gatk_sv_gd_src must contain a top-level gatk_sv_gd/ directory" >&2; exit 1; }
      export PYTHONPATH="$PWD/gd_pkg${PYTHONPATH:+:$PYTHONPATH}"
      GD_PY=(python -m gatk_sv_gd.plot_mosaic)
    elif [ -n "$GD_SCRIPT" ]; then
      GD_PY=(python "$GD_SCRIPT")
    else
      GD_PY=(python -m gatk_sv_gd.plot_mosaic)
    fi

    "${GD_PY[@]}" plot \
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
      ~{if defined(plot_title) then "--title '" + sub(select_first([plot_title, ""]), "'", "'\\''") + "'" else ""} \
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
