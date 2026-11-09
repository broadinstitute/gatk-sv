version 1.0

import "CollectCoverage.wdl" as cov
import "Structs.wdl"

# Nominates potential mosaic carriers of genomic disorders (GDs) from read depth.
#
# Inputs are per-batch coverage matrices, parallel to batches: low_res_rd_matrices
# (used as is) and/or high_res_rd_matrices (bgzipped, with a .tbi next to each, as
# for CallGenomicDisorderCNVs). When low_res_rd_matrices is not given, the
# high-res matrices are condensed with CondenseReadCounts (reference_dict,
# rebinned_interval_size); the low-res matrices are then workflow outputs. When
# high_res_rd_matrices are given, they replace the low-res bins for intervals
# shorter than 10x the low-res median bin size, as plot.py does. At least one of
# the two is required.
#
# For every autosomal interval (distinct chr/start/end) in gd_table, each batch is
# processed independently (one task per batch, so hundreds of batches scale out):
#   1. Normalization: each sample's counts are scaled so its autosomal median is
#      CN 2 (the plot.py raw-depth procedure).
#   2. Abnormal depth outside the interval: a sample is dropped for that interval
#      when its depth deviates (more than max_flank_deviation from the median of
#      the batch) in at least flank_min_deviant_bins of any flank_window_bins
#      consecutive bins, in the same direction, in the left/right flank (the
#      cluster's body size clipped to 50 kb - 1 Mb) or in any other interval of
#      the interval's cluster. Bins overlapping a segmental duplication, and bins
#      where more than flank_max_population_fraction of the batch deviates, are
#      ignored (the PlotGDMosaicDepth procedure).
#   3. Candidate test: the median normalized depth over the interval, excluding
#      bins that overlap a segmental duplication by at least 1 bp, must lie in
#      [del_min, del_max] (default 1.2-1.8, mosaic DEL) or [dup_min, dup_max]
#      (default 2.2-2.8, mosaic DUP), and at least min_interval_bins bins must
#      remain.
# chrX and chrY intervals are skipped. DEL and DUP rows of one interval share one
# depth; a candidate is reported under the GD_ID whose svtype matches the
# direction of its depth (the interval's first GD_ID if none does).
#
# Outputs: candidate_table (gd_id, sample, batch, median_depth, plus mosaic_type,
# interval coordinates, svtype and n_bins), and an interval_summary table (per
# batch and interval: samples, bins, samples dropped for abnormal depth, and
# candidates). nominate_mosaic_script is nominate_mosaic.py from gatk-sv-gd; the
# gatk_sv_gd package must be installed in the image.
workflow NominateMosaicGDCarriers {
  input {
    Array[String] batches
    Array[File]? low_res_rd_matrices
    Array[File]? high_res_rd_matrices
    File? reference_dict
    Int rebinned_interval_size = 10000
    File gd_table
    File segdup_bed

    Float del_min = 1.2
    Float del_max = 1.8
    Float dup_min = 2.2
    Float dup_max = 2.8
    Int min_interval_bins = 1
    Boolean nahr_only = true
    Float max_flank_deviation = 0.4
    Int flank_window_bins = 5
    Int flank_min_deviant_bins = 3
    Float flank_max_population_fraction = 0.05

    File nominate_mosaic_script

    String output_prefix = "gd_mosaic_candidates"
    String sv_pipeline_docker
    String? gatk_docker
    RuntimeAttr? runtime_attr_condense
    RuntimeAttr? runtime_attr_nominate
    RuntimeAttr? runtime_attr_merge
  }

  Boolean use_low_res = defined(low_res_rd_matrices)
  Boolean use_high_res = defined(high_res_rd_matrices)
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

    if (use_high_res) {
      File hr_counts = high_res_inputs[i]
      File hr_index = high_res_inputs[i] + ".tbi"
    }

    File low_res_counts = if use_low_res then low_res_inputs[i] else select_first([CondenseReadCounts.out])

    call NominateBatch {
      input:
        batch = batches[i],
        low_res_counts = low_res_counts,
        high_res_counts = hr_counts,
        high_res_counts_index = hr_index,
        gd_table = gd_table,
        segdup_bed = segdup_bed,
        del_min = del_min,
        del_max = del_max,
        dup_min = dup_min,
        dup_max = dup_max,
        min_interval_bins = min_interval_bins,
        nahr_only = nahr_only,
        max_flank_deviation = max_flank_deviation,
        flank_window_bins = flank_window_bins,
        flank_min_deviant_bins = flank_min_deviant_bins,
        flank_max_population_fraction = flank_max_population_fraction,
        nominate_mosaic_script = nominate_mosaic_script,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_nominate
    }
  }

  call MergeTables as MergeCandidates {
    input:
      tables = NominateBatch.candidates,
      out_name = output_prefix + ".candidates.tsv.gz",
      columns = "gd_id,sample,batch,median_depth,mosaic_type,chr,start,end,svtype,n_bins",
      nominate_mosaic_script = nominate_mosaic_script,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_merge
  }

  call MergeTables as MergeSummaries {
    input:
      tables = NominateBatch.interval_summary,
      out_name = output_prefix + ".interval_summary.tsv.gz",
      columns = "batch,gd_id,chr,start,end,n_samples,n_bins,n_abnormal,n_candidates",
      nominate_mosaic_script = nominate_mosaic_script,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_merge
  }

  output {
    File candidate_table = MergeCandidates.merged
    File interval_summary = MergeSummaries.merged
    Array[File] batch_candidate_tables = NominateBatch.candidates
    Array[File] condensed_low_res_matrices = if use_low_res then low_res_inputs else select_all(CondenseReadCounts.out)
  }
}

# Nominates candidate mosaic events for one batch.
task NominateBatch {
  input {
    String batch
    File low_res_counts
    File? high_res_counts
    File? high_res_counts_index
    File gd_table
    File segdup_bed
    Float del_min
    Float del_max
    Float dup_min
    Float dup_max
    Int min_interval_bins
    Boolean nahr_only
    Float max_flank_deviation
    Int flank_window_bins
    Int flank_min_deviant_bins
    Float flank_max_population_fraction
    File nominate_mosaic_script
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Int disk_gb = ceil(20 + 3 * size(low_res_counts, "GiB") + size(high_res_counts, "GiB") + size(segdup_bed, "GiB"))

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

    python ~{nominate_mosaic_script} nominate \
      --batch '~{batch}' \
      --gd-table ~{gd_table} \
      --low-res-counts ~{low_res_counts} \
      ~{if defined(high_res_counts) then "--high-res-counts " + high_res_counts else ""} \
      --segdup-bed ~{segdup_bed} \
      --del-min ~{del_min} --del-max ~{del_max} \
      --dup-min ~{dup_min} --dup-max ~{dup_max} \
      --min-bins ~{min_interval_bins} \
      ~{if nahr_only then "--nahr-only" else "--no-nahr-only"} \
      --max-deviation ~{max_flank_deviation} \
      --window-bins ~{flank_window_bins} \
      --min-deviant-bins ~{flank_min_deviant_bins} \
      --max-population-fraction ~{flank_max_population_fraction} \
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
    File candidates = batch + ".candidates.tsv.gz"
    File interval_summary = batch + ".interval_summary.tsv.gz"
  }
}

# Concatenates per-batch tables. Paths are passed in a list file, not on the
# command line, so hundreds of batches do not exceed argument limits; tables are
# streamed one at a time.
task MergeTables {
  input {
    Array[File] tables
    String out_name
    String columns
    File nominate_mosaic_script
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Int disk_gb = ceil(10 + 2 * size(tables, "GiB"))

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 4.0,
    disk_gb: disk_gb,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail

    python ~{nominate_mosaic_script} merge \
      --input-list ~{write_lines(tables)} \
      --columns '~{columns}' \
      --out '~{out_name}'
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
    File merged = out_name
  }
}
