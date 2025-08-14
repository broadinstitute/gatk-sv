version 1.0

import "Structs.wdl"
import "MakeBincovMatrix.wdl" as mbm
import "PloidyEstimation.wdl" as pe
import "RawVcfQC.wdl" as vcfqc
import "WGD.wdl" as wgd
import "MedianCov.wdl" as mc

# Runs single sample tasks on collected evidence:
#   - Ploidy determination
#   - Dosage scoring
#   - QC for raw SV calls (optional)

workflow EvidenceQC {
  input {
    # Batch info
    String batch
    Array[String] samples

    # Optional QC tasks
    Boolean run_vcf_qc

    # Global files
    File genome_file

    # Coverage files
    Array[File] counts

    # SV tool calls
    Array[File]? dragen_vcfs       # Dragen VCF
    Array[File]? manta_vcfs        # Manta VCF
    Array[File]? melt_vcfs         # Melt VCF
    Array[File]? wham_vcfs         # Wham VCF
    Array[File]? scramble_vcfs     # Scramble VCF

    # WGD files
    File wgd_scoring_mask

    # Runtime parameters
    String sv_base_mini_docker
    String sv_base_docker
    String sv_pipeline_docker
    String sv_pipeline_qc_docker

    Int? disk_overhead_bincov_gb

    Boolean run_ploidy = true

    Array[Float]? melt_insert_size

    RuntimeAttr? runtime_attr_qc
    RuntimeAttr? runtime_attr_qc_outlier
    RuntimeAttr? runtime_attr_qc_counts
    RuntimeAttr? ploidy_score_runtime_attr
    RuntimeAttr? ploidy_build_runtime_attr

    RuntimeAttr? wgd_build_runtime_attr
    RuntimeAttr? wgd_score_runtime_attr
    RuntimeAttr? runtime_attr_bincov_attr
    RuntimeAttr? runtime_attr_make_qc_table
    RuntimeAttr? runtime_attr_variant_count_plots

    RuntimeAttr? runtime_attr_mediancov       # Memory ignored, use median_cov_mem_gb_per_sample
    Float? median_cov_mem_gb_per_sample

    # Do not use
    File? NONE_FILE_
  }

  call mbm.MakeBincovMatrix as MakeBincovMatrix {
    input:
      samples = samples,
      count_files = counts,
      batch = batch,
      disk_overhead_gb = disk_overhead_bincov_gb,
      sv_base_mini_docker = sv_base_mini_docker,
      sv_base_docker = sv_base_docker,
      runtime_attr_override = runtime_attr_bincov_attr
  }

  Float median_cov_mem_gb = select_first([median_cov_mem_gb_per_sample, 0.5]) * length(samples) + 7.5
  call mc.MedianCov as MedianCov {
    input:
      bincov_matrix = MakeBincovMatrix.merged_bincov,
      cohort_id = batch,
      sv_pipeline_qc_docker = sv_pipeline_qc_docker,
      runtime_attr = runtime_attr_mediancov,
      mem_gb_override = median_cov_mem_gb
  }

  if (run_ploidy) {
    call pe.Ploidy as Ploidy {
      input:
        bincov_matrix = MakeBincovMatrix.merged_bincov,
        batch = batch,
        sv_base_mini_docker = sv_base_mini_docker,
        sv_pipeline_qc_docker = sv_pipeline_qc_docker,
        runtime_attr_score = ploidy_score_runtime_attr,
        runtime_attr_build = ploidy_build_runtime_attr
    }
  }

  call wgd.WGD as WGD {
    input:
      batch = batch,
      wgd_scoring_mask = wgd_scoring_mask,
      bincov_matrix = MakeBincovMatrix.merged_bincov,
      sv_pipeline_qc_docker = sv_pipeline_qc_docker,
      runtime_attr_build = wgd_build_runtime_attr,
      runtime_attr_score = wgd_score_runtime_attr
  }

  if (run_vcf_qc) {
    if (defined(dragen_vcfs) && (length(select_first([dragen_vcfs])) > 0)) {
      call vcfqc.RawVcfQC as RawVcfQC_Dragen {
        input:
          vcfs = select_first([dragen_vcfs]),
          prefix = batch,
          caller = "Dragen",
          runtime_attr_qc = runtime_attr_qc,
          sv_pipeline_docker = sv_pipeline_docker,
          runtime_attr_outlier = runtime_attr_qc_outlier
      }
    }
    if (defined(manta_vcfs) && (length(select_first([manta_vcfs])) > 0)) {
      call vcfqc.RawVcfQC as RawVcfQC_Manta {
        input:
          vcfs = select_first([manta_vcfs]),
          prefix = batch,
          caller = "Manta",
          sv_pipeline_docker = sv_pipeline_docker,
          runtime_attr_qc = runtime_attr_qc,
          runtime_attr_outlier = runtime_attr_qc_outlier,
          runtime_attr_counts = runtime_attr_qc_counts
      }
    }
    if (defined(melt_vcfs) && (length(select_first([melt_vcfs])) > 0)) {
      call vcfqc.RawVcfQC as RawVcfQC_Melt {
        input:
          vcfs = select_first([melt_vcfs]),
          prefix = batch,
          caller = "Melt",
          sv_pipeline_docker = sv_pipeline_docker,
          runtime_attr_qc = runtime_attr_qc,
          runtime_attr_outlier = runtime_attr_qc_outlier,
          runtime_attr_counts = runtime_attr_qc_counts
      }
    }
    if (defined(wham_vcfs) && (length(select_first([wham_vcfs])) > 0)) {
      call vcfqc.RawVcfQC as RawVcfQC_Wham {
        input:
          vcfs = select_first([wham_vcfs]),
          prefix = batch,
          caller = "Wham",
          sv_pipeline_docker = sv_pipeline_docker,
          runtime_attr_qc = runtime_attr_qc,
          runtime_attr_outlier = runtime_attr_qc_outlier,
          runtime_attr_counts = runtime_attr_qc_counts
      }
    }
    if (defined(scramble_vcfs) && (length(select_first([scramble_vcfs])) > 0)) {
      call vcfqc.RawVcfQC as RawVcfQC_Scramble {
        input:
          vcfs = select_first([scramble_vcfs]),
          prefix = batch,
          caller = "Scramble",
          sv_pipeline_docker = sv_pipeline_docker,
          runtime_attr_qc = runtime_attr_qc,
          runtime_attr_outlier = runtime_attr_qc_outlier,
          runtime_attr_counts = runtime_attr_qc_counts
      }
    }
  }

  if (run_ploidy) {
    if (run_vcf_qc) {
      Array[File] variant_count_files = select_all([
        RawVcfQC_Dragen.variant_counts,
        RawVcfQC_Manta.variant_counts,
        RawVcfQC_Melt.variant_counts,
        RawVcfQC_Wham.variant_counts,
        RawVcfQC_Scramble.variant_counts
      ])

      call CreateVariantCountPlots {
        input:
          variant_count_files = variant_count_files,
          ploidy_plots_tarball = select_first([Ploidy.ploidy_plots]),
          output_prefix = batch,
          sv_pipeline_docker = sv_pipeline_docker,
          runtime_attr_override = runtime_attr_variant_count_plots
      }
    }

    call MakeQcTable {
      input:
        output_prefix = batch,
        samples = samples,
        ploidy_plots = select_first([CreateVariantCountPlots.ploidy_plots, Ploidy.ploidy_plots]),
        bincov_median = MedianCov.medianCov,
        WGD_scores = WGD.WGD_scores,
        melt_insert_size = select_first([melt_insert_size, []]),

        dragen_qc_low = RawVcfQC_Dragen.low,
        dragen_qc_high = RawVcfQC_Dragen.high,
        dragen_variant_counts = RawVcfQC_Dragen.variant_counts,

        manta_qc_low = RawVcfQC_Manta.low,
        manta_qc_high = RawVcfQC_Manta.high,
        manta_variant_counts = RawVcfQC_Manta.variant_counts,

        melt_qc_low = RawVcfQC_Melt.low,
        melt_qc_high = RawVcfQC_Melt.high,
        melt_variant_counts = RawVcfQC_Melt.variant_counts,

        wham_qc_low = RawVcfQC_Wham.low,
        wham_qc_high = RawVcfQC_Wham.high,
        wham_variant_counts = RawVcfQC_Wham.variant_counts,

        scramble_qc_low = RawVcfQC_Scramble.low,
        scramble_qc_high = RawVcfQC_Scramble.high,
        scramble_variant_counts = RawVcfQC_Scramble.variant_counts,

        sv_pipeline_docker = sv_pipeline_docker,
        runtime_override = runtime_attr_make_qc_table
    }
  }

  output {
    File? dragen_qc_low = RawVcfQC_Dragen.low
    File? dragen_qc_high = RawVcfQC_Dragen.high
    File? manta_qc_low = RawVcfQC_Manta.low
    File? manta_qc_high = RawVcfQC_Manta.high
    File? melt_qc_low = RawVcfQC_Melt.low
    File? melt_qc_high = RawVcfQC_Melt.high
    File? wham_qc_low = RawVcfQC_Wham.low
    File? wham_qc_high = RawVcfQC_Wham.high
    File? scramble_qc_low = RawVcfQC_Scramble.low
    File? scramble_qc_high = RawVcfQC_Scramble.high
    
    File? dragen_variant_counts = RawVcfQC_Dragen.variant_counts
    File? manta_variant_counts = RawVcfQC_Manta.variant_counts
    File? melt_variant_counts = RawVcfQC_Melt.variant_counts
    File? wham_variant_counts = RawVcfQC_Wham.variant_counts
    File? scramble_variant_counts = RawVcfQC_Scramble.variant_counts

    File? ploidy_matrix = Ploidy.ploidy_matrix
    File? ploidy_plots = if run_ploidy then select_first([CreateVariantCountPlots.ploidy_plots, Ploidy.ploidy_plots]) else NONE_FILE_

    File WGD_dist = WGD.WGD_dist
    File WGD_matrix = WGD.WGD_matrix
    File WGD_scores = WGD.WGD_scores

    File bincov_matrix = MakeBincovMatrix.merged_bincov
    File bincov_matrix_index = MakeBincovMatrix.merged_bincov_idx
    File bincov_median = MedianCov.medianCov

    File? qc_table = MakeQcTable.qc_table
  }
}

task CreateVariantCountPlots {
  input {
    Array[File] variant_count_files
    File ploidy_plots_tarball
    String output_prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr runtime_default = object {
    cpu_cores: 1,
    mem_gb: 4,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1,
    disk_gb: 50
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, runtime_default])

  output {
    File ploidy_plots = "${output_prefix}.ploidy_est.tar.gz"
  }

  command <<<
    set -euo pipefail

    tar -xzf ~{ploidy_plots_tarball}

    mkdir -p ./ploidy_est/variant_count_plots
    
    for file in ~{sep=' ' variant_count_files}; do
      if [[ -f "$file" ]]; then
        caller=$(basename "$file" | cut -d'.' -f2)
        echo "Processing $caller variant counts from $file"
        
        cd ./ploidy_est/variant_count_plots
        Rscript /opt/sv-pipeline/scripts/plot_variant_counts.R "$file" "$caller"
        cd ../..
      fi
    done

    tar -czf ~{output_prefix}.ploidy_est.tar.gz ./ploidy_est/
  >>>

  runtime {
    docker: sv_pipeline_docker
    cpu: select_first([runtime_attr.cpu_cores, runtime_default.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, runtime_default.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, runtime_default.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, runtime_default.boot_disk_gb])
    preemptible: select_first([runtime_attr.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, runtime_default.max_retries])
  }
}

task MakeQcTable {
  input {

    File ploidy_plots
    File bincov_median
    File WGD_scores
    Array[Float] melt_insert_size
    Array[String] samples


    File? dragen_qc_low
    File? dragen_qc_high
    File? dragen_variant_counts
    
    File? manta_qc_low
    File? manta_qc_high
    File? manta_variant_counts

    File? melt_qc_low
    File? melt_qc_high
    File? melt_variant_counts

    File? wham_qc_low
    File? wham_qc_high
    File? wham_variant_counts

    File? scramble_qc_low
    File? scramble_qc_high
    File? scramble_variant_counts

    String sv_pipeline_docker
    String output_prefix
    RuntimeAttr? runtime_override
  }

  output {
    File qc_table = "${output_prefix}.evidence_qc_table.tsv"
  }

  command <<<
    set -euo pipefail

    if ~{length(melt_insert_size) > 0} ; then
      echo -e "sample_ID\tmean_insert_size" > mean_insert_size.tsv
      paste ~{write_lines(samples)} ~{write_lines(melt_insert_size)} >> mean_insert_size.tsv
    fi

    tar -xvf ~{ploidy_plots}

    python /opt/sv-pipeline/scripts/make_evidence_qc_table.py \
      ~{"--estimated-copy-number-filename " + "./ploidy_est/estimated_copy_numbers.txt.gz"} \
      ~{"--sex-assignments-filename " + "./ploidy_est/sample_sex_assignments.txt.gz"} \
      ~{"--median-cov-filename " + bincov_median} \
      ~{"--wgd-scores-filename " + WGD_scores} \
      ~{"--binwise-cnv-qvalues-filename " + "./ploidy_est/binwise_CNV_qValues.bed.gz"} \
      ~{"--dragen-qc-outlier-high-filename " + dragen_qc_high} \
      ~{"--manta-qc-outlier-high-filename " + manta_qc_high} \
      ~{"--melt-qc-outlier-high-filename " + melt_qc_high} \
      ~{"--wham-qc-outlier-high-filename " + wham_qc_high} \
      ~{"--scramble-qc-outlier-high-filename " + scramble_qc_high} \
      ~{"--dragen-qc-outlier-low-filename " + dragen_qc_low} \
      ~{"--manta-qc-outlier-low-filename " + manta_qc_low} \
      ~{"--melt-qc-outlier-low-filename " + melt_qc_low} \
      ~{"--wham-qc-outlier-low-filename " + wham_qc_low} \
      ~{"--scramble-qc-outlier-low-filename " + scramble_qc_low} \
      ~{"--dragen-variant-counts-filename " + dragen_variant_counts} \
      ~{"--manta-variant-counts-filename " + manta_variant_counts} \
      ~{"--melt-variant-counts-filename " + melt_variant_counts} \
      ~{"--wham-variant-counts-filename " + wham_variant_counts} \
      ~{"--scramble-variant-counts-filename " + scramble_variant_counts} \
      ~{if (length(melt_insert_size) > 0) then "--melt-insert-size mean_insert_size.tsv" else ""} \
      ~{"--output-prefix " + output_prefix}
  >>>
  
  RuntimeAttr runtime_default = object {
    cpu_cores: 1,
    mem_gb: 4,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1,
    disk_gb: 100
  }

  RuntimeAttr runtime_attr = select_first([runtime_override, runtime_default])

    runtime {
        docker: sv_pipeline_docker
        cpu: select_first([runtime_attr.cpu_cores, runtime_default.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, runtime_default.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, runtime_default.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, runtime_default.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, runtime_default.max_retries])
    }
}
