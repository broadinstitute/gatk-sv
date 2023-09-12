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
    RuntimeAttr? ploidy_score_runtime_attr
    RuntimeAttr? ploidy_build_runtime_attr

    RuntimeAttr? wgd_build_runtime_attr
    RuntimeAttr? wgd_score_runtime_attr
    RuntimeAttr? runtime_attr_bincov_attr
    RuntimeAttr? runtime_attr_make_qc_table

    RuntimeAttr? runtime_attr_mediancov       # Memory ignored, use median_cov_mem_gb_per_sample
    Float? median_cov_mem_gb_per_sample
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
    if (defined(manta_vcfs) && (length(select_first([manta_vcfs])) > 0)) {
      call vcfqc.RawVcfQC as RawVcfQC_Manta {
        input:
          vcfs = select_first([manta_vcfs]),
          prefix = batch,
          caller = "Manta",
          runtime_attr_qc = runtime_attr_qc,
          sv_pipeline_docker = sv_pipeline_docker,
          runtime_attr_outlier = runtime_attr_qc_outlier
      }
    }
    if (defined(melt_vcfs) && (length(select_first([melt_vcfs])) > 0)) {
      call vcfqc.RawVcfQC as RawVcfQC_Melt {
        input:
          vcfs = select_first([melt_vcfs]),
          prefix = batch,
          caller = "Melt",
          runtime_attr_qc = runtime_attr_qc,
          sv_pipeline_docker = sv_pipeline_docker,
          runtime_attr_outlier = runtime_attr_qc_outlier
      }
    }
    if (defined(wham_vcfs) && (length(select_first([wham_vcfs])) > 0)) {
      call vcfqc.RawVcfQC as RawVcfQC_Wham {
        input:
          vcfs = select_first([wham_vcfs]),
          prefix = batch,
          caller = "Wham",
          runtime_attr_qc = runtime_attr_qc,
          sv_pipeline_docker = sv_pipeline_docker,
          runtime_attr_outlier = runtime_attr_qc_outlier
      }
    }
    if (defined(scramble_vcfs) && (length(select_first([scramble_vcfs])) > 0)) {
      call vcfqc.RawVcfQC as RawVcfQC_Scramble {
        input:
          vcfs = select_first([scramble_vcfs]),
          prefix = batch,
          caller = "Scramble",
          runtime_attr_qc = runtime_attr_qc,
          sv_pipeline_docker = sv_pipeline_docker,
          runtime_attr_outlier = runtime_attr_qc_outlier
      }
    }
  }
  if (run_ploidy) {
      call MakeQcTable {
        input:
          output_prefix = batch,
          samples = samples,
          ploidy_plots = select_first([Ploidy.ploidy_plots]),
          bincov_median = MedianCov.medianCov,
          WGD_scores = WGD.WGD_scores,
          melt_insert_size = select_first([melt_insert_size, []]),

          manta_qc_low = RawVcfQC_Manta.low,
          manta_qc_high = RawVcfQC_Manta.high,
          melt_qc_low = RawVcfQC_Melt.low,
          melt_qc_high = RawVcfQC_Melt.high,
          wham_qc_low = RawVcfQC_Wham.low,
          wham_qc_high = RawVcfQC_Wham.high,
          scramble_qc_low = RawVcfQC_Scramble.low,
          scramble_qc_high = RawVcfQC_Scramble.high,

          sv_pipeline_docker = sv_pipeline_docker,
          runtime_override = runtime_attr_make_qc_table
    }
  }

  output {
    File? manta_qc_low = RawVcfQC_Manta.low
    File? manta_qc_high = RawVcfQC_Manta.high
    File? melt_qc_low = RawVcfQC_Melt.low
    File? melt_qc_high = RawVcfQC_Melt.high
    File? wham_qc_low = RawVcfQC_Wham.low
    File? wham_qc_high = RawVcfQC_Wham.high
    File? scramble_qc_low = RawVcfQC_Scramble.low
    File? scramble_qc_high = RawVcfQC_Scramble.high

    File? ploidy_matrix = Ploidy.ploidy_matrix
    File? ploidy_plots = Ploidy.ploidy_plots

    File WGD_dist = WGD.WGD_dist
    File WGD_matrix = WGD.WGD_matrix
    File WGD_scores = WGD.WGD_scores

    File bincov_matrix = MakeBincovMatrix.merged_bincov
    File bincov_matrix_index = MakeBincovMatrix.merged_bincov_idx
    File bincov_median = MedianCov.medianCov

    File? qc_table = MakeQcTable.qc_table
  }
}

task MakeQcTable {
  input {

    File ploidy_plots
    File bincov_median
    File WGD_scores
    Array[Float] melt_insert_size
    Array[String] samples


    File? manta_qc_low
    File? manta_qc_high
    File? melt_qc_low
    File? melt_qc_high
    File? wham_qc_low
    File? wham_qc_high
    File? scramble_qc_low
    File? scramble_qc_high

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
      ~{"--median-cov-filename " + bincov_median} \
      ~{"--wgd-scores-filename " + WGD_scores} \
      ~{"--binwise-cnv-qvalues-filename " + "./ploidy_est/binwise_CNV_qValues.bed.gz"} \
      ~{"--manta-qc-outlier-high-filename " + manta_qc_high} \
      ~{"--melt-qc-outlier-high-filename " + melt_qc_high} \
      ~{"--wham-qc-outlier-high-filename " + wham_qc_high} \
      ~{"--manta-qc-outlier-low-filename " + manta_qc_low} \
      ~{"--melt-qc-outlier-low-filename " + melt_qc_low} \
      ~{"--wham-qc-outlier-low-filename " + wham_qc_low} \
      ~{if (length(melt_insert_size) > 0) then "--melt-insert-size mean_insert_size.tsv" else ""} \
      ~{"--output-prefix " + output_prefix}
  >>>
#
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
