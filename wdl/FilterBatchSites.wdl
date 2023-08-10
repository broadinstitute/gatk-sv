version 1.0

import "Structs.wdl"
import "PlotSVCountsPerSample.wdl" as sv_counts

workflow FilterBatchSites {
  input {
    String batch
    File? manta_vcf
    File? wham_vcf
    File? melt_vcf
    File? scramble_vcf
    File? depth_vcf
    File evidence_metrics
    File evidence_metrics_common
    String sv_pipeline_docker

    # PlotSVCountsPerSample metrics
    Int N_IQR_cutoff_plotting = 6

    RuntimeAttr? runtime_attr_adjudicate
    RuntimeAttr? runtime_attr_rewrite_scores
    RuntimeAttr? runtime_attr_filter_annotate_vcf
    RuntimeAttr? runtime_attr_merge_pesr_vcfs
    RuntimeAttr? runtime_attr_count_svs
    RuntimeAttr? runtime_attr_plot_svcounts
    RuntimeAttr? runtime_attr_cat_outliers_preview

  }

  Array[String] algorithms = ["manta", "wham", "melt", "scramble", "depth"]
  Array[File?] vcfs_array = [manta_vcf, wham_vcf, melt_vcf, scramble_vcf, depth_vcf]
  Int num_algorithms = length(algorithms)

  call AdjudicateSV {
    input:
      metrics = evidence_metrics,
      batch = batch,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_adjudicate
  }

  call RewriteScores {
    input:
      metrics = evidence_metrics_common,
      cutoffs = AdjudicateSV.cutoffs,
      scores = AdjudicateSV.scores,
      batch = batch,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_rewrite_scores
  }

  scatter (i in range(num_algorithms)) {
    if (defined(vcfs_array[i])) {
      call FilterAnnotateVcf {
        input:
          vcf = select_first([vcfs_array[i]]),
          metrics = evidence_metrics,
          prefix = "${batch}.${algorithms[i]}",
          scores = RewriteScores.updated_scores,
          cutoffs = AdjudicateSV.cutoffs,
          sv_pipeline_docker = sv_pipeline_docker,
          runtime_attr_override = runtime_attr_filter_annotate_vcf
      }
    }
  }

  call sv_counts.PlotSVCountsPerSample {
    input:
      prefix = batch,
      vcfs=[FilterAnnotateVcf.annotated_vcf[0], FilterAnnotateVcf.annotated_vcf[1], FilterAnnotateVcf.annotated_vcf[2], FilterAnnotateVcf.annotated_vcf[3], FilterAnnotateVcf.annotated_vcf[4]],
      N_IQR_cutoff = N_IQR_cutoff_plotting,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_count_svs = runtime_attr_count_svs,
      runtime_attr_plot_svcounts = runtime_attr_plot_svcounts,
      runtime_attr_cat_outliers_preview = runtime_attr_cat_outliers_preview
  }

  output {
    File? sites_filtered_manta_vcf = FilterAnnotateVcf.annotated_vcf[0]
    File? sites_filtered_wham_vcf = FilterAnnotateVcf.annotated_vcf[1]
    File? sites_filtered_melt_vcf = FilterAnnotateVcf.annotated_vcf[2]
    File? sites_filtered_scramble_vcf = FilterAnnotateVcf.annotated_vcf[3]
    File? sites_filtered_depth_vcf = FilterAnnotateVcf.annotated_vcf[4]
    File cutoffs = AdjudicateSV.cutoffs
    File scores = RewriteScores.updated_scores
    File RF_intermediate_files = AdjudicateSV.RF_intermediate_files
    Array[File] sites_filtered_sv_counts = PlotSVCountsPerSample.sv_counts
    Array[File] sites_filtered_sv_count_plots = PlotSVCountsPerSample.sv_count_plots
    File sites_filtered_outlier_samples_preview = PlotSVCountsPerSample.outlier_samples_preview
    File sites_filtered_outlier_samples_with_reason = PlotSVCountsPerSample.outlier_samples_with_reason
    Int sites_filtered_num_outlier_samples = PlotSVCountsPerSample.num_outlier_samples
  }

}

task AdjudicateSV {
  input {
    File metrics
    String batch
    String sv_pipeline_docker
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
    File scores = "${batch}.scores"
    File cutoffs = "${batch}.cutoffs"
    File RF_intermediate_files = "${batch}.RF_intermediate_files.tar.gz"
  }
  command <<<

    set -euo pipefail
    svtk adjudicate ~{metrics} ~{batch}.scores ~{batch}.cutoffs
    mkdir ~{batch}.RF_intermediate_files
    mv *_trainable.txt ~{batch}.RF_intermediate_files/
    mv *_testable.txt ~{batch}.RF_intermediate_files/
    tar -czvf ~{batch}.RF_intermediate_files.tar.gz ~{batch}.RF_intermediate_files

  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

task RewriteScores {
  input {
    File metrics
    File cutoffs
    File scores
    String batch
    String sv_pipeline_docker
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
    File updated_scores = "${batch}.updated_scores"
  }
  command <<<

    set -euo pipefail
    Rscript /opt/sv-pipeline/03_variant_filtering/scripts/modify_cutoffs.R \
      -c ~{cutoffs} \
      -m ~{metrics} \
      -s ~{scores}  \
      -o ~{batch}.updated_scores

  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

task FilterAnnotateVcf {
  input {
    File vcf
    File metrics
    File scores
    File cutoffs
    String prefix
    String sv_pipeline_docker
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
    File annotated_vcf = "${prefix}.with_evidence.vcf.gz"
  }
  command <<<

    set -euo pipefail
    cat \
    <(sed -e '1d' ~{scores} | fgrep -e DEL -e DUP | awk '($3!="NA" && $3>=0.5)' | cut -f1 | fgrep -w -f - <(zcat ~{vcf})) \
    <(sed -e '1d' ~{scores} | fgrep -e INV -e BND -e INS | awk '($3!="NA" && $3>=0.5)' | cut -f1 | fgrep -w -f - <(zcat ~{vcf}) | sed -e 's/SVTYPE=DEL/SVTYPE=BND/' -e 's/SVTYPE=DUP/SVTYPE=BND/' -e 's/<DEL>/<BND>/' -e 's/<DUP>/<BND>/') \
      | cat <(sed -n -e '/^#/p' <(zcat ~{vcf})) - \
      | vcf-sort -c \
      | bgzip -c \
      > filtered.vcf.gz

    /opt/sv-pipeline/03_variant_filtering/scripts/rewrite_SR_coords.py filtered.vcf.gz ~{metrics} ~{cutoffs} stdout \
      | vcf-sort -c \
      | bgzip -c \
      > filtered.corrected_coords.vcf.gz

    /opt/sv-pipeline/03_variant_filtering/scripts/annotate_RF_evidence.py filtered.corrected_coords.vcf.gz ~{scores} ~{prefix}.with_evidence.vcf
    bgzip ~{prefix}.with_evidence.vcf

  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

}