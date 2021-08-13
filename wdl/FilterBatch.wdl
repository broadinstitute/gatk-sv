version 1.0

import "FilterOutliers.wdl" as filter_outliers
import "Utils.wdl" as util
import "FilterBatchMetrics.wdl" as metrics

workflow FilterBatch {
  input {
    String batch
    File? manta_vcf
    File? delly_vcf
    File? wham_vcf
    File? melt_vcf
    File? depth_vcf
    File evidence_metrics
    File evidence_metrics_common

    Int outlier_cutoff_nIQR
    File? outlier_cutoff_table

    # Module metrics parameters
    # Run module metrics workflow at the end - on by default
    Boolean? run_module_metrics
    String? sv_pipeline_base_docker  # required if run_module_metrics = true
    File? primary_contigs_list  # required if run_module_metrics = true
    File? ped_file  # required if run_module_metrics = true
    File? baseline_filtered_depth_vcf  # baseline files are optional for metrics workflow
    File? baseline_filtered_pesr_vcf

    String sv_pipeline_docker
    String sv_base_mini_docker
    String linux_docker

    RuntimeAttr? runtime_attr_adjudicate
    RuntimeAttr? runtime_attr_rewrite_scores
    RuntimeAttr? runtime_attr_filter_annotate_vcf
    RuntimeAttr? runtime_attr_ids_from_vcf
    RuntimeAttr? runtime_attr_merge_pesr_vcfs

    RuntimeAttr? runtime_attr_identify_outliers
    RuntimeAttr? runtime_attr_exclude_outliers
    RuntimeAttr? runtime_attr_cat_outliers
    RuntimeAttr? runtime_attr_filter_samples
  }

  Array[String] algorithms = ["manta", "delly", "wham", "melt", "depth"]
  Array[File?] vcfs_array = [manta_vcf, delly_vcf, wham_vcf, melt_vcf, depth_vcf]
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

  call util.GetSampleIdsFromVcf {
    input:
      vcf = select_first(vcfs_array),
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_ids_from_vcf
  }

  call filter_outliers.FilterOutlierSamples as FilterOutlierSamples {
    input:
      batch = batch,
      algorithms = algorithms,
      outlier_cutoff_table = outlier_cutoff_table,
      N_IQR_cutoff = outlier_cutoff_nIQR,
      vcfs = FilterAnnotateVcf.annotated_vcf,
      samples = GetSampleIdsFromVcf.out_array,
      linux_docker = linux_docker,
      sv_pipeline_docker = sv_pipeline_docker,
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_identify_outliers = runtime_attr_identify_outliers,
      runtime_attr_exclude_outliers= runtime_attr_exclude_outliers,
      runtime_attr_cat_outliers = runtime_attr_cat_outliers,
      runtime_attr_filter_samples = runtime_attr_filter_samples
  }
  
  call MergePesrVcfs {
    input:
      manta_vcf = FilterOutlierSamples.vcfs_noOutliers[0],
      delly_vcf = FilterOutlierSamples.vcfs_noOutliers[1],
      wham_vcf = FilterOutlierSamples.vcfs_noOutliers[2],
      melt_vcf = FilterOutlierSamples.vcfs_noOutliers[3],
      batch = batch,
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_merge_pesr_vcfs
  }

  Boolean run_module_metrics_ = if defined(run_module_metrics) then select_first([run_module_metrics]) else true
  if (run_module_metrics_) {
    call metrics.FilterBatchMetrics {
      input:
        name = batch,
        samples = GetSampleIdsFromVcf.out_array,
        filtered_pesr_vcf = MergePesrVcfs.merged_pesr_vcf,
        filtered_depth_vcf = select_first([FilterOutlierSamples.vcfs_noOutliers[4]]),
        cutoffs = AdjudicateSV.cutoffs,
        outlier_list = FilterOutlierSamples.outlier_samples_excluded_file,
        ped_file = select_first([ped_file]),
        samples_post_filtering_file = FilterOutlierSamples.filtered_batch_samples_file,
        baseline_filtered_pesr_vcf = baseline_filtered_pesr_vcf,
        baseline_filtered_depth_vcf = baseline_filtered_depth_vcf,
        contig_list = select_first([primary_contigs_list]),
        linux_docker = linux_docker,
        sv_pipeline_base_docker = select_first([sv_pipeline_base_docker]),
        sv_base_mini_docker = sv_base_mini_docker
    }
  }

  output {
    File? filtered_manta_vcf = FilterOutlierSamples.vcfs_noOutliers[0]
    File? filtered_delly_vcf = FilterOutlierSamples.vcfs_noOutliers[1]
    File? filtered_wham_vcf = FilterOutlierSamples.vcfs_noOutliers[2]
    File? filtered_melt_vcf = FilterOutlierSamples.vcfs_noOutliers[3]
    File? filtered_depth_vcf = FilterOutlierSamples.vcfs_noOutliers[4]
    File filtered_pesr_vcf = MergePesrVcfs.merged_pesr_vcf
    File cutoffs = AdjudicateSV.cutoffs
    File scores = RewriteScores.updated_scores
    File RF_intermediate_files = AdjudicateSV.RF_intermediate_files
    Array[String] outlier_samples_excluded = FilterOutlierSamples.outlier_samples_excluded
    Array[String] batch_samples_postOutlierExclusion = FilterOutlierSamples.filtered_batch_samples_list
    File outlier_samples_excluded_file = FilterOutlierSamples.outlier_samples_excluded_file
    File batch_samples_postOutlierExclusion_file = FilterOutlierSamples.filtered_batch_samples_file

    File? metrics_file_filterbatch = FilterBatchMetrics.metrics_file
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

task UpdatePedFile {
  input {
    File ped_file
    Array[String] excluded_samples
    String batch
    String linux_docker
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
    File filtered_ped_file = "${batch}.outlier_samples_removed.fam"
  }
  command <<<

    fgrep -wvf ~{write_lines(excluded_samples)} ~{ped_file} > ~{batch}.outlier_samples_removed.fam
  
  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: linux_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

}

task MergePesrVcfs {
  input {
    File? manta_vcf
    File? delly_vcf
    File? wham_vcf
    File? melt_vcf
    String batch
    String sv_base_mini_docker
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

  Array[File] vcfs_array = select_all([manta_vcf, delly_vcf, wham_vcf, melt_vcf])

  output {
    File merged_pesr_vcf = "${batch}.filtered_pesr_merged.vcf.gz"
  }
  command <<<

    set -euo pipefail
    vcf-concat ~{sep=" " vcfs_array} \
      | vcf-sort -c \
      | bgzip -c \
      > ~{batch}.filtered_pesr_merged.vcf.gz
  
  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_base_mini_docker
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
      <(sed -e '1d' ~{scores} | fgrep -e DEL -e DUP | awk '($3>=0.5)' | cut -f1 | fgrep -w -f - <(zcat ~{vcf})) \
      <(sed -e '1d' ~{scores} | fgrep -e INV -e BND -e INS | awk '($3>=0.5)' | cut -f1 | fgrep -w -f - <(zcat ~{vcf}) | sed -e 's/SVTYPE=DEL/SVTYPE=BND/' -e 's/SVTYPE=DUP/SVTYPE=BND/' -e 's/<DEL>/<BND>/' -e 's/<DUP>/<BND>/') \
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

