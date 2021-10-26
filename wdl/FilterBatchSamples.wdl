version 1.0

import "FilterOutlierSamples.wdl" as filter_outliers
import "Structs.wdl"
import "Utils.wdl" as util
import "IdentifyOutlierSamples.wdl" as identify_outliers

# Workflow to identify & filter outliers from VCFs as part of FilterBatch after FilterBatchSites & PlotSVCountsPerSample
workflow FilterBatchSamples {
  input {
    String batch
    File? manta_vcf
    File? delly_vcf
    File? wham_vcf
    File? melt_vcf
    File? scramble_vcf
    File? depth_vcf
    File? manta_counts  # SV counts files from PlotSVCountsPerSample. If not provided, SV counts will be calculated as part of this workflow
    File? delly_counts
    File? wham_counts
    File? melt_counts
    File? scramble_counts
    File? depth_counts
    Int N_IQR_cutoff
    File? outlier_cutoff_table
    String sv_pipeline_docker
    String sv_base_mini_docker
    String linux_docker
    RuntimeAttr? runtime_attr_identify_outliers
    RuntimeAttr? runtime_attr_exclude_outliers
    RuntimeAttr? runtime_attr_cat_outliers
    RuntimeAttr? runtime_attr_filter_samples
    RuntimeAttr? runtime_attr_ids_from_vcf
    RuntimeAttr? runtime_attr_count_svs
    RuntimeAttr? runtime_attr_merge_pesr_vcfs
  }

  Array[File?] vcfs = [manta_vcf, delly_vcf, wham_vcf, melt_vcf, scramble_vcf, depth_vcf]
  Array[String] algorithms = ["manta", "delly", "wham", "melt", "scramble", "depth"]  # fixed algorithms to enable File outputs to be determined
  Int num_algorithms = length(algorithms)
  Array[File?] sv_counts_ = [manta_counts, delly_counts, wham_counts, melt_counts, scramble_counts, depth_counts]

  scatter (i in range(num_algorithms)) {
    if (defined(vcfs[i])) {
      call identify_outliers.IdentifyOutlierSamples {
        input:
          vcf = select_first([vcfs[i]]),
          name = batch,
          sv_counts = sv_counts_[i],
          N_IQR_cutoff = N_IQR_cutoff,
          outlier_cutoff_table = outlier_cutoff_table,
          vcf_identifier = algorithms[i],
          sv_pipeline_docker = sv_pipeline_docker,
          linux_docker = linux_docker,
          runtime_attr_identify_outliers = runtime_attr_identify_outliers,
          runtime_attr_cat_outliers = runtime_attr_cat_outliers,
          runtime_attr_count_svs = runtime_attr_count_svs
      }
    }
  }

  # Merge list of outliers from all algorithms
  call identify_outliers.CatOutliers {
    input:
      outliers = select_all(IdentifyOutlierSamples.outlier_samples_file),
      batch = batch,
      linux_docker = linux_docker,
      runtime_attr_override = runtime_attr_cat_outliers
  }

  scatter (i in range(num_algorithms)) {
    if (defined(vcfs[i])) {
      call filter_outliers.ExcludeOutliers {
        input:
          vcf = select_first([vcfs[i]]),
          outliers_list = CatOutliers.outliers_list,
          outfile = "${batch}.${algorithms[i]}.outliers_removed.vcf.gz",
          sv_base_mini_docker = sv_base_mini_docker,
          runtime_attr_override = runtime_attr_exclude_outliers
      }
    }
  }
  
  call util.GetSampleIdsFromVcf {
    input:
      vcf = select_first(vcfs),
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_ids_from_vcf
  }

  # Write new list of samples without outliers
  call filter_outliers.FilterSampleList {
    input:
      original_samples = GetSampleIdsFromVcf.out_array,
      outlier_samples = CatOutliers.outliers_list,
      batch = batch,
      linux_docker = linux_docker,
      runtime_attr_override = runtime_attr_filter_samples
  }

  call MergePesrVcfs {
    input:
      manta_vcf = ExcludeOutliers.vcf_no_outliers[0],
      delly_vcf = ExcludeOutliers.vcf_no_outliers[1],
      wham_vcf = ExcludeOutliers.vcf_no_outliers[2],
      melt_vcf = ExcludeOutliers.vcf_no_outliers[3],
      scramble_vcf = ExcludeOutliers.vcf_no_outliers[4],
      batch = batch,
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_merge_pesr_vcfs
  }

  output {
    File? outlier_filtered_manta_vcf = ExcludeOutliers.vcf_no_outliers[0]
    File? outlier_filtered_delly_vcf = ExcludeOutliers.vcf_no_outliers[1]
    File? outlier_filtered_wham_vcf = ExcludeOutliers.vcf_no_outliers[2]
    File? outlier_filtered_melt_vcf = ExcludeOutliers.vcf_no_outliers[3]
    File? outlier_filtered_scramble_vcf = ExcludeOutliers.vcf_no_outliers[4]
    File? outlier_filtered_depth_vcf = ExcludeOutliers.vcf_no_outliers[5]
    File? outlier_filtered_pesr_vcf = MergePesrVcfs.merged_pesr_vcf
    Array[String] filtered_batch_samples_list = FilterSampleList.filtered_samples_list
    File filtered_batch_samples_file = FilterSampleList.filtered_samples_file
    Array[String] outlier_samples_excluded = CatOutliers.outliers_list
    File outlier_samples_excluded_file = CatOutliers.outliers_file
  }
}

task MergePesrVcfs {
  input {
    File? manta_vcf
    File? delly_vcf
    File? wham_vcf
    File? melt_vcf
    File? scramble_vcf
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

  Array[File] vcfs_array = select_all([manta_vcf, delly_vcf, wham_vcf, melt_vcf, scramble_vcf])

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

