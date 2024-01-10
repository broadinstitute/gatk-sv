version 1.0

import "FilterOutlierSamples.wdl" as filter_outliers
import "Structs.wdl"
import "Utils.wdl" as util
import "IdentifyOutlierSamples.wdl" as identify_outliers
import "TasksMakeCohortVcf.wdl" as tasks_mcv

# Workflow to identify & filter outliers from VCFs as part of FilterBatch after FilterBatchSites & PlotSVCountsPerSample
workflow FilterBatchSamples {
  input {
    String batch
    File? manta_vcf
    File? wham_vcf
    File? melt_vcf
    File? scramble_vcf
    File? depth_vcf
    File? manta_counts  # SV counts files from PlotSVCountsPerSample. If not provided, SV counts will be calculated as part of this workflow
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
    RuntimeAttr? runtime_attr_subset_vcf
    RuntimeAttr? runtime_attr_cat_outliers
    RuntimeAttr? runtime_attr_filter_samples
    RuntimeAttr? runtime_attr_ids_from_vcf
    RuntimeAttr? runtime_attr_count_svs
    RuntimeAttr? runtime_attr_merge_pesr_vcfs
  }

  Array[File?] vcfs = [manta_vcf, wham_vcf, melt_vcf, scramble_vcf, depth_vcf]
  Array[String] algorithms = ["manta", "wham", "melt", "scramble", "depth"]  # fixed algorithms to enable File outputs to be determined
  Int num_algorithms = length(algorithms)
  Array[File?] sv_counts_ = [manta_counts, wham_counts, melt_counts, scramble_counts, depth_counts]

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
      call util.SubsetVcfBySamplesList {
        input:
          vcf = select_first([vcfs[i]]),
          list_of_samples = CatOutliers.outliers_file,
          outfile_name = "${batch}.${algorithms[i]}.outliers_removed.vcf.gz",
          remove_samples = true,
          sv_base_mini_docker = sv_base_mini_docker,
          runtime_attr_override = runtime_attr_subset_vcf
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

  Array[File] pesr_vcfs_no_outliers = select_all([SubsetVcfBySamplesList.vcf_subset[0], SubsetVcfBySamplesList.vcf_subset[1], SubsetVcfBySamplesList.vcf_subset[2], SubsetVcfBySamplesList.vcf_subset[3]])
  Array[File] pesr_vcfs_no_outliers_index = select_all([SubsetVcfBySamplesList.vcf_subset_index[0], SubsetVcfBySamplesList.vcf_subset_index[1], SubsetVcfBySamplesList.vcf_subset_index[2], SubsetVcfBySamplesList.vcf_subset_index[3]])
  call tasks_mcv.ConcatVcfs as MergePesrVcfs {
    input:
      vcfs=pesr_vcfs_no_outliers,
      vcfs_idx=pesr_vcfs_no_outliers_index,
      allow_overlaps=true,
      outfile_prefix="~{batch}.filtered_pesr_merged",
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_attr_merge_pesr_vcfs
  }

  output {
    File? outlier_filtered_manta_vcf = SubsetVcfBySamplesList.vcf_subset[0]
    File? outlier_filtered_wham_vcf = SubsetVcfBySamplesList.vcf_subset[1]
    File? outlier_filtered_melt_vcf = SubsetVcfBySamplesList.vcf_subset[2]
    File? outlier_filtered_scramble_vcf = SubsetVcfBySamplesList.vcf_subset[3]
    File? outlier_filtered_depth_vcf = SubsetVcfBySamplesList.vcf_subset[4]

    File? outlier_filtered_manta_vcf_index = SubsetVcfBySamplesList.vcf_subset_index[0]
    File? outlier_filtered_wham_vcf_index = SubsetVcfBySamplesList.vcf_subset_index[1]
    File? outlier_filtered_melt_vcf_index = SubsetVcfBySamplesList.vcf_subset_index[2]
    File? outlier_filtered_scramble_vcf_index = SubsetVcfBySamplesList.vcf_subset_index[3]
    File? outlier_filtered_depth_vcf_index = SubsetVcfBySamplesList.vcf_subset_index[4]

    File outlier_filtered_pesr_vcf = MergePesrVcfs.concat_vcf
    File outlier_filtered_pesr_vcf_index = MergePesrVcfs.concat_vcf_idx
    Array[String] filtered_batch_samples_list = FilterSampleList.filtered_samples_list
    File filtered_batch_samples_file = FilterSampleList.filtered_samples_file
    Array[String] outlier_samples_excluded = CatOutliers.outliers_list
    File outlier_samples_excluded_file = CatOutliers.outliers_file
  }
}
