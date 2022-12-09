version 1.0

import "Structs.wdl"
import "Utils.wdl" as util
import "IdentifyOutlierSamples.wdl" as identify_outliers

# Filter outlier samples by IQR or cutoff table for a single VCF. Recommended to run PlotSVCountsPerSample first to choose cutoff
workflow FilterOutlierSamples {
  input {
    String name  # batch or cohort
    File vcf
    File? sv_counts  # SV counts file from PlotSVCountsPerSample - if not provided, will create
    Int N_IQR_cutoff
    File? outlier_cutoff_table
    String? vcf_identifier  # required (enter algorithm here) if providing outlier_cutoff_table, otherwise used in some file prefixes
    String sv_pipeline_docker
    String sv_base_mini_docker
    String linux_docker
    RuntimeAttr? runtime_attr_identify_outliers
    RuntimeAttr? runtime_attr_subset_vcf
    RuntimeAttr? runtime_attr_cat_outliers
    RuntimeAttr? runtime_attr_filter_samples
    RuntimeAttr? runtime_attr_ids_from_vcf
    RuntimeAttr? runtime_attr_count_svs
  }

  call identify_outliers.IdentifyOutlierSamples {
    input:
      vcf = vcf,
      name = name,
      sv_counts = sv_counts, 
      N_IQR_cutoff = N_IQR_cutoff,
      outlier_cutoff_table = outlier_cutoff_table,
      vcf_identifier = vcf_identifier,
      sv_pipeline_docker = sv_pipeline_docker,
      linux_docker = linux_docker,
      runtime_attr_identify_outliers = runtime_attr_identify_outliers,
      runtime_attr_cat_outliers = runtime_attr_cat_outliers,
      runtime_attr_count_svs = runtime_attr_count_svs
  }

  call util.SubsetVcfBySamplesList {
    input:
      vcf = vcf,
      list_of_samples = IdentifyOutlierSamples.outlier_samples_file,
      outfile_name = "${name}.outliers_removed.vcf.gz",
      remove_samples = true,
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_subset_vcf
  }
  
  call util.GetSampleIdsFromVcf {
    input:
      vcf = vcf,
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_ids_from_vcf
  }

  # Write new list of samples without outliers
  call FilterSampleList {
    input:
      original_samples = GetSampleIdsFromVcf.out_array,
      outlier_samples = IdentifyOutlierSamples.outlier_samples_list,
      batch = name,
      linux_docker = linux_docker,
      runtime_attr_override = runtime_attr_filter_samples
  }

  output {
    File outlier_filtered_vcf = SubsetVcfBySamplesList.vcf_subset
    Array[String] filtered_samples_list = FilterSampleList.filtered_samples_list
    File filtered_samples_file = FilterSampleList.filtered_samples_file
    Array[String] outlier_samples_excluded = IdentifyOutlierSamples.outlier_samples_list
    File outlier_samples_excluded_file = IdentifyOutlierSamples.outlier_samples_file
    File sv_counts_file = IdentifyOutlierSamples.sv_counts_file
  }
}


# Write new list of samples per batch after outlier filtering
task FilterSampleList {
  input {
    Array[String] original_samples
    Array[String] outlier_samples
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
    Array[String] filtered_samples_list = read_lines("${batch}.outliers_excluded.samples.list")
    File filtered_samples_file = "${batch}.outliers_excluded.samples.list"
  }
  command <<<

    fgrep -wvf ~{write_lines(outlier_samples)} ~{write_lines(original_samples)} > ~{batch}.outliers_excluded.samples.list
  
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

