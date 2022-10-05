version 1.0

import "Structs.wdl"
import "Utils.wdl" as util

workflow SubsetVcfBySamplesList {
  input {
    File vcf
    File? vcf_index
    File list_of_samples_to_keep
    String? subset_name

    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_subset_by_samples
  }

  call util.SubsetVcfBySamplesList {
    input:
      vcf = vcf,
      vcf_idx = vcf_index,
      list_of_samples_to_keep = list_of_samples_to_keep,
      subset_name = subset_name,
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_subset_by_samples
  }

  output {
    File vcf_subset = SubsetVcfBySamplesList.vcf_subset
    File vcf_subset_index = SubsetVcfBySamplesList.vcf_subset_idx
  }
}
