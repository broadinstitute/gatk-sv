version 1.0

import "Structs.wdl"
import "Utils.wdl" as util

workflow SubsetVcfBySamples {
  input {
    File vcf
    File? vcf_index
    File list_of_samples  # List of samples to keep (default, remove_samples = false) or remove (remove_samples = true)
    String? outfile_name
    Boolean? remove_samples  # If false (default), keep samples in provided list. If true, remove them.
    Boolean? remove_private_sites  # If true (default), remove sites that are private to excluded samples. If false, keep sites even if no remaining samples are non-ref.

    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_subset_by_samples
  }

  call util.SubsetVcfBySamplesList {
    input:
      vcf = vcf,
      vcf_idx = vcf_index,
      list_of_samples = list_of_samples,
      outfile_name = outfile_name,
      remove_samples = remove_samples,
      remove_private_sites = remove_private_sites,
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_subset_by_samples
  }

  output {
    File vcf_subset = SubsetVcfBySamplesList.vcf_subset
    File vcf_subset_index = SubsetVcfBySamplesList.vcf_subset_index
  }
}
