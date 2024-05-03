version 1.0

import "Structs.wdl"
import "SubsetVcfBySamples.wdl" as subset

workflow SubsetVcfBySamplesScattered {
  input {
    Array[File] vcfs
    Array[File]? vcf_indices
    File list_of_samples  # List of samples to keep (default, remove_samples = false) or remove (remove_samples = true)
    String? outfile_name
    Boolean? remove_samples  # If false (default), keep samples in provided list. If true, remove them.
    Boolean? remove_private_sites  # If true (default), remove sites that are private to excluded samples. If false, keep sites even if no remaining samples are non-ref.

    # Do not use
    File? NONE_FILE_

    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_subset_by_samples
  }

  scatter (i in range(length(vcfs))) {
    call subset.SubsetVcfBySamples {
      input:
        vcf = vcfs[i],
        vcf_index = if defined(vcf_indices) then select_first([vcf_indices])[i] else NONE_FILE_,
        list_of_samples = list_of_samples,
        outfile_name = outfile_name,
        remove_samples = remove_samples,
        remove_private_sites = remove_private_sites,
        sv_base_mini_docker = sv_base_mini_docker,
        runtime_attr_subset_by_samples = runtime_attr_subset_by_samples
    }
  }

  output {
    Array[File] vcf_subset = SubsetVcfBySamples.vcf_subset
    Array[File] vcf_subset_index = SubsetVcfBySamples.vcf_subset_index
  }
}
