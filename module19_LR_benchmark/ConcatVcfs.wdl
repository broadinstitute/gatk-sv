version 1.0

import "Structs.wdl"
import "MergeVcfs.wdl" as MergeVcfs
workflow ConcatVcfs {
  input {
    Array[File] input_vcfs     
    Array[File]? input_vcfs_idx
    String output_prefix 
    String sv_base_mini_docker
    String sv_pipeline_base_docker

    RuntimeAttr? runtime_attr_merge_vcfs
    RuntimeAttr? runtime_attr_concat_vcfs
    RuntimeAttr? runtime_attr_extract_chromosome_vcf

  }



  if (!defined(input_vcfs_idx)) {
    scatter (idx in range(length(input_vcfs))) {
      call MergeVcfs.IndexVcf as IndexVcf{
        input:
          vcf = input_vcfs[idx],
          sv_base_mini_docker = sv_base_mini_docker
      }
    }
  }

  Array[File] vcfs_idx = select_first([IndexVcf.indexed_vcf_idx,input_vcfs_idx])

  call MergeVcfs.ConcatVcfs {
    input:
      input_vcfs = input_vcfs,
      input_vcfs_idx = vcfs_idx,
      output_name = "${output_prefix}.vcf.gz",
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_concat_vcfs
  }

  output {
    File final_merged_vcf = ConcatVcfs.output_vcf
    File final_merged_vcf_index = ConcatVcfs.output_vcf_idx
  }
}
