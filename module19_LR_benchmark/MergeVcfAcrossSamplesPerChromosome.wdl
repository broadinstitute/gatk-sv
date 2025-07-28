version 1.0

import "Structs.wdl"
import  "LongReadGenotypeTasks.wdl" as LongReadGenotypeTasks

workflow MergeVcfAcrossSamplesPerChromosome {
  input {

    Array[File] input_vcfs
    Array[File] input_vcf_idxes
    String contig

    File? monitoring_script
    
    String sv_base_mini_docker
    String sv_pipeline_base_docker

    RuntimeAttr? runtime_attr_merge_vcfs
  }

  scatter (i in range(length(input_vcfs))){
    call LongReadGenotypeTasks.ExtractChromosomeVcf {
        input:
          input_vcf = input_vcfs[i],
          input_vcf_idx = input_vcf_idxes[i],
          chromosome = contig,
          docker_image = sv_pipeline_base_docker,
          monitoring_script = monitoring_script
    }
  }

  call LongReadGenotypeTasks.MergeVcfs {
    input:
      vcfs = ExtractChromosomeVcf.output_vcf,
      vcfs_idx = ExtractChromosomeVcf.output_vcf_idx,
      output_name = "~{contig}.merged.vcf.gz"
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_merge_vcfs
  }

  output {
    File merged_vcf = MergeVcfs.output_merged_vcf
    File merged_vcf_idx = MergeVcfs.output_merged_vcf_idx
  }
}
