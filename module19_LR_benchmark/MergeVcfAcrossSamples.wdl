version 1.0

import "Structs.wdl"
import  "LongReadGenotypeTasks.wdl" as LongReadGenotypeTasks
import "MergeVcfAcrossSamplesPerChromosome.wdl" as MergeVcfAcrossSamplesPerChromosome

workflow MergeVcfAcrossSamples {
  input {

    Array[File] input_vcfs
    Array[File] input_vcf_idxes
    Array[String] contigs

    String output_prefix        # Prefix for final merged output

    File? monitoring_script

    String sv_base_mini_docker
    String sv_pipeline_base_docker

    RuntimeAttr? runtime_attr_merge_vcfs
    RuntimeAttr? runtime_attr_concat_vcfs
  }

  
  scatter(contig in contigs){
    call MergeVcfAcrossSamplesPerChromosome.MergeVcfAcrossSamplesPerChromosome as MergeVcfAcrossSamplesPerChromosome{
        input:
          input_vcfs = input_vcfs,
          input_vcf_idxes = input_vcf_idxes,
          contig = contig,
          monitoring_script = monitoring_script,
          sv_base_mini_docker = sv_base_mini_docker,
          sv_pipeline_base_docker = sv_pipeline_base_docker,
          runtime_attr_merge_vcfs = runtime_attr_merge_vcfs

    }
  }

  call LongReadGenotypeTasks.ConcatVcfs{
      input:
        vcfs = MergeVcfAcrossSamplesPerChromosome.merged_vcf,
        vcfs_idx = MergeVcfAcrossSamplesPerChromosome.merged_vcf_idx,
        remove_dup = false,
        outfile_prefix = output_prefix, 
        sv_base_mini_docker = sv_base_mini_docker,
        runtime_attr_override = runtime_attr_concat_vcfs
  }

  output {
    File output_vcf = ConcatVcfs.concat_vcf
    File output_vcf_idx = ConcatVcfs.concat_vcf_idx
  }
}

