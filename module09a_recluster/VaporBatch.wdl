version 1.0

import "Structs.wdl"
import "Vapor.wdl" as vapor_bed

workflow VaporBatch {
  input {
    Array[String] samples
    Array[File] bam_or_cram_files
    Array[File] bam_or_cram_indexes
    File bed_file  # Multi-sample bed file, generated with svtk vcf2bed (e.g. as in MainVcfQc)

    File ref_fasta
    File ref_fai
    File ref_dict
    File contigs

    String vapor_docker
    String sv_base_mini_docker
    String sv_pipeline_docker

    RuntimeAttr? runtime_attr_vapor
    RuntimeAttr? runtime_attr_bcf2vcf
    RuntimeAttr? runtime_attr_vcf2bed
    RuntimeAttr? runtime_attr_split_vcf
    RuntimeAttr? runtime_attr_concat_beds
  }

  scatter (i in range(length(bam_or_cram_files))) {
    call vapor_bed.Vapor {
      input:
        prefix = samples[i],
        bam_or_cram_file = bam_or_cram_files[i],
        bam_or_cram_index = bam_or_cram_indexes[i],
        bed_file = bed_file,
        sample_id = samples[i],
        ref_fasta = ref_fasta,
        ref_fai = ref_fai,
        ref_dict = ref_dict,
        contigs = contigs,
        vapor_docker = vapor_docker,
        sv_base_mini_docker = sv_base_mini_docker,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_vapor = runtime_attr_vapor,
        runtime_attr_bcf2vcf = runtime_attr_bcf2vcf,
        runtime_attr_vcf2bed = runtime_attr_vcf2bed,
        runtime_attr_split_vcf = runtime_attr_split_vcf,
        runtime_attr_concat_beds = runtime_attr_concat_beds
    }
  }
  output {
    Array[File] vapor_batch_beds = Vapor.vapor_bed
    Array[File] vapor_batch_plots = Vapor.vapor_plots
  }
}


