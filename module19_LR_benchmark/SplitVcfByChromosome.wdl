version 1.0

import "Structs.wdl"
import "LongReadGenotypeTasks.wdl" as LongReadGenotypeTasks

workflow SplitVcfByChromosome {
  input {
    Array[String] chromosome_list
    File input_vcf
    File input_vcfs_idx
    String sv_pipeline_base_docker
  }


  scatter (chrom in chromosome_list) {

    call LongReadGenotypeTasks.ExtractChromosomeVcf {
        input:
          input_vcf = input_vcf,
          input_vcf_idx = input_vcfs_idx,
          chromosome = chrom,
          docker_image = sv_base_mini_docker
      }
    }

  output {
    Array[File] per_chr_vcf = ExtractChromosomeVcf.output_vcf
    Array[File] per_chr_vcf_idx = ExtractChromosomeVcf.output_vcf_idx
  }
}


