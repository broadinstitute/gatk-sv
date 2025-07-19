version 1.0

import "Structs.wdl"
import "MergeVcfsByChromosome.wdl" as MergeVcfsByChromosome

workflow ConvertBubblesToBiallelicByChromosome {
    input{
        File input_vcf
        File input_vcf_idx

        Array[File] panel_biallelic_vcf_list
        Array[File] panel_biallelic_vcf_idx_list

        Array[String] chromosomes

        File convert_to_biallelic_script

        String sv_base_mini_docker
        String sv_pipeline_base_docker

        RuntimeAttr? runtime_attr_extract_chrom_vcf
        RuntimeAttr? runtime_attr_convert_bubbles_to_biallelic
        RuntimeAttr? runtime_attr_concat_biallelic_vcf
        
    }

  scatter (i in range(length(chromosomes))) {

    call LongReadGenotypeTasks.ExtractChromosomeVcf {
        input:
          input_vcf = input_vcf,
          input_vcf_idx = input_vcf_idx,
          chromosome = chromosomes[i],
          docker_image = sv_base_mini_docker,
          runtime_attr_override = runtime_attr_extract_chrom_vcf
        }


    call LongReadGenotypeTasks.ConvertBubblesToBiallelic{
        input:
            input_vcf = ExtractChromosomeVcf.output_vcf,
            input_vcf_idx = ExtractChromosomeVcf.output_vcf_idx,

            panel_biallelic_vcf = panel_biallelic_vcf_list[i],
            panel_biallelic_vcf_idx = panel_biallelic_vcf_idx_list[i],

            convert_to_biallelic_script = convert_to_biallelic_script,
            docker_image = sv_pipeline_base_docker,
            runtime_attr_override = runtime_attr_convert_bubbles_to_biallelic
        }
    }

    String prefix = basename(input_vcf, ".vcf.gz")

    call LongReadGenotypeTasks.ConcatVcfs {
        input:
            vcfs = ConvertBubblesToBiallelic.biallelic_vcf,
            vcfs_idx = ConvertBubblesToBiallelic.biallelic_vcf_idx,
            remove_dup = false,

            outfile_prefix = "~{prefix}.biallelic",
            sv_base_mini_docker = sv_base_mini_docker,
            runtime_attr_override = runtime_attr_concat_biallelic_vcf
    }


    output{
        File converted_biallelic_vcf = ConcatVcfs.concat_vcf
        File converted_biallelic_vcf_idx = ConcatVcfs.concat_vcf_idx
    }
}









