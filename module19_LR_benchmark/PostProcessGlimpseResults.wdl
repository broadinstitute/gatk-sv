version 1.0

import "Structs.wdl"
import "PostProcessGlimpseResultsPerContig.wdl" as PostProcessGlimpseResultsPerContig
import "LongReadGenotypeTasks.wdl" as LongReadGenotypeTasks

workflow PostProcessGlimpseResults {
    input {
        File ref_panel_vcf
        File ref_panel_vcf_tbi

        File genotype_output_vcf
        File genotype_output_vcf_tbi

        Array[File] panel_biallelic_vcf_list
        Array[File] panel_biallelic_vcf_idx_list

        Array[String] chromosomes
        String output_prefix        

        File convert_to_biallelic_script
        File? monitoring_script

        String sv_base_mini_docker
        String sv_pipeline_base_docker

        RuntimeAttr? runtime_attr_extract_chrom_vcf_ref_panel
        RuntimeAttr? runtime_attr_split_ref_panel
        RuntimeAttr? runtime_attr_add_id_to_info_column
        RuntimeAttr? runtime_attr_convert_bubbles_to_biallelic
        RuntimeAttr? runtime_attr_concat_vcf

    }


    scatter (i in range(length(chromosomes))){

        call LongReadGenotypeTasks.ExtractChromosomeVcf as extract_chrom_vcf_ref_panel {
            input:
              input_vcf = ref_panel_vcf,
              input_vcf_idx = ref_panel_vcf_tbi,
              chromosome = chromosomes[i],
              docker_image = sv_base_mini_docker,
              runtime_attr_override = runtime_attr_extract_chrom_vcf_ref_panel
        }

        call LongReadGenotypeTasks.ExtractChromosomeVcf as extract_chrom_vcf_genotype_output {
            input:
              input_vcf = genotype_output_vcf,
              input_vcf_idx = genotype_output_vcf_tbi,
              chromosome = chromosomes[i],
              docker_image = sv_base_mini_docker,
              runtime_attr_override = runtime_attr_extract_chrom_vcf_ref_panel
        }


        call PostProcessGlimpseResultsPerContig.PostProcessGlimpseResultsPerContig as post_process_glimpse_output_per_contig{
            input:
                ref_panel_vcf = extract_chrom_vcf_ref_panel.output_vcf,
                ref_panel_vcf_tbi = extract_chrom_vcf_ref_panel.output_vcf_idx,
                genotype_output_vcf = extract_chrom_vcf_genotype_output.output_vcf,
                genotype_output_vcf_tbi = extract_chrom_vcf_genotype_output.output_vcf_idx,
                panel_biallelic_vcf = panel_biallelic_vcf_list[i],
                panel_biallelic_vcf_idx = panel_biallelic_vcf_idx_list[i],
                convert_to_biallelic_script = convert_to_biallelic_script,
                sv_pipeline_base_docker = sv_pipeline_base_docker,
                monitoring_script = monitoring_script,
                runtime_attr_split_ref_panel = runtime_attr_split_ref_panel,
                runtime_attr_add_id_to_info_column = runtime_attr_add_id_to_info_column,
                runtime_attr_convert_bubbles_to_biallelic = runtime_attr_convert_bubbles_to_biallelic
        }
    }

    call LongReadGenotypeTasks.ConcatVcfs {
        input:
            vcfs = post_process_glimpse_output_per_contig.output_biallelic_vcf,
            vcfs_idx = post_process_glimpse_output_per_contig.output_biallelic_vcf_idx,
            remove_dup = false,

            outfile_prefix = "~{output_prefix}.biallelic",
            sv_base_mini_docker = sv_base_mini_docker,
            runtime_attr_override = runtime_attr_concat_vcf
    }

    output{
        File output_vcf = ConcatVcfs.concat_vcf,
        File output_vcf_idx = ConcatVcfs.concat_vcf_idx
    }
}






