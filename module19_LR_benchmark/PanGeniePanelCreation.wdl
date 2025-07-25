version 1.0

import "PanGeniePanelCreationPerContig.wdl" as PanGeniePanelCreationPerContig
import "LongReadGenotypeTasks.wdl"

workflow PanGeniePanelCreation {
    input {

        File phased_vcf
        File phased_vcf_idx

        File reference_fasta
        
        File prepare_vcf_script
        File add_ids_script
        File merge_vcfs_script
        Float frac_missing = 0.2
        
        String output_prefix

        Array[String] contig_list

        String docker_image
        String sv_pipeline_base_docker
        String sv_base_mini_docker
        File? monitoring_script
    }

    scatter(chrom in contig_list){

        call LongReadGenotypeTasks.ExtractChromosomeVariants{
            input:
                input_vcf = phased_vcf,
                input_vcf_index = phased_vcf_idx,
                chromosome = chrom,
                output_name = "~{chrom}.vcf.gz",
                docker_image = sv_pipeline_base_docker  
            }



        call PanGeniePanelCreationPerContig.PanGeniePanelCreationPerContig as PanGeniePanelCreationPerContig{
            input:
                phased_bcf = ExtractChromosomeVariants.chr_vcf,
                phased_bcf_idx = ExtractChromosomeVariants.chr_vcf_idx,
                reference_fasta = reference_fasta,
                prepare_vcf_script = prepare_vcf_script,
                add_ids_script = add_ids_script,
                merge_vcfs_script = merge_vcfs_script,
                frac_missing = frac_missing,
                output_prefix = "~{chrom}.panel",
                docker_image =docker_image,
                monitoring_script = monitoring_script
        }
    }

    call LongReadGenotypeTasks.ConcatVcfs as concat_panel_vcf{
        input:
            vcfs = PanGeniePanelCreationPerContig.panel_vcf_gz,
            vcfs_idx = PanGeniePanelCreationPerContig.panel_vcf_gz_tbi,
            merge_sort = false,
            remove_dup = false,
            outfile_prefix = "~{output_prefix}.prepare.id.split.mergehap",
            sv_base_mini_docker = sv_base_mini_docker
    }

    call LongReadGenotypeTasks.ConcatVcfs as concat_panel_id_split_vcf{
        input:
            vcfs = PanGeniePanelCreationPerContig.panel_id_split_vcf_gz,
            vcfs_idx = PanGeniePanelCreationPerContig.panel_id_split_vcf_gz_tbi,
            merge_sort = false,
            remove_dup = false,
            outfile_prefix = "~{output_prefix}.prepare.id.split",
            sv_base_mini_docker = sv_base_mini_docker

    }


    output {
        File panel_vcf_gz = concat_panel_vcf.concat_vcf
        File panel_vcf_gz_tbi = concat_panel_vcf.concat_vcf_idx
        File panel_id_split_vcf_gz = concat_panel_id_split_vcf.concat_vcf
        File panel_id_split_vcf_gz_tbi = concat_panel_id_split_vcf.concat_vcf_idx
    }
}











