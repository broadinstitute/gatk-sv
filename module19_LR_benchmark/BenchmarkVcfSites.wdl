version 1.0

import "Structs.wdl"
import "ExtractFileByIndex.wdl" as ExtractFileByIndex
import "LongReadGenotypeTasks.wdl" as LongReadGenotypeTasks
import "BenchmarkVcfSitesPerContig.wdl" as BenchmarkVcfSitesPerContig

workflow BenchmarkVcfSites{
    input{
        File query_vcf
        File query_vcf_idx  
        File ref_vcf
        File ref_vcf_idx
        Array[String] chromosomes

        File ref_dict

        File anno_script_bash
        File anno_script_helper_R
        File benchmark_script_bash
        File banchmark_script_helper_R

        File repeat_mask
        File simple_repeats
        File segmental_duplicates


        Boolean short_read_benchmark = false

        String? truvari_params
        String sv_base_mini_docker
        String sv_pipeline_base_docker
    }

    String query_prefix = basename(query_vcf, ".vcf.gz")
    String ref_prefix = basename(ref_vcf, ".vcf.gz")

    scatter (index in range(length(chromosomes))){

        call LongReadGenotypeTasks.ExtractChromosomeVariants as extract_chrom_variants_query{
            input:
                input_vcf = query_vcf,
                input_vcf_index = query_vcf_idx,
                chromosome = chromosomes[index],
                output_name = "~{query_prefix}.~{chromosomes[index]}.vcf.gz",
                docker_image = sv_pipeline_base_docker  
            }

        call LongReadGenotypeTasks.ExtractChromosomeVariants as extract_chrom_variants_ref{
            input:
                input_vcf = ref_vcf,
                input_vcf_index = ref_vcf_idx,
                chromosome = chromosomes[index],
                output_name = "~{ref_prefix}.~{chromosomes[index]}.vcf.gz",
                docker_image = sv_pipeline_base_docker  
            }

        call BenchmarkVcfSitesPerContig.BenchmarkVcfSitesPerContig{
            input:
                query_vcf = extract_chrom_variants_query.chr_vcf,
                query_vcf_idx = extract_chrom_variants_query.chr_vcf_idx,
                ref_vcf = extract_chrom_variants_ref.chr_vcf,
                ref_vcf_idx = extract_chrom_variants_ref.chr_vcf_idx,
                chromosome = chromosomes[index],
                ref_dict = ref_dict,
                short_read_benchmark = false,

                anno_script_bash = anno_script_bash,
                anno_script_helper_R = anno_script_helper_R,
                benchmark_script_bash = benchmark_script_bash,
                banchmark_script_helper_R = banchmark_script_helper_R,

                repeat_mask = repeat_mask,
                simple_repeats = simple_repeats,
                segmental_duplicates = segmental_duplicates,

                sv_base_mini_docker = sv_base_mini_docker,
                sv_pipeline_base_docker = sv_pipeline_base_docker
          }
    }

    call LongReadGenotypeTasks.ConcatVcfs as combine_vcfs_tp_ref{
      input:
        vcfs = BenchmarkVcfSitesPerContig.tp_ref,
        outfile_prefix = "~{query_prefix}.vs.~{ref_prefix}.tp_ref",
        sv_base_mini_docker = sv_base_mini_docker
    }

    call LongReadGenotypeTasks.ConcatVcfs as combine_vcfs_tp_query{
      input:
        vcfs = BenchmarkVcfSitesPerContig.tp_query,
        outfile_prefix = "~{query_prefix}.vs.~{ref_prefix}.tp_query",
        sv_base_mini_docker = sv_base_mini_docker
    }

    call LongReadGenotypeTasks.ConcatVcfs as combine_vcfs_fp_ref{
      input:
        vcfs = BenchmarkVcfSitesPerContig.fp_ref,
        outfile_prefix = "~{query_prefix}.vs.~{ref_prefix}.fp_ref",
        sv_base_mini_docker = sv_base_mini_docker
    }

    call LongReadGenotypeTasks.ConcatVcfs as combine_vcfs_fp_query{
      input:
        vcfs = BenchmarkVcfSitesPerContig.fp_query,
        outfile_prefix = "~{query_prefix}.vs.~{ref_prefix}.fp_query",
        sv_base_mini_docker = sv_base_mini_docker
    }

    call LongReadGenotypeTasks.CalcuCompStat as calcu_comp_stat{
      input:
        tp_query = combine_vcfs_tp_query.concat_vcf,
        tp_ref = combine_vcfs_tp_ref.concat_vcf,
        fp_query = combine_vcfs_fp_query.concat_vcf,
        fp_ref = combine_vcfs_fp_ref.concat_vcf,
        docker_image = sv_pipeline_base_docker
      }

    call LongReadGenotypeTasks.PlotCompResults as plot_comp_results{
        input:
            comp_stat = calcu_comp_stat.comp_stat,
            docker_image = sv_pipeline_base_docker
        }

  output {
    File tp_query = combine_vcfs_tp_query.concat_vcf  
    File tp_ref = combine_vcfs_tp_ref.concat_vcf
    File fp_query = combine_vcfs_fp_query.concat_vcf  
    File fp_ref = combine_vcfs_fp_ref.concat_vcf
    File benchmark_stat = calcu_comp_stat.comp_stat
    File benchmark_plot = plot_comp_results.figure
  }

}

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}
