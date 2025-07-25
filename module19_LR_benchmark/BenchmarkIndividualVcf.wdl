version 1.0

import "Structs.wdl"
import "ExtractFileByIndex.wdl" as ExtractFileByIndex
import "LongReadGenotypeTasks.wdl" as LongReadGenotypeTasks
import "BenchmarkIndividualVcfPerContig.wdl" as BenchmarkIndividualVcfPerContig


workflow BenchmarkIndividualVcf{
    input{
        File ref_vcf
        File ref_vcf_idx
        File? ref_filter_vcf
        File? ref_filter_vcf_idx
        File query_vcf
        File query_vcf_idx  
        File? query_filter_vcf
        File? query_filter_vcf_idx

        Array[String] chromosomes

        File anno_script_bash
        File anno_script_Rscript

        File anno_script_bash
        File anno_script_helper_R
        File benchmark_script_bash
        File banchmark_script_helper_R


        File repeat_mask
        File simple_repeats
        File segmental_duplicates

        Array[String] sample_ids
        File ref_dict

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

        if (defined(query_filter_vcf) && defined(ref_filter_vcf)){

          call LongReadGenotypeTasks.ExtractChromosomeVariants as extract_chrom_variants_query_filter_a{
              input:
                  input_vcf = query_filter_vcf,
                  input_vcf_index = query_filter_vcf_idx,
                  chromosome = chromosomes[index],
                  output_name = "~{query_prefix}.~{chromosomes[index]}.filter.vcf.gz",
                  docker_image = sv_pipeline_base_docker  
            }

          call LongReadGenotypeTasks.ExtractChromosomeVariants as extract_chrom_variants_ref_filter_a{
            input:
                input_vcf = ref_filter_vcf,
                input_vcf_index = ref_filter_vcf_idx,
                chromosome = chromosomes[index],
                output_name = "~{ref_prefix}.~{chromosomes[index]}.filter.vcf.gz",
                docker_image = sv_pipeline_base_docker  
            }  
        
          call BenchmarkIndividualVcfPerContig.BenchmarkIndividualVcfPerContig as benchmark_individual_vcf_per_contig_filter_query_ref{
              input:
                query_vcf = extract_chrom_variants_query.chr_vcf,
                query_vcf_idx = extract_chrom_variants_query.chr_vcf_idx, 

                ref_vcf = extract_chrom_variants_ref.chr_vcf,
                ref_vcf_idx = extract_chrom_variants_ref.chr_vcf_idx,

                query_filter_vcf = extract_chrom_variants_query_filter_a.chr_vcf, 
                query_filter_vcf_idx = extract_chrom_variants_query_filter_a.chr_vcf_idx,

                ref_filter_vcf = extract_chrom_variants_ref_filter_a.chr_vcf, 
                ref_filter_vcf_idx = extract_chrom_variants_ref_filter_a.chr_vcf_idx, 

                chromosome = chromosomes[index],

                anno_script_bash = anno_script_bash,
                anno_script_helper_R = anno_script_helper_R,
                benchmark_script_bash = benchmark_script_bash,
                banchmark_script_helper_R = banchmark_script_helper_R,

                repeat_mask = repeat_mask,
                simple_repeats = simple_repeats,
                segmental_duplicates = segmental_duplicates,

                short_read_benchmark = short_read_benchmark,

                sample_ids = sample_ids,
                ref_dict = ref_dict,
                truvari_params = truvari_params,
                sv_base_mini_docker = sv_base_mini_docker,
                sv_pipeline_base_docker = sv_pipeline_base_docker
          }
        
        }

        if (defined(query_filter_vcf) && !defined(ref_filter_vcf)){
          call LongReadGenotypeTasks.ExtractChromosomeVariants as extract_chrom_variants_query_filter_b{
              input:
                  input_vcf = query_filter_vcf,
                  input_vcf_index = query_filter_vcf_idx,
                  chromosome = chromosomes[index],
                  output_name = "~{query_prefix}.~{chromosomes[index]}.filter.vcf.gz",
                  docker_image = sv_pipeline_base_docker  
            }

          call BenchmarkIndividualVcfPerContig.BenchmarkIndividualVcfPerContig as benchmark_individual_vcf_per_contig_filter_query{
              input:
                query_vcf = extract_chrom_variants_query.chr_vcf,
                query_vcf_idx = extract_chrom_variants_query.chr_vcf_idx, 

                ref_vcf = extract_chrom_variants_ref.chr_vcf,
                ref_vcf_idx = extract_chrom_variants_ref.chr_vcf_idx,

                query_filter_vcf = extract_chrom_variants_query_filter_b.chr_vcf, 
                query_filter_vcf_idx = extract_chrom_variants_query_filter_b.chr_vcf_idx,

                chromosome = chromosomes[index],

                anno_script_bash = anno_script_bash,
                anno_script_helper_R = anno_script_helper_R,
                benchmark_script_bash = benchmark_script_bash,
                banchmark_script_helper_R = banchmark_script_helper_R,

                repeat_mask = repeat_mask,
                simple_repeats = simple_repeats,
                segmental_duplicates = segmental_duplicates,

                short_read_benchmark = short_read_benchmark,

                sample_ids = sample_ids,
                ref_dict = ref_dict,
                truvari_params = truvari_params,
                sv_base_mini_docker = sv_base_mini_docker,
                sv_pipeline_base_docker = sv_pipeline_base_docker
          }

        }

        if (!defined(query_filter_vcf) && defined(ref_filter_vcf)){

          call LongReadGenotypeTasks.ExtractChromosomeVariants as extract_chrom_variants_ref_filter_b{
            input:
                input_vcf = ref_filter_vcf,
                input_vcf_index = ref_filter_vcf_idx,
                chromosome = chromosomes[index],
                output_name = "~{ref_prefix}.~{chromosomes[index]}.filter.vcf.gz",
                docker_image = sv_pipeline_base_docker  
            }  
        
          call BenchmarkIndividualVcfPerContig.BenchmarkIndividualVcfPerContig as benchmark_individual_vcf_per_contig_filter_ref{
              input:
                query_vcf = extract_chrom_variants_query.chr_vcf,
                query_vcf_idx = extract_chrom_variants_query.chr_vcf_idx, 

                ref_vcf = extract_chrom_variants_ref.chr_vcf,
                ref_vcf_idx = extract_chrom_variants_ref.chr_vcf_idx,

                ref_filter_vcf = extract_chrom_variants_ref_filter_b.chr_vcf, 
                ref_filter_vcf_idx = extract_chrom_variants_ref_filter_b.chr_vcf_idx, 

                chromosome = chromosomes[index],

                anno_script_bash = anno_script_bash,
                anno_script_helper_R = anno_script_helper_R,
                benchmark_script_bash = benchmark_script_bash,
                banchmark_script_helper_R = banchmark_script_helper_R,

                repeat_mask = repeat_mask,
                simple_repeats = simple_repeats,
                segmental_duplicates = segmental_duplicates,

                short_read_benchmark = short_read_benchmark,

                sample_ids = sample_ids,
                ref_dict = ref_dict,
                truvari_params = truvari_params,
                sv_base_mini_docker = sv_base_mini_docker,
                sv_pipeline_base_docker = sv_pipeline_base_docker
          }
        
        }

        if (!defined(query_filter_vcf) && !defined(ref_filter_vcf)){
        
          call BenchmarkIndividualVcfPerContig.BenchmarkIndividualVcfPerContig as benchmark_individual_vcf_per_contig{
              input:
                query_vcf = extract_chrom_variants_query.chr_vcf,
                query_vcf_idx = extract_chrom_variants_query.chr_vcf_idx, 

                ref_vcf = extract_chrom_variants_ref.chr_vcf,
                ref_vcf_idx = extract_chrom_variants_ref.chr_vcf_idx,

                chromosome = chromosomes[index],

                anno_script_bash = anno_script_bash,
                anno_script_helper_R = anno_script_helper_R,
                benchmark_script_bash = benchmark_script_bash,
                banchmark_script_helper_R = banchmark_script_helper_R,

                repeat_mask = repeat_mask,
                simple_repeats = simple_repeats,
                segmental_duplicates = segmental_duplicates,

                short_read_benchmark = short_read_benchmark,

                sample_ids = sample_ids,
                ref_dict = ref_dict,
                truvari_params = truvari_params,
                sv_base_mini_docker = sv_base_mini_docker,
                sv_pipeline_base_docker = sv_pipeline_base_docker
          }

        }
      }







    scatter (i in range(length(sample_ids))){
        call ExtractFileByIndex.ExtractFileByIndex as extract_tp_query{
          input:
            index = i,
            nested_files = select_first([benchmark_individual_vcf_per_contig_filter_query_ref.tp_query,benchmark_individual_vcf_per_contig_filter_query.tp_query,benchmark_individual_vcf_per_contig_filter_ref.tp_query,benchmark_individual_vcf_per_contig.tp_query])
        }

        call LongReadGenotypeTasks.ConcatVcfs as combine_vcfs_tp_query{
          input:
            vcfs = extract_tp_query.processed_files,
            outfile_prefix = "~{sample_ids[i]}.~{query_prefix}.vs.~{ref_prefix}.tp_query",
            sv_base_mini_docker = sv_base_mini_docker
        }

        call ExtractFileByIndex.ExtractFileByIndex as extract_tp_ref{
          input:
            index = i,
            nested_files = select_first([benchmark_individual_vcf_per_contig_filter_query_ref.tp_ref,benchmark_individual_vcf_per_contig_filter_query.tp_ref,benchmark_individual_vcf_per_contig_filter_ref.tp_ref,benchmark_individual_vcf_per_contig.tp_ref])
        }

        call LongReadGenotypeTasks.ConcatVcfs as combine_vcfs_tp_ref{
          input:
            vcfs = extract_tp_ref.processed_files,
            outfile_prefix = "~{sample_ids[i]}.~{query_prefix}.vs.~{ref_prefix}.tp_ref",
            sv_base_mini_docker = sv_base_mini_docker
        }

        call ExtractFileByIndex.ExtractFileByIndex as extract_fp_query{
          input:
            index = i,
            nested_files = select_first([benchmark_individual_vcf_per_contig_filter_query_ref.fp_query,benchmark_individual_vcf_per_contig_filter_query.fp_query,benchmark_individual_vcf_per_contig_filter_ref.fp_query,benchmark_individual_vcf_per_contig.fp_query])
        }

        call LongReadGenotypeTasks.ConcatVcfs as combine_vcfs_fp_query{
          input:
            vcfs = extract_fp_query.processed_files,
            outfile_prefix = "~{sample_ids[i]}.~{query_prefix}.vs.~{ref_prefix}.fp_query",
            sv_base_mini_docker = sv_base_mini_docker
        }

        call ExtractFileByIndex.ExtractFileByIndex as extract_fp_ref{
          input:
            index = i,
            nested_files = select_first([benchmark_individual_vcf_per_contig_filter_query_ref.fp_ref,benchmark_individual_vcf_per_contig_filter_query.fp_ref,benchmark_individual_vcf_per_contig_filter_ref.fp_ref,benchmark_individual_vcf_per_contig.fp_ref])
        }

        call LongReadGenotypeTasks.ConcatVcfs as combine_vcfs_fp_ref{
          input:
            vcfs = extract_fp_ref.processed_files,
            outfile_prefix = "~{sample_ids[i]}.~{query_prefix}.vs.~{ref_prefix}.fp_ref",
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
    }

  output {
    Array[File] tp_query = combine_vcfs_tp_query.concat_vcf  
    Array[File] tp_ref = combine_vcfs_tp_ref.concat_vcf
    Array[File] fp_query = combine_vcfs_fp_query.concat_vcf  
    Array[File] fp_ref = combine_vcfs_fp_ref.concat_vcf
    Array[File] benchmark_stat = calcu_comp_stat.comp_stat
    Array[File] benchmark_figure = plot_comp_results.figure
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
