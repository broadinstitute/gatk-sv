version 1.0

import "ExtractFileByIndex.wdl" as ExtractFileByIndex
import "LongReadGenotypeTasks.wdl" as LongReadGenotypeTasks
import "BenchmarkIndividualVcfPerContig.wdl" as BenchmarkIndividualVcfPerContig


workflow BenchmarkIndividualVcf{
    input{
        File query_vcf
        File query_vcf_idx  
        File ref_vcf
        File ref_vcf_idx
        Array[String] chromosomes

        File anno_script_bash
        File anno_script_Rscript
        File repeat_mask
        File simple_repeats
        File segmental_duplicates

        Array[String] sample_ids
        File ref_dict

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

        call BenchmarkIndividualVcfPerContig.BenchmarkIndividualVcfPerContig{
            input:
                query_vcf = extract_chrom_variants_query.chr_vcf,
                ref_vcf = extract_chrom_variants_ref.chr_vcf,
                chromosome = chromosomes[index],
                anno_script_bash = anno_script_bash,
                anno_script_Rscript = anno_script_Rscript,
                repeat_mask = repeat_mask,
                simple_repeats = simple_repeats,
                segmental_duplicates = segmental_duplicates,
                sample_ids = sample_ids,
                ref_dict = ref_dict,
                truvari_params = truvari_params,
                sv_pipeline_base_docker = sv_pipeline_base_docker
        }

    }

    scatter (i in range(length(sample_ids))){
        call ExtractFileByIndex.ExtractFileByIndex as extract_tp_query{
          input:
            index = i,
            nested_files = BenchmarkIndividualVcfPerContig.tp_query
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
            nested_files = BenchmarkIndividualVcfPerContig.tp_ref
        }

        call LongReadGenotypeTasks.ConcatVcfs as combine_vcfs_tp_ref{
          input:
            vcfs = extract_tp_ref.processed_files,
            outfile_prefix = "~{sample_ids[i]}.~{ref_prefix}.vs.~{ref_prefix}.tp_ref",
            sv_base_mini_docker = sv_base_mini_docker
        }

        call ExtractFileByIndex.ExtractFileByIndex as extract_fp_query{
          input:
            index = i,
            nested_files = BenchmarkIndividualVcfPerContig.fp_query
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
            nested_files = BenchmarkIndividualVcfPerContig.fp_ref
        }

        call LongReadGenotypeTasks.ConcatVcfs as combine_vcfs_fp_ref{
          input:
            vcfs = extract_fp_ref.processed_files,
            outfile_prefix = "~{sample_ids[i]}.~{ref_prefix}.vs.~{ref_prefix}.fp_ref",
            sv_base_mini_docker = sv_base_mini_docker
        }
    }

  output {
    Array[File] tp_query = combine_vcfs_tp_query.concat_vcf  
    Array[File] tp_ref = combine_vcfs_tp_ref.concat_vcf
    Array[File] fp_query = combine_vcfs_fp_query.concat_vcf  
    Array[File] fp_ref = combine_vcfs_fp_ref.concat_vcf
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
