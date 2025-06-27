version 1.0

import BenchmarkIndividualVcfPerContig.wdl as BenchmarkIndividualVcfPerContig

import ExtractFileByIndex.wdl as ExtractFileByIndex

import LongReadGenotypeTasks.wdl as LongReadGenotypeTasks


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
}#script to benchmark SVs in matched genome using truvari

version 1.0

import "https://api.firecloud.org/ga4gh/v1/tools/LR_genotype:ExtractIndividualFromVCF/versions/2/plain-WDL/descriptor" as ExtractIndividualFromVCF
import "https://api.firecloud.org/ga4gh/v1/tools/LR_genotype:ExtractTriosFromVCFByGenomicContext.wdl/versions/5/plain-WDL/descriptor"  as ExtractTriosFromVCFByGenomicContext
import "https://api.firecloud.org/ga4gh/v1/tools/LR_genotype:TruvariBench/versions/1/plain-WDL/descriptor" as TruvariBench

workflow BenchmarkIndividualVcf{
  input{
    File query_vcf
    File ref_vcf

    File anno_script_Rscript
    File anno_script_bash
    File repeat_mask
    File segmental_duplicates
    File simple_repeats

    Array[String] sample_ids
    Array[String] chromosomes
    File ref_dict

    String? truvari_params
    String sv_pipeline_base_docker
  }


  call ExtractIndividualFromVCF.ExtractIndividualFromVCF as extract_individual_query{
    input:
      vcf_file = query_vcf,
      sample_ids = sample_ids,
      sv_pipeline_base_docker = sv_pipeline_base_docker
  }

  call ExtractIndividualFromVCF.ExtractIndividualFromVCF as extract_individual_ref{
    input:
      vcf_file = ref_vcf,
      sample_ids = sample_ids,
      sv_pipeline_base_docker = sv_pipeline_base_docker
  }

  String prefix_query = basename(query_vcf,'.vcf.gz')
  String prefix_ref = basename(ref_vcf,'.vcf.gz')

  scatter (index in range(length(sample_ids))){

    call TruvariBench.CallTruvariBench as truvari_bench_snvs{
      input:
          chromosomes = chromosomes,
          truth_vcf  = extract_individual_query.all_snv_vcfs[index],
          comp_vcf   = extract_individual_ref.all_snv_vcfs[index],
          prefix = "~{sample_ids[index]}.SNVs",
          ref_dict = ref_dict, 
          truvari_params = truvari_params
    }

    call TruvariBench.CallTruvariBench as truvari_bench_indels_sm{
      input:
          chromosomes = chromosomes,
          truth_vcf  = extract_individual_query.all_indel_1_30[index],
          comp_vcf   = extract_individual_ref.all_indel_1_50[index],
          prefix = "~{sample_ids[index]}.indels_sm",
          ref_dict = ref_dict, 
          truvari_params = truvari_params
    }

    call TruvariBench.CallTruvariBench as truvari_bench_indels_lg_1{
      input:
          chromosomes = chromosomes,
          truth_vcf  = extract_individual_query.all_indel_31_50[index],
          comp_vcf   = extract_individual_ref.all_indel_1_50[index],
          prefix = "~{sample_ids[index]}.indels_lg_1",
          ref_dict = ref_dict, 
          truvari_params = truvari_params
    }

    call TruvariBench.CallTruvariBench as truvari_bench_indels_lg_2{
      input:
          chromosomes = chromosomes,
          truth_vcf  = extract_individual_query.all_indel_31_50[index],
          comp_vcf   = extract_individual_ref.all_sv_over30[index],
          prefix = "~{sample_ids[index]}.indels_lg_2",
          ref_dict = ref_dict, 
          truvari_params = truvari_params
    }

    call TruvariBench.CallTruvariBench as truvari_bench_sv{
      input:
          chromosomes = chromosomes,
          truth_vcf  = extract_individual_query.all_sv_over50[index],
          comp_vcf   = extract_individual_ref.all_sv_over30[index],
          prefix = "~{sample_ids[index]}.sv",
          ref_dict = ref_dict, 
          truvari_params = truvari_params
    }

    call MergeVcfs as merge_tp_query{
      input:
        input_vcfs = [truvari_bench_snvs.tp_comp_vcf, truvari_bench_indels_sm.tp_comp_vcf, truvari_bench_indels_lg_1.tp_comp_vcf, truvari_bench_indels_lg_2.tp_comp_vcf, truvari_bench_sv.tp_comp_vcf],
        output_name  = "~{prefix_query}.~{sample_ids[index]}.query_tp",
        docker_image = sv_pipeline_base_docker  
    }


    call MergeVcfs as merge_tp_ref{
      input:
        input_vcfs = [truvari_bench_snvs.tp_base_vcf, truvari_bench_indels_sm.tp_base_vcf, truvari_bench_indels_lg_1.tp_base_vcf, truvari_bench_indels_lg_2.tp_base_vcf, truvari_bench_sv.tp_base_vcf],
        output_name  = "~{prefix_query}.~{sample_ids[index]}.ref_tp",
        docker_image = sv_pipeline_base_docker  
    }

  }

  output {
    Array[File] tp_query = merge_tp_query.merged_vcf  
    Array[File] tp_ref = merge_tp_ref.merged_vcf
  }
}



task MergeVcfs {
  input {
    Array[File] input_vcfs
    String output_name
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  command <<<
    set -e

    # Merge the input VCFs
    bcftools merge \
      ~{sep=' ' input_vcfs} \
      -Oz -o merged.tmp.vcf.gz

    # Index the temp merged VCF
    tabix -p vcf merged.tmp.vcf.gz

    # Remove duplicates
    bcftools norm -d all \
      -Oz -o ~{output_name} merged.tmp.vcf.gz

    # Index the final deduplicated VCF
    tabix -p vcf ~{output_name}
  >>>

  output {
    File merged_vcf = "~{output_name}"
    File merged_vcf_index = "~{output_name}.tbi"
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 15,
    disk_gb: 20,
    boot_disk_gb: 10,
    preemptible_tries: 1,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: docker_image
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
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
