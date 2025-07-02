version 1.0

import "Structs.wdl"
import "TruvariBench.wdl" as TruvariBench
import "LongReadGenotypeTasks.wdl" as LongReadGenotypeTasks
import "ExtractIndividualFromVCF.wdl" as ExtractIndividualFromVCF
import "ExtractTriosFromVCFByGenomicContext.wdl" as ExtractTriosFromVCFByGenomicContext


workflow BenchmarkIndividualVcfPerContig{
  input{
    File query_vcf
    File ref_vcf

    File anno_script_bash
    File anno_script_Rscript
    File repeat_mask
    File simple_repeats
    File segmental_duplicates

    Array[String] sample_ids
    String chromosome
    File ref_dict

    String? truvari_params
    String sv_base_mini_docker
    String sv_pipeline_base_docker
  }

  call LongReadGenotypeTasks.ExtractVariantSites as extract_variant_sites_query{
    input:
      input_vcf = query_vcf,
      docker_image = sv_pipeline_base_docker
  }

  call LongReadGenotypeTasks.ExtractVariantSites as extract_variant_sites_ref{
    input:
      input_vcf = ref_vcf,
      docker_image = sv_pipeline_base_docker
  }

  call ExtractIndividualFromVCF.ExtractIndividualFromVCF as extract_individual_query{
    input:
      vcf_file = extract_variant_sites_query.updated_vcf,
      sample_ids = sample_ids,
      sv_base_mini_docker = sv_base_mini_docker,
      sv_pipeline_base_docker = sv_pipeline_base_docker
  }

  call ExtractIndividualFromVCF.ExtractIndividualFromVCF as extract_individual_ref{
    input:
      vcf_file = extract_variant_sites_ref.updated_vcf,
      sample_ids = sample_ids,
      sv_base_mini_docker = sv_base_mini_docker,
      sv_pipeline_base_docker = sv_pipeline_base_docker
  }

  call LongReadGenotypeTasks.AnnotateGenomicContext as annotate_genomic_context_query{
    input:
      variant_sites = extract_variant_sites_query.variant_sites,
      anno_script_bash = anno_script_bash,
      anno_script_Rscript = anno_script_Rscript,
      repeat_mask = repeat_mask,
      simple_repeats = simple_repeats,
      segmental_duplicates = segmental_duplicates,
      docker_image = sv_pipeline_base_docker
  }

  call LongReadGenotypeTasks.AnnotateGenomicContext as annotate_genomic_context_ref{
    input:
      variant_sites = extract_variant_sites_ref.variant_sites,
      anno_script_bash = anno_script_bash,
      anno_script_Rscript = anno_script_Rscript,
      repeat_mask = repeat_mask,
      simple_repeats = simple_repeats,
      segmental_duplicates = segmental_duplicates,
      docker_image = sv_pipeline_base_docker
  }

  String prefix_query = basename(query_vcf,'.vcf.gz')
  String prefix_ref = basename(ref_vcf,'.vcf.gz')

  scatter (index in range(length(sample_ids))){

    call LongReadGenotypeTasks.BenchmarkSNVs as truvari_bench_snvs{
      input:
          comp_vcf = extract_individual_query.all_snv_vcfs[index],
          base_vcf = extract_individual_ref.all_snv_vcfs[index],
          docker_image = sv_pipeline_base_docker
    }

    call TruvariBench.CallTruvariBench as truvari_bench_indels_sm{
      input:
          chromosomes = [chromosome],
          comp_vcf  = extract_individual_query.all_indel_1_30[index],
          truth_vcf   = extract_individual_ref.all_indel_1_50[index],
          prefix = "~{sample_ids[index]}.indels_sm",
          ref_dict = ref_dict, 
          truvari_params = "-s 1 -S 1"
    }

    call TruvariBench.CallTruvariBench as truvari_bench_indels_lg_1{
      input:
          chromosomes = [chromosome],
          comp_vcf  = extract_individual_query.all_indel_31_50[index],
          truth_vcf   = extract_individual_ref.all_indel_1_50[index],
          prefix = "~{sample_ids[index]}.indels_lg_1",
          ref_dict = ref_dict, 
          truvari_params = "-s 1 -S 1"
    }

    call TruvariBench.CallTruvariBench as truvari_bench_indels_lg_2{
      input:
          chromosomes = [chromosome],
          comp_vcf  = extract_individual_query.all_indel_31_50[index],
          truth_vcf   = extract_individual_ref.all_sv_over30[index],
          prefix = "~{sample_ids[index]}.indels_lg_2",
          ref_dict = ref_dict, 
          truvari_params = "-s 1 -S 1"
    }

    call TruvariBench.CallTruvariBench as truvari_bench_sv{
      input:
          chromosomes = [chromosome],
          comp_vcf  = extract_individual_query.all_sv_over50[index],
          truth_vcf   = extract_individual_ref.all_sv_over30[index],
          prefix = "~{sample_ids[index]}.sv",
          ref_dict = ref_dict, 
          truvari_params = "-s 1 -S 1"
    }


    call LongReadGenotypeTasks.ConcatVcfs as merge_tp_query{
      input:
        vcfs = [truvari_bench_snvs.tp_comp_vcf, 
                truvari_bench_indels_sm.tp_comp_vcf, 
                truvari_bench_indels_lg_1.tp_comp_vcf, 
                truvari_bench_indels_lg_2.tp_comp_vcf, 
                truvari_bench_sv.tp_comp_vcf],
        outfile_prefix  = "~{prefix_query}.~{sample_ids[index]}.query_tp",
        sv_base_mini_docker = sv_base_mini_docker
    }

    call LongReadGenotypeTasks.ConcatVcfs as merge_tp_ref{
      input:
        vcfs = [truvari_bench_snvs.tp_base_vcf, 
                truvari_bench_indels_sm.tp_base_vcf, 
                truvari_bench_indels_lg_1.tp_base_vcf, 
                truvari_bench_indels_lg_2.tp_base_vcf, 
                truvari_bench_sv.tp_base_vcf],
        outfile_prefix  = "~{prefix_query}.~{sample_ids[index]}.ref_tp",
        sv_base_mini_docker = sv_base_mini_docker
    }

    call LongReadGenotypeTasks.ConcatVcfs as merge_fp_query{
      input:
        vcfs = [truvari_bench_snvs.fp_vcf, 
                      truvari_bench_indels_sm.fp_vcf, 
                      truvari_bench_indels_lg_1.fp_vcf, 
                      truvari_bench_indels_lg_2.fp_vcf, 
                      truvari_bench_sv.fp_vcf],
        outfile_prefix  = "~{prefix_query}.~{sample_ids[index]}.query_fp",
        sv_base_mini_docker = sv_base_mini_docker
    }

    call LongReadGenotypeTasks.ConcatVcfs as merge_fp_ref{
      input:
        vcfs = [truvari_bench_snvs.fn_vcf, 
                      truvari_bench_indels_sm.fn_vcf, 
                      truvari_bench_indels_lg_1.fn_vcf, 
                      truvari_bench_indels_lg_2.fn_vcf, 
                      truvari_bench_sv.fn_vcf],
        outfile_prefix  = "~{prefix_query}.~{sample_ids[index]}.ref_fp",
        sv_base_mini_docker = sv_base_mini_docker
    }

    call LongReadGenotypeTasks.AddGenomicContextToVcfR as add_genomic_context_query_tp{
      input:
          vcf_file = merge_tp_query.concat_vcf,
          svid_annotation = annotate_genomic_context_query.variant_anno,
          docker_image = sv_pipeline_base_docker
    }

    call LongReadGenotypeTasks.AddGenomicContextToVcfR as add_genomic_context_ref_tp{
      input:
          vcf_file = merge_tp_ref.concat_vcf,
          svid_annotation = annotate_genomic_context_ref.variant_anno,
          docker_image = sv_pipeline_base_docker
    }

    call LongReadGenotypeTasks.AddGenomicContextToVcfR as add_genomic_context_query_fp{
      input:
          vcf_file = merge_fp_query.concat_vcf,
          svid_annotation = annotate_genomic_context_query.variant_anno,
          docker_image = sv_pipeline_base_docker
    }

    call LongReadGenotypeTasks.AddGenomicContextToVcfR as add_genomic_context_ref_fp{
      input:
          vcf_file = merge_fp_ref.concat_vcf,
          svid_annotation = annotate_genomic_context_ref.variant_anno,
          docker_image = sv_pipeline_base_docker
    }
  }

  output {
    File updated_query_vcf = extract_variant_sites_query.updated_vcf,
    File updated_ref_vcf = extract_variant_sites_ref.updated_vcf,
    Array[File] tp_query = add_genomic_context_query_tp.annotated_vcf  
    Array[File] tp_ref = add_genomic_context_ref_tp.annotated_vcf
    Array[File] fp_query = add_genomic_context_query_fp.annotated_vcf  
    Array[File] fp_ref = add_genomic_context_ref_fp.annotated_vcf
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

