version 1.0

import "Structs.wdl"
import "TruvariBench.wdl" as TruvariBench
import "LongReadGenotypeTasks.wdl" as LongReadGenotypeTasks
import "ExtractIndividualFromVCF.wdl" as ExtractIndividualFromVCF
import "ExtractTriosFromVCFByGenomicContext.wdl" as ExtractTriosFromVCFByGenomicContext


workflow BenchmarkIndividualVcfPerContig{
  input{
    File ref_vcf
    File ref_vcf_idx
    File query_vcf
    File query_vcf_idx  
    File? ref_filter_vcf
    File? ref_filter_vcf_idx
    File? query_filter_vcf
    File? query_filter_vcf_idx

    File anno_script_bash
    File anno_script_helper_R
    File benchmark_script_bash
    File banchmark_script_helper_R

    File repeat_mask
    File simple_repeats
    File segmental_duplicates

    Array[String] sample_ids
    String chromosome
    File ref_dict

    Boolean short_read_benchmark = false

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


  if (defined(query_filter_vcf)){
    call LongReadGenotypeTasks.FilterVcfByAnotherVcf as filter_query_vcf_by_ovr {
      input:
        vcf_file = extract_variant_sites_query.updated_vcf,
        vcf_idx = extract_variant_sites_query.updated_vcf_idx,
        vcf_file_b = query_filter_vcf,
        vcf_file_b_idx = query_filter_vcf_idx,
        docker_image = sv_pipeline_base_docker
    }
  }

  File updated_query_vcf = select_first([filter_query_vcf_by_ovr.filtered_vcf, extract_variant_sites_query.updated_vcf])
  File updated_query_vcf_idx = select_first([filter_query_vcf_by_ovr.filtered_vcf_idx, extract_variant_sites_query.updated_vcf_idx])


  if (defined(ref_filter_vcf)){
    call LongReadGenotypeTasks.FilterVcfByAnotherVcf as filter_ref_vcf_by_ovr {
    input:
      vcf_file = extract_variant_sites_ref.updated_vcf,
      vcf_idx = extract_variant_sites_ref.updated_vcf_idx,
      vcf_file_b = ref_filter_vcf,
      vcf_file_b_idx = ref_filter_vcf_idx,
      docker_image = sv_pipeline_base_docker
    }
  }

  File updated_ref_vcf = select_first([filter_ref_vcf_by_ovr.filtered_vcf, extract_variant_sites_ref.updated_vcf])
  File updated_ref_vcf_idx = select_first([filter_ref_vcf_by_ovr.filtered_vcf_idx, extract_variant_sites_ref.updated_vcf_idx])

  call ExtractIndividualFromVCF.ExtractIndividualFromVCF as extract_individual_query{
    input:
      vcf_file = updated_query_vcf,
      sample_ids = sample_ids,
      sv_base_mini_docker = sv_base_mini_docker,
      sv_pipeline_base_docker = sv_pipeline_base_docker
  }

  call ExtractIndividualFromVCF.ExtractIndividualFromVCF as extract_individual_ref{
    input:
      vcf_file = updated_ref_vcf,
      sample_ids = sample_ids,
      sv_base_mini_docker = sv_base_mini_docker,
      sv_pipeline_base_docker = sv_pipeline_base_docker
  }

  call LongReadGenotypeTasks.AnnotateGenomicContext as annotate_genomic_context_query{
    input:
      variant_sites = extract_variant_sites_query.variant_sites,
      anno_script_bash = anno_script_bash,
      anno_script_Rscript = anno_script_helper_R,
      repeat_mask = repeat_mask,
      simple_repeats = simple_repeats,
      segmental_duplicates = segmental_duplicates,
      docker_image = sv_pipeline_base_docker
  }

  call LongReadGenotypeTasks.AnnotateGenomicContext as annotate_genomic_context_ref{
    input:
      variant_sites = extract_variant_sites_ref.variant_sites,
      anno_script_bash = anno_script_bash,
      anno_script_Rscript = anno_script_helper_R,
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
          prefix = "~{sample_ids[index]}.SNVs",
          docker_image = sv_pipeline_base_docker
    }

    call TruvariBench.CallTruvariBench as truvari_bench_indels{
      input:
          chromosomes = [chromosome],
          comp_vcf  = extract_individual_query.all_indel_1_50[index],
          truth_vcf = extract_individual_ref.all_indel_1_50[index],
          prefix = "~{sample_ids[index]}.indels",
          ref_dict = ref_dict, 
          truvari_params = "-s 1 -S 1"
    }

    if (short_read_benchmark){
      call LongReadGenotypeTasks.BenchmarkSVs as truvari_bench_sv_sr{
        input:
          comp_vcf  = extract_individual_query.all_sv_over30[index],
          base_vcf  = extract_individual_ref.all_sv_over30[index],
          prefix = "~{sample_ids[index]}.SVs",
          benchmark_bash = benchmark_script_bash,
          banchmark_helper_R = banchmark_script_helper_R,
          docker_image = sv_pipeline_base_docker
          }
      }

    if (!short_read_benchmark){
      call TruvariBench.CallTruvariBench as truvari_bench_sv_lr{
        input:
            chromosomes = [chromosome],
            comp_vcf  = extract_individual_query.all_sv_over30[index],
            truth_vcf = extract_individual_ref.all_sv_over30[index],
            prefix = "~{sample_ids[index]}.sv",
            ref_dict = ref_dict, 
            truvari_params = "-s 1 -S 1"
      }    
    }

    File SV_tp_comp_vcf =  select_first([truvari_bench_sv_sr.tp_comp_vcf, truvari_bench_sv_lr.tp_comp_vcf])
    File SV_tp_base_vcf =  select_first([truvari_bench_sv_sr.tp_base_vcf, truvari_bench_sv_lr.tp_base_vcf])
    File SV_fp_vcf =  select_first([truvari_bench_sv_sr.fp_vcf, truvari_bench_sv_lr.fp_vcf])
    File SV_fn_vcf =  select_first([truvari_bench_sv_sr.fn_vcf, truvari_bench_sv_lr.fn_vcf])

    call LongReadGenotypeTasks.ConcatVcfs as merge_tp_query{
      input:
        vcfs = [truvari_bench_snvs.tp_comp_vcf, truvari_bench_indels.tp_comp_vcf, SV_tp_comp_vcf],
        outfile_prefix  = "~{prefix_query}.~{sample_ids[index]}.query_tp",
        sv_base_mini_docker = sv_base_mini_docker
    }

    call LongReadGenotypeTasks.ConcatVcfs as merge_tp_ref{
      input:
        vcfs = [truvari_bench_snvs.tp_base_vcf,  truvari_bench_indels.tp_base_vcf,  SV_tp_base_vcf],
        outfile_prefix  = "~{prefix_query}.~{sample_ids[index]}.ref_tp",
        sv_base_mini_docker = sv_base_mini_docker
    }

    call LongReadGenotypeTasks.ConcatVcfs as merge_fp_query{
      input:
        vcfs = [truvari_bench_snvs.fp_vcf, truvari_bench_indels.fp_vcf, SV_fp_vcf],
        outfile_prefix  = "~{prefix_query}.~{sample_ids[index]}.query_fp",
        sv_base_mini_docker = sv_base_mini_docker
    }

    call LongReadGenotypeTasks.ConcatVcfs as merge_fp_ref{
      input:
        vcfs = [truvari_bench_snvs.fn_vcf,  truvari_bench_indels.fn_vcf, SV_fn_vcf],
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

