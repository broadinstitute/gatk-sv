version 1.0

import "Structs.wdl"
import "TruvariBench.wdl" as TruvariBench
import "LongReadGenotypeTasks.wdl" as LongReadGenotypeTasks
import "ExtractIndividualFromVCF.wdl" as ExtractIndividualFromVCF
import "ExtractTriosFromVCFByGenomicContext.wdl" as ExtractTriosFromVCFByGenomicContext


workflow BenchmarkVcfSitesPerContig{
  input{
    File query_vcf
    File query_vcf_idx
    File ref_vcf
    File ref_vcf_idx

    String chromosome
    File ref_dict

    File anno_script_bash
    File anno_script_helper_R
    File benchmark_script_bash
    File banchmark_script_helper_R

    File repeat_mask
    File simple_repeats
    File segmental_duplicates

    Boolean short_read_benchmark = false
    Boolean simplify_comp_vcf_lt_20bp = false
    Boolean simplify_base_vcf_lt_20bp = false


    String? truvari_params
    String sv_base_mini_docker
    String sv_pipeline_base_docker

    RuntimeAttr? runtime_attr_benchmark_SNVs
    RuntimeAttr? runtime_attr_split_ref
    RuntimeAttr? runtime_attr_split_query
    RuntimeAttr? runtime_attr_add_dummy_gt_ref
    RuntimeAttr? runtime_attr_add_dummy_gt_query
    RuntimeAttr? runtime_attr_extract_variant_sites_ref
    RuntimeAttr? runtime_attr_extract_variant_sites_query
  }

  String prefix_query = basename(query_vcf,'.vcf.gz')
  String prefix_ref = basename(ref_vcf,'.vcf.gz')

  call LongReadGenotypeTasks.SplitVcfToSites as split_query_vcf_into_sites{
    input:
      vcf_file = query_vcf,
      vcf_idx = query_vcf_idx,
      docker_image = sv_base_mini_docker
  }

  call LongReadGenotypeTasks.SplitVcfToSites as split_ref_vcf_into_sites{
    input:
      vcf_file = ref_vcf,
      vcf_idx = ref_vcf_idx,
      docker_image = sv_base_mini_docker
  }

  call LongReadGenotypeTasks.AddDummyGT as add_dummy_gt_query{
    input:
      sites_file = split_query_vcf_into_sites.vcf_sites,
      sites_idx  = split_query_vcf_into_sites.vcf_sites_idx,
      docker_image   = sv_pipeline_base_docker,
      runtime_attr_override = runtime_attr_add_dummy_gt_query
  }

  call LongReadGenotypeTasks.AddDummyGT as add_dummy_gt_ref{
    input:
      sites_file = split_ref_vcf_into_sites.vcf_sites,
      sites_idx  = split_ref_vcf_into_sites.vcf_sites_idx,
      docker_image   = sv_pipeline_base_docker,
      runtime_attr_override = runtime_attr_add_dummy_gt_ref
  }

  call LongReadGenotypeTasks.ExtractVariantSites as extract_variant_sites_query{
    input:
      input_vcf = add_dummy_gt_query.vcf_file,
      docker_image = sv_pipeline_base_docker,
      runtime_attr_override = runtime_attr_extract_variant_sites_query
  }

  call LongReadGenotypeTasks.ExtractVariantSites as extract_variant_sites_ref{
    input:
      input_vcf = add_dummy_gt_ref.vcf_file,
      docker_image = sv_pipeline_base_docker,
      runtime_attr_override = runtime_attr_extract_variant_sites_ref
  }

  call LongReadGenotypeTasks.SplitVariantsBySizeAt20bp as split_query{
    input:
      input_vcf      = extract_variant_sites_query.updated_vcf,
      input_vcf_idx  = extract_variant_sites_query.updated_vcf_idx, 
      docker_image   = sv_pipeline_base_docker,
      runtime_attr_override = runtime_attr_split_query
  }

  call LongReadGenotypeTasks.SplitVariantsBySizeAt20bp as split_ref{
    input:
      input_vcf      = extract_variant_sites_ref.updated_vcf,
      input_vcf_idx  = extract_variant_sites_ref.updated_vcf_idx, 
      docker_image   = sv_pipeline_base_docker,
      runtime_attr_override = runtime_attr_split_ref
    }
    

  call LongReadGenotypeTasks.BenchmarkSNVs as truvari_bench_lt_20bp{
    input:
        comp_vcf = split_query.indels_lt_20_vcf,
        base_vcf = split_ref.indels_lt_20_vcf,
        prefix = "~{chromosome}.lt_20bp",
        simplify_comp_vcf = simplify_comp_vcf_lt_20bp,
        simplify_base_vcf = simplify_base_vcf_lt_20bp,
        docker_image = sv_pipeline_base_docker,
        runtime_attr_override = runtime_attr_benchmark_SNVs
  }


  if (short_read_benchmark){
    call LongReadGenotypeTasks.BenchmarkSVsV2 as truvari_bench_gt_20bp_sr{
      input:
        comp_vcf  = split_query.svs_gt_20_vcf,
        base_vcf  = split_ref.svs_gt_20_vcf,
        prefix = "~{chromosome}.gt_20bp",
        benchmark_bash = benchmark_script_bash,
        banchmark_helper_R = banchmark_script_helper_R,
        docker_image = sv_pipeline_base_docker
      }
  }

  if (!short_read_benchmark){
    call TruvariBench.CallTruvariBench as truvari_bench_gt_20bp_lr{
      input:
          chromosomes = [chromosome],
          comp_vcf  = split_query.svs_gt_20_vcf,
          truth_vcf = split_ref.svs_gt_20_vcf,
          prefix = "~{chromosome}.gt_20bp",
          ref_dict = ref_dict, 
          truvari_params = "-s 10 -S 10"
    }
  }

  File SV_tp_comp_vcf =  select_first([truvari_bench_gt_20bp_sr.tp_comp_vcf, truvari_bench_gt_20bp_lr.tp_comp_vcf])
  File SV_tp_base_vcf =  select_first([truvari_bench_gt_20bp_sr.tp_base_vcf, truvari_bench_gt_20bp_lr.tp_base_vcf])
  File SV_fp_vcf      =  select_first([truvari_bench_gt_20bp_sr.fp_vcf,      truvari_bench_gt_20bp_lr.fp_vcf])
  File SV_fn_vcf      =  select_first([truvari_bench_gt_20bp_sr.fn_vcf,      truvari_bench_gt_20bp_lr.fn_vcf])

  call LongReadGenotypeTasks.ConcatVcfs as merge_tp_query{
    input:
      vcfs = [truvari_bench_lt_20bp.tp_comp_vcf, SV_tp_comp_vcf],
      outfile_prefix  = "~{prefix_query}.query_tp",
      sv_base_mini_docker = sv_base_mini_docker
  }

  call LongReadGenotypeTasks.ConcatVcfs as merge_tp_ref{
    input:
      vcfs = [truvari_bench_lt_20bp.tp_base_vcf,  SV_tp_base_vcf],
      outfile_prefix  = "~{prefix_query}.ref_tp",
      sv_base_mini_docker = sv_base_mini_docker
  }

  call LongReadGenotypeTasks.ConcatVcfs as merge_fp_query{
    input:
      vcfs = [truvari_bench_lt_20bp.fp_vcf, SV_fp_vcf],
      outfile_prefix  = "~{prefix_query}.query_fp",
      sv_base_mini_docker = sv_base_mini_docker
  }

  call LongReadGenotypeTasks.ConcatVcfs as merge_fp_ref{
    input:
      vcfs = [truvari_bench_lt_20bp.fn_vcf,  SV_fn_vcf],
      outfile_prefix  = "~{prefix_query}.ref_fp",
      sv_base_mini_docker = sv_base_mini_docker
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

  output {
    File tp_query = add_genomic_context_query_tp.annotated_vcf  
    File tp_ref = add_genomic_context_ref_tp.annotated_vcf
    File fp_query = add_genomic_context_query_fp.annotated_vcf  
    File fp_ref = add_genomic_context_ref_fp.annotated_vcf
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

