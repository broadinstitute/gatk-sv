version 1.0

import "Structs.wdl"
import "ExtractTriosFromVCF.wdl" as ExtractTriosFromVCF
import "LongReadGenotypeTasks.wdl" as LongReadGenotypeTasks

workflow ExtractTriosFromVCFByGenomicContext {
  input {
    File input_vcf
    File input_vcf_index

    Array[Array[String]] families  # List of trios: [[fa1, mo1, ch1], [fa2, mo2, ch2], ...]
    File inheri_table
    String prefix
    Int variants_per_shard = 200000

    File anno_script_bash
    File anno_script_Rscript
    File repeat_mask
    File simple_repeats
    File segmental_duplicates

    String sv_base_mini_docker 
    String sv_pipeline_base_docker 

    RuntimeAttr? runtime_attr_override
    RuntimeAttr? runtime_attr_annotate_genomic_context
    RuntimeAttr? runtime_attr_calcu_inheri_table_snv
    RuntimeAttr? runtime_attr_calcu_inheri_table_sv
    RuntimeAttr? runtime_attr_calcu_inheri_table_indel_sm
    RuntimeAttr? runtime_attr_calcu_inheri_table_indel_lg
    RuntimeAttr? runtime_attr_evaluate_inheri_by_gq_snv
    RuntimeAttr? runtime_attr_evaluate_inheri_by_gq_indel_sm
    RuntimeAttr? runtime_attr_evaluate_inheri_by_gq_indel_lg
    RuntimeAttr? runtime_attr_evaluate_inheri_by_gq_sv

  }

  call LongReadGenotypeTasks.SplitVcfIntoShards as SplitVcfIntoShards {
    input:
        input_vcf = input_vcf,
        input_vcf_index = input_vcf_index,
        variants_per_shard = 200000,
        output_prefix = prefix,
        docker_image = sv_pipeline_base_docker
    }
  
  scatter (index in range(length(SplitVcfIntoShards.split_vcfs))) {

    call LongReadGenotypeTasks.ExtractVariantSites as ExtractVariantSites{
        input:
          input_vcf = SplitVcfIntoShards.split_vcfs[index],
          docker_image = sv_pipeline_base_docker
        }

    call LongReadGenotypeTasks.AnnotateGenomicContext{
        input:
          variant_sites = ExtractVariantSites.variant_sites,
          anno_script_bash = anno_script_bash,
          anno_script_Rscript = anno_script_Rscript,
          repeat_mask = repeat_mask,
          segmental_duplicates = segmental_duplicates,
          simple_repeats = simple_repeats,
          docker_image = sv_pipeline_base_docker,
          runtime_attr_override = runtime_attr_annotate_genomic_context
        }

    call LongReadGenotypeTasks.SplitVcfByAnnotationR as split_vcf_by_annotation_sr {
        input:
          vcf_file = ExtractVariantSites.updated_vcf,
          vcf_idx = ExtractVariantSites.updated_vcf_idx,
          svid_annotation = AnnotateGenomicContext.variant_anno, 
          anno_list = "SR",
          midfix = "SR",
          docker_image = sv_pipeline_base_docker
      }

    call LongReadGenotypeTasks.SplitVcfByAnnotationR as split_vcf_by_annotation_sd {
        input:
          vcf_file = ExtractVariantSites.updated_vcf,
          vcf_idx = ExtractVariantSites.updated_vcf_idx,
          svid_annotation = AnnotateGenomicContext.variant_anno, 
          anno_list = "SD",
          midfix = "SD",
          docker_image = sv_pipeline_base_docker
      }

    call LongReadGenotypeTasks.SplitVcfByAnnotationR as split_vcf_by_annotation_rm {
        input:
          vcf_file = ExtractVariantSites.updated_vcf,
          vcf_idx = ExtractVariantSites.updated_vcf_idx,
          svid_annotation = AnnotateGenomicContext.variant_anno, 
          anno_list = "RM",
          midfix = "RM",
          docker_image = sv_pipeline_base_docker
      }

    call LongReadGenotypeTasks.SplitVcfByAnnotationR as split_vcf_by_annotation_us {
        input:
          vcf_file = ExtractVariantSites.updated_vcf,
          vcf_idx = ExtractVariantSites.updated_vcf_idx,
          svid_annotation = AnnotateGenomicContext.variant_anno, 
          anno_list = "US",
          midfix = "US",
          docker_image = sv_pipeline_base_docker
      }

    call LongReadGenotypeTasks.SplitVcfByAnnotationR as split_vcf_by_annotation_clean {
        input:
          vcf_file = ExtractVariantSites.updated_vcf,
          vcf_idx = ExtractVariantSites.updated_vcf_idx,
          svid_annotation = AnnotateGenomicContext.variant_anno, 
          anno_list = "US,RM",
          midfix = "US_RM",
          docker_image = sv_pipeline_base_docker
      }
  }

  call LongReadGenotypeTasks.ConcatVcfs as concat_vcf_by_annotation_sr{
      input:
        vcfs = split_vcf_by_annotation_sr.split_vcf,
        vcfs_idx = split_vcf_by_annotation_sr.split_vcf_idx,
        outfile_prefix = "SR",
        sv_base_mini_docker = sv_base_mini_docker
  }

  call LongReadGenotypeTasks.ConcatVcfs as concat_vcf_by_annotation_sd{
      input:
        vcfs = split_vcf_by_annotation_sd.split_vcf,
        vcfs_idx = split_vcf_by_annotation_sd.split_vcf_idx,
        outfile_prefix = "SD",
        sv_base_mini_docker = sv_base_mini_docker
  }

  call LongReadGenotypeTasks.ConcatVcfs as concat_vcf_by_annotation_rm{
      input:
        vcfs = split_vcf_by_annotation_rm.split_vcf,
        vcfs_idx = split_vcf_by_annotation_rm.split_vcf_idx,
        outfile_prefix = "RM",
        sv_base_mini_docker = sv_base_mini_docker
  }

  call LongReadGenotypeTasks.ConcatVcfs as concat_vcf_by_annotation_us{
      input:
        vcfs = split_vcf_by_annotation_us.split_vcf,
        vcfs_idx = split_vcf_by_annotation_us.split_vcf_idx,
        outfile_prefix = "US",
        sv_base_mini_docker = sv_base_mini_docker
  }

  call LongReadGenotypeTasks.ConcatVcfs as concat_vcf_by_annotation_clean{
      input:
        vcfs = split_vcf_by_annotation_clean.split_vcf,
        vcfs_idx = split_vcf_by_annotation_clean.split_vcf_idx,
        outfile_prefix = "US_RM",
        sv_base_mini_docker = sv_base_mini_docker
  }

  call ExtractTriosFromVCF.ExtractTriosFromVCF as extract_trios_from_SRs{
    input:
      input_vcf = concat_vcf_by_annotation_sr.concat_vcf,
      families = families,
      inheri_table = inheri_table,
      prefix = "~{prefix}.SR",
      sv_pipeline_base_docker = sv_pipeline_base_docker,
      runtime_attr_override = runtime_attr_override,
      runtime_attr_calcu_inheri_table_snv = runtime_attr_calcu_inheri_table_snv,
      runtime_attr_calcu_inheri_table_sv = runtime_attr_calcu_inheri_table_sv,
      runtime_attr_calcu_inheri_table_indel_sm = runtime_attr_calcu_inheri_table_indel_sm,
      runtime_attr_calcu_inheri_table_indel_lg = runtime_attr_calcu_inheri_table_indel_lg,
      runtime_attr_evaluate_inheri_by_gq_snv = runtime_attr_evaluate_inheri_by_gq_snv,
      runtime_attr_evaluate_inheri_by_gq_sv = runtime_attr_evaluate_inheri_by_gq_sv,
      runtime_attr_evaluate_inheri_by_gq_indel_sm = runtime_attr_evaluate_inheri_by_gq_indel_sm,
      runtime_attr_evaluate_inheri_by_gq_indel_lg = runtime_attr_evaluate_inheri_by_gq_indel_lg
    }

  call ExtractTriosFromVCF.ExtractTriosFromVCF as extract_trios_from_SDs{
    input:
      input_vcf = concat_vcf_by_annotation_sd.concat_vcf,
      families = families,
      inheri_table = inheri_table,
      prefix = "~{prefix}.SD",
      sv_pipeline_base_docker = sv_pipeline_base_docker,
      runtime_attr_override = runtime_attr_override,
      runtime_attr_calcu_inheri_table_snv = runtime_attr_calcu_inheri_table_snv,
      runtime_attr_calcu_inheri_table_sv = runtime_attr_calcu_inheri_table_sv,
      runtime_attr_calcu_inheri_table_indel_sm = runtime_attr_calcu_inheri_table_indel_sm,
      runtime_attr_calcu_inheri_table_indel_lg = runtime_attr_calcu_inheri_table_indel_lg,
      runtime_attr_evaluate_inheri_by_gq_snv = runtime_attr_evaluate_inheri_by_gq_snv,
      runtime_attr_evaluate_inheri_by_gq_sv = runtime_attr_evaluate_inheri_by_gq_sv,
      runtime_attr_evaluate_inheri_by_gq_indel_sm = runtime_attr_evaluate_inheri_by_gq_indel_sm,
      runtime_attr_evaluate_inheri_by_gq_indel_lg = runtime_attr_evaluate_inheri_by_gq_indel_lg
    }

  call ExtractTriosFromVCF.ExtractTriosFromVCF as extract_trios_from_RMs{
    input:
      input_vcf = concat_vcf_by_annotation_rm.concat_vcf,
      families = families,
      inheri_table = inheri_table,
      prefix = "~{prefix}.RM",
      sv_pipeline_base_docker = sv_pipeline_base_docker,
      runtime_attr_override = runtime_attr_override,
      runtime_attr_calcu_inheri_table_snv = runtime_attr_calcu_inheri_table_snv,
      runtime_attr_calcu_inheri_table_sv = runtime_attr_calcu_inheri_table_sv,
      runtime_attr_calcu_inheri_table_indel_sm = runtime_attr_calcu_inheri_table_indel_sm,
      runtime_attr_calcu_inheri_table_indel_lg = runtime_attr_calcu_inheri_table_indel_lg,
      runtime_attr_evaluate_inheri_by_gq_snv = runtime_attr_evaluate_inheri_by_gq_snv,
      runtime_attr_evaluate_inheri_by_gq_sv = runtime_attr_evaluate_inheri_by_gq_sv,
      runtime_attr_evaluate_inheri_by_gq_indel_sm = runtime_attr_evaluate_inheri_by_gq_indel_sm,
      runtime_attr_evaluate_inheri_by_gq_indel_lg = runtime_attr_evaluate_inheri_by_gq_indel_lg
    }

  call ExtractTriosFromVCF.ExtractTriosFromVCF as extract_trios_from_USs{
    input:
      input_vcf = concat_vcf_by_annotation_us.concat_vcf,
      families = families,
      inheri_table = inheri_table,
      prefix = "~{prefix}.US",
      sv_pipeline_base_docker = sv_pipeline_base_docker,
      runtime_attr_override = runtime_attr_override,
      runtime_attr_calcu_inheri_table_snv = runtime_attr_calcu_inheri_table_snv,
      runtime_attr_calcu_inheri_table_sv = runtime_attr_calcu_inheri_table_sv,
      runtime_attr_calcu_inheri_table_indel_sm = runtime_attr_calcu_inheri_table_indel_sm,
      runtime_attr_calcu_inheri_table_indel_lg = runtime_attr_calcu_inheri_table_indel_lg,
      runtime_attr_evaluate_inheri_by_gq_snv = runtime_attr_evaluate_inheri_by_gq_snv,
      runtime_attr_evaluate_inheri_by_gq_sv = runtime_attr_evaluate_inheri_by_gq_sv,
      runtime_attr_evaluate_inheri_by_gq_indel_sm = runtime_attr_evaluate_inheri_by_gq_indel_sm,
      runtime_attr_evaluate_inheri_by_gq_indel_lg = runtime_attr_evaluate_inheri_by_gq_indel_lg
    }

  call ExtractTriosFromVCF.ExtractTriosFromVCF as extract_trios_from_USRMs{
    input:
      input_vcf = concat_vcf_by_annotation_clean.concat_vcf,
      families = families,
      inheri_table = inheri_table,
      prefix = "~{prefix}.US_RM",
      sv_pipeline_base_docker = sv_pipeline_base_docker,
      runtime_attr_override = runtime_attr_override,
      runtime_attr_calcu_inheri_table_snv = runtime_attr_calcu_inheri_table_snv,
      runtime_attr_calcu_inheri_table_sv = runtime_attr_calcu_inheri_table_sv,
      runtime_attr_calcu_inheri_table_indel_sm = runtime_attr_calcu_inheri_table_indel_sm,
      runtime_attr_calcu_inheri_table_indel_lg = runtime_attr_calcu_inheri_table_indel_lg,
      runtime_attr_evaluate_inheri_by_gq_snv = runtime_attr_evaluate_inheri_by_gq_snv,
      runtime_attr_evaluate_inheri_by_gq_sv = runtime_attr_evaluate_inheri_by_gq_sv,
      runtime_attr_evaluate_inheri_by_gq_indel_sm = runtime_attr_evaluate_inheri_by_gq_indel_sm,
      runtime_attr_evaluate_inheri_by_gq_indel_lg = runtime_attr_evaluate_inheri_by_gq_indel_lg
    }

  call ExtractTriosFromVCF.ExtractTriosFromVCF as extract_trios_from_all{
    input:
      input_vcf = input_vcf,
      families = families,
      inheri_table = inheri_table,
      prefix = "~{prefix}.all",
      sv_pipeline_base_docker = sv_pipeline_base_docker,
      runtime_attr_override = runtime_attr_override,
      runtime_attr_calcu_inheri_table_snv = runtime_attr_calcu_inheri_table_snv,
      runtime_attr_calcu_inheri_table_sv = runtime_attr_calcu_inheri_table_sv,
      runtime_attr_calcu_inheri_table_indel_sm = runtime_attr_calcu_inheri_table_indel_sm,
      runtime_attr_calcu_inheri_table_indel_lg = runtime_attr_calcu_inheri_table_indel_lg,
      runtime_attr_evaluate_inheri_by_gq_snv = runtime_attr_evaluate_inheri_by_gq_snv,
      runtime_attr_evaluate_inheri_by_gq_sv = runtime_attr_evaluate_inheri_by_gq_sv,
      runtime_attr_evaluate_inheri_by_gq_indel_sm = runtime_attr_evaluate_inheri_by_gq_indel_sm,
      runtime_attr_evaluate_inheri_by_gq_indel_lg = runtime_attr_evaluate_inheri_by_gq_indel_lg
    }
  
  output{
    Array[File] inheritance_table_inte_SR = extract_trios_from_SRs.inheritance_table_inte
    Array[File] inheritance_table_inte_SD = extract_trios_from_SDs.inheritance_table_inte
    Array[File] inheritance_table_inte_RM = extract_trios_from_RMs.inheritance_table_inte
    Array[File] inheritance_table_inte_US = extract_trios_from_USs.inheritance_table_inte
    Array[File] inheritance_table_inte_US_RM = extract_trios_from_USRMs.inheritance_table_inte
    Array[File] inheritance_table_inte_all = extract_trios_from_all.inheritance_table_inte

    Array[File] inheritance_by_gq_table_SR = extract_trios_from_SRs.inheritance_by_gq_table_inte
    Array[File] inheritance_by_gq_table_SD = extract_trios_from_SDs.inheritance_by_gq_table_inte
    Array[File] inheritance_by_gq_table_RM = extract_trios_from_RMs.inheritance_by_gq_table_inte
    Array[File] inheritance_by_gq_table_US = extract_trios_from_USs.inheritance_by_gq_table_inte
    Array[File] inheritance_by_gq_table_US_RM = extract_trios_from_USRMs.inheritance_by_gq_table_inte
    Array[File] inheritance_by_gq_table_all = extract_trios_from_all.inheritance_by_gq_table_inte


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


