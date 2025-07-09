version 1.0

import "Structs.wdl"
import  "LongReadGenotypeTasks.wdl" as LongReadGenotypeTasks

workflow ExtractTriosFromVCF {

  input {
    File input_vcf
    Array[Array[String]] families  # List of trios: [[fa1, mo1, ch1], [fa2, mo2, ch2], ...]
    String sv_pipeline_base_docker 
    File inheri_table
    String prefix

    RuntimeAttr? runtime_attr_override
    RuntimeAttr? runtime_attr_ovr_calcu_inheri_table_snv
    RuntimeAttr? runtime_attr_ovr_calcu_inheri_table_sv
    RuntimeAttr? runtime_attr_ovr_calcu_inheri_table_indel_lg
    RuntimeAttr? runtime_attr_ovr_calcu_inheri_table_indel_sm
    RuntimeAttr? runtime_attr_ovr_calcu_inheri_table_snv
    RuntimeAttr? runtime_attr_ovr_calcu_inheri_table_indel_lg
    RuntimeAttr? runtime_attr_ovr_calcu_inheri_table_indel_sm
    RuntimeAttr? runtime_attr_ovr_calcu_inheri_table_sv
  }


  scatter (family in families) {
    call LongReadGenotypeTasks.WriteTrioSampleFile as WriteTrioSampleFile {
      input:
        family = family,
        docker_image = sv_pipeline_base_docker
    }

    call LongReadGenotypeTasks.ExtractTrioVCF as ExtractTrioVCF {
      input:
        input_vcf = input_vcf,
        sample_file = WriteTrioSampleFile.sample_file,
        family_id = WriteTrioSampleFile.family_id,
        docker_image = sv_pipeline_base_docker
    }

    call LongReadGenotypeTasks.SplitVariantsBySize as SplitVariantsBySize {
      input:
        input_vcf = ExtractTrioVCF.output_vcf,
        docker_image = sv_pipeline_base_docker
    }

    call LongReadGenotypeTasks.CalculateInheritanceTable as calcu_inheri_table_snv{
      input:
        input_vcf = SplitVariantsBySize.snv_vcf,
        input_vcf_idx = SplitVariantsBySize.snv_vcf_idx,
        docker_image = sv_pipeline_base_docker,
        runtime_attr_override = runtime_attr_ovr_calcu_inheri_table_snv
    }

    call LongReadGenotypeTasks.CalculateInheritanceTable as calcu_inheri_table_indel_sm{
      input:
        input_vcf = SplitVariantsBySize.indel_1_30_vcf,
        input_vcf_idx = SplitVariantsBySize.indel_1_30_vcf_idx,
        docker_image = sv_pipeline_base_docker,
        runtime_attr_override = runtime_attr_ovr_calcu_inheri_table_indel_sm
    }

    call LongReadGenotypeTasks.CalculateInheritanceTable as calcu_inheri_table_indel_lg{
      input:
        input_vcf = SplitVariantsBySize.indel_31_50_vcf,
        input_vcf_idx = SplitVariantsBySize.indel_31_50_vcf_idx,
        docker_image = sv_pipeline_base_docker,
        runtime_attr_override = runtime_attr_ovr_calcu_inheri_table_indel_lg
    }
    
    call LongReadGenotypeTasks.CalculateInheritanceTable as calcu_inheri_table_sv{
      input:
        input_vcf = SplitVariantsBySize.sv_vcf,
        input_vcf_idx = SplitVariantsBySize.sv_vcf_idx,
        docker_image = sv_pipeline_base_docker,
        runtime_attr_override = runtime_attr_ovr_calcu_inheri_table_sv
    }

    call LongReadGenotypeTasks.EvaluateInheriByGQ as evaluate_inheri_by_gq_snv{
      input:
        vcf_file = SplitVariantsBySize.snv_vcf,
        vcf_idx_file = SplitVariantsBySize.snv_vcf_idx,
        inheri_stat = inheri_table,
        docker_image   = sv_pipeline_base_docker,
        runtime_attr_override = runtime_attr_evaluate_inheri_by_gq_snv
      }

    call LongReadGenotypeTasks.EvaluateInheriByGQ as evaluate_inheri_by_gq_indel_sm{
      input:
        vcf_file = SplitVariantsBySize.indel_1_30_vcf,
        vcf_idx_file = SplitVariantsBySize.indel_1_30_vcf_idx,
        inheri_stat = inheri_table,
        docker_image   = sv_pipeline_base_docker,
        runtime_attr_override = runtime_attr_evaluate_inheri_by_gq_indel_sm
      }

    call LongReadGenotypeTasks.EvaluateInheriByGQ as evaluate_inheri_by_gq_indel_lg{
      input:
        vcf_file = SplitVariantsBySize.indel_31_50_vcf,
        vcf_idx_file = SplitVariantsBySize.indel_31_50_vcf_idx,
        inheri_stat = inheri_table,
        docker_image   = sv_pipeline_base_docker,
        runtime_attr_override = runtime_attr_evaluate_inheri_by_gq_indel_lg
      }

    call LongReadGenotypeTasks.EvaluateInheriByGQ as evaluate_inheri_by_gq_sv{
      input:
        vcf_file = SplitVariantsBySize.sv_vcf,
        vcf_idx_file = SplitVariantsBySize.sv_vcf_idx,
        inheri_stat = inheri_table,
        docker_image   = sv_pipeline_base_docker,
        runtime_attr_override = runtime_attr_evaluate_inheri_by_gq_sv
      }

    call LongReadGenotypeTasks.IntegrateInheritanceTable as IntegrateInheritanceTable{
      input:
        inheri_table_snv = calcu_inheri_table_snv.inheri_stat,
        inheri_table_indel_sm = calcu_inheri_table_indel_sm.inheri_stat,
        inheri_table_indel_lg = calcu_inheri_table_indel_lg.inheri_stat,
        inheri_table_sv = calcu_inheri_table_sv.inheri_stat,
        inheri_table_ref = inheri_table,
        family_id = WriteTrioSampleFile.family_id,
        prefix   = prefix,
        docker_image = sv_pipeline_base_docker
    }    

    call LongReadGenotypeTasks.IntegrateInheriByGQTable as IntegrateInheriByGQTable{
      input:
        inheri_gq_table_snv = evaluate_inheri_by_gq_snv.inheri_by_GQ_stat ,
        inheri_gq_table_indel_sm = evaluate_inheri_by_gq_indel_sm.inheri_by_GQ_stat,
        inheri_gq_table_indel_lg = evaluate_inheri_by_gq_indel_lg.inheri_by_GQ_stat,
        inheri_gq_table_sv = evaluate_inheri_by_gq_sv.inheri_by_GQ_stat,
        family_id = WriteTrioSampleFile.family_id,
        prefix   = prefix,
        docker_image = sv_pipeline_base_docker
    }

  }

  output{
    Array[File] inheritance_table_inte = IntegrateInheritanceTable.integrated_inheri_stat
    Array[File] inheritance_by_gq_table_inte = IntegrateInheriByGQTable.integrated_inheri_by_gq_stat
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


