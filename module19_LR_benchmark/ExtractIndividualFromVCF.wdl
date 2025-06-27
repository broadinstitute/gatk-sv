version 1.0

import "LongReadGenotypeTasks.wdl" as LongReadGenotypeTasks

workflow ExtractIndividualFromVCF {
  input {
    File vcf_file
    Array[String] sample_ids
    String sv_base_mini_docker
    String sv_pipeline_base_docker
  }

  scatter (sample_id in sample_ids) {
    call LongReadGenotypeTasks.ExtractVariantIndividualGenome {
      input:
        vcf_file = vcf_file,
        sample_id = sample_id,
        docker_image = sv_pipeline_base_docker
    }

    call LongReadGenotypeTasks.SplitVariantsBySize{
      input:
        input_vcf = ExtractVariantIndividualGenome.non_ref_vcf,
        docker_image = sv_pipeline_base_docker
    }

    String prefix = basename(ExtractVariantIndividualGenome.non_ref_vcf, ".vcf.gz")

    call LongReadGenotypeTasks.ConcatVcfs as concat_vcf_1to50bp{
      input:
        vcfs = [SplitVariantsBySize.indel_1_30_vcf, SplitVariantsBySize.indel_31_50_vcf],
        vcfs_idx = [SplitVariantsBySize.indel_1_30_vcf_idx, SplitVariantsBySize.indel_31_50_vcf_idx],
        outfile_prefix = "~{prefix}.indel_1_50",
        sv_base_mini_docker = sv_base_mini_docker
    }

    call LongReadGenotypeTasks.ConcatVcfs as concat_vcf_over30bp{
      input:
        vcfs = [SplitVariantsBySize.indel_31_50_vcf, SplitVariantsBySize.sv_vcf],
        vcfs_idx = [SplitVariantsBySize.indel_31_50_vcf_idx, SplitVariantsBySize.sv_vcf_idx],
        outfile_prefix = "~{prefix}.sv_over30",
        sv_base_mini_docker = sv_base_mini_docker
    }
  }


  output {
    Array[File] all_snv_vcfs = SplitVariantsBySize.snv_vcf  
    Array[File] all_indel_1_30 = SplitVariantsBySize.indel_1_30_vcf
    Array[File] all_indel_31_50 = SplitVariantsBySize.indel_31_50_vcf
    Array[File] all_sv_over50 = SplitVariantsBySize.sv_vcf
    Array[File] all_sv_over30 = concat_vcf_over30bp.concat_vcf
    Array[File] all_indel_1_50 = concat_vcf_1to50bp.concat_vcf
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
