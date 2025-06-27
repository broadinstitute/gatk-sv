version 1.0

import AnnotateGenomicContext.wdl as AnnotateGenomicContext
workflow ExtractIndividualFromVCF {
  input {
    File vcf_file
    Array[String] sample_ids
    String sv_pipeline_base_docker
  }

  scatter (sample_id in sample_ids) {
    call AnnotateGenomicContext.ExtractVariantIndividualGenome {
      input:
        vcf_file = vcf_file,
        sample_id = sample_id,
        docker_image = sv_pipeline_base_docker
    }

    call AnnotateGenomicContext.SplitVariantsBySize{
      input:
        input_vcf = ExtractVariantIndividualGenome.non_ref_vcf,
        docker_image = sv_pipeline_base_docker
    }

    call AnnotateGenomicContext.ConcatVcfs as concat_vcf_1to50bp{
      input:
        vcf1 = SplitVariantsBySize.indel_1_30_vcf,
        vcf2 = SplitVariantsBySize.indel_31_50_vcf,
        idx1 = SplitVariantsBySize.indel_1_30_vcf_idx,
        idx2 = SplitVariantsBySize.indel_31_50_vcf_idx,
        input_vcf = ExtractVariantIndividualGenome.non_ref_vcf,
        appendix = 'indel_1_50',
        docker_image = sv_pipeline_base_docker
    }

    call AnnotateGenomicContext.ConcatVcfs as concat_vcf_over30bp{
      input:
        vcf1 = SplitVariantsBySize.indel_31_50_vcf,
        vcf2 = SplitVariantsBySize.sv_vcf,
        idx1 = SplitVariantsBySize.indel_31_50_vcf_idx,
        idx2 = SplitVariantsBySize.sv_vcf_idx,
        input_vcf = ExtractVariantIndividualGenome.non_ref_vcf,
        appendix = 'sv_over30',
        docker_image = sv_pipeline_base_docker
    }
  }


  output {
    Array[File] all_snv_vcfs = SplitVariantsBySize.snv_vcf  
    Array[File] all_indel_1_30 = SplitVariantsBySize.indel_1_30_vcf
    Array[File] all_indel_31_50 = SplitVariantsBySize.indel_31_50_vcf
    Array[File] all_sv_over50 = SplitVariantsBySize.sv_vcf
    Array[File] all_sv_over30 = concat_vcf_over30bp.combined_vcf
    Array[File] all_indel_1_50 = concat_vcf_1to50bp.combined_vcf
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
