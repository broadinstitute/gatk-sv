version 1.0

import "Structs.wdl"
import "SVConcordance.wdl" as concordance
import "FormatVcfForGatk.wdl" as format
import "RecalibrateGq.wdl" as recalibrate
import "TasksMakeCohortVcf.wdl" as tasks
import "Utils.wdl" as utils

workflow RecalibrateCleanedVcf {
  input {
    File vcf
    String prefix

    String samtools_cloud_docker
    String gatk_docker
    String sv_base_mini_docker
    String sv_pipeline_docker
    String sv_utils_docker

    RuntimeAttr? runtime_attr_shard
  }

  call format.FormatVcfForGatk {
    input:
      prefix=prefix,
      vcf=vcf,
      sv_pipeline_docker=sv_pipeline_docker,
      sv_base_mini_docker=sv_base_mini_docker
  }

  call concordance.SVConcordance {
    input:
      output_prefix=prefix,
      eval_vcf=FormatVcfForGatk.gatk_formatted_vcf,
      sv_base_mini_docker=sv_base_mini_docker,
      gatk_docker=gatk_docker
  }

  call recalibrate.RecalibrateGq {
    input:
      vcf=SVConcordance.concordance_vcf,
      vcf_index=SVConcordance.concordance_vcf_index,
      sv_base_mini_docker=sv_base_mini_docker,
      samtools_cloud_docker=samtools_cloud_docker,
      gatk_docker=gatk_docker,
      sv_base_mini_docker=sv_base_mini_docker,
      sv_pipeline_docker=sv_pipeline_docker,
      sv_utils_docker=sv_utils_docker
  }

  output {
    File recalibrated_vcf = RecalibrateGq.filtered_vcf
    File recalibrated_vcf_index = RecalibrateGq.filtered_vcf_index

    File concordance_vcf = SVConcordance.concordance_vcf
    File concordance_vcf_vcf_index = SVConcordance.concordance_vcf_index
  }
}
