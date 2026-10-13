version 1.0

import "Structs.wdl"
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

    call LongReadGenotypeTasks.SplitMultiAllelicToBiAllelicSingleSample{
      input:
        vcf_file = ExtractVariantIndividualGenome.non_ref_vcf,
        vcf_idx = ExtractVariantIndividualGenome.non_ref_vcf_idx,
        docker_image = sv_pipeline_base_docker
    }

    # SNV / indel 1-20 / indel 21-49 / SV >= 50 bp; each variant in exactly one bin.
    call LongReadGenotypeTasks.SplitVariantsBySizeBenchmark{
      input:
        input_vcf = SplitMultiAllelicToBiAllelicSingleSample.biallelic_vcf,
        docker_image = sv_pipeline_base_docker
    }

    call LongReadGenotypeTasks.ExtractVariantSitesSimple{
      input:
        input_vcf = SplitMultiAllelicToBiAllelicSingleSample.biallelic_vcf,
        docker_image = sv_pipeline_base_docker
    }
  }

  output {
    Array[File] all_variant_sites = ExtractVariantSitesSimple.variant_sites
    Array[File] all_snv_vcfs = SplitVariantsBySizeBenchmark.snv_vcf
    Array[File] all_indel_1_20 = SplitVariantsBySizeBenchmark.indel_1_20_vcf
    Array[File] all_indel_1_20_idx = SplitVariantsBySizeBenchmark.indel_1_20_vcf_idx
    Array[File] all_indel_21_49 = SplitVariantsBySizeBenchmark.indel_21_49_vcf
    Array[File] all_indel_21_49_idx = SplitVariantsBySizeBenchmark.indel_21_49_vcf_idx
    Array[File] all_sv_ge50 = SplitVariantsBySizeBenchmark.sv_ge50_vcf
    Array[File] all_sv_ge50_idx = SplitVariantsBySizeBenchmark.sv_ge50_vcf_idx
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
