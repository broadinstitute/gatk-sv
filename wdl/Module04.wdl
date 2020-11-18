##########################################################################################

## Base script:   https://portal.firecloud.org/#methods/Talkowski-SV/04_v2_genotype_batch/57/wdl

## Github commit: talkowski-lab/gatk-sv-v1:<ENTER HASH HERE IN FIRECLOUD>

##########################################################################################

version 1.0

import "GenotypePESRPart1.wdl" as gp1
import "GenotypePESRPart2.wdl" as gp2
import "GenotypeDepthPart1.wdl" as gd1
import "GenotypeDepthPart2.wdl" as gd2
import "Tasks04.wdl" as tasks04
import "Utils.wdl" as util

workflow Module04 {
  input {
    File batch_pesr_vcf
    File batch_depth_vcf
    File cohort_pesr_vcf
    File cohort_depth_vcf
    String batch

    Int n_per_split
    File coveragefile       # batch coverage file
    File medianfile         # post-exclusion batch median file
    File famfile            # post-exclusion batch famfile
    File? rf_cutoffs         # Random forest cutoffs; required unless skipping training
    File? seed_cutoffs      # Required unless skipping training
    Int n_RD_genotype_bins  # number of RdTest bins
    File discfile
    File? pesr_exclude_list    # Required unless skipping training
    File splitfile
    String? reference_build  #hg19 or hg38, Required unless skipping training
    File bin_exclude
    File ref_dict
    # If all specified, training will be skipped (for single sample pipeline)
    File? genotype_pesr_pesr_sepcutoff
    File? genotype_pesr_depth_sepcutoff
    File? genotype_depth_pesr_sepcutoff
    File? genotype_depth_depth_sepcutoff
    File? SR_metrics
    File? PE_metrics

    String sv_base_mini_docker
    String sv_pipeline_docker
    String sv_pipeline_rdtest_docker
    String linux_docker

    # Common
    RuntimeAttr? runtime_attr_split_vcf
    RuntimeAttr? runtime_attr_merge_counts
    RuntimeAttr? runtime_attr_split_variants
    RuntimeAttr? runtime_attr_make_subset_vcf
    RuntimeAttr? runtime_attr_rdtest_genotype
    RuntimeAttr? runtime_attr_add_genotypes
    RuntimeAttr? runtime_attr_concat_vcfs

    # Master
    RuntimeAttr? runtime_attr_add_batch
    RuntimeAttr? runtime_attr_index_vcf
    RuntimeAttr? runtime_attr_ids_from_vcf

    # PE train
    RuntimeAttr? runtime_attr_make_batch_bed
    RuntimeAttr? runtime_attr_count_pe
    RuntimeAttr? runtime_attr_pe_genotype

    # SR train
    RuntimeAttr? runtime_attr_count_sr
    RuntimeAttr? runtime_attr_sr_genotype

    # RD train
    RuntimeAttr? runtime_attr_training_bed
    RuntimeAttr? runtime_attr_genotype_train
    RuntimeAttr? runtime_attr_generate_cutoff
    RuntimeAttr? runtime_attr_update_cutoff
    RuntimeAttr? runtime_attr_split_variants
    RuntimeAttr? runtime_attr_merge_genotypes

    # PESR part 2
    RuntimeAttr? runtime_attr_count_pe
    RuntimeAttr? runtime_attr_genotype_pe
    RuntimeAttr? runtime_attr_count_sr
    RuntimeAttr? runtime_attr_genotype_sr
    RuntimeAttr? runtime_attr_integrate_gq
    RuntimeAttr? runtime_attr_integrate_pesr_gq
    RuntimeAttr? runtime_attr_triple_stream_cat

    # Depth part 2
    RuntimeAttr? runtime_attr_integrate_depth_gq
    RuntimeAttr? runtime_attr_concat_vcfs

  }

  Boolean single_sample_mode = defined(genotype_pesr_pesr_sepcutoff) && defined(genotype_pesr_depth_sepcutoff) && defined(genotype_depth_depth_sepcutoff) && defined(SR_metrics) && defined(PE_metrics)
  call tasks04.AddBatchSamples as AddBatchSamplesPESR {
    input:
      batch_vcf = batch_pesr_vcf,
      cohort_vcf = cohort_pesr_vcf,
      prefix = "${batch}.pesr",
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_add_batch
  }

  call tasks04.AddBatchSamples as AddBatchSamplesDepth {
    input:
      batch_vcf = batch_depth_vcf,
      cohort_vcf = cohort_depth_vcf,
      prefix = "${batch}.depth",
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_add_batch
  }

  call util.GetSampleIdsFromVcf {
    input:
      vcf = batch_pesr_vcf,
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_ids_from_vcf
  }

  if (!single_sample_mode) {
    call gp1.GenotypePESRPart1 as GenotypePESRPart1 {
      input:
        bin_exclude=bin_exclude,
        samples = GetSampleIdsFromVcf.out_array,
        pesr_exclude_list = select_first([pesr_exclude_list]),
        discfile = discfile,
        n_RD_genotype_bins = n_RD_genotype_bins,
        batch_vcf = batch_pesr_vcf,
        seed_cutoffs = select_first([seed_cutoffs]),
        medianfile = medianfile,
        batch = batch,
        rf_cutoffs = select_first([rf_cutoffs]),
        coveragefile = coveragefile,
        reference_build = select_first([reference_build]),
        n_per_PE_split = n_per_split,
        famfile = famfile,
        splitfile = splitfile,
        n_per_RD_split = n_per_split,
        n_per_SR_split = n_per_split,
        ref_dict = ref_dict,
        sv_base_mini_docker = sv_base_mini_docker,
        sv_pipeline_docker = sv_pipeline_docker,
        sv_pipeline_rdtest_docker = sv_pipeline_rdtest_docker,
        runtime_attr_split_vcf = runtime_attr_split_vcf,
        runtime_attr_merge_counts = runtime_attr_merge_counts,
        runtime_attr_make_batch_bed = runtime_attr_make_batch_bed,
        runtime_attr_count_pe = runtime_attr_count_pe,
        runtime_attr_pe_genotype = runtime_attr_pe_genotype,
        runtime_attr_count_sr = runtime_attr_count_sr,
        runtime_attr_sr_genotype = runtime_attr_sr_genotype,
        runtime_attr_training_bed = runtime_attr_training_bed,
        runtime_attr_genotype_train = runtime_attr_genotype_train,
        runtime_attr_generate_cutoff = runtime_attr_generate_cutoff,
        runtime_attr_update_cutoff = runtime_attr_update_cutoff,
        runtime_attr_split_variants = runtime_attr_split_variants,
        runtime_attr_rdtest_genotype = runtime_attr_rdtest_genotype,
        runtime_attr_merge_genotypes = runtime_attr_merge_genotypes
    }
  }

  call gp2.GenotypePESRPart2 as GenotypePESRPart2 {
    input:
      bin_exclude=bin_exclude,
      samples = GetSampleIdsFromVcf.out_array,
      discfile = discfile,
      PE_metrics = select_first([PE_metrics, GenotypePESRPart1.PE_metrics]),
      n_RdTest_bins = n_RD_genotype_bins,
      medianfile = medianfile,
      cohort_vcf = AddBatchSamplesPESR.updated_vcf,
      batch = batch,
      RD_depth_sepcutoff = select_first([genotype_pesr_depth_sepcutoff, GenotypePESRPart1.RD_depth_sepcutoff]),
      RD_pesr_sepcutoff = select_first([genotype_pesr_pesr_sepcutoff, GenotypePESRPart1.RD_pesr_sepcutoff]),
      coveragefile = coveragefile,
      SR_metrics = select_first([SR_metrics, GenotypePESRPart1.SR_metrics]),
      n_per_split = n_per_split,
      famfile = famfile,
      splitfile = splitfile,
      ref_dict = ref_dict,
      sv_base_mini_docker = sv_base_mini_docker,
      sv_pipeline_docker = sv_pipeline_docker,
      sv_pipeline_rdtest_docker = sv_pipeline_rdtest_docker,
      linux_docker = linux_docker,
      runtime_attr_split_variants = runtime_attr_split_variants,
      runtime_attr_make_subset_vcf = runtime_attr_make_subset_vcf,
      runtime_attr_count_pe = runtime_attr_count_pe,
      runtime_attr_genotype_pe = runtime_attr_genotype_pe,
      runtime_attr_count_sr = runtime_attr_count_sr,
      runtime_attr_genotype_sr = runtime_attr_genotype_sr,
      runtime_attr_rdtest_genotype = runtime_attr_rdtest_genotype,
      runtime_attr_integrate_gq = runtime_attr_integrate_gq,
      runtime_attr_integrate_pesr_gq = runtime_attr_integrate_pesr_gq,
      runtime_attr_add_genotypes = runtime_attr_add_genotypes,
      runtime_attr_triple_stream_cat = runtime_attr_triple_stream_cat,
      runtime_attr_concat_vcfs = runtime_attr_concat_vcfs
  }

  if (!single_sample_mode) {
    call gd1.GenotypeDepthPart1 as GenotypeDepthPart1 {
      input:
        bin_exclude=bin_exclude,
        samples = GetSampleIdsFromVcf.out_array,
        n_RD_genotype_bins = n_RD_genotype_bins,
        batch_vcf = batch_depth_vcf,
        seed_cutoffs = select_first([seed_cutoffs]),
        medianfile = medianfile,
        batch = batch,
        rf_cutoffs = select_first([rf_cutoffs]),
        coveragefile = coveragefile,
        reference_build = select_first([reference_build]),
        famfile = famfile,
        n_per_RD_split = n_per_split,
        ref_dict = ref_dict,
        sv_base_mini_docker = sv_base_mini_docker,
        sv_pipeline_docker = sv_pipeline_docker,
        sv_pipeline_rdtest_docker = sv_pipeline_rdtest_docker,
        runtime_attr_training_bed = runtime_attr_training_bed,
        runtime_attr_genotype_train = runtime_attr_genotype_train,
        runtime_attr_generate_cutoff = runtime_attr_generate_cutoff,
        runtime_attr_update_cutoff = runtime_attr_update_cutoff,
        runtime_attr_split_variants = runtime_attr_split_variants,
        runtime_attr_rdtest_genotype = runtime_attr_rdtest_genotype,
        runtime_attr_merge_genotypes = runtime_attr_merge_genotypes
    }
  }
  call gd2.GenotypeDepthPart2 as GenotypeDepthPart2 {
    input:
      bin_exclude=bin_exclude,
      samples = GetSampleIdsFromVcf.out_array,
      n_RdTest_bins = n_RD_genotype_bins,
      medianfile = medianfile,
      cohort_vcf = AddBatchSamplesDepth.updated_vcf,
      batch = batch,
      RD_depth_sepcutoff = select_first([genotype_depth_depth_sepcutoff, GenotypeDepthPart1.RD_depth_sepcutoff]),
      RD_pesr_sepcutoff = select_first([genotype_depth_pesr_sepcutoff, GenotypeDepthPart1.RD_pesr_sepcutoff]),
      coveragefile = coveragefile,
      n_per_split = n_per_split,
      famfile = famfile,
      ref_dict = ref_dict,
      sv_base_mini_docker = sv_base_mini_docker,
      sv_pipeline_docker = sv_pipeline_docker,
      sv_pipeline_rdtest_docker = sv_pipeline_rdtest_docker,
      runtime_attr_split_variants = runtime_attr_split_variants,
      runtime_attr_rdtest_genotype = runtime_attr_rdtest_genotype,
      runtime_attr_make_subset_vcf = runtime_attr_make_subset_vcf,
      runtime_attr_integrate_depth_gq = runtime_attr_integrate_depth_gq,
      runtime_attr_add_genotypes = runtime_attr_add_genotypes,
      runtime_attr_concat_vcfs = runtime_attr_concat_vcfs
  }
  output {
    File sr_bothside_pass = GenotypePESRPart2.bothside_pass
    File sr_background_fail = GenotypePESRPart2.background_fail

    File? trained_PE_metrics = GenotypePESRPart1.PE_metrics
    File? trained_SR_metrics = GenotypePESRPart1.SR_metrics

    File? trained_genotype_pesr_pesr_sepcutoff = GenotypePESRPart1.RD_pesr_sepcutoff
    File? trained_genotype_pesr_depth_sepcutoff = GenotypePESRPart1.RD_depth_sepcutoff
    File? trained_genotype_depth_pesr_sepcutoff = GenotypeDepthPart1.RD_pesr_sepcutoff
    File? trained_genotype_depth_depth_sepcutoff = GenotypeDepthPart1.RD_depth_sepcutoff
    
    File genotyped_depth_vcf = GenotypeDepthPart2.genotyped_vcf
    File genotyped_depth_vcf_index = GenotypeDepthPart2.genotyped_vcf_index
    File genotyped_pesr_vcf = GenotypePESRPart2.genotyped_vcf
    File genotyped_pesr_vcf_index = GenotypePESRPart2.genotyped_vcf_index
    Array[File] regeno_coverage_medians = GenotypeDepthPart2.regeno_coverage_medians
  }
}


