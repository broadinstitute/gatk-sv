version 1.0

import "TrainRDGenotyping.wdl" as rd_train
import "TrainPEGenotyping.wdl" as pe_train
import "TrainSRGenotyping.wdl" as sr_train

workflow GenotypePESRPart1 {
  input {
    File bin_exclude
    File batch_vcf
    String batch
    File coveragefile        # batch coverage file
    File? coveragefile_index # batch coverage index file
    File medianfile          # batch median file
    File rf_cutoffs          # Random forest cutoffs
    File seed_cutoffs
    Array[String] samples    # List of samples in batch
    Int n_RD_genotype_bins   # number of RdTest bins
    Int n_per_RD_split       # number of variants per RdTest split
    Int n_per_PE_split
    File discfile
    File? discfile_index
    File pesr_exclude_list
    File splitfile
    File? splitfile_index
    Int n_per_SR_split
    String reference_build   #hg19 or hg38
    File ref_dict

    String sv_base_mini_docker
    String sv_pipeline_docker
    String sv_pipeline_rdtest_docker

    # Runtime attributes
    RuntimeAttr? runtime_attr_split_vcf
    RuntimeAttr? runtime_attr_merge_counts

    # PE
    RuntimeAttr? runtime_attr_make_batch_bed
    RuntimeAttr? runtime_attr_count_pe
    RuntimeAttr? runtime_attr_pe_genotype

    # SR
    RuntimeAttr? runtime_attr_count_sr
    RuntimeAttr? runtime_attr_sr_genotype

    # RD
    RuntimeAttr? runtime_attr_training_bed
    RuntimeAttr? runtime_attr_genotype_train
    RuntimeAttr? runtime_attr_generate_cutoff
    RuntimeAttr? runtime_attr_update_cutoff
    RuntimeAttr? runtime_attr_split_variants
    RuntimeAttr? runtime_attr_rdtest_genotype
    RuntimeAttr? runtime_attr_merge_genotypes
  }

  call rd_train.TrainRDGenotyping as TrainRDGenotyping {
    input:
      bin_exclude=bin_exclude,
      rf_cutoffs = rf_cutoffs,
      medianfile = medianfile,
      n_bins = n_RD_genotype_bins,
      vcf = batch_vcf,
      coveragefile = coveragefile,
      coveragefile_index = coveragefile_index,
      n_per_split = n_per_RD_split,
      prefix = "~{batch}.pesr",
      seed_cutoffs = seed_cutoffs,
      reference_build = reference_build,
      samples = samples,
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

  call pe_train.TrainPEGenotyping as TrainPEGenotyping {
    input:
      RD_melted_genotypes = TrainRDGenotyping.melted_genotypes,
      batch_vcf = batch_vcf,
      medianfile = medianfile,
      RF_cutoffs = rf_cutoffs,
      RD_genotypes = TrainRDGenotyping.genotypes,
      exclude_list = pesr_exclude_list,
      n_per_split = n_per_PE_split,
      discfile = discfile,
      discfile_index = discfile_index,
      samples = samples,
      batch_ID = batch,
      ref_dict = ref_dict,
      sv_base_mini_docker = sv_base_mini_docker,
      sv_pipeline_docker = sv_pipeline_docker,
      sv_pipeline_rdtest_docker = sv_pipeline_rdtest_docker,
      runtime_attr_split_vcf = runtime_attr_split_vcf,
      runtime_attr_make_batch_bed = runtime_attr_make_batch_bed,
      runtime_attr_merge_counts = runtime_attr_merge_counts,
      runtime_attr_count_pe = runtime_attr_count_pe,
      runtime_attr_genotype = runtime_attr_pe_genotype
  }

  call sr_train.TrainSRGenotyping as TrainSRGenotyping {
    input:
      splitfile = splitfile,
      splitfile_index = splitfile_index,
      RD_melted_genotypes = TrainRDGenotyping.melted_genotypes,
      batch_vcf = batch_vcf,
      medianfile = medianfile,
      RF_cutoffs = rf_cutoffs,
      n_per_split = n_per_SR_split,
      PE_train = TrainPEGenotyping.PE_train,
      PE_genotypes = TrainPEGenotyping.PE_genotypes,
      samples = samples,
      batch_ID = batch,
      ref_dict = ref_dict,
      sv_base_mini_docker = sv_base_mini_docker,
      sv_pipeline_docker = sv_pipeline_docker,
      sv_pipeline_rdtest_docker = sv_pipeline_rdtest_docker,
      runtime_attr_split_vcf = runtime_attr_split_vcf,
      runtime_attr_merge_counts = runtime_attr_merge_counts,
      runtime_attr_count_sr = runtime_attr_count_sr,
      runtime_attr_genotype = runtime_attr_sr_genotype
  }

  output {
    File RD_depth_sepcutoff = TrainRDGenotyping.depth_sepcutoff
    File SR_metrics = TrainSRGenotyping.SR_metrics
    File PE_metrics = TrainPEGenotyping.PE_metrics
    File RD_pesr_sepcutoff = TrainRDGenotyping.pesr_sepcutoff
  }
}


