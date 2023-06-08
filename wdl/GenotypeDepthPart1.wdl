version 1.0

import "TrainRDGenotyping.wdl" as rd_train

workflow GenotypeDepthPart1 {
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
    String reference_build   #hg19 or hg38
    File ref_dict

    String sv_base_mini_docker
    String sv_pipeline_rdtest_docker
    String sv_pipeline_docker
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
      seed_cutoffs = seed_cutoffs,
      medianfile = medianfile,
      coveragefile = coveragefile,
      coveragefile_index = coveragefile_index,
      prefix = "~{batch}.depth",
      n_bins = n_RD_genotype_bins,
      reference_build = reference_build,
      samples = samples,
      n_per_split = n_per_RD_split,
      vcf = batch_vcf,
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

  output {
    File RD_pesr_sepcutoff = TrainRDGenotyping.pesr_sepcutoff
    File RD_depth_sepcutoff = TrainRDGenotyping.depth_sepcutoff
  }
}


