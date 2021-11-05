version 1.0

import "TasksGenerateBatchMetrics.wdl" as tasksbatchmetrics
import "SRTestChromosome.wdl" as src

workflow SRTest {
  input {
    File autosome_contigs
    File allosome_contigs
    File splitfile
    Int split_size
    File medianfile
    File ped_file
    String batch
    File vcf
    String algorithm
    File male_samples
    File female_samples
    File male_only_variant_ids
    File samples
    Boolean run_common
    Int? common_cnv_size_cutoff  # Required if run_common is true
    File ref_dict

    String sv_pipeline_docker
    String linux_docker
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_sex_list
    RuntimeAttr? runtime_attr_split_vcf
    RuntimeAttr? runtime_attr_srtest
    RuntimeAttr? runtime_attr_merge_allo
    RuntimeAttr? runtime_attr_merge_stats
  }

  Array[Array[String]] autosomes = read_tsv(autosome_contigs)
  Array[Array[String]] allosomes = read_tsv(allosome_contigs)

  # Run srtest on each autosome
  scatter (autosome in autosomes) {
    call src.SRTestChromosome as SRTestAutosome {
      input:
        splitfile = splitfile,
        split_size = split_size,
        ped_file = ped_file,
        batch = batch,
        medianfile = medianfile,
        chrom = autosome[0],
        vcf = vcf,
        algorithm = algorithm,
        samples = samples,
        male_samples = male_samples,
        female_samples = female_samples,
        male_only_variant_ids = male_only_variant_ids,
        allosome = false,
        run_common = run_common,
        common_cnv_size_cutoff = common_cnv_size_cutoff,
        ref_dict = ref_dict,
        sv_base_mini_docker = sv_base_mini_docker,
        linux_docker = linux_docker,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_split_vcf = runtime_attr_split_vcf,
        runtime_attr_srtest = runtime_attr_srtest,
        runtime_attr_merge_allo = runtime_attr_merge_allo,
        runtime_attr_merge_stats = runtime_attr_merge_stats
    }
  }

  # Run srtest on each allosome
  scatter (allosome in allosomes) {
    call src.SRTestChromosome as SRTestAllosome {
      input:
        splitfile = splitfile,
        split_size = split_size,
        ped_file = ped_file,
        batch = batch,
        medianfile = medianfile,
        chrom = allosome[0],
        vcf = vcf,
        algorithm = algorithm,
        samples = samples,
        male_samples = male_samples,
        female_samples = female_samples,
        male_only_variant_ids = male_only_variant_ids,
        allosome = true,
        run_common = run_common,
        common_cnv_size_cutoff = common_cnv_size_cutoff,
        ref_dict = ref_dict,
        sv_base_mini_docker = sv_base_mini_docker,
        linux_docker = linux_docker,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_split_vcf = runtime_attr_split_vcf,
        runtime_attr_srtest = runtime_attr_srtest,
        runtime_attr_merge_allo = runtime_attr_merge_allo,
        runtime_attr_merge_stats = runtime_attr_merge_stats
    }
  }

  # Combine srtest results into single file
  call tasksbatchmetrics.MergeStats as MergeStats {
    input:
      stats = flatten([SRTestAutosome.stats, SRTestAllosome.stats]),
      prefix = "${batch}.${algorithm}",
      linux_docker = linux_docker,
      runtime_attr_override = runtime_attr_merge_stats
  }

  if (run_common) {
    call tasksbatchmetrics.MergeStats as MergeStatsCommon {
      input:
        stats = select_all(SRTestAutosome.stats_common),
        prefix = "${batch}.${algorithm}.common",
        linux_docker = linux_docker,
        runtime_attr_override = runtime_attr_merge_stats
    }
  }

  output {
    File srtest = MergeStats.merged_stats
    File? srtest_common = MergeStatsCommon.merged_stats
  }
}


