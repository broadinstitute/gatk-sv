##########################################################################################

## Base script:   https://portal.firecloud.org/#methods/Talkowski-SV/02_petest/11/wdl

## Github commit: talkowski-lab/gatk-sv-v1:<ENTER HASH HERE IN FIRECLOUD>

##########################################################################################

version 1.0

import "Tasks02.wdl" as tasks02
import "PETestChromosome.wdl" as pec

workflow PETest {
  input {
    File vcf
    String algorithm
    String batch
    File discfile
    File medianfile
    Int split_size
    File allosome_contigs
    File autosome_contigs
    File ped_file
    File male_samples
    File female_samples
    File samples
    Int common_cnv_size_cutoff
    Int tabix_retries

    String sv_base_mini_docker
    String linux_docker
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_sex_list
    RuntimeAttr? runtime_attr_split_vcf
    RuntimeAttr? runtime_attr_petest
    RuntimeAttr? runtime_attr_merge_allo
    RuntimeAttr? runtime_attr_merge_stats
  }

  Array[Array[String]] autosomes = read_tsv(autosome_contigs)
  Array[Array[String]] allosomes = read_tsv(allosome_contigs)

  scatter (autosome in autosomes) {
    call pec.PETestChromosome as PETestAutosome {
      input:
        medianfile = medianfile,
        algorithm = algorithm,
        vcf = vcf,
        chrom = autosome[0],
        ped_file = ped_file,
        split_size = split_size,
        batch = batch,
        discfile = discfile,
        samples = samples,
        male_samples = male_samples,
        female_samples = female_samples,
        allosome = false,
        common_cnv_size_cutoff = common_cnv_size_cutoff,
        tabix_retries = tabix_retries,
        sv_base_mini_docker = sv_base_mini_docker,
        linux_docker = linux_docker,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_split_vcf = runtime_attr_split_vcf,
        runtime_attr_petest = runtime_attr_petest,
        runtime_attr_merge_allo = runtime_attr_merge_allo,
        runtime_attr_merge_stats = runtime_attr_merge_stats
    }

  }

  scatter (allosome in allosomes) {
    call pec.PETestChromosome as PETestAllosome {
      input:
        medianfile = medianfile,
        algorithm = algorithm,
        vcf = vcf,
        chrom = allosome[0],
        split_size = split_size,
        batch = batch,
        ped_file = ped_file,
        discfile = discfile,
        samples = samples,
        male_samples = male_samples,
        female_samples = female_samples,
        allosome = true,
        common_cnv_size_cutoff = common_cnv_size_cutoff,
        tabix_retries = tabix_retries,
        sv_base_mini_docker = sv_base_mini_docker,
        linux_docker = linux_docker,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_split_vcf = runtime_attr_split_vcf,
        runtime_attr_petest = runtime_attr_petest,
        runtime_attr_merge_allo = runtime_attr_merge_allo,
        runtime_attr_merge_stats = runtime_attr_merge_stats
    }
  }

  call tasks02.MergeStats as MergeStats {
    input:
      stats = flatten([PETestAutosome.stats, PETestAllosome.stats]),
      prefix = "${batch}.${algorithm}",
      linux_docker = linux_docker,
      runtime_attr_override = runtime_attr_merge_stats
  }

  call tasks02.MergeStats as MergeStatsCommon {
    input:
      stats = select_all(PETestAutosome.stats_common),
      prefix = "${batch}.${algorithm}.common",
      linux_docker = linux_docker,
      runtime_attr_override = runtime_attr_merge_stats
  }

  output {
    File petest = MergeStats.merged_stats
    File petest_common = MergeStatsCommon.merged_stats
  }
}
