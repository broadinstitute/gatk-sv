##########################################################################################

## Base script:   https://portal.firecloud.org/#methods/Talkowski-SV/02_rdtest/10/wdl

## Github commit: talkowski-lab/gatk-sv-v1:<ENTER HASH HERE IN FIRECLOUD>

##########################################################################################

version 1.0

import "Tasks02.wdl" as tasks02
import "RDTestChromosome.wdl" as rdc

workflow RDTest {
  input {
    String batch
    File medianfile
    File autosome_contigs
    File coveragefile
    String flags
    File vcf
    Int split_size
    String algorithm
    File allosome_contigs
    File ped_file
    File male_samples
    File female_samples
    File samples
    Int tabix_retries

    String sv_pipeline_docker
    String sv_pipeline_rdtest_docker
    String linux_docker
    RuntimeAttr? runtime_attr_rdtest
    RuntimeAttr? runtime_attr_split_rd_vcf
    RuntimeAttr? runtime_attr_merge_allo
    RuntimeAttr? runtime_attr_merge_stats
  }

  Array[Array[String]] autosomes = read_tsv(autosome_contigs)
  Array[Array[String]] allosomes = read_tsv(allosome_contigs)

  scatter (autosome in autosomes) {
    call rdc.RDTestChromosome as RDTestAutosome {
      input:
        batch = batch,
        medianfile = medianfile,
        coveragefile = coveragefile,
        flags = flags,
        vcf = vcf,
        chrom = autosome[0],
        split_size = split_size,
        algorithm = algorithm,
        ped_file = ped_file,
        samples = samples,
        male_samples = male_samples,
        female_samples = female_samples,
        allosome = false,
        tabix_retries = tabix_retries,
        sv_pipeline_docker = sv_pipeline_docker,
        sv_pipeline_rdtest_docker = sv_pipeline_rdtest_docker,
        linux_docker = linux_docker,
        runtime_attr_rdtest = runtime_attr_rdtest,
        runtime_attr_split_rd_vcf = runtime_attr_split_rd_vcf,
        runtime_attr_merge_allo = runtime_attr_merge_allo,
        runtime_attr_merge_stats = runtime_attr_merge_stats
    }
  }

  scatter (allosome in allosomes) {
    call rdc.RDTestChromosome as RDTestAllosome {
      input:
        batch = batch,
        medianfile = medianfile,
        coveragefile = coveragefile,
        flags = flags,
        vcf = vcf,
        chrom = allosome[0],
        split_size = split_size,
        algorithm = algorithm,
        ped_file = ped_file,
        samples = samples,
        male_samples = male_samples,
        female_samples = female_samples,
        allosome = true,
        tabix_retries = tabix_retries,
        sv_pipeline_docker = sv_pipeline_docker,
        sv_pipeline_rdtest_docker = sv_pipeline_rdtest_docker,
        linux_docker = linux_docker,
        runtime_attr_rdtest = runtime_attr_rdtest,
        runtime_attr_split_rd_vcf = runtime_attr_split_rd_vcf,
        runtime_attr_merge_allo = runtime_attr_merge_allo,
        runtime_attr_merge_stats = runtime_attr_merge_stats
    }
  }

  call tasks02.MergeStats as MergeStats {
    input:
      stats = flatten([RDTestAutosome.stats, RDTestAllosome.stats]),
      prefix = "${batch}.${algorithm}",
      linux_docker = linux_docker,
      runtime_attr_override = runtime_attr_merge_stats
  }

  output {
    File rdtest = MergeStats.merged_stats
  }
}
