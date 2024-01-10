version 1.0

import "TasksGenerateBatchMetrics.wdl" as tasksbatchmetrics
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
    File male_only_variant_ids
    File samples
    File ref_dict

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
        male_only_variant_ids = male_only_variant_ids,
        allosome = false,
        ref_dict = ref_dict,
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
        male_only_variant_ids = male_only_variant_ids,
        allosome = true,
        ref_dict = ref_dict,
        sv_pipeline_docker = sv_pipeline_docker,
        sv_pipeline_rdtest_docker = sv_pipeline_rdtest_docker,
        linux_docker = linux_docker,
        runtime_attr_rdtest = runtime_attr_rdtest,
        runtime_attr_split_rd_vcf = runtime_attr_split_rd_vcf,
        runtime_attr_merge_allo = runtime_attr_merge_allo,
        runtime_attr_merge_stats = runtime_attr_merge_stats
    }
  }

  call tasksbatchmetrics.MergeStats as MergeStats {
    input:
      stats = flatten([RDTestAutosome.out_stats, RDTestAllosome.out_stats]),
      prefix = "${batch}.${algorithm}",
      linux_docker = linux_docker,
      runtime_attr_override = runtime_attr_merge_stats
  }

  output {
    File rdtest = MergeStats.merged_stats
  }
}
