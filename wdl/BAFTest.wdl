##########################################################################################

## Base script:   https://portal.firecloud.org/#methods/Talkowski-SV/02_baftest/13/wdl

## Github commit: talkowski-lab/gatk-sv-v1:<ENTER HASH HERE IN FIRECLOUD>

##########################################################################################

version 1.0

import "Tasks02.wdl" as tasks02
import "BAFTestChromosome.wdl" as bafc

workflow BAFTest {

  input {
    File vcf
    File baf_metrics
    String batch
    Array[String] samples
    String algorithm
    Int split_size
    File autosome_contigs
    Int tabix_retries

    String linux_docker
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_baftest
    RuntimeAttr? runtime_attr_split_baf_vcf
    RuntimeAttr? runtime_attr_merge_baf
    RuntimeAttr? runtime_attr_merge_stats
  }

  Array[Array[String]] autosomes = read_tsv(autosome_contigs)

  scatter (autosome in autosomes) {
    call bafc.BAFTestChromosome as BAFTestAutosome {
      input:
        batch = batch,
        samples = samples,
        baf_metrics = baf_metrics,
        algorithm = algorithm,
        vcf = vcf,
        split_size = split_size,
        chrom = autosome[0],
        tabix_retries = tabix_retries,
        linux_docker = linux_docker,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_baftest = runtime_attr_baftest,
        runtime_attr_split_baf_vcf = runtime_attr_split_baf_vcf,
        runtime_attr_merge_baf = runtime_attr_merge_baf
    }
  }

  call tasks02.MergeStats as MergeStats {
    input:
      stats = BAFTestAutosome.stats,
      prefix = "${batch}.${algorithm}",
      linux_docker = linux_docker,
      runtime_attr_override = runtime_attr_merge_stats
  }

  output {
    File baftest = MergeStats.merged_stats
  }
}

