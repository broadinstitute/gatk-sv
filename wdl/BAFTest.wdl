version 1.0

import "TasksGenerateBatchMetrics.wdl" as tasksbatchmetrics
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
    File ref_dict

    String linux_docker
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_baftest
    RuntimeAttr? runtime_attr_split_baf_vcf
    RuntimeAttr? runtime_attr_merge_baf
    RuntimeAttr? runtime_attr_merge_stats
  }

  Array[Array[String]] autosomes = read_tsv(autosome_contigs)

  scatter (autosome in autosomes) {
    if (autosome[0] != "") {
      call bafc.BAFTestChromosome as BAFTestAutosome {
        input:
          batch = batch,
          samples = samples,
          baf_metrics = baf_metrics,
          algorithm = algorithm,
          vcf = vcf,
          split_size = split_size,
          chrom = autosome[0],
          ref_dict = ref_dict,
          linux_docker = linux_docker,
          sv_pipeline_docker = sv_pipeline_docker,
          runtime_attr_baftest = runtime_attr_baftest,
          runtime_attr_split_baf_vcf = runtime_attr_split_baf_vcf,
          runtime_attr_merge_baf = runtime_attr_merge_baf
      }
    }
  }

  Array[File] baf_stats = select_all(BAFTestAutosome.stats)
  if (length(baf_stats) > 0) {
    call tasksbatchmetrics.MergeStats as MergeStats {
      input:
        stats = baf_stats,
        prefix = "${batch}.${algorithm}",
        linux_docker = linux_docker,
        runtime_attr_override = runtime_attr_merge_stats
    }
  }
  

  output {
    File? baftest = MergeStats.merged_stats
  }
}

