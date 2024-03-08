version 1.0

import "Module07XfBatchEffect.wdl" as batch_effect

workflow BatchEffectAcrossContigs {
  input {
    Array[File] vcfs
    Array[File] vcf_idxs
    File sample_batch_assignments
    File batches_list
    File sample_pop_assignments
    File famfile
    File contiglist
    File? par_bed
    Int? onevsall_cutoff=2
    String prefix

    String sv_pipeline_docker
    String sv_base_mini_docker
  }

  Array[String] contigs = read_lines(contiglist)

  scatter ( i in range(length(vcfs)) ) {
    call batch_effect.XfBatchEffect {
      input:
        vcf=vcfs[i],
        vcf_idx=vcf_idxs[i],
        sample_batch_assignments=sample_batch_assignments,
        batches_list=batches_list,
        sample_pop_assignments=sample_pop_assignments,
        famfile=famfile,
        par_bed=par_bed,
        onevsall_cutoff=onevsall_cutoff,
        prefix="~{prefix}.{contigs[i]}",
        sv_base_mini_docker=sv_base_mini_docker,
        sv_pipeline_docker=sv_pipeline_docker
    }
  }

  output {
    Array[File] batch_effects_labeled_vcfs = XfBatchEffect.labeled_vcf
    Array[File] batch_effects_labeled_vcf_indexes = XfBatchEffect.labeled_vcf_idx
  }
}
