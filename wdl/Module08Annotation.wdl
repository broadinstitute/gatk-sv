version 1.0

import "AnnotateVcf.wdl" as ann
import "PruneAndAddVafs.wdl" as pav

workflow Module07 {
  
  input {
    String vcf
    File   vcf_idx
    String prefix
    File   contig_list

    File protein_coding_gtf
    File linc_rna_gtf
    File promoter_bed
    File noncoding_bed

    File? sample_pop_assignments  # Two-column file with sample ID & pop assignment. "." for pop will ignore sample
    File? prune_list              # List of samples to be excluded from the output vcf
    File? ped_file                # Used for M/F AF calculations
    Int   sv_per_shard

    String sv_base_mini_docker
    String sv_pipeline_docker

    RuntimeAttr? runtime_attr_annotate_intervals
    RuntimeAttr? runtime_attr_merge_annotations
    RuntimeAttr? runtime_attr_subset_vcf
    RuntimeAttr? runtime_attr_concat_vcfs
    RuntimeAttr? runtime_attr_prune_vcf
    RuntimeAttr? runtime_attr_shard_vcf
    RuntimeAttr? runtime_attr_compute_AFs
    RuntimeAttr? runtime_attr_combine_vcfs
  }

  call ann.AnnotateVcf as AnnotateVcf {
    input:
      vcf                = vcf,
      vcf_idx            = vcf_idx,
      prefix             = prefix,
      contig_list        = contig_list,
      protein_coding_gtf = protein_coding_gtf,
      linc_rna_gtf       = linc_rna_gtf,
      promoter_bed       = promoter_bed,
      noncoding_bed      = noncoding_bed,
      sv_base_mini_docker     = sv_base_mini_docker,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_annotate_intervals = runtime_attr_annotate_intervals,
      runtime_attr_merge_annotations  = runtime_attr_merge_annotations,
      runtime_attr_subset_vcf         = runtime_attr_subset_vcf,
      runtime_attr_concat_vcfs        = runtime_attr_concat_vcfs
  }

  call pav.PruneAndAddVafs as PruneAndAddVafs {
    input:
      vcf                    = AnnotateVcf.annotated_vcf,
      vcf_idx                = AnnotateVcf.annotated_vcf_idx,
      prefix                 = prefix,
      sample_pop_assignments = sample_pop_assignments,
      prune_list             = prune_list,
      ped_file               = ped_file,
      sv_per_shard           = sv_per_shard,
      contig_list            = contig_list,
      sv_base_mini_docker     = sv_base_mini_docker,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_prune_vcf    = runtime_attr_prune_vcf,
      runtime_attr_shard_vcf    = runtime_attr_shard_vcf,
      runtime_attr_compute_AFs  = runtime_attr_compute_AFs,
      runtime_attr_combine_vcfs = runtime_attr_combine_vcfs,
      runtime_attr_concat_vcfs  = runtime_attr_concat_vcfs
  }

  output {
    File output_vcf     = PruneAndAddVafs.output_vcf
    File output_vcf_idx = PruneAndAddVafs.output_vcf_idx
  }
}
