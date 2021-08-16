version 1.0

import "ScatterAnnotateVcfByChrom.wdl" as ann
import "PruneAndAddVafs.wdl" as pav
import "AnnotateExternalAF.wdl" as eaf

workflow AnnotateVcf {

  input {
    File vcf
    File vcf_idx
    File contig_list
    String prefix

    File protein_coding_gtf
    File linc_rna_gtf
    File promoter_bed
    File noncoding_bed

    Int max_shards_per_chrom_step1
    Int min_records_per_shard_step1

    File? sample_pop_assignments  # Two-column file with sample ID & pop assignment. "." for pop will ignore sample
    File? prune_list              # List of samples to be excluded from the output vcf
    File? ped_file                # Used for M/F AF calculations
    Int   sv_per_shard

    File? ref_bed              # File with external allele frequencies
    String? ref_prefix         # prefix name for external AF call set (required if ref_bed set)
    Array[String]? population  # populations to annotate external AF for (required if ref_bed set)

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
    RuntimeAttr? runtime_attr_modify_vcf
    RuntimeAttr? runtime_override_combine_vcfs
    RuntimeAttr? runtime_override_split_vcf
  }

  call ann.ScatterAnnotateVcfByChrom as ScatterAnnotateVcfByChrom {
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
      vcf                    = ScatterAnnotateVcfByChrom.annotated_vcf,
      vcf_idx                = ScatterAnnotateVcfByChrom.annotated_vcf_idx,
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

  if (defined(ref_bed)) {
    call eaf.AnnotateExternalAF as AnnotateExternalAF {
      input:
        vcf     = PruneAndAddVafs.output_vcf,
        vcf_idx = PruneAndAddVafs.output_vcf_idx,
        ref_bed = select_first([ref_bed]),
        population = select_first([population]),
        ref_prefix = select_first([ref_prefix]),
        prefix = prefix,
        contigs = read_lines(contig_list),
        max_shards_per_chrom_step1 = max_shards_per_chrom_step1,
        min_records_per_shard_step1 = min_records_per_shard_step1,
        sv_base_mini_docker = sv_base_mini_docker,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_modify_vcf = runtime_attr_modify_vcf,
        runtime_override_split_vcf = runtime_override_split_vcf,
        runtime_override_combine_vcfs = runtime_override_combine_vcfs
    }
  }

  output {
    File output_vcf     = select_first([AnnotateExternalAF.annotated_vcf, PruneAndAddVafs.output_vcf])
    File output_vcf_idx = select_first([AnnotateExternalAF.annotated_vcf_tbi, PruneAndAddVafs.output_vcf_idx])
  }
}
