version 1.0

import "AnnotateFunctionalConsequences.wdl" as func
import "PruneAndAddVafs.wdl" as pav
import "AnnotateExternalAF.wdl" as eaf

workflow AnnotateVcf {

  input {
    File vcf
    File vcf_idx
    File contig_list
    String prefix

    File protein_coding_gtf
    File? noncoding_bed
    Int? promoter_window
    Int? max_breakend_as_cnv_length
    String? svannotate_additional_args

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
    String gatk_docker

    RuntimeAttr? runtime_attr_svannotate
    RuntimeAttr? runtime_attr_concat_vcfs
    RuntimeAttr? runtime_attr_prune_vcf
    RuntimeAttr? runtime_attr_shard_vcf
    RuntimeAttr? runtime_attr_compute_AFs
    RuntimeAttr? runtime_attr_combine_vcfs
    RuntimeAttr? runtime_attr_modify_vcf
    RuntimeAttr? runtime_override_combine_vcfs
    RuntimeAttr? runtime_override_split_vcf
    RuntimeAttr? runtime_attr_split_ref_bed
    RuntimeAttr? runtime_attr_split_query_vcf
    RuntimeAttr? runtime_attr_bedtools_closest
    RuntimeAttr? runtime_attr_select_matched_svs
  }

  call func.AnnotateFunctionalConsequences {
    input:
      vcf = vcf,
      vcf_index = vcf_idx,
      prefix = prefix,
      protein_coding_gtf = protein_coding_gtf,
      noncoding_bed = noncoding_bed,
      promoter_window = promoter_window,
      max_breakend_as_cnv_length = max_breakend_as_cnv_length,
      additional_args = svannotate_additional_args,
      gatk_docker = gatk_docker,
      runtime_attr_svannotate = runtime_attr_svannotate
  }

  call pav.PruneAndAddVafs as PruneAndAddVafs {
    input:
      vcf                    = AnnotateFunctionalConsequences.annotated_vcf,
      vcf_idx                = AnnotateFunctionalConsequences.annotated_vcf_index,
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
        runtime_override_combine_vcfs = runtime_override_combine_vcfs,
        runtime_attr_split_ref_bed = runtime_attr_split_ref_bed,
        runtime_attr_split_query_vcf = runtime_attr_split_query_vcf,
        runtime_attr_bedtools_closest = runtime_attr_bedtools_closest,
        runtime_attr_select_matched_svs = runtime_attr_select_matched_svs
    }
  }

  output {
    File output_vcf     = select_first([AnnotateExternalAF.annotated_vcf, PruneAndAddVafs.output_vcf])
    File output_vcf_idx = select_first([AnnotateExternalAF.annotated_vcf_tbi, PruneAndAddVafs.output_vcf_idx])
  }
}
