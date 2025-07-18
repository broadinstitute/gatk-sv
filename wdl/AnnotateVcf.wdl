version 1.0

import "Structs.wdl"
import "ShardedAnnotateVcf.wdl" as sharded_annotate_vcf
import "TasksMakeCohortVcf.wdl" as MiniTasks

workflow AnnotateVcf {

  input {
    Array[File] vcfs  # GATK-SV VCF for annotation. Index .tbi must be located at the same path
    File contig_list  # Ordered list of contigs to annotate that are present in the input VCF
    String prefix

    File protein_coding_gtf
    File? noncoding_bed
    Int? promoter_window
    Int? max_breakend_as_cnv_length
    String? svannotate_additional_args

    File? sample_pop_assignments  # Two-column file with sample ID & pop assignment. "." for pop will ignore sample
    File? sample_keep_list              # List of samples to be retained from the output vcf
    File? ped_file                # Used for M/F AF calculations
    File? par_bed
    File? allosomes_list
    Int   sv_per_shard

    File? ref_bed              # File with external allele frequencies
    String? ref_prefix         # prefix name for external AF call set (required if ref_bed set)
    Array[String]? population  # populations to annotate external AF for (required if ref_bed set)

    String sv_pipeline_docker
    String sv_base_mini_docker
    String gatk_docker

    RuntimeAttr? runtime_attr_svannotate
    RuntimeAttr? runtime_attr_scatter_vcf
    RuntimeAttr? runtime_attr_subset_vcf_by_samples_list
    RuntimeAttr? runtime_attr_compute_AFs
    RuntimeAttr? runtime_attr_modify_vcf
    RuntimeAttr? runtime_attr_split_ref_bed
    RuntimeAttr? runtime_attr_split_query_vcf
    RuntimeAttr? runtime_attr_bedtools_closest
    RuntimeAttr? runtime_attr_select_matched_svs
    RuntimeAttr? runtime_attr_concat
  }

  Array[String] contigs = read_lines(contig_list)

  scatter (i in range(length(contigs))) {
    call sharded_annotate_vcf.ShardedAnnotateVcf {
      input:
        vcf = vcfs[i],
        vcf_idx = vcfs[i] + ".tbi",
        contig = contigs[i],
        prefix = prefix,
        protein_coding_gtf = protein_coding_gtf,
        noncoding_bed = noncoding_bed,
        promoter_window = promoter_window,
        svannotate_additional_args = svannotate_additional_args,
        max_breakend_as_cnv_length = max_breakend_as_cnv_length,

        sample_pop_assignments = sample_pop_assignments,
        sample_keep_list = sample_keep_list,
        ped_file = ped_file,
        par_bed = par_bed,
        sv_per_shard = sv_per_shard,
        allosomes_list = allosomes_list,

        ref_bed = ref_bed,
        ref_prefix = ref_prefix,
        population = population,

        gatk_docker = gatk_docker,
        sv_pipeline_docker = sv_pipeline_docker,
        sv_base_mini_docker = sv_base_mini_docker,

        runtime_attr_svannotate = runtime_attr_svannotate,
        runtime_attr_scatter_vcf = runtime_attr_scatter_vcf,
        runtime_attr_subset_vcf_by_samples_list = runtime_attr_subset_vcf_by_samples_list,
        runtime_attr_compute_AFs  = runtime_attr_compute_AFs,
        runtime_attr_modify_vcf = runtime_attr_modify_vcf,
        runtime_attr_split_ref_bed  = runtime_attr_split_ref_bed,
        runtime_attr_split_query_vcf  = runtime_attr_split_query_vcf,
        runtime_attr_bedtools_closest = runtime_attr_bedtools_closest,
        runtime_attr_select_matched_svs = runtime_attr_select_matched_svs,
        runtime_attr_concat = runtime_attr_concat
    }
  }

  output {
    Array[File] annotated_vcfs = ShardedAnnotateVcf.annotated_vcf
    Array[File] annotated_vcf_indexes = ShardedAnnotateVcf.annotated_vcf_idx
  }
}
