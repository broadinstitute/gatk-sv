version 1.0

import "Structs.wdl"
import "ShardedAnnotateVcf.wdl" as sharded_annotate_vcf

workflow AnnotateVcf {

  input {
    Array[File] vcf_list
    Array[File] vcf_idx_list
    File contig_list
    Array[String] prefix_list

    File protein_coding_gtf
    File? noncoding_bed
    Int? promoter_window
    Int? max_breakend_as_cnv_length
    String? svannotate_additional_args

    Int max_shards_per_chrom_step1
    Int min_records_per_shard_step1

    File? sample_pop_assignments  # Two-column file with sample ID & pop assignment. "." for pop will ignore sample
    File? sample_keep_list              # List of samples to be retained from the output vcf
    File? ped_file                # Used for M/F AF calculations
    File? par_bed
    File? allosomes_list
    Int   sv_per_shard

    File? ref_bed              # File with external allele frequencies
    String? ref_prefix         # prefix name for external AF call set (required if ref_bed set)
    Array[String]? population  # populations to annotate external AF for (required if ref_bed set)

    Boolean use_hail
    String? gcs_project

    String sv_pipeline_docker
    String? sv_pipeline_hail_docker
    String sv_base_mini_docker
    String gatk_docker

    RuntimeAttr? runtime_attr_svannotate
    RuntimeAttr? runtime_attr_concat_vcfs
    RuntimeAttr? runtime_attr_shard_vcf
    RuntimeAttr? runtime_attr_compute_AFs
    RuntimeAttr? runtime_attr_combine_vcfs
    RuntimeAttr? runtime_attr_modify_vcf
    RuntimeAttr? runtime_attr_combine_vcfs
    RuntimeAttr? runtime_attr_split_vcf
    RuntimeAttr? runtime_attr_split_ref_bed
    RuntimeAttr? runtime_attr_split_query_vcf
    RuntimeAttr? runtime_attr_bedtools_closest
    RuntimeAttr? runtime_attr_select_matched_svs
    RuntimeAttr? runtime_attr_concat_sharded_cluster
    RuntimeAttr? runtime_attr_preconcat_sharded_cluster
    RuntimeAttr? runtime_attr_hail_merge_sharded_cluster
    RuntimeAttr? runtime_attr_fix_header_sharded_cluster
    RuntimeAttr? runtime_attr_get_vcf_header_with_members_info_line
  }

  Array[String] contigs = read_lines(contig_list)

  scatter (i in range(length(vcf_list))) {
    call sharded_annotate_vcf.ShardedAnnotateVcf as ShardedAnnotateVcf{
      input:
        vcf = vcf_list[i],
        vcf_idx = vcf_idx_list[i],
        contig = contigs[i],
        prefix = prefix_list[i],
        protein_coding_gtf = protein_coding_gtf,
        noncoding_bed = noncoding_bed,
        promoter_window = promoter_window,
        svannotate_additional_args = svannotate_additional_args,
        max_breakend_as_cnv_length = max_breakend_as_cnv_length,

        max_shards_per_chrom_step1 = max_shards_per_chrom_step1,
        min_records_per_shard_step1 = min_records_per_shard_step1,
        sample_pop_assignments = sample_pop_assignments,
        sample_keep_list = sample_keep_list,
        ped_file = ped_file,
        par_bed = par_bed,
        sv_per_shard = sv_per_shard,
        allosomes_list = allosomes_list,

        ref_bed = ref_bed,
        ref_prefix = ref_prefix,
        population = population,

        use_hail = use_hail,
        gcs_project = gcs_project,

        gatk_docker = gatk_docker,
        sv_pipeline_docker = sv_pipeline_docker,
        sv_base_mini_docker = sv_base_mini_docker,
        sv_pipeline_hail_docker = sv_pipeline_hail_docker,

        runtime_attr_svannotate = runtime_attr_svannotate ,
        runtime_attr_concat_vcfs  = runtime_attr_concat_vcfs  ,
        runtime_attr_shard_vcf  = runtime_attr_shard_vcf  ,
        runtime_attr_compute_AFs  = runtime_attr_compute_AFs  ,
        runtime_attr_combine_vcfs = runtime_attr_combine_vcfs ,
        runtime_attr_modify_vcf = runtime_attr_modify_vcf ,
        runtime_attr_combine_vcfs = runtime_attr_combine_vcfs ,
        runtime_attr_split_vcf  = runtime_attr_split_vcf  ,
        runtime_attr_split_ref_bed  = runtime_attr_split_ref_bed  ,
        runtime_attr_split_query_vcf  = runtime_attr_split_query_vcf  ,
        runtime_attr_bedtools_closest = runtime_attr_bedtools_closest ,
        runtime_attr_select_matched_svs = runtime_attr_select_matched_svs ,
        runtime_attr_concat_sharded_cluster = runtime_attr_concat_sharded_cluster ,
        runtime_attr_preconcat_sharded_cluster  = runtime_attr_preconcat_sharded_cluster  ,
        runtime_attr_hail_merge_sharded_cluster = runtime_attr_hail_merge_sharded_cluster ,
        runtime_attr_fix_header_sharded_cluster = runtime_attr_fix_header_sharded_cluster ,
        runtime_attr_get_vcf_header_with_members_info_line  = runtime_attr_get_vcf_header_with_members_info_line
    }
  }

  output {
    Array[File] output_vcf_list     = ShardedAnnotateVcf.output_vcf
    Array[File] output_vcf_idx_list = ShardedAnnotateVcf.output_vcf_idx
  }
}