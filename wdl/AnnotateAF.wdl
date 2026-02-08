version 1.0

import "Structs.wdl"
import "ShardedAnnotateVcf.wdl" as sharded_annotate_vcf
import "TasksMakeCohortVcf.wdl" as MiniTasks

workflow AnnotateAF {
    input {
        File vcf
        File vcf_index
        Array[String] contigs
        String prefix

        File sample_pop_assignments
        File ped_file
        File par_bed
        File lps_tsv
        Array[String]? strip_info_fields

        Int records_per_shard
        String annotate_af_docker
        String sv_base_mini_docker
        String gatk_docker

        RuntimeAttr? runtime_attr_split_ref_bed
        RuntimeAttr? runtime_attr_scatter_vcf
        RuntimeAttr? runtime_attr_strip_info_fields
        RuntimeAttr? runtime_attr_compute_AFs
        RuntimeAttr? runtime_attr_concat
    }

    scatter (contig in contigs) {
        call sharded_annotate_vcf.ShardedAnnotateVcf {
            input:
                vcf = vcf,
                vcf_idx = vcf_index,
                contig = contig,
                prefix = "~{prefix}.~{contig}",
                sample_pop_assignments = sample_pop_assignments,
                ped_file = ped_file,
                par_bed = par_bed,
                lps_tsv = lps_tsv,
                strip_info_fields = strip_info_fields,
                records_per_shard = records_per_shard,
                gatk_docker = gatk_docker,
                sv_pipeline_docker = annotate_af_docker,
                sv_base_mini_docker = sv_base_mini_docker,
                runtime_attr_scatter_vcf = runtime_attr_scatter_vcf,
                runtime_attr_strip_info_fields = runtime_attr_strip_info_fields,
                runtime_attr_compute_AFs = runtime_attr_compute_AFs
        }
    }

    call MiniTasks.ConcatVcfs {
        input:
            vcfs = flatten(ShardedAnnotateVcf.sharded_annotated_vcf),
            vcfs_idx = flatten(ShardedAnnotateVcf.sharded_annotated_vcf_idx),
            allow_overlaps = true,
            outfile_prefix = "~{prefix}.annotated",
            sv_base_mini_docker = sv_base_mini_docker,
            runtime_attr_override = runtime_attr_concat
    }

    output {
        File af_annotated_vcf = ConcatVcfs.concat_vcf
        File af_annotated_vcf_idx = ConcatVcfs.concat_vcf_idx
    }
}
