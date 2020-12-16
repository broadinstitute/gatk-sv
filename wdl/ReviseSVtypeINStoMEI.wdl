version 1.0

import "Structs.wdl"
import "TasksMakeCohortVcf.wdl" as MiniTasks
import "ReviseSVtypeINStoMEIperContig.wdl" as ReviseSVtypePerContig

workflow ReviseSVtypeINStoMEI {
  input{
    File vcf
    File vcf_idx
    String prefix
    File contiglist

    Int max_shards_per_chrom_step1
    Int min_records_per_shard_step1

    String sv_base_mini_docker
    String sv_pipeline_docker

    RuntimeAttr? runtime_override_split_vcf_to_clean
    RuntimeAttr? runtime_attr_ReviseSVtypeMEI
    RuntimeAttr? runtime_override_combine_step_1_vcfs
  }

  Array[String] contigs = transpose(read_tsv(contiglist))[0]
  scatter ( contig in contigs ) {
    call ReviseSVtypePerContig.ReviseSVtypeINStoMEIperContig as ReviseSVtypeINStoMEIperContig{
      input:
        vcf = vcf,
        vcf_idx = vcf_idx,
        prefix = prefix,
        contig = contig,
        max_shards_per_chrom_step1 = max_shards_per_chrom_step1,
        min_records_per_shard_step1 = min_records_per_shard_step1,
        sv_base_mini_docker = sv_base_mini_docker,
        runtime_override_split_vcf_to_clean = runtime_override_split_vcf_to_clean,
        runtime_attr_ReviseSVtypeMEI = runtime_attr_ReviseSVtypeMEI,
        runtime_override_combine_step_1_vcfs = runtime_override_combine_step_1_vcfs
    }
  }

  call MiniTasks.ConcatVcfs as CombineStep2Vcfs {
      input:
        vcfs = ReviseSVtypeINStoMEIperContig.updated_vcf,
        vcfs_idx = ReviseSVtypeINStoMEIperContig.updated_vcf_idx,
        naive = true,
        outfile_prefix = "~{prefix}.SVtypeRevisedINStoMEI",
        sv_base_mini_docker = sv_base_mini_docker,
        runtime_attr_override = runtime_override_combine_step_1_vcfs
  }

  output{
      File updated_vcf = CombineStep2Vcfs.concat_vcf
      File updated_vcf_idx = CombineStep2Vcfs.concat_vcf_idx
  }
}




