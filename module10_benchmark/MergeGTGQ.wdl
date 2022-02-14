version 1.0

import "TasksBenchmark.wdl" as mini_tasks


workflow MergeGTGQ {

  input {
    Array[Array[File]] shards_chr_sample_gtgq
    Array[Array[File]] shards_chr_sample_beds
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_concat_gtgqs
    RuntimeAttr? runtime_attr_concat_refs
  }

  Array[Array[File]] gtgq_shards_chr_sample = transpose(shards_chr_sample_gtgq)
  Array[Array[File]] beds_shards_chr_sample = transpose(shards_chr_sample_beds)

  call mini_tasks.ConcatGTGQs as ConcatGTGQs {
      input:
        files = gtgq_shards_chr_sample,
        sv_base_mini_docker = sv_base_mini_docker,
        runtime_attr_override = runtime_attr_concat_gtgqs
  }

  call mini_tasks.ConcatREFs as ConcatREFs {
      input:
        files = beds_shards_chr_sample,
        sv_base_mini_docker = sv_base_mini_docker,
        runtime_attr_override = runtime_attr_concat_refs

  }


  output {
    Array[File] gtgq = ConcatGTGQs.merged_gtgq_files
    Array[File] refs = ConcatREFs.merged_ref_files
  }
  
}


