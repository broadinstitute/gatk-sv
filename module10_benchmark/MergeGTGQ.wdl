version 1.0

import "TasksBenchmark.wdl" as mini_tasks


workflow MergeGTGQ {

  input {
    Array[Array[File]] shards_chr_sample
    String sv_base_mini_docker
  }

  Array[Array[File]] shards_sample_chr = transpose(shards_chr_sample)

  scatter (i in range(length(shards_sample_chr))) {
    call mini_tasks.ConcatGTGQs as ConcatGTGQs {
      input:
        files = shards_sample_chr[i],
        sv_base_mini_docker = sv_base_mini_docker
    }
  }

  output {
    Array[File] gtgq = ConcatGTGQs.merged_gtgq_file
  }
}


