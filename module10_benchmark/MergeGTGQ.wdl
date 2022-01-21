version 1.0

import "TasksBenchmark.wdl" as mini_tasks


workflow MergeGTGQ {

  input {
    Array[Array[File]] gtgq_shards_chr_sample
    Array[Array[File]] bed_shards_chr_sample
    String sv_base_mini_docker
  }

  Array[Array[File]] gtgq_shards_chr_sample = transpose(shards_chr_sample)

  scatter (i in range(length(gtgq_shards_chr_sample))) {
    call mini_tasks.ConcatGTGQs as ConcatGTGQs {
      input:
        files = gtgq_shards_chr_sample[i],
        sv_base_mini_docker = sv_base_mini_docker
    }
  }

  scatter (i in range(length(bed_shards_chr_sample))) {
    call mini_tasks.ConcatGTGQs as ConcatBEDs {
      input:
        files = bed_shards_chr_sample[i],
        sv_base_mini_docker = sv_base_mini_docker
    }
  }



  output {
    Array[File] gtgq = ConcatGTGQs.merged_gtgq_file
    Array[File] bed = ConcatBEDs.merged_gtgq_file
  }
}


