version 1.0

import "TasksMakeCohortVcf.wdl" as tasks

workflow ConcatVcfFiles {
  input {
    Array[File] vcfs
    Array[File]? vcfs_idx
    String outfile_prefix
    String sv_base_mini_docker
  }

  call tasks.ConcatVcfs {
    input:
      vcfs = vcfs,
      vcfs_idx = vcfs_idx,
      outfile_prefix = outfile_prefix,
      sv_base_mini_docker = sv_base_mini_docker
  }

  output {
    File concat_vcf = ConcatVcfs.concat_vcf
    File concat_vcf_index = ConcatVcfs.concat_vcf_idx
  }
}
