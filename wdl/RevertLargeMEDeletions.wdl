version 1.0

import "RevertLargeMEDeletionsPerContig.wdl" as revert

workflow RevertLargeMEDeletions {
  input {
    Array[File] vcfs

    Int records_per_shard = 20000
    String sv_base_mini_docker
    String sv_pipeline_docker
  }

  scatter ( vcf in vcfs ) {
    call revert.RevertLargeMEDeletionsPerContig {
      input:
        vcf = vcf,
        records_per_shard = records_per_shard,
        sv_pipeline_docker = sv_pipeline_docker,
        sv_base_mini_docker = sv_base_mini_docker
    }
  }

  output {
    Array[File] vcfs_large_me_dels_reverted = RevertLargeMEDeletionsPerContig.vcf_large_me_dels_reverted
    Array[File] vcfs_large_me_dels_reverted_idxs = RevertLargeMEDeletionsPerContig.vcf_large_me_dels_reverted_index
  }
}
