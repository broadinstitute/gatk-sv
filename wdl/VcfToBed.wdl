version 1.0

import "CollectQcVcfWide.wdl" as collect
import "TasksMakeCohortVcf.wdl" as MiniTasks

workflow VcfToBed {
  input {
    File vcf
    Int records_per_shard = 20000
    String prefix
    String? flags

    String sv_pipeline_docker
    String sv_base_mini_docker

    RuntimeAttr? runtime_attr_svtk_vcf2bed
  }

  call MiniTasks.ScatterVcf {
    input:
      vcf = vcf,
      vcf_index = vcf + ".tbi",
      records_per_shard = records_per_shard,
      prefix = prefix,
      sv_pipeline_docker = sv_pipeline_docker
  }

  scatter (i in range(length(ScatterVcf.shards))) {
    call collect.SvtkVcf2bed {
      input:
        vcf=ScatterVcf.shards[i],
        flags=flags,
        prefix="~{prefix}.shard_~{i}",
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_attr_svtk_vcf2bed
    }
  }

  call MiniTasks.ConcatBeds {
    input:
      shard_bed_files = SvtkVcf2bed.vcf2bed_subworkflow_out,
      prefix=prefix,
      index_output=false,
      sv_base_mini_docker=sv_base_mini_docker
  }

  output {
    File bed = ConcatBeds.merged_bed_file
  }
}
