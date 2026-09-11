version 1.0

import "CollectQcVcfWide.wdl" as collect
import "TasksMakeCohortVcf.wdl" as MiniTasks

workflow VcfToBed {
  input {
    File vcf
    Int records_per_shard = 20000
    String prefix
    String? flags
    Boolean pass_multiallelic_only = false

    String sv_pipeline_docker
    String sv_base_mini_docker

    RuntimeAttr? runtime_attr_svtk_vcf2bed
    RuntimeAttr? runtime_attr_concat_bed
    RuntimeAttr? runtime_attr_scatter_vcf
    RuntimeAttr? runtime_attr_filter_bed
  }

  call MiniTasks.ScatterVcf {
    input:
      vcf = vcf,
      vcf_index = vcf + ".tbi",
      records_per_shard = records_per_shard,
      prefix = prefix,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override=runtime_attr_scatter_vcf
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
      index_output=true,
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_attr_concat_bed
  }

  if (pass_multiallelic_only) {
    call FilterBed {
      input:
        bed=ConcatBeds.merged_bed_file,
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_attr_override=runtime_attr_filter_bed
    }
  }

  output {
    File bed = select_first([FilterBed.out, ConcatBeds.merged_bed_file])
    File bed_index = select_first([FilterBed.out_index, ConcatBeds.merged_bed_idx])
  }
}


task FilterBed {
  input {
    File bed

    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  String outfile_name = basename(bed, ".bed.gz") + ".filtered.bed.gz"

  RuntimeAttr runtime_default = object {
                                  mem_gb: 3.75,
                                  disk_gb: ceil(10.0 +  size(bed, "GB") * 3),
                                  cpu_cores: 1,
                                  preemptible_tries: 3,
                                  max_retries: 1,
                                  boot_disk_gb: 10
                                }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

  Float runtime_mem_gb = select_first([runtime_override.mem_gb, runtime_default.mem_gb])

  command <<<
    set -euo pipefail
    zcat ~{bed} | awk -F'\t' '
    NR==1 && /^#/ {
        for(i=1; i<=NF; i++) {
            if($i == "FILTER") col=i
        }
        print $0
        next
    }
    col && ($col == "PASS" || $col == "MULTIALLELIC") {
        print $0
    }' | bgzip -c > ~{outfile_name}
    tabix -p bed ~{outfile_name}
  >>>

  output {
    File out = outfile_name
    File out_index = outfile_name + ".tbi"
  }
  runtime {
    memory: runtime_mem_gb + " GiB"
    disks: "local-disk " + select_first([runtime_override.disk_gb, runtime_default.disk_gb]) + " HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_base_mini_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }
}
