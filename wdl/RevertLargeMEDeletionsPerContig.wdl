version 1.0

import "Structs.wdl"
import "TasksMakeCohortVcf.wdl" as tasks

workflow RevertLargeMEDeletionsPerContig {
  input {
    File vcf
    Int records_per_shard
    String sv_pipeline_docker
    String sv_base_mini_docker
  }

  String prefix = basename(vcf, ".vcf.gz")

  call tasks.ScatterVcf {
    input:
      vcf = vcf,
      records_per_shard = records_per_shard,
      prefix = "~{prefix}.scatter",
      sv_pipeline_docker = sv_pipeline_docker
  }

  scatter (i in range(length(ScatterVcf.shards))) {
    call RevertLargeMEDels {
      input:
        vcf = ScatterVcf.shards[i],
        prefix = "~{prefix}.shard_~{i}.ME_DELs_gt10kb_reverted",
        sv_pipeline_docker = sv_pipeline_docker
    }
  }

  call tasks.ConcatVcfs {
    input:
      vcfs=RevertLargeMEDels.out,
      naive=true,
      outfile_prefix=prefix,
      sv_base_mini_docker=sv_base_mini_docker
  }

  output {
    File vcf_large_me_dels_reverted = ConcatVcfs.concat_vcf
    File vcf_large_me_dels_reverted_index = ConcatVcfs.concat_vcf_idx
  }
}


task RevertLargeMEDels {
  input {
    File vcf
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 3.75,
                               disk_gb: ceil(100 + size(vcf, "GB") * 2),
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File out = "~{prefix}.vcf.gz"
    File out_index = "~{prefix}.vcf.gz.tbi"
  }

  command <<<
  set -euo pipefail

  python <<CODE
import pysam
vcf = pysam.VariantFile("~{vcf}", 'r')
out = pysam.VariantFile("~{prefix}.vcf.gz", 'w', header=vcf.header)
for record in vcf:
  if "DEL:ME" in record.alts[0] and record.info['SVLEN'] > 10000:
    record.alts = ("<BND>",)
    record.info['SVTYPE'] = "BND"
    record.filter.clear()
    record.filter.add("PASS")
    record.info['CHR2'] = record.chrom
    record.info['END2'] = record.pos + record.info['SVLEN']
    record.info.pop('SVLEN')
    record.stop = record.pos + 1
  out.write(record)
vcf.close()
out.close()
CODE

  tabix ~{prefix}.vcf.gz
  >>>

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

