version 1.0

import "Structs.wdl"
import "TasksMakeCohortVcf.wdl" as MiniTasks

workflow AnnotateBiallelicCNCounts {

  input {
    Array[File] vcfs
    String prefix

    String sv_pipeline_docker
    String sv_base_mini_docker

    RuntimeAttr? runtime_attr_annotate
    RuntimeAttr? runtime_attr_concat
  }

  scatter (vcf in vcfs) {
    call AnnotateCNCounts {
      input:
        vcf = vcf,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_annotate
    }
  }

  call MiniTasks.ConcatVcfs {
    input:
      vcfs = AnnotateCNCounts.annotated_vcf,
      naive = true,
      sites_only = true,
      outfile_prefix = "~{prefix}.cn_counts",
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_concat
  }

  output {
    File annotated_vcf = ConcatVcfs.concat_vcf
    File annotated_vcf_idx = ConcatVcfs.concat_vcf_idx
  }
}

task AnnotateCNCounts {
  input {
    File vcf
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  String output_file = basename(vcf, ".vcf.gz") + ".cn_counts.vcf.gz"

  Float input_size = size(vcf, "GB")
  RuntimeAttr runtime_default = object {
    mem_gb: 3.75,
    disk_gb: ceil(10.0 + input_size * 2.0),
    cpu_cores: 1,
    preemptible_tries: 3,
    max_retries: 1,
    boot_disk_gb: 10
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, runtime_default])

  runtime {
    memory: select_first([runtime_attr.mem_gb, runtime_default.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, runtime_default.disk_gb]) + " HDD"
    cpu: select_first([runtime_attr.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_attr.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euo pipefail

    python3 <<CODE
import pysam
from collections import defaultdict

with pysam.VariantFile("~{vcf}") as vcf_in:
    out_header = vcf_in.header
    out_header.add_line('##INFO=<ID=CN_COUNT,Number=.,Type=Integer,Description="Number of samples observed at each copy state, starting from CN=0 (multiallelic CNVs only)">')

    with pysam.VariantFile("~{output_file}", 'w', header=out_header) as vcf_out:
        for record in vcf_in:
            svtype = record.info.get("SVTYPE", None)
            if svtype in ("DEL", "DUP"):
                cn_counter = defaultdict(int)
                for sample in record.samples:
                    rd_cn = record.samples[sample].get("RD_CN", None)
                    if rd_cn is not None:
                        cn_counter[rd_cn] += 1
                if cn_counter:
                    max_cn = max(cn_counter.keys())
                    cn_count = tuple(cn_counter.get(i, 0) for i in range(max_cn + 1))
                    record.info["CN_COUNT"] = cn_count
            vcf_out.write(record)

CODE

    tabix -p vcf "~{output_file}"
  >>>

  output {
    File annotated_vcf = output_file
    File annotated_vcf_idx = output_file + ".tbi"
  }
}
