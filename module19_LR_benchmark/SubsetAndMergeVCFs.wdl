version 1.0

import "Structs.wdl"
import  "LongReadGenotypeTasks.wdl" as LongReadGenotypeTasks

workflow SubsetAndMergeVCFs {
  input {
    Array[File] input_vcfs      # List of .vcf.gz files
    File sample_list            # Text file with samples to include
    String output_prefix        # Prefix for final merged output
    String sv_base_mini_docker
  }

  scatter (vcf in input_vcfs) {
    call SubsetVCF {
      input:
        input_vcf = vcf,
        sample_list = sample_list,
        sv_base_mini_docker = sv_base_mini_docker
    }

    call RemoveSV{
      input:
        input_vcf = SubsetVCF.output_vcf,
        sv_base_mini_docker = sv_base_mini_docker
    }
  }

  call LongReadGenotypeTasks.ConcatVcfs {
    input:
      vcfs = RemoveSV.output_vcf,
      vcfs_idx = RemoveSV.output_vcf_idx,
      outfile_prefix = output_prefix,
      sv_base_mini_docker = sv_base_mini_docker
  }

  output {
    File merged_vcf = ConcatVcfs.concat_vcf
    File merged_vcf_idx = ConcatVcfs.concat_vcf_idx
  }
}

task SubsetVCF {
  input {
    File input_vcf
    File sample_list
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  command <<<
    set -euo pipefail

    bcftools view -S ~{sample_list} -c 1 ~{input_vcf} -Oz -o subset.vcf.gz
    tabix -p vcf subset.vcf.gz
  >>>

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 2,
    disk_gb: ceil(size(input_vcf, "GB") * 2),
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File output_vcf = "subset.vcf.gz"
    File output_vcf_idx = "subset.vcf.gz.tbi"
  }

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_base_mini_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

task RemoveSV {
  input {
    File input_vcf
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  command <<<
    set -euo pipefail

    zcat ~{input_vcf} | grep -v "SV" | bgzip > subset.non_SV.vcf.gz
    tabix -p vcf subset.non_SV.vcf.gz
  >>>

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 2,
    disk_gb: ceil(size(input_vcf, "GB") * 2),
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File output_vcf = "subset.non_SV.vcf.gz"
    File output_vcf_idx = "subset.non_SV.vcf.gz.tbi"
  }

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_base_mini_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

