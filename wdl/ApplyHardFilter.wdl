version 1.0

import "Structs.wdl"

workflow ApplyHardFilter {
  input {
    String prefix
    File vcf
    String? filter_name

    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_hard_filter_wham
  }

  call HardFilterWham {
    input:
      prefix = prefix,
      vcf = vcf,
      filter_name = filter_name,
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_hard_filter_wham
  }

  output {
    File hard_filtered_vcf = HardFilterWham.hard_filtered_vcf
    File hard_filtered_vcf_index = HardFilterWham.hard_filtered_vcf_index
  }
}


# Removes Wham-only, SR-only deletions <1kb
task HardFilterWham {
  input {
    String prefix
    File vcf
    String? filter_name

    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_hard_filter_wham
  }

  String filter_name_ = select_first([filter_name, "filter_wham_del_sr_1kb"])
  String hard_filtered_vcf_name = "~{prefix}.~{filter_name_}.vcf.gz"

  # Disk must be scaled proportionally to the size of the VCF
  Float input_size = size(vcf, "GiB")
  RuntimeAttr default_attr = object {
    mem_gb: 3.75,
    disk_gb: ceil(10.0 + (input_size * 1.5)),
    cpu_cores: 1,
    preemptible_tries: 3,
    max_retries: 1,
    boot_disk_gb: 10
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<

    set -euo pipefail

    bcftools view -e 'SVTYPE=="DEL" && COUNT(ALGORITHMS)==1 && ALGORITHMS=="wham" && SVLEN<1000 && COUNT(EVIDENCE)==1 && EVIDENCE=="SR"' ~{vcf} -Oz -o "~{hard_filtered_vcf_name}"

    tabix -f -p vcf "~{hard_filtered_vcf_name}"

  >>>

  output {
    File hard_filtered_vcf = "~{hard_filtered_vcf_name}"
    File hard_filtered_vcf_index = "~{hard_filtered_vcf_name}.tbi"
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
