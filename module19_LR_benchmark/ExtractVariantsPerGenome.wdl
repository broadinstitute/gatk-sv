version 1.0
import "Structs.wdl"


workflow ExtractVariantsPerGenome {
  input {
    File input_vcf
    File input_vcf_idx
    Array[String] sample_list
    String sv_base_pipeline_docker
    RuntimeAttr? runtime_attr_extract_sample

  }

  # Scatter across each sample in the list

  scatter (sample in sample_list) {
    call ExtractSampleNonRef {
      input:
        vcf = input_vcf,
        vcf_idx = input_vcf_idx,
        sample = sample,
        docker_image   = sv_base_pipeline_docker,
        runtime_attr_override = runtime_attr_extract_sample
    }
  }

  output {
    Array[File] extracted_vcfs = ExtractSampleNonRef.output_vcf
  }
}

task ExtractSampleNonRef {
  input {
    File vcf
    File vcf_idx
    String sample_name
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  command <<<
    set -euo pipefail

    bcftools view \
      --samples ~{sample_name} \
      --min-ac 1 \
      -O z \
      -o ~{sample_name}.nonref.vcf.gz \
      ~{vcf}

    tabix -p vcf ~{sample_name}.nonref.vcf.gz
  >>>

  output {
    File output_vcf = "~{sample_name}.nonref.vcf.gz"
    File output_vcf_idx = "~{sample_name}.nonref.vcf.gz.tbi"
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 10 + ceil(size(vcf, "GiB")*1.5),
    disk_gb: 15 + ceil(size(vcf, "GiB")*1.5),
    boot_disk_gb: 10,
    preemptible_tries: 1,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: docker_image
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}
