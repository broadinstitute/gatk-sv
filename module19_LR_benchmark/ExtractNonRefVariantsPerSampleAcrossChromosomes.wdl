version 1.0

import "Structs.wdl"

workflow ExtractNonRefVariantsPerSampleAcrossChromosomes {
  input {
    Array[File] input_vcfs
    Array[File]? input_vcfs_idx
    File sample_list
    String output_prefix
    String sv_pipeline_base_docker

    RuntimeAttr? runtime_attr_index_vcf
    RuntimeAttr? runtime_attr_extract_sample
    RuntimeAttr? runtime_attr_concat_vcfs
  }

  if (!defined(input_vcfs_idx)) {
    scatter (vcf in input_vcfs) {
      call IndexVcf {
        input:
          vcf = vcf,
          docker_image = sv_pipeline_base_docker,
          runtime_attr_override = runtime_attr_index_vcf
      }
    }
  }

  Array[File] vcf_indices = select_first([input_vcfs_idx, IndexVcf.indexed_vcf_idx])
  Array[String] sample_ids = read_lines(sample_list)

  scatter (sample_id in sample_ids) {
    scatter (idx in range(length(input_vcfs))) {
      call ExtractSampleNonRef {
        input:
          vcf = input_vcfs[idx],
          vcf_index = vcf_indices[idx],
          sample_id = sample_id,
          output_prefix = output_prefix,
          docker_image = sv_pipeline_base_docker,
          runtime_attr_override = runtime_attr_extract_sample
      }
    }

    call ConcatVcfs as ConcatPerSampleVcfs {
      input:
        input_vcfs = ExtractSampleNonRef.non_ref_vcf,
        input_vcfs_idx = ExtractSampleNonRef.non_ref_vcf_idx,
        output_name = output_prefix + "." + sample_id + ".nonref.vcf.gz",
        docker_image = sv_pipeline_base_docker,
        runtime_attr_override = runtime_attr_concat_vcfs
    }
  }

  output {
    Array[String] output_sample_ids = sample_ids
    Array[File] per_sample_non_ref_vcfs = ConcatPerSampleVcfs.output_vcf
    Array[File] per_sample_non_ref_vcf_indices = ConcatPerSampleVcfs.output_vcf_idx
  }
}

task IndexVcf {
  input {
    File vcf
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 2,
    disk_gb: ceil(size(vcf, "GB") * 2 + 10),
    boot_disk_gb: 10,
    preemptible_tries: 1,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  String vcf_name = basename(vcf)

  command <<<
    set -euo pipefail
    tabix -p vcf ~{vcf}
  >>>

  output {
    File indexed_vcf_idx = "~{vcf_name}.tbi"
  }

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

task ExtractSampleNonRef {
  input {
    File vcf
    File vcf_index
    String sample_id
    String output_prefix
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 4,
    disk_gb: ceil(size(vcf, "GB") * 2 + 10),
    boot_disk_gb: 10,
    preemptible_tries: 1,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  String vcf_base = basename(vcf, ".vcf.gz")
  String out_name = output_prefix + "." + sample_id + "." + vcf_base + ".nonref.vcf.gz"

  command <<<
    set -euo pipefail
    bcftools view \
      -s ~{sample_id} \
      -c 1 \
      -Oz \
      -o ~{out_name} \
      ~{vcf}
    tabix -p vcf ~{out_name}
  >>>

  output {
    File non_ref_vcf = out_name
    File non_ref_vcf_idx = out_name + ".tbi"
  }

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

task ConcatVcfs {
  input {
    Array[File] input_vcfs
    Array[File] input_vcfs_idx
    String output_name
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 8,
    disk_gb: ceil(size(input_vcfs, "GB") * 2 + 10),
    boot_disk_gb: 10,
    preemptible_tries: 1,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail
    bcftools concat ~{sep=' ' input_vcfs} -Oz -o ~{output_name}
    tabix -p vcf ~{output_name}
  >>>

  output {
    File output_vcf = output_name
    File output_vcf_idx = output_name + ".tbi"
  }

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