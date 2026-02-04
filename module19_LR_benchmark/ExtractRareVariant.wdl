version 1.0

import "Structs.wdl"

workflow ExtractRareVariant {
  input {
    File vcf_gz
    File vcf_tbi
    String docker_image

    RuntimeAttr? runtime_attr_split
    RuntimeAttr? runtime_attr_biallelic
    RuntimeAttr? runtime_attr_filter
    RuntimeAttr? runtime_attr_concat
  }

  call DetectContigs {
    input:
      vcf_gz = vcf_gz,
      docker_image = docker_image
  }

  scatter (c in DetectContigs.contigs) {
    call SplitByContig {
      input:
        vcf_gz = vcf_gz,
        vcf_tbi = vcf_tbi,
        contig = c,
        docker_image = docker_image,
        runtime_attr_override = runtime_attr_split
    }

    call SplitToBiallelic {
      input:
        vcf_gz = SplitByContig.out_vcf,
        docker_image = docker_image,
        runtime_attr_override = runtime_attr_biallelic
    }

    call FilterByAF_001 {
      input:
        vcf_gz = SplitToBiallelic.out_vcf,
        af_threshold = 0.01,
        docker_image = docker_image,
        runtime_attr_override = runtime_attr_filter
    }

    call FilterByAF_0001 {
      input:
        vcf_gz = SplitToBiallelic.out_vcf,
        af_threshold = 0.001,
        docker_image = docker_image,
        runtime_attr_override = runtime_attr_filter
    }
  }

  call ConcatVCFs {
    input:
      vcfs = flatten(FilterByAF_001.out_vcf),
      output_name = "AF_lt_0.01.vcf.gz",
      docker_image = docker_image,
      runtime_attr_override = runtime_attr_concat
  }

  call ConcatVCFs as ConcatAF0001 {
    input:
      vcfs = flatten(FilterByAF_0001.out_vcf),
      output_name = "AF_lt_0.001.vcf.gz",
      docker_image = docker_image,
      runtime_attr_override = runtime_attr_concat
  }

  output {
    File af_lt_001_vcf  = ConcatVCFs.out_vcf
    File af_lt_0001_vcf = ConcatAF0001.out_vcf
  }
}

# =================================================
# Detect contigs
# =================================================
task DetectContigs {
  input {
    File vcf_gz
    String docker_image
  }

  command <<<
    bcftools index -n ~{vcf_gz} > contigs.txt
  >>>

  output {
    Array[String] contigs = read_lines("contigs.txt")
  }

  runtime {
    docker: docker_image
    cpu: 1
    memory: "2 GiB"
  }
}

# =================================================
# Split by contig
# =================================================
task SplitByContig {
  input {
    File vcf_gz
    File vcf_tbi
    String contig
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: ceil(2 + size(vcf_gz, "GB")),
    disk_gb: ceil(10 + size(vcf_gz, "GB") * 2),
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    bcftools view -r ~{contig} ~{vcf_gz} -Oz -o ~{contig}.vcf.gz
    bcftools index ~{contig}.vcf.gz
  >>>

  output {
    File out_vcf = "~{contig}.vcf.gz"
  }

  runtime {
    docker: docker_image
    cpu: runtime_attr.cpu_cores
    memory: runtime_attr.mem_gb + " GiB"
    disks: "local-disk " + runtime_attr.disk_gb + " HDD"
  }
}

# =================================================
# Split multiallelic to biallelic
# =================================================
task SplitToBiallelic {
  input {
    File vcf_gz
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: ceil(3 + size(vcf_gz, "GB")),
    disk_gb: ceil(10 + size(vcf_gz, "GB") * 2),
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    bcftools norm -m -any ~{vcf_gz} -Oz -o biallelic.vcf.gz
    bcftools index biallelic.vcf.gz
  >>>

  output {
    File out_vcf = "biallelic.vcf.gz"
  }

  runtime {
    docker: docker_image
    cpu: runtime_attr.cpu_cores
    memory: runtime_attr.mem_gb + " GiB"
    disks: "local-disk " + runtime_attr.disk_gb + " HDD"
  }
}

# =================================================
# Filter by AF
# =================================================
task FilterByAF_001 {
  input {
    File vcf_gz
    Float af_threshold
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: ceil(2 + size(vcf_gz, "GB")),
    disk_gb: ceil(10 + size(vcf_gz, "GB") * 2),
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    bcftools view -i "INFO/AF < ~{af_threshold}" ~{vcf_gz} -Oz -o filtered.vcf.gz
    bcftools index filtered.vcf.gz
  >>>

  output {
    File out_vcf = "filtered.vcf.gz"
  }

  runtime {
    docker: docker_image
    cpu: runtime_attr.cpu_cores
    memory: runtime_attr.mem_gb + " GiB"
    disks: "local-disk " + runtime_attr.disk_gb + " HDD"
  }
}

# =================================================
# Concat VCFs
# =================================================
task ConcatVCFs {
  input {
    Array[File] vcfs
    String output_name
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 4,
    disk_gb: 20,
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    bcftools concat -a -Oz ~{sep=' ' vcfs} -o ~{output_name}
    bcftools index ~{output_name}
  >>>

  output {
    File out_vcf = "~{output_name}"
  }

  runtime {
    docker: docker_image
    cpu: runtime_attr.cpu_cores
    memory: runtime_attr.mem_gb + " GiB"
  }
}