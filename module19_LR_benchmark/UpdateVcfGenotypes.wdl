version 1.0

import "Structs.wdl"

workflow UpdateVcfGenotypes {
  input {
    File vcf_A
    File vcf_A_idx
    File vcf_B
    File vcf_B_idx
    File patch_script          # patch_vcf.py
    String output_prefix
    Array[String] contigs
    String docker_image

    RuntimeAttr? runtime_attr_concat_vcfs
    RuntimeAttr? runtime_attr_patch_genotypes
    RuntimeAttr? runtime_attr_extract_contig_vcf
    RuntimeAttr? runtime_attr_extract_variant_sites

  }

  scatter (ctg in contigs) {
    call ExtractContigVCF {
      input:
        vcf = vcf_A,
        vcf_idx = vcf_A_idx,
        contig = ctg,
        docker_image = docker_image,
        runtime_attr_override = runtime_attr_extract_contig_vcf
    }

    call ExtractVariantSites {
      input:
        vcf = ExtractContigVCF.out_vcf,
        docker_image = docker_image,
        runtime_attr_override = runtime_attr_extract_variant_sites
    }

    call PatchGenotypes {
      input:
        contig = ctg,
        vcf_B = vcf_B,
        vcf_B_idx = vcf_B_idx,
        vcf_A_contig = ExtractVariantSites.out_sites,
        patch_script = patch_script,
        docker_image = docker_image,
        runtime_attr_override = runtime_attr_patch_genotypes
    }
  }

  call ConcatPatchedVCFs{
    input:
      vcfs = PatchGenotypes.out_vcf,
      output_prefix = output_prefix, 
      docker_image = docker_image,
      runtime_attr_override = runtime_attr_concat_vcfs
  }

  output {
    File updated_vcf = ConcatPatchedVCFs.out_vcf
    File updated_vcf_idx = ConcatPatchedVCFs.out_vcf_idx
  }
}

task ExtractContigVCF {
  input {
    File vcf
    File vcf_idx
    String contig
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  command <<<
    set -euo pipefail

    bcftools view -r ~{contig} -Oz -o ~{contig}.A.vcf.gz ~{vcf}
    tabix -p vcf ~{contig}.A.vcf.gz
  >>>

  output {
    File out_vcf = "~{contig}.A.vcf.gz"
    File out_vcf_idx = "~{contig}.A.vcf.gz.tbi"
  }

  RuntimeAttr default_attr = object {
      cpu_cores: 1,
      mem_gb: 10,
      disk_gb: 15 + ceil(size(vcf, "GiB") *3),
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

task ExtractVariantSites {
  input {
    File vcf
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  String prefix = basename(vcf, ".vcf.gz")
  command <<<
    set -euo pipefail

    bcftools view ~{vcf} | cut -f1-8 | bgzip > ~{prefix}.sites.gz
    tabix -p vcf ~{prefix}.sites.gz
  >>>

  output {
    File out_sites = "~{prefix}.sites.gz"
    File out_sites_idx = "~{prefix}.sites.gz.tbi"
  }

  RuntimeAttr default_attr = object {
      cpu_cores: 1,
      mem_gb: 10,
      disk_gb: 15 + ceil(size(vcf, "GiB") *3),
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


task PatchGenotypes {
  input {
    File vcf_A_contig
    File vcf_B
    File vcf_B_idx
    File patch_script
    String contig
    String docker_image
    RuntimeAttr? runtime_attr_override
  }


  String prefix = basename(vcf_B, ".vcf.gz")
  command <<<
    set -euo pipefail

    bcftools view ~{vcf_B} ~{contig} |bgzip > "~{prefix}.~{contig}.target.vcf.gz"
    bcftools sort "~{prefix}.~{contig}.target.vcf.gz" -Oz -o "~{prefix}.~{contig}.sorted.vcf.gz"
    tabix -p vcf "~{prefix}.~{contig}.sorted.vcf.gz"

    python3 ~{patch_script} \
      ~{vcf_A_contig} \
      "~{prefix}.~{contig}.sorted.vcf.gz" \
      ~{contig}.patched.vcf.gz

    tabix -p vcf ~{contig}.patched.vcf.gz
  >>>

  output {
    File out_vcf = "~{contig}.patched.vcf.gz"
    File out_vcf_idx = "~{contig}.patched.vcf.gz.tbi"
  }

  RuntimeAttr default_attr = object {
      cpu_cores: 1,
      mem_gb: 10,
      disk_gb: 15 + ceil(size(vcf_B, "GiB") *2),
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

task ConcatPatchedVCFs {
  input {
    Array[File] vcfs
    String output_prefix
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  command <<<
    set -euo pipefail

    # Ensure deterministic order
    printf "%s\n" ~{sep='\n' vcfs} > vcf.list

    bcftools concat \
      -f vcf.list \
      -Oz \
      -o ~{output_prefix}.vcf.gz

    tabix -p vcf ~{output_prefix}.vcf.gz
  >>>

  output {
    File out_vcf = "~{output_prefix}.vcf.gz"
    File out_vcf_idx = "~{output_prefix}.vcf.gz.tbi"
  }

  RuntimeAttr default_attr = object {
      cpu_cores: 1,
      mem_gb: 10,
      disk_gb: 15 + ceil(size(vcfs, "GiB") *3),
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
