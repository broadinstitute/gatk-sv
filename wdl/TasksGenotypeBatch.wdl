version 1.0

import "Structs.wdl"

task AddGenotypes {
  input {
    File vcf
    File genotypes
    File varGQ
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 3.75,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File genotyped_vcf = "~{prefix}.genotyped.vcf.gz"
    File genotyped_vcf_index = "~{prefix}.genotyped.vcf.gz.tbi"
  }
  command <<<

    set -euxo pipefail

    # in some cases a vargq cannot be computed and is returned as '.'. Remove these from the final vcf.
    gzip -cd ~{varGQ} | awk '$5 == "." {print $1}' > bad.vargq.list
    gzip -cd ~{vcf} | { grep -wvf bad.vargq.list || [[ $? == 1 ]]; } | bgzip > clean.vcf.gz
    gzip -cd ~{genotypes} | { grep -wvf bad.vargq.list || [[ $? == 1 ]]; } | bgzip > clean.genotypes.txt.gz
    gzip -cd ~{varGQ} | { grep -wvf bad.vargq.list || [[ $? == 1 ]]; } | bgzip > clean.vargq.txt.gz

    /opt/sv-pipeline/04_variant_resolution/scripts/add_genotypes.py \
      clean.vcf.gz \
      clean.genotypes.txt.gz \
      clean.vargq.txt.gz \
      ~{prefix}.genotyped.vcf

    mkdir tmp
    bcftools sort -T tmp/ ~{prefix}.genotyped.vcf -Oz -o ~{prefix}.genotyped.vcf.gz
    tabix ~{prefix}.genotyped.vcf.gz

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

task MakeSubsetVcf {
  input {
    File vcf
    File bed
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 3.75,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  String prefix = basename(bed, ".bed")

  output {
    File subset_vcf = "${prefix}.vcf.gz"
  }
  command <<<

    set -euxo pipefail
    zcat ~{vcf} | fgrep -e "#" > ~{prefix}.vcf
    zcat ~{vcf} | { fgrep -w -f <(cut -f4 ~{bed}) || true; } >> ~{prefix}.vcf
    bgzip ~{prefix}.vcf

  >>>
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

task AddBatchSamples {
  input {
    File batch_vcf
    File cohort_vcf
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 3.75,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File updated_vcf = "${prefix}.vcf.gz"
  }
  command <<<

    set -euo pipefail
    /opt/sv-pipeline/04_variant_resolution/scripts/add_batch_samples.py ~{batch_vcf} ~{cohort_vcf} ~{prefix}.vcf
    bgzip ~{prefix}.vcf
  
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
task IntegrateDepthGq {
  input {
    File vcf
    File RD_melted_genotypes
    File RD_vargq
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 3.75,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File genotypes = "genotype.indiv.depth.txt.gz"
    File varGQ = "genotype.variant.depth.txt.gz"
  }
  command <<<

    /opt/sv-pipeline/04_variant_resolution/scripts/IntegrateGQ_depthonly.sh \
      ~{vcf} \
      ~{RD_melted_genotypes} \
      ~{RD_vargq}
  
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
