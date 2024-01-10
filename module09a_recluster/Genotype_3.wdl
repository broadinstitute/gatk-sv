version 1.0
import "Structs.wdl"
workflow Regeno{
  input{
  File depth_vcf
  String batch
  File regeno_vcf
  File regeno_variants
  String sv_pipeline_docker
  }
  call ConcatRegenotypedVcfs{
    input:
    batch=batch,
    depth_vcf=depth_vcf,
    regeno_vcf=regeno_vcf,
    sv_pipeline_docker=sv_pipeline_docker,
    regeno_variants=regeno_variants
  }
  output {
    File genotyped_vcf = ConcatRegenotypedVcfs.genotyped_vcf
    File genotyped_vcf_idx = ConcatRegenotypedVcfs.genotyped_vcf_idx
  }
}

task ConcatRegenotypedVcfs {
    input{
    String batch
    File regeno_variants
    File depth_vcf
    File regeno_vcf
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 16,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  command <<<
    set -euo pipefail
    zcat ~{regeno_vcf} |fgrep "#" > head.txt
    zcat ~{regeno_vcf} |fgrep -f ~{regeno_variants} >body.txt
    cat head.txt body.txt|bgzip -c > regeno.vcf.gz
    zcat ~{depth_vcf} |fgrep -f ~{regeno_variants} -v |bgzip -c > no_variant.vcf.gz
    vcf-concat regeno.vcf.gz no_variant.vcf.gz \
      | vcf-sort -c \
      | bgzip -c > ~{batch}.depth.regeno_final.vcf.gz
    tabix ~{batch}.depth.regeno_final.vcf.gz
       >>>
  output {
    File genotyped_vcf = "~{batch}.depth.regeno_final.vcf.gz"
    File genotyped_vcf_idx = "~{batch}.depth.regeno_final.vcf.gz.tbi"
  }
  
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


