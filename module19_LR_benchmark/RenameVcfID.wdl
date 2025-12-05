version 1.0

import "Structs.wdl"
import "LongReadGenotypeTasks.wdl" as LongReadGenotypeTasks

workflow RenameVcfID {
  input {
    File vcf
    File? vcf_idx
    String chrom
    Int pos
    Int end
    String mid_fix
    String gatk_docker
    String sv_base_mini_docker
    String sv_pipeline_base_docker
    RuntimeAttr? runtime_attr_extract_target_variants
    RuntimeAttr? runtime_attr_extract_variant_sites
    RuntimeAttr? runtime_attr_split_to_biallelic
    RuntimeAttr? runtime_attr_index_vcf
  }

  if (!defined(vcf_idx)) {

    call IndexVcf{
      vcf = vcf,
      docker_image = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_index_vcf
    }

  }

  vcf_index = select_first([IndexVcf.vcf_idx, vcf_idx])
 
 

  call ExtractRegionFromVCF{
    input:
      vcf = vcf,
      vcf_idx = vcf_index,
      chrom = chrom,
      pos = pos,
      end = end,
      mid_fix = mid_fix,
      docker_image = sv_pipeline_base_docker,
      runtime_attr_override = runtime_attr_extract_target_variants
  }

  call LongReadGenotypeTasks.ExtractVariantSites {
    input:
      input_vcf = ExtractRegionFromVCF.region_vcf,
      docker_image = sv_pipeline_base_docker,
      runtime_attr_override = runtime_attr_extract_variant_sites
  }

  call LongReadGenotypeTasks.SplitMultiAllelicToBiAllelicSingleSample{
    input:
      vcf_file = ExtractVariantSites.updated_vcf,
      vcf_idx  = ExtractVariantSites.updated_vcf_idx,
      docker_image = sv_pipeline_base_docker,
      runtime_attr_override = runtime_attr_split_to_biallelic
  }

  output {
    File output_vcf = SplitMultiAllelicToBiAllelicSingleSample.biallelic_vcf
    File output_vcf_idx = SplitMultiAllelicToBiAllelicSingleSample.biallelic_vcf_idx
  }
}

task IndexVcf {
    input {
        File vcf                                  # .vcf.gz or .bcf
        String docker_image
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail
        tabix -p vcf 
    >>>

    output {
        File vcf_idx = "~{vcf}.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 10 + ceil(size(vcf, "GiB")*2),
        disk_gb: 15 + ceil(size(vcf, "GiB")*2),
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

task ExtractRegionFromVCF {
    input {
        File vcf        # .vcf.gz or .bcf
        File vcf_idx  # .tbi or .csi
        String chrom
        Int pos
        Int end
        String mid_fix
        String docker_image
        RuntimeAttr? runtime_attr_override
    }

    String prefix = basename(vcf, ".vcf.gz")

    command <<<
        bcftools view \
            -r ~{chrom}:~{pos}-~{end} \
            -O z \
            -o "~{prefix}.~{mid_fix}.vcf.gz" \
            ~{vcf}
        
        tabix -p vcf "~{prefix}.~{mid_fix}.vcf.gz"
    >>>

    output {
        File region_vcf = "~{prefix}.~{mid_fix}.vcf.gz"
        File region_vcf_idx = "~{prefix}.~{mid_fix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 10 + ceil(size(vcf, "GiB")*2),
        disk_gb: 15 + ceil(size(vcf, "GiB")*2),
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

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}
