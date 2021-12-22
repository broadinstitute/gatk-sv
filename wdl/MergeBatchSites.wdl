version 1.0

import "Structs.wdl"

workflow MergeBatchSites {
  input {
    Array[File] depth_vcfs    # Filtered depth VCFs across batches
    Array[File] pesr_vcfs     # Filtered PESR VCFs across batches
    String cohort             # Cohort name or project prefix for all cohort-level outputs
    String sv_pipeline_docker
    String contig
    File original_cohort_depth_vcf
    File original_cohort_pesr_vcf
    RuntimeAttr? runtime_attr_merge_pesr
    RuntimeAttr? runtime_attr_merge_depth
    RuntimeAttr? runtime_attr_subset_depth
    RuntimeAttr? runtime_attr_subset_pesr
    RuntimeAttr? runtime_attr_update_pesr
    RuntimeAttr? runtime_attr_update_depth
  }

  scatter (i in range(length(depth_vcfs))) {
    call SubsetVcfToContig as SubsetDepthVcf {
      input:
        vcf = depth_vcfs[i],
        contig = contig,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_subset_depth
    }
    call SubsetVcfToContig as SubsetPESRVcf {
      input:
        vcf = pesr_vcfs[i],
        contig = contig,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_subset_pesr
    }
  }

  call MergeVcfs as MergePESRVcfs {
    input:
      vcfs = SubsetPESRVcf.subset_vcf,
      prefix = cohort + ".chrX.all_batches.pesr",
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_merge_pesr
  }

  call MergeDepthVcfs {
    input:
      vcfs = SubsetDepthVcf.subset_vcf,
      cohort = cohort,
      prefix = cohort + ".chrX.all_batches.depth",
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_merge_depth
  }

  call UpdateChromosomeX as UpdateCohortDepthVcf {
    input:
      original_vcf = original_cohort_depth_vcf,
      chrx_vcf = MergeDepthVcfs.merged_vcf,
      prefix = cohort + ".all_batches.depth",
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_update_depth
  }

  output {
    File subset_cohort_pesr_vcf = MergePESRVcfs.merged_vcf
    File subset_cohort_depth_vcf = MergeDepthVcfs.merged_vcf
    File cohort_depth_vcf = UpdateCohortDepthVcf.updated_vcf
    Array[File] subset_depth_vcfs = SubsetDepthVcf.subset_vcf
    Array[File] subset_pesr_vcfs = SubsetPESRVcf.subset_vcf
  }
}

task UpdateChromosomeX {
  input {
    File original_vcf
    File? original_vcf_index
    File chrx_vcf
    String prefix
    String sv_pipeline_docker
    Boolean? create_index
    RuntimeAttr? runtime_attr_override
  }

  Boolean create_index_ = if (defined(create_index)) then create_index else false

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
    File updated_vcf = "~{prefix}.vcf.gz"
    File updated_vcf_index = "~{prefix}.vcf.gz.tbi"
  }
  command <<<
    set -euxo pipefail
    if [ -z ~{original_vcf_index} ]; then
      tabix -p vcf ~{original_vcf}
    fi
    bcftools view -r chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22 ~{original_vcf} | bgzip -c > before.vcf.gz
    bcftools view -r chrY ~{original_vcf} | bgzip -c > after.vcf.gz
    bcftools concat --naive before.vcf.gz ~{chrx_vcf} after.vcf.gz > "~{prefix}.vcf.gz"
    if [ ~{create_index_} == "true" ]; then
      tabix -p vcf ~{prefix}.vcf.gz
    else
      touch ~{prefix}.vcf.gz.tbi
    fi
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


task SubsetVcfToContig {
  input {
    File vcf
    String contig
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  String vcf_basename = basename(vcf, ".vcf.gz")

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
    File subset_vcf = "~{vcf_basename}.~{contig}.vcf.gz"
  }
  command <<<
    set -euo pipefail
    bcftools view -t ~{contig} ~{vcf} | bgzip -c > "~{vcf_basename}.~{contig}.vcf.gz"
    
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

task MergeVcfs {
  input {
    Array[File] vcfs
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
    File merged_vcf = "~{prefix}.vcf.gz"
  }
  command <<<

    set -euo pipefail
    /opt/sv-pipeline/04_variant_resolution/scripts/merge_vcfs.py ~{write_lines(vcfs)} ~{prefix}.vcf
    rm ~{sep=' ' vcfs}
    vcf-sort -c ~{prefix}.vcf | bgzip -c > ~{prefix}.vcf.gz
  
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

task MergeDepthVcfs {
  input {
    Array[File] vcfs
    String cohort
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
    File merged_vcf = "~{prefix}.vcf.gz"
  }
  command <<<
    set -euxo pipefail
    /opt/sv-pipeline/04_variant_resolution/scripts/merge_vcfs.py ~{write_lines(vcfs)} ~{prefix}.vcf
    vcf-sort -c ~{prefix}.vcf | bgzip -c > ~{prefix}.vcf.gz
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
