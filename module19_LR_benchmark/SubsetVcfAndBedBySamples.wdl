version 1.0

import "Structs.wdl"

workflow SubsetVcfAndBedBySamples {
  input {
    File vcf_gz
    File vcf_index
    File sample_annotation
    String output_prefix
    String mid_fix  
    String sv_pipeline_base_docker
    RuntimeAttr? runtime_attr_subset_vcf_by_samples
    RuntimeAttr? runtime_attr_vcf2bed
  }

  call subset_vcf_by_samples {
    input:
      vcf_gz = vcf_gz,
      vcf_index = vcf_index,
      sample_annotation = sample_annotation,
      output_prefix = output_prefix,
      mid_fix = mid_fix,
      docker_image = sv_pipeline_base_docker,
      runtime_attr_override = runtime_attr_subset_vcf_by_samples
  }

  call vcf2bed {
    input: 
      vcf_gz  = subset_vcf_by_samples.subset_vcf,
      vcf_index = subset_vcf_by_samples.subset_vcf_index, 
      docker_image = sv_pipeline_base_docker,
      runtime_attr_override = runtime_attr_vcf2bed
  }


  output {
    File subset_vcf = subset_vcf_by_samples.subset_vcf
    File subset_vcf_index = subset_vcf_by_samples.subset_vcf_index
    File subset_bed = vcf2bed.subset_bed
  }
}


task subset_vcf_by_samples {
  input {
    File vcf_gz          # Input VCF (bgzipped)
    File vcf_index       # .tbi index
    File sample_annotation     # Text file: one sample ID per line
    String output_prefix # Output prefix
    String mid_fix 
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  String prefix = basename(vcf_gz, ".vcf.gz")
  command <<<
    set -euo pipefail

    # Subset VCF
    bcftools view -f PASS \
    -S <(cut -f1 ~{sample_annotation} | tail -n+2 ) \
    ~{vcf_gz} \
    -Oz -o ~{prefix}.subset.vcf.gz

    tabix -p vcf  ~{prefix}.subset.vcf.gz


    #replace sample name
    bcftools reheader \
    -s <(tail -n+2 ~{sample_annotation}) \
    -o ~{output_prefix}.~{mid_fix}.vcf.gz \
    ~{prefix}.subset.vcf.gz

    tabix -p vcf ~{output_prefix}.~{mid_fix}.vcf.gz
  >>>

  output {
    File subset_vcf = "~{output_prefix}.~{mid_fix}.vcf.gz"
    File subset_vcf_index = "~{output_prefix}.~{mid_fix}.vcf.gz.tbi"
  }

  RuntimeAttr default_attr = object {
      cpu_cores: 4,
      mem_gb: 10 + ceil(size(vcf_gz,"GiB"))*2,
      disk_gb: 20 + ceil(size(vcf_gz,"GiB"))*2,
      boot_disk_gb: 20,
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

task vcf2bed {
  input {
    File vcf_gz          # Input VCF (bgzipped)
    File vcf_index       # .tbi index
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  String prefix = basename(vcf_gz, ".vcf.gz")
  command <<<
    set -euo pipefail
    bcftools view -c 1 ~{vcf_gz} ~{prefix}.non_ref.vcf.gz
    tabix -p vcf ~{prefix}.non_ref.vcf.gz 
    svtk vcf2bed -i SVTYPE -i SVLEN -i AC -i AN --include-filters ~{prefix}.non_ref.vcf.gz ~{prefix}.bed
    bgzip ~{prefix}.bed
  >>>

  output {
    File subset_bed = "~{prefix}.bed.gz"
  }

  RuntimeAttr default_attr = object {
      cpu_cores: 4,
      mem_gb: 10 + ceil(size(vcf_gz,"GiB"))*2,
      disk_gb: 20 + ceil(size(vcf_gz,"GiB"))*2,
      boot_disk_gb: 20,
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


