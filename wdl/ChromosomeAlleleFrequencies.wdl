# Helper workflow to calculate basic AF statistics for a single chromosome on an input VCF

version 1.0

import "Structs.wdl"

# Add VAF annotation
workflow ChromosomeAlleleFrequencies {

  input {

    File   vcf
    File   vcf_idx
    String contig
    String prefix

    File? sample_pop_assignments   # Two-column file with sample ID & pop assignment. "." for pop will ignore sample
    File? ped_file                 # Used for M/F AF calculations
    File? par_bed
    File? allosomes_list


    String sv_pipeline_docker
    String sv_base_mini_docker

    RuntimeAttr? runtime_attr_shard_vcf
    RuntimeAttr? runtime_attr_compute_AFs
    RuntimeAttr? runtime_attr_combine_vcfs
  }

  call ComputeShardAFs {
      input:
        vcf = vcf,
        prefix = "${prefix}.${contig}",
        sample_pop_assignments = sample_pop_assignments,
        ped_file = ped_file,
        par_bed  = par_bed,
        allosomes_list = allosomes_list,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_compute_AFs
  }


  # Final output
  output {
    File vcf_wAFs = ComputeShardAFs.shard_wAFs
    File vcf_wAFs_idx = ComputeShardAFs.shard_wAFs_idx
  }
}


task ComputeShardAFs {
  input {
    File vcf
    String prefix
    String sv_pipeline_docker
    File? sample_pop_assignments
    File? ped_file
    File? par_bed
    File? allosomes_list
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 1.5,
    disk_gb: ceil(20 + size(vcf, "GB") * 2),
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail
    /opt/sv-pipeline/05_annotation/scripts/compute_AFs.py "~{vcf}" stdout \
      ~{"-p " + sample_pop_assignments} \
      ~{"-f " + ped_file} \
      ~{"-par " + par_bed} \
      ~{"--allosomes-list " + allosomes_list} \
    | bgzip -c \
    > "~{prefix}.wAFs.vcf.gz"

    tabix -p vcf "~{prefix}.wAFs.vcf.gz"
  >>>

  output {
    File shard_wAFs = "~{prefix}.wAFs.vcf.gz"
    File shard_wAFs_idx = "~{prefix}.wAFs.vcf.gz.tbi"
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
