# Helper workflow to calculate basic AF statistics for a single chromosome on an input VCF

version 1.0

import "Structs.wdl"

# Add VAF annotation
workflow ChromosomeAlleleFrequencies {

  input {

    File   vcf
    File   vcf_idx
    Int    sv_per_shard
    String contig
    String prefix

    File? sample_pop_assignments   # Two-column file with sample ID & pop assignment. "." for pop will ignore sample
    File? ped_file                 # Used for M/F AF calculations

    String sv_pipeline_docker
    String sv_base_mini_docker

    RuntimeAttr? runtime_attr_shard_vcf
    RuntimeAttr? runtime_attr_compute_AFs
    RuntimeAttr? runtime_attr_combine_vcfs
  }

  # Tabix to chromosome of interest, and shard input VCF for stats collection
  call ShardVcf {
    input:
      vcf          = vcf,
      vcf_idx      = vcf_idx,
      contig       = contig,
      sv_per_shard = sv_per_shard,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_shard_vcf
  }

  # Scatter over VCF shards
  scatter ( shard in ShardVcf.shard_vcfs ) {
    # Collect AF summary stats
    call ComputeShardAlleleFrequencies {
      input:
        vcf                    = shard,
        prefix                 = "${prefix}.${contig}",
        sample_pop_assignments = sample_pop_assignments,
        ped_file               = ped_file,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_compute_AFs
    }
  }

  # Merge shards into single VCF
  call CombineShardedVcfs {
    input:
      vcfs   = ComputeShardAlleleFrequencies.shard_wAFs,
      prefix = "${prefix}.${contig}",
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_combine_vcfs
  }

  # Final output
  output {
    File vcf_wAFs = CombineShardedVcfs.vcf_out
    File vcf_wAFs_idx = CombineShardedVcfs.vcf_out_idx
  }
}

# Shard VCF into fixed size chunks
task ShardVcf {

  input {
    File   vcf
    File   vcf_idx
    Int    sv_per_shard
    String contig

    String sv_pipeline_docker

    RuntimeAttr? runtime_attr_override
  }
  
  output {
    Array[File] shard_vcfs = glob("vcf.shard.*.vcf.gz")
  }

  command <<<

    set -euo pipefail

    # Tabix chromosome of interest
    tabix -h ~{vcf} ~{contig} | bgzip -c > ~{contig}.vcf.gz
    
    # Then shard VCF
    /opt/sv-pipeline/scripts/shard_VCF.sh \
      ~{contig}.vcf.gz \
      ~{sv_per_shard} \
      "vcf.shard."

    # if there were no shards created just make an empty one
    if [ ! -e vcf.shard.000000.vcf.gz ]; then
      cp ~{contig}.vcf.gz vcf.shard.000000.vcf.gz
    fi
  >>>
  
  #########################
  RuntimeAttr default_attr = object {
    cpu_cores:          1, 
    mem_gb:             3.75, 
    disk_gb:            250,
    boot_disk_gb:       10,
    preemptible_tries:  3,
    max_retries:        0
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  runtime {
    cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
    memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
    disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
    preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
    docker:                 sv_pipeline_docker
  }
}

# Subset a vcf to a single chromosome, and add global AF information (no subpop)
task ComputeShardAlleleFrequencies {

  input {

    File   vcf
    String prefix
    
    File? sample_pop_assignments
    File? ped_file
    
    String sv_pipeline_docker

    RuntimeAttr? runtime_attr_override
  }
  
  output {
    File shard_wAFs = "${prefix}.wAFs.vcf.gz"
  }

  command <<<

    set -euo pipefail
    
    optionals=" "
    if ~{defined(sample_pop_assignments)}; then
      optionals="$( echo "$optionals" ) -p ~{sample_pop_assignments}"
    fi
    
    if ~{defined(ped_file)}; then
      optionals="$( echo "$optionals" ) -f ~{ped_file}"
    fi
    
    echo -e "OPTIONALS INTERPRETED AS: $optionals"
    echo -e "NOW RUNNING: /opt/sv-pipeline/05_annotation/scripts/compute_AFs.py $( echo "$optionals" ) ~{vcf} stdout"
    # Tabix chromosome of interest & compute AN, AC, and AF
    /opt/sv-pipeline/05_annotation/scripts/compute_AFs.py $optionals "~{vcf}" stdout \
      | bgzip -c \
      > "~{prefix}.wAFs.vcf.gz"
  
  >>>
  
  RuntimeAttr default_attr = object {
    cpu_cores:          1, 
    mem_gb:             3.75, 
    disk_gb:            20,
    boot_disk_gb:       10,
    preemptible_tries:  3,
    max_retries:        0
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  runtime {
    cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
    memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
    disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
    preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
    docker:                 sv_pipeline_docker
  }
}

# Merge VCF shards
task CombineShardedVcfs {

  input {
    
    Array[File] vcfs
    String      prefix

    String sv_base_mini_docker
    
    RuntimeAttr? runtime_attr_override
  }

  
  output {
    File vcf_out     = "${prefix}.wAFs.vcf.gz"
    File vcf_out_idx = "${prefix}.wAFs.vcf.gz.tbi"
  }

  command <<<

    set -euo pipefail
    vcf-concat ~{sep=" "  vcfs} \
      | vcf-sort \
      | bgzip -c \
      > "~{prefix}.wAFs.vcf.gz";
    tabix -p vcf "~{prefix}.wAFs.vcf.gz"
  
  >>>
 
  #########################
  RuntimeAttr default_attr = object {
    cpu_cores:          1, 
    mem_gb:             3.75, 
    disk_gb:            50,
    boot_disk_gb:       10,
    preemptible_tries:  3,
    max_retries:        0
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr]) 
  runtime {
    cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
    memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
    disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
    preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
    docker:                 sv_base_mini_docker
  }
}
