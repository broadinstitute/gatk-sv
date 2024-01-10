version 1.0

import "Structs.wdl"
import "CleanVcf5.wdl" as cleanvcf5
import "TasksMakeCohortVcf.wdl" as tmc

workflow CalcAF {
  input {
    File vcf
    File vcf_idx
    Int sv_per_shard
    String prefix
    String sv_pipeline_docker
    String sv_pipeline_updates_docker
    File? sample_pop_assignments  #Two-column file with sample ID & pop assignment. "." for pop will ignore sample
    File? famfile                 #Used for M/F AF calculations
    File? par_bed                 #Used for marking hemizygous males on X & Y
    File? allosomes_list          #allosomes .fai used to override default sex chromosome assignments
    String? drop_empty_records

    RuntimeAttr? runtime_attr_compute_shard_af
    RuntimeAttr? runtime_attr_scatter_vcf
    RuntimeAttr? runtime_attr_combine_sharded_vcfs
  }


  # Tabix to chromosome of interest, and shard input VCF for stats collection
  call tmc.ScatterVcf {
    input:
      vcf=vcf,
      prefix=prefix,
      sv_pipeline_docker=sv_pipeline_updates_docker,
      records_per_shard=sv_per_shard,
      runtime_attr_override = runtime_attr_scatter_vcf
  }

  # Scatter over VCF shards
  scatter ( shard in ScatterVcf.shards ) {
    # Collect AF summary stats
    call ComputeShardAFs {
      input:
        vcf=shard,
        sv_pipeline_docker=sv_pipeline_docker,
        prefix=prefix,
        sample_pop_assignments=sample_pop_assignments,
        famfile=famfile,
        par_bed=par_bed,
        allosomes_list=allosomes_list,
        runtime_attr_override = runtime_attr_compute_shard_af
      }
  	}

  # Merge shards into single VCF
  call CombineShardedVcfs {
    input:
      vcfs=ComputeShardAFs.shard_wAFs,
      sv_pipeline_docker=sv_pipeline_docker,
      prefix=prefix,
      drop_empty_records=drop_empty_records,
      runtime_attr_override = runtime_attr_combine_sharded_vcfs
  }

  # Final output
  output {
    File vcf_wAFs = CombineShardedVcfs.vcf_out
    File vcf_wAFs_idx = CombineShardedVcfs.vcf_out_idx
  }
}

# Subset a vcf to a single chromosome, and add global AF information (no subpop)
task ComputeShardAFs {
  input {
    File vcf
    String prefix
    String sv_pipeline_docker
    File? sample_pop_assignments
    File? famfile
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
    optionals=" "
    if [ ~{default="SKIP" sample_pop_assignments} != "SKIP" ]; then
      optionals="$( echo "$optionals" ) -p ~{sample_pop_assignments}"
    fi
    if [ ~{default="SKIP" famfile} != "SKIP" ]; then
      optionals="$( echo "$optionals" ) -f ~{famfile}"
    fi
    if [ ~{default="SKIP" par_bed} != "SKIP" ]; then
      optionals="$( echo "$optionals" ) --par ~{par_bed}"
    fi
    if [ ~{default="SKIP" allosomes_list} != "SKIP" ]; then
      optionals="$( echo "$optionals" ) --allosomes-list ~{allosomes_list}"
    fi
    echo -e "OPTIONALS INTERPRETED AS: $optionals"
    echo -e "NOW RUNNING: /opt/sv-pipeline/05_annotation/scripts/compute_AFs.py $( echo "$optionals" ) ~{vcf} stdout"
    #Tabix chromosome of interest & compute AN, AC, and AF
    /opt/sv-pipeline/05_annotation/scripts/compute_AFs.py $optionals "~{vcf}" stdout \
    | bgzip -c \
    > "~{prefix}.wAFs.vcf.gz"
  >>>

  output {
    File shard_wAFs = "~{prefix}.wAFs.vcf.gz"
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


# Merge VCF shards & drop records with zero remaining non-ref alleles
task CombineShardedVcfs {
  input {
    Array[File] vcfs
    String prefix
    String sv_pipeline_docker
    String? drop_empty_records
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 4,
    disk_gb: 50,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command {
    set -euo pipefail
    vcf-concat -f ~{write_lines(vcfs)} \
    | vcf-sort \
    > merged.vcf
    if [ ~{default="TRUE" drop_empty_records} == "TRUE" ]; then
      /opt/sv-pipeline/05_annotation/scripts/prune_allref_records.py \
        merged.vcf stdout \
      | bgzip -c \
      > "~{prefix}.wAFs.vcf.gz"
    else
      cat merged.vcf | bgzip -c > "~{prefix}.wAFs.vcf.gz"
    fi
    tabix -p vcf "~{prefix}.wAFs.vcf.gz"
  }


  output {
    File vcf_out = "~{prefix}.wAFs.vcf.gz"
    File vcf_out_idx = "~{prefix}.wAFs.vcf.gz.tbi"
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

