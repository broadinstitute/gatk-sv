version 1.0

import "Structs.wdl"
import "GetShardInputs.wdl"
import "Utils.wdl" as util

# Does prelim translocation resolve from raw manta calls
workflow TinyResolve {
  input {
    Array[String] samples         # Sample ID
    File vcf_tar                  # Tarballed VCFs
    File cytoband
    Array[File] discfile
    File mei_bed
    Int samples_per_shard = 25
    String sv_pipeline_docker
    String linux_docker
    RuntimeAttr? runtime_attr_resolve
    RuntimeAttr? runtime_attr_untar
  }

  scatter (disc in discfile) {
    File discfile_idx = disc + ".tbi"
  }
  File cytoband_idx = cytoband + ".tbi"

  Int num_samples = length(samples)
  Float num_samples_float = num_samples
  Int num_shards = ceil(num_samples_float / samples_per_shard)

  call util.UntarFiles {
    input:
      tar = vcf_tar,
      glob_suffix = ".vcf.gz",
      linux_docker = linux_docker,
      runtime_attr_override = runtime_attr_untar
  }

  scatter (i in range(num_shards)) {
    call GetShardInputs.GetShardInputs as GetShardSamples {
      input:
        items_per_shard = samples_per_shard,
        shard_number = i,
        num_items = num_samples,
        all_items = samples
    }

    call GetShardInputs.GetShardInputs as GetShardDiscfiles {
      input:
        items_per_shard = samples_per_shard,
        shard_number = i,
        num_items = num_samples,
        all_items = discfile
    }

    call GetShardInputs.GetShardInputs as GetShardDiscfileIndexes {
      input:
        items_per_shard = samples_per_shard,
        shard_number = i,
        num_items = num_samples,
        all_items = discfile_idx
    }

    call GetShardInputs.GetShardInputs as GetShardVcfs {
      input:
        items_per_shard = samples_per_shard,
        shard_number = i,
        num_items = num_samples,
        all_items = UntarFiles.out
    }

    call ResolveManta {
      input:
        raw_vcfs=GetShardVcfs.shard_items,
        samples=GetShardSamples.shard_items,
        sv_pipeline_docker = sv_pipeline_docker,
        cytoband=cytoband,
        cytoband_idx=cytoband_idx,
        discfile=GetShardDiscfiles.shard_items,
        discfile_idx=GetShardDiscfileIndexes.shard_items,
        mei_bed=mei_bed,
        runtime_attr_override=runtime_attr_resolve
    }
  }

  output {
    Array[File] tloc_manta_vcf = flatten(ResolveManta.tloc_vcf)
  }
}


task ResolveManta {
  input {
    Array[File] raw_vcfs
    Array[String] samples
    File cytoband_idx
    Array[File] discfile
    Array[File] discfile_idx
    File cytoband
    File mei_bed
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Int num_samples = length(samples)
  Float input_size = size(discfile,"GiB")
  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 3.75, 
    disk_gb: ceil(10+input_size),
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail
    vcfs=(~{sep=" " raw_vcfs})
    sample_ids=(~{sep=" " samples})
    discfiles=(~{sep=" " discfile})
    for (( i=0; i<~{num_samples}; i++ ));
    do
      vcf=${vcfs[$i]}
      tabix -p vcf $vcf
      sample_id=${sample_ids[$i]}
      pe=${discfiles[$i]}
      sample_no=`printf %03d $i`
      bash /opt/sv-pipeline/00_preprocessing/scripts/mantatloccheck.sh $vcf $pe ${sample_id} ~{mei_bed} ~{cytoband}
      mv ${sample_id}.manta.complex.vcf.gz tloc_${sample_no}.${sample_id}.manta.complex.vcf.gz
    done
  >>>

  output {
    Array[File] tloc_vcf = glob("tloc_*.vcf.gz")
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
