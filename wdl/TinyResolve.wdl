version 1.0

import "Structs.wdl"
# Does perlim translocation resolve from raw manta calls
workflow TinyResolve {
  input {
    Array[String] samples         # Sample ID
    Array[File] manta_vcfs        # Manta VCF
    File cytoband
    Array[File] discfile
    File mei_bed
    Int samples_per_shard = 25
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr
  }

  scatter (disc in discfile) {
    File discfile_idx = disc + ".tbi"
  }
  File cytoband_idx = cytoband + ".tbi"

  Int num_samples = length(samples)
  Float num_samples_float = num_samples
  Int num_shards = ceil(num_samples_float / samples_per_shard)

  scatter (i in range(num_shards)) {
    scatter (j in range(samples_per_shard)) {
      Int idx_ = i * samples_per_shard + j
      if (idx_ < num_samples) {
        String shard_samples_ = samples[idx_]
        File shard_vcfs_ = manta_vcfs[idx_]
        File shard_discfiles_ = discfile[idx_]
        File shard_discfiles_idx_ = discfile_idx[idx_]
      }
    }
    Array[String] shard_samples = select_all(shard_samples_)
    Array[File] shard_vcfs = select_all(shard_vcfs_)
    Array[File] shard_discfiles = select_all(shard_discfiles_)
    Array[File] shard_discfiles_idx = select_all(shard_discfiles_idx_)

    call ResolveManta {
      input:
        raw_vcfs=shard_vcfs,
        samples=shard_samples,
        sv_pipeline_docker = sv_pipeline_docker,
        cytoband=cytoband,
        cytoband_idx=cytoband_idx,
        discfile=shard_discfiles,
        discfile_idx=shard_discfiles_idx,
        mei_bed=mei_bed,
        runtime_attr_override=runtime_attr
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
