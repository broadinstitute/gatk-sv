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
    String? python_docker
    RuntimeAttr? runtime_attr
    RuntimeAttr? runtime_attr_shard
  }

  scatter (disc in discfile) {
    File discfile_idx = disc + ".tbi"
  }
  File cytoband_idx = cytoband + ".tbi"

  Int num_samples = length(samples)
  Float num_samples_float = num_samples
  Int num_shards = ceil(num_samples_float / samples_per_shard)

  scatter (i in range(num_shards)) {
    call GetShardInputs {
      input:
        samples_per_shard = samples_per_shard,
        shard_number = i,
        num_samples = num_samples,
        all_samples = samples,
        all_discfiles = discfile,
        all_discfile_idxs = discfile_idx,
        all_vcfs = manta_vcfs,
        python_docker = python_docker,
        runtime_attr_override = runtime_attr_shard
    }

    call ResolveManta {
      input:
        raw_vcfs=GetShardInputs.shard_vcfs,
        samples=GetShardInputs.shard_samples,
        sv_pipeline_docker = sv_pipeline_docker,
        cytoband=cytoband,
        cytoband_idx=cytoband_idx,
        discfile=GetShardInputs.shard_discfiles,
        discfile_idx=GetShardInputs.shard_discfiles_idx,
        mei_bed=mei_bed,
        runtime_attr_override=runtime_attr
    }
  }

  output {
    Array[File] tloc_manta_vcf = flatten(ResolveManta.tloc_vcf)
  }
}

task GetShardInputs {
  input {
    Int samples_per_shard
    Int shard_number
    Int num_samples
    Array[String] all_samples
    Array[String] all_discfiles
    Array[String] all_discfile_idxs
    Array[String] all_vcfs
    String? python_docker
    RuntimeAttr? runtime_attr_override
  }
  String python_docker_ = select_first([python_docker, "python:3.7.12-bullseye"])
  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 3.75, 
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  command <<<
    set -euo pipefail
    python3 <<CODE
    all_samples = ['~{sep = "','" all_samples}']
    all_discfiles = ['~{sep = "','" all_discfiles}']
    all_discfile_idxs = ['~{sep = "','" all_discfile_idxs}']
    all_vcfs = ['~{sep = "','" all_vcfs}']
    start = ~{shard_number} * ~{samples_per_shard}
    end = min(start + ~{samples_per_shard}, ~{num_samples})
    with open("shard_discfiles.txt", 'w') as disc, open("shard_samples.txt", 'w') as samp, open("shard_disc_idx.txt", 'w') as idxs, open("shard_vcfs.txt", 'w') as vcfs:
      disc.write("\n".join(all_discfiles[start:end]))
      idxs.write("\n".join(all_discfile_idxs[start:end]))
      samp.write("\n".join(all_samples[start:end]))
      vcfs.write("\n".join(all_vcfs[start:end]))
    CODE
  >>>
  output {
    Array[String] shard_discfiles = read_lines("shard_discfiles.txt")
    Array[String] shard_discfiles_idx = read_lines("shard_disc_idx.txt")
    Array[String] shard_vcfs = read_lines("shard_vcfs.txt")
    Array[String] shard_samples = read_lines("shard_samples.txt")
  }
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: python_docker_
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
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
