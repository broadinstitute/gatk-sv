version 1.0

# Derived from gatk-sv wdl/{Structs,Utils,TinyResolve}.wdl -- see build_tloc_wdl.py
# Only the tloc path of GatherBatchEvidence; unrelated siblings (gCNV, CNMOPS,
# CondenseReadCounts, matrices, metrics) are not called.

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow TlocResolveOnly {
  input {
    Array[String] samples         # Sample IDs, production order
    File vcf_tar                  # standardized VCF tar, production PreprocessPESR output
    File cytoband
    Array[File] discfile
    File mei_bed
    Int samples_per_shard = 25
    String sv_pipeline_docker
    String linux_docker
    RuntimeAttr? runtime_attr_resolve
    RuntimeAttr? runtime_attr_untar
  }

  # --- identical to TinyResolve ---
  scatter (disc in discfile) {
    File discfile_idx = disc + ".tbi"
  }
  File cytoband_idx = cytoband + ".tbi"

  Int num_samples = length(samples)
  Float num_samples_float = num_samples
  Int num_shards = ceil(num_samples_float / samples_per_shard)

  call UntarFiles {
    input:
      tar = vcf_tar,
      glob_suffix = ".vcf.gz",
      linux_docker = linux_docker,
      runtime_attr_override = runtime_attr_untar
  }

  # GetShardInputs upstream is a pure-WDL workflow (scatter + guarded index + select_all) and
  # therefore costs no VMs. Inlined with the same semantics rather than rewritten as tasks, so
  # sharding stays free and the slices are the same slices production computed.
  scatter (i in range(num_shards)) {
    scatter (j in range(samples_per_shard)) {
      Int idx = i * samples_per_shard + j
      if (idx < num_samples) {
        String s_item = samples[idx]
      }
      if (idx < num_samples) {
        File d_item = discfile[idx]
      }
      if (idx < num_samples) {
        File di_item = discfile_idx[idx]
      }
      if (idx < num_samples) {
        File v_item = UntarFiles.out[idx]
      }
    }
    Array[String] shard_samples = select_all(s_item)
    Array[File] shard_discfile = select_all(d_item)
    Array[File] shard_discfile_idx = select_all(di_item)
    Array[File] shard_vcfs = select_all(v_item)

    call ResolveManta {
      input:
        raw_vcfs = shard_vcfs,
        samples = shard_samples,
        sv_pipeline_docker = sv_pipeline_docker,
        cytoband = cytoband,
        cytoband_idx = cytoband_idx,
        discfile = shard_discfile,
        discfile_idx = shard_discfile_idx,
        mei_bed = mei_bed,
        runtime_attr_override = runtime_attr_resolve
    }
  }

  output {
    Array[File] tloc_manta_vcf = flatten(ResolveManta.tloc_vcf)
  }
}

task UntarFiles {
  input {
    File tar
    String? glob_suffix
    String linux_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 1.0,
                               disk_gb: ceil(10 + 2 * size(tar, "GB")),
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  String glob_arg = "out/*" + select_first([glob_suffix, ""])

  command <<<
    set -euo pipefail
    mkdir out
    tar xzf ~{tar} -C out/
  >>>

  output {
    Array[File] out = glob("~{glob_arg}")
  }

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: linux_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    noAddress: true
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
    noAddress: true
  }
}
