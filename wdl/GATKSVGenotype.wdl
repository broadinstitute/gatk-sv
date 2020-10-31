version 1.0

import "Structs.wdl"
import "GATKSVGenotypeInnerScatter.wdl" as InnerScatter

workflow GATKSVGenotype {
  input {
    File vcf
    File sample_mean_depth_file
    String batch

    Int records_per_shard = 2000
    Int discrete_samples = 1000

    # Model parameters
    Float? eps_pe
    Float? eps_sr1
    Float? eps_sr2

    String genotyping_gatk_docker
    String sharding_gatk_docker
    String sv_base_mini_docker

    # cpu or cuda
    String device_train = "cpu"
    String device_genotype = "cpu"

    String gpu_type = "nvidia-tesla-k80"
    String nvidia_driver_version = "450.80.02"

    RuntimeAttr? runtime_attr_split
    RuntimeAttr? runtime_attr_shard
    RuntimeAttr? runtime_attr_train
    RuntimeAttr? runtime_attr_infer
    RuntimeAttr? runtime_attr_concat
  }

  call SplitDepthCalls {
    input:
      vcf = vcf,
      vcf_index = vcf + ".tbi",
      output_name = batch,
      gatk_docker = sharding_gatk_docker,
      runtime_attr_override = runtime_attr_split
  }

  # Note INV records were converted to BND by SVCluster
  Array[String] svtypes = ["DEL", "DUP", "INS", "BND"]

  scatter (svtype in svtypes) {
    String shard_name = batch + "." + svtype
    call ShardVcf {
      input:
        vcf = SplitDepthCalls.out_pesr,
        vcf_index = SplitDepthCalls.out_pesr_index,
        records_per_shard = records_per_shard,
        svtype = svtype,
        basename = shard_name,
        gatk_docker = sharding_gatk_docker,
        runtime_attr_override = runtime_attr_shard
    }
    call InnerScatter.GATKSVGenotypeInnerScatter {
      input:
        vcfs = ShardVcf.out,
        sample_mean_depth_file = sample_mean_depth_file,
        records_per_shard = records_per_shard,
        discrete_samples = discrete_samples,
        eps_pe = eps_pe,
        eps_sr1 = eps_sr1,
        eps_sr2 = eps_sr2,
        genotyping_gatk_docker = genotyping_gatk_docker,
        device_train = device_train,
        device_genotype = device_genotype,
        gpu_type = gpu_type,
        nvidia_driver_version = nvidia_driver_version,
        runtime_attr_train = runtime_attr_train,
        runtime_attr_infer = runtime_attr_infer
    }
  }

  Array[File] genotyped_vcf_shards = flatten([[SplitDepthCalls.out_depth], flatten(GATKSVGenotypeInnerScatter.out)])
  Array[File] genotyped_vcf_shard_indexes = flatten([[SplitDepthCalls.out_depth_index], flatten(GATKSVGenotypeInnerScatter.out_index)])

  call ConcatVcfs {
    input:
      vcfs = genotyped_vcf_shards,
      vcfs_idx = genotyped_vcf_shard_indexes,
      outfile_prefix = "~{batch}.final",
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_concat
  }

  output {
    File genotyped_vcf = ConcatVcfs.out
    File genotyped_vcf_index = ConcatVcfs.out_index
  }
}

task SplitDepthCalls {
  input {
    File vcf
    File vcf_index
    String output_name
    String gatk_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 1.0,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  Float mem_gb = select_first([runtime_attr.mem_gb, default_attr.mem_gb])
  Int java_mem_mb = ceil(mem_gb * 1000 * 0.8)

  output {
    File out_depth = "~{output_name}.depth.vcf.gz"
    File out_depth_index = "~{output_name}.depth.vcf.gz.tbi"
    File out_pesr = "~{output_name}.pesr.vcf.gz"
    File out_pesr_index = "~{output_name}.pesr.vcf.gz.tbi"
  }
  command <<<
    set -euo pipefail

    gatk --java-options -Xmx~{java_mem_mb}M SelectVariants \
      -V ~{vcf} \
      -O ~{output_name}.depth.vcf.gz \
      -select "ALGORITHMS == 'depth'"

    gatk --java-options -Xmx~{java_mem_mb}M SelectVariants \
      -V ~{vcf} \
      -O ~{output_name}.pesr.vcf.gz \
      -select "ALGORITHMS == 'depth'" \
      --invertSelect
  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: mem_gb + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: gatk_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

# Combine multiple sorted VCFs
task ConcatVcfs {
  input {
    Array[File] vcfs
    Array[File]? vcfs_idx
    String? outfile_prefix
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  String outfile_name = outfile_prefix + ".vcf.gz"

  # when filtering/sorting/etc, memory usage will likely go up (much of the data will have to
  # be held in memory or disk while working, potentially in a form that takes up more space)
  Float input_size = size(vcfs, "GB")
  Float compression_factor = 5.0
  Float base_disk_gb = 5.0
  Float base_mem_gb = 2.0
  RuntimeAttr runtime_default = object {
    mem_gb: base_mem_gb + compression_factor * input_size,
    disk_gb: ceil(base_disk_gb + input_size * (2.0 + compression_factor)),
    cpu_cores: 1,
    preemptible_tries: 3,
    max_retries: 1,
    boot_disk_gb: 10
  }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_base_mini_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euo pipefail
    VCFS="~{write_lines(vcfs)}"
    if ~{!defined(vcfs_idx)}; then
      cat ${VCFS} | xargs -n1 tabix
    fi
    bcftools concat -a --output-type z --file-list ${VCFS} --output "~{outfile_name}"
    tabix -p vcf -f "~{outfile_name}"
  >>>

  output {
    File out = outfile_name
    File out_index = outfile_name + ".tbi"
  }
}

task ShardVcf {
  input {
    File vcf
    File vcf_index
    String? svtype
    Int records_per_shard
    String basename
    String gatk_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 0.9,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  Float mem_gb = select_first([runtime_attr.mem_gb, default_attr.mem_gb])
  Int java_mem_mb = ceil(mem_gb * 1000 * 0.8)

  output {
    Array[File] out = glob("*.vcf.gz")
  }
  command <<<

    set -euo pipefail
    gatk --java-options -Xmx~{java_mem_mb}M SelectVariants \
      -V ~{vcf} \
      -O ~{basename} \
      --max-variants-per-shard ~{records_per_shard} \
      ~{"-select \"SVTYPE == '" + svtype + "'\""} \

  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: gatk_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}