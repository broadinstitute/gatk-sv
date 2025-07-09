version 1.0

# this wdl script would merge vcfs across individuals by size
# for SNVs and small indels 1-19bp, it'll require perfect match to merge variants
# for large indels 20bp and above as well as SVs, it'll use truvari to merge variants with flexibility in matching variants
# it should not generate multi-allelic variants in the output vcf

import "Structs.wdl"
import "MergeVcfsByChromosome.wdl" as MergeVcfsByChromosome
import "TruvariIntersample.wdl" as TruvariIntersample
import "LongReadGenotypeTasks.wdl" as LongReadGenotypeTasks

workflow MergeVcfsBySize {
  input {
    Array[File] input_vcfs     
    Array[File]? input_vcfs_idx
    Array[String] sample_list
    Array[String] chromosomes    

    File reference_fai
    File density_counter_py

    Boolean convert_to_biallelic = false
    String output_prefix 
    String sv_base_mini_docker
    String sv_pipeline_base_docker
  }

  if (!defined(input_vcfs_idx)) {
    scatter (idx in range(length(input_vcfs))) {
      call IndexVcf{
        input:
          vcf = input_vcfs[idx],
          sv_base_mini_docker = sv_base_mini_docker
      }
    }
  }

  scatter (idx in range(length(input_vcfs))) {
    call LongReadGenotypeTasks.SplitVariantsBySizeAt20bp {
      input:
        input_vcfs = input_vcfs[idx],
        docker_image = sv_pipeline_base_docker
    }
  }

  Array[File] sm_vcfs = SplitVariantsBySizeAt20bp.indels_lt_20_vcf
  Array[File] sm_vcf_idxes = SplitVariantsBySizeAt20bp.indels_lt_20_vcf_idx
  Array[File] lg_vcfs = SplitVariantsBySizeAt20bp.svs_gt_20_vcf
  Array[File] lg_vcf_idxes = SplitVariantsBySizeAt20bp.svs_gt_20_vcf_idx

  scatter (chrom in chromosomes) {
    call MergeVcfsByChromosome.MergeVcfsByChromosome as merge_vcfs_by_chrom_sm {
        input:
          chrom = chrom,
          input_vcfs = sm_vcfs,
          input_vcfs_idx = sm_vcf_idxes,
          sample_list = sample_list,
          convert_to_biallelic = convert_to_biallelic,
          sv_base_mini_docker = sv_base_mini_docker,
          sv_pipeline_base_docker = sv_pipeline_base_docker
    }

    call  TruvariIntersample.TruvariIntersample as merge_vcfs_by_chrom_lg {
      input:
        input_vcfs = lg_vcfs,
        input_vcfs_idx = lg_vcf_idxes,
        sample_ids = sample_list,
        chrom = chrom,
        reference_fai = reference_fai,
        density_counter_py = density_counter_py
    }

  }

  call ConcatVcfs {
    input:
      input_vcfs = merge_vcfs_by_chrom_sm.merged_vcf,
      input_vcfs_idx = merge_vcfs_by_chrom_sm.merged_vcf_idx,
      output_name = "${output_prefix}.vcf.gz",
      sv_base_mini_docker = sv_base_mini_docker
  }

  output {
    File final_merged_vcf = ConcatVcfs.output_vcf
    File final_merged_vcf_index = ConcatVcfs.output_vcf_idx
  }
}

# Task 1: Extract a chromosome from a VCF
task ExtractChromosomeVcf {
  input {
    File input_vcf
    String chromosome
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 2,
    disk_gb: ceil(size(input_vcf, "GB") * 2),
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])


  command <<<
    set -e
    bcftools view -r ~{chromosome} ~{input_vcf} -Oz -o ~{chromosome}.vcf.gz
    tabix -p vcf ~{chromosome}.vcf.gz
  >>>

  output {
    File output_vcf = "~{chromosome}.vcf.gz"
    File output_vcf_idx = "~{chromosome}.vcf.gz.tbi"
  }

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_base_mini_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

# Task 2: Merge multiple VCFs
task MergeVcfs {
  input {
    Array[File] input_vcfs
    String output_name
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 10,
    disk_gb: ceil(10 + size(input_vcfs, "GB") * 2),
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -e
    bcftools merge ~{sep=' ' input_vcfs} -Oz -o ~{output_name}
    tabix -p vcf ~{output_name}
  >>>

  output {
    File output_merged_vcf = output_name
    File output_merged_vcf_idx = "${output_name}.tbi"
  }

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_base_mini_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

# Task 3: Concatenate per-chromosome VCFs
task ConcatVcfs {
  input {
    Array[File] input_vcfs
    Array[File] input_vcfs_idx
    String output_name
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 10,
    disk_gb: ceil(10 + size(input_vcfs, "GB") * 2),
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])


  command <<<
    set -e
    bcftools concat ~{sep=' ' input_vcfs} -Oz -o ~{output_name}
    tabix -p vcf ~{output_name}
  >>>

  output {
    File output_vcf = output_name
    File output_vcf_idx = "${output_name}.tbi"
  }

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_base_mini_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

# Task 4: index VCFs
task IndexVcf {
  input {
    File vcf                # input VCF (.vcf.gz)
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 2,
    disk_gb: ceil(size(vcf, "GB") * 2),
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  String vcf_name = basename(vcf)
  command <<<
    set -e
      tabix -p vcf ~{vcf}
  >>>

  output {
    File indexed_vcf_idx = "~{vcf_name}.tbi"
  }

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_base_mini_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}
