version 1.0

import "Structs.wdl"
import "TasksMakeCohortVcf.wdl" as tasks_cohort

workflow RefineComplexVariants {
  input {
    File vcf
    String prefix

    Array[File] batch_sample_lists
    Array[File] pe_files
    Array[File] depth_del_beds
    Array[File] depth_dup_beds

    Int n_per_split = 10000
    Int min_pe_cpx = 3
    Int min_pe_ctx = 3

    String gatk_docker
    String sv_base_mini_docker

    RuntimeAttr? runtime_attr_scatter_vcf
    RuntimeAttr? runtime_attr_refine
    RuntimeAttr? runtime_attr_concat
  }

  # Create tarballs of depth BED files with ordered manifests
  call TarFilesWithOrderedManifest as TarDelBeds {
    input:
      files = depth_del_beds,
      tarball_name = "~{prefix}.depth_del_beds",
      folder_name = "depth_del_beds",
      sv_base_mini_docker = sv_base_mini_docker
  }

  call TarFilesWithOrderedManifest as TarDupBeds {
    input:
      files = depth_dup_beds,
      tarball_name = "~{prefix}.depth_dup_beds",
      folder_name = "depth_dup_beds",
      sv_base_mini_docker = sv_base_mini_docker
  }

  # Scatter the VCF
  call tasks_cohort.ScatterVcf {
    input:
      vcf = vcf,
      prefix = prefix,
      records_per_shard = n_per_split,
      sv_pipeline_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_scatter_vcf
  }

  # Run RefineComplexVariants on each shard
  scatter (i in range(length(ScatterVcf.shards))) {
    call RunRefineComplexVariants {
      input:
        vcf = ScatterVcf.shards[i],
        output_prefix = "~{prefix}.shard_~{i}",
        batch_sample_lists = batch_sample_lists,
        pe_files = pe_files,
        depth_del_beds_tarball = TarDelBeds.tarball,
        depth_dup_beds_tarball = TarDupBeds.tarball,
        min_pe_cpx = min_pe_cpx,
        min_pe_ctx = min_pe_ctx,
        gatk_docker = gatk_docker,
        runtime_attr_override = runtime_attr_refine
    }
  }

  # Gather refined VCFs and concatenate
  call tasks_cohort.ConcatVcfs {
    input:
      vcfs = RunRefineComplexVariants.refined_vcf,
      vcfs_idx = RunRefineComplexVariants.refined_vcf_index,
      naive = true,
      outfile_prefix = "~{prefix}.refined",
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_concat
  }

  output {
    File cpx_refined_vcf = ConcatVcfs.concat_vcf
    File cpx_refined_vcf_index = ConcatVcfs.concat_vcf_idx
  }
}

task TarFilesWithOrderedManifest {
  input {
    Array[File] files
    String tarball_name
    String folder_name
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  String output_tarball = "~{tarball_name}.tar.gz"
  String ordered_manifest_path = "~{folder_name}/manifest/ordered_manifest.txt"

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 2,
    disk_gb: ceil(10 + size(files, "GB") * 2),
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail
    mkdir "~{folder_name}"
    mkdir -p "~{folder_name}/manifest"
    : > "~{ordered_manifest_path}"

    # Preserve declared array order in a manifest that is embedded in the tarball.
    while IFS= read -r file; do
      base_name=$(basename "$file")
      target_path="~{folder_name}/${base_name}"
      if [[ -e "$target_path" ]]; then
        echo "Duplicate filename detected while preserving original names: ${base_name}" >&2
        exit 1
      fi

      cp "$file" "$target_path"
      echo "$target_path" >> "~{ordered_manifest_path}"
    done < ~{write_lines(files)}

    tar -czf ~{output_tarball} "~{folder_name}"
  >>>

  output {
    File tarball = output_tarball
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

task RunRefineComplexVariants {
  input {
    File vcf
    String output_prefix
    Array[File] batch_sample_lists
    Array[File] pe_files
    File depth_del_beds_tarball
    File depth_dup_beds_tarball
    Int min_pe_cpx
    Int min_pe_ctx
    String gatk_docker
    RuntimeAttr? runtime_attr_override
  }

  parameter_meta {
    pe_files: {
      localization_optional: true
    }
  }

  String output_vcf = "~{output_prefix}.vcf.gz"

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 16,
    disk_gb: ceil(20 + size([vcf, depth_del_beds_tarball, depth_dup_beds_tarball], "GB") * 3),
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail

    tar -xzf ~{depth_del_beds_tarball}
    tar -xzf ~{depth_dup_beds_tarball}

    gatk RefineComplexVariants \
      -V ~{vcf} \
      --min-pe-cpx ~{min_pe_cpx} \
      --min-pe-ctx ~{min_pe_ctx} \
      -O ~{output_vcf} \
      --batch-sample-lists ~{sep=" --batch-sample-lists " batch_sample_lists} \
      --discordant-pairs-files ~{sep=" --discordant-pairs-files " pe_files} \
      --depth-del-beds ~{sep=" --depth-del-beds " read_lines("depth_del_beds/manifest/ordered_manifest.txt")} \
      --depth-dup-beds ~{sep=" --depth-dup-beds " read_lines("depth_dup_beds/manifest/ordered_manifest.txt")}
  >>>

  output {
    File refined_vcf = output_vcf
    File refined_vcf_index = output_vcf + ".tbi"
  }

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
