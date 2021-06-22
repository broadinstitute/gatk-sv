version 1.0

# Author: Ryan Collins <rlcollins@g.harvard.edu>

import "TasksMakeCohortVcf.wdl" as MiniTasks

# Workflow to gather lists of variant IDs per sample from an SV VCF
workflow CollectQcPerSample {
  input {
    File vcf
    File samples_list
    String prefix
    Int samples_per_shard

    String sv_base_mini_docker
    String sv_pipeline_docker

    # overrides for local tasks
    RuntimeAttr? runtime_override_collect_vids_per_sample

    # overrides for mini tasks
    RuntimeAttr? runtime_override_split_samples_list
    RuntimeAttr? runtime_override_tar_shard_vid_lists
  }

  # Shard sample list
  call MiniTasks.SplitUncompressed as SplitSamplesList {
    input:
      whole_file=samples_list,
      lines_per_shard=samples_per_shard,
      shard_prefix=prefix + ".list_shard.",
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_override_split_samples_list
  }

  # Collect VCF-wide summary stats per sample list
  scatter (sublist in SplitSamplesList.shards) {
    call CollectVidsPerSample {
      input:
        vcf=vcf,
        samples_list=sublist,
        prefix=prefix,
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_override_collect_vids_per_sample
    }
  }

  # Merge all VID lists into single output directory and tar it
  call TarShardVidLists {
    input:
      in_tarballs=CollectVidsPerSample.vid_lists_tarball,
      folder_name=prefix + "_perSample_VIDs_merged",
      tarball_prefix=prefix + "_perSample_VIDs",
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_override_tar_shard_vid_lists
  }

  # Final output
  output {
    File vid_lists = TarShardVidLists.tarball
  }
}


# Task to collect list of VIDs per sample
task CollectVidsPerSample {
  input {
    File vcf
    File samples_list
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  String outdirprefix = prefix + "_perSample_VIDs"
  
  # Must scale disk proportionally to size of input VCF
  Float input_size = size([vcf, samples_list], "GiB")
  Float disk_scaling_factor = 1.5
  Float base_disk_gb = 10.0
  Float base_mem_gb = 3.75
  RuntimeAttr runtime_default = object {
    mem_gb: 3.75,
    disk_gb: ceil(base_disk_gb + (input_size * disk_scaling_factor)),
    cpu_cores: 1,
    preemptible_tries: 3,
    max_retries: 1,
    boot_disk_gb: 10
  }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GiB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<

    set -eu -o pipefail

    # Make output directory
    mkdir -p ~{outdirprefix}

    # Filter VCF to list of samples of interest, split into list of genotypes per 
    # sample, and write one .tsv file per sample to output directory
    bcftools view -S ~{samples_list} ~{vcf} \
    | bcftools view --min-ac 1 \
    | bcftools query -f '[%SAMPLE\t%ID\t%ALT\t%GT\t%GQ\t%CN\t%CNQ\n]' \
    | awk '{OFS="\t"; gt = $4; gq = $5; if ($3 == "<CNV>") { gq = $7; if ($6 == 2) { gt = "0/0" } else if ($6 == 1 || $6 == 3) { gt = "0/1" } else { gt = "1/1"} }; print $1, $2, gt, gq}' \
    | awk -v outprefix="~{outdirprefix}" '$3 != "0/0" && $3 != "./." {OFS="\t"; print $2, $3, $4 >> outprefix"/"$1".VIDs_genotypes.txt" }'

    # Gzip all output lists
    for FILE in ~{outdirprefix}/*.VIDs_genotypes.txt; do
      gzip -f "$FILE"
      rm -f "$FILE"
    done

    # Bundle all files as a tarball (to make it easier on call caching for large cohorts)
    cd ~{outdirprefix} && \
    tar -czvf ../~{outdirprefix}.tar.gz *.VIDs_genotypes.txt.gz && \
    cd -
  >>>

  output {
    File vid_lists_tarball = "~{outdirprefix}.tar.gz"
  }
}


# Task to merge VID lists across shards
task TarShardVidLists {
  input {
    Array[File] in_tarballs
    String? folder_name
    String? tarball_prefix
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  String tar_folder_name = select_first([folder_name, "merged"])
  String outfile_name = select_first([tarball_prefix, tar_folder_name]) + ".tar.gz"

  # Since the input files are often/always compressed themselves, assume compression factor for tarring is 1.0
  Float input_size = size(in_tarballs, "GB")
  Float base_disk_gb = 5.0
  Float base_mem_gb = 2.0
  RuntimeAttr runtime_default = object {
    mem_gb: base_mem_gb,
    disk_gb: ceil(base_disk_gb + input_size * 2.0),
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
    # Create final output directory
    mkdir "~{tar_folder_name}"

    while read tarball_path; do
      tar -xzvf "$tarball_path" --directory ~{tar_folder_name}/
    done < ~{write_lines(in_tarballs)}

    # Compress final output directory
    tar -czvf "~{outfile_name}" "~{tar_folder_name}"
  >>>

  output {
    File tarball = outfile_name
  }
}
