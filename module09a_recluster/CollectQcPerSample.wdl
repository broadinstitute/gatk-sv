version 1.0

import "TasksMakeCohortVcf.wdl" as MiniTasks

# Workflow to gather lists of variant IDs per sample from one or more SV VCFs
workflow CollectQcPerSample {
  input {
    Array[File] vcfs
    Boolean vcf_format_has_cn = true
    File samples_list
    String prefix

    String sv_base_mini_docker
    String sv_pipeline_docker

    # overrides for local tasks
    RuntimeAttr? runtime_override_collect_vids_per_sample

    # overrides for mini tasks
    RuntimeAttr? runtime_override_split_samples_list
    RuntimeAttr? runtime_override_merge_sharded_per_sample_vid_lists
  }

  String output_prefix = "~{prefix}.per_sample_qc"

  # Collect VCF-wide summary stats per sample list per VCF
  scatter ( vcf in vcfs ) {
    call CollectVidsPerSample {
      input:
        vcf=vcf,
        vcf_format_has_cn=vcf_format_has_cn,
        samples_list=samples_list,
        prefix=output_prefix,
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_override_collect_vids_per_sample
    }
  }

  # Merge all VID lists into single output directory and tar it
  call MergeShardedPerSampleVidLists {
    input:
      tarballs=CollectVidsPerSample.vid_lists_tarball,
      samples_list=samples_list,
      prefix=output_prefix,
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_override_merge_sharded_per_sample_vid_lists
  }

  # Final output
  output {
    File vid_lists = MergeShardedPerSampleVidLists.merged_tarball
  }
}


# Task to collect list of VIDs per sample
task CollectVidsPerSample {
  input {
    File vcf
    Boolean vcf_format_has_cn = true
    File samples_list
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  String outdirprefix = "~{prefix}_perSample_VIDs"
  
  # Must scale disk proportionally to size of input VCF
  Float input_size = size([vcf, samples_list], "GiB")
  RuntimeAttr runtime_default = object {
    mem_gb: 1.5,
    disk_gb: ceil(10.0 + input_size),
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
    if [ ~{vcf_format_has_cn} == "true" ]; then
      bcftools view -S ~{samples_list} ~{vcf} \
      | bcftools +fill-tags -- -t AC \
      | bcftools view -i 'SVTYPE=="CNV" || AC>0' \
      | bcftools query -f '[%SAMPLE\t%ID\t%ALT\t%GT\t%GQ\t%CN\t%CNQ\n]' \
      | awk '{OFS="\t"; gt = $4; gq = $5; if ($3 == "<CNV>") { gq = $7; if ($6 == 2) { gt = "0/0" } else if ($6 == 1 || $6 == 3) { gt = "0/1" } else { gt = "1/1"} }; print $1, $2, gt, gq}' \
      | awk -v outprefix="~{outdirprefix}" '$3 != "0/0" && $3 != "./." {OFS="\t"; print $2, $3, $4 >> outprefix"/"$1".VIDs_genotypes.txt" }'
    else
      bcftools view -S ~{samples_list} ~{vcf} \
      | bcftools +fill-tags -- -t AC \
      | bcftools view -i 'SVTYPE=="CNV" || AC>0' \
      | bcftools query -f '[%SAMPLE\t%ID\t%ALT\t%GT\t%GQ\n]' \
      | awk '{OFS="\t"; gt = $4; gq = $5; if ($3 ~ /CN0/) { if ($4 == "0/2") { gt = "0/0" } else if ($4 == "0/1" || $4 == "0/3") { gt = "0/1" } else { gt = "1/1"} }; print $1, $2, gt, gq}' \
      | awk -v outprefix="~{outdirprefix}" '$3 != "0/0" && $3 != "./." {OFS="\t"; print $2, $3, $4 >> outprefix"/"$1".VIDs_genotypes.txt" }'
    fi

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


# Merge multiple tarballs of per-sample VID lists
task MergeShardedPerSampleVidLists {
  input {
    Array[File] tarballs
    File samples_list
    String prefix
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr runtime_default = object {
    mem_gb: 3.75,
    disk_gb: 20,
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
    docker: sv_base_mini_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -eu -o pipefail

    # Create final output directory
    mkdir "~{prefix}_perSample_VID_lists"

    # Extract each tarball into its own unique directory
    mkdir shards
    while read i tarball_path; do
      mkdir "shards/shard_$i"
      tar -xzvf "$tarball_path" --directory "shards/shard_$i"/
    done < <( awk -v OFS="\t" '{ print NR, $1 }' ~{write_lines(tarballs)} )

    # Merge all shards per sample and write to final directory
    while read sample; do
      find shards/ -name "$sample.VIDs_genotypes.txt.gz" \
      | xargs -I {} zcat {} \
      | sort -Vk1,1 -k2,2n -k3,3n \
      | gzip -c \
      > "~{prefix}_perSample_VID_lists/$sample.VIDs_genotypes.txt.gz"
    done < ~{samples_list}

    # Compress final output directory
    tar -czvf "~{prefix}_perSample_VID_lists.tar.gz" "~{prefix}_perSample_VID_lists"
  >>>

  output {
    File merged_tarball = "~{prefix}_perSample_VID_lists.tar.gz"
  }
}
