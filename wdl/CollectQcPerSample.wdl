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
  call MiniTasks.FilesToTarredFolder as TarShardVidLists {
    input:
      in_files=flatten(CollectVidsPerSample.vid_lists),
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
  
  # when filtering/sorting/etc, memory usage will likely go up (much of the data will have to
  # be held in memory or disk while working, potentially in a form that takes up more space)
  Float input_size = size([vcf, samples_list], "GiB")
  Float compression_factor = 5.0
  Float base_disk_gb = 5.0
  Float base_mem_gb = 2.0
  RuntimeAttr runtime_default = object {
    mem_gb: base_mem_gb + compression_factor * input_size,
    disk_gb: ceil(base_disk_gb + input_size * (2.0 + 2.0 * compression_factor)),
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
    
    # For purposes of memory, cut vcf to samples of interest
    zcat ~{vcf} > uncompressed.vcf
    rm ~{vcf}
    grep -B9999999999 -m1 -Ev "^#" uncompressed.vcf  | sed '$ d' > header.vcf \
      || cp uncompressed.vcf header.vcf
    N_HEADER=$(wc -l < header.vcf)
    IDXS=$( grep -Ev "^##" header.vcf \
            | sed 's/\t/\n/g' \
            | awk -v OFS="\t" '{ print $1, NR }' \
            | fgrep -wf ~{samples_list} \
            | cut -f2 \
            | sort -nk1,1 \
            | uniq \
            | paste -s -d, \
            | awk '{ print "1-9,"$1 }' \
            || printf "")

    if [ -z "$IDXS" ]; then
      # nothing to find, make empty dir for output glob to look at
      mkdir -p "~{prefix}_perSample_VIDs"
    else
      cut -f"$IDXS" uncompressed.vcf \
        | vcftools --vcf - --stdout --non-ref-ac-any 1 --recode --recode-INFO-all \
        | grep -v "^#" \
        | cut -f3 \
        > VIDs_to_keep.list

      # Gather list of VIDs and genotypes per sample
      {
        grep -Ev "^##" header.vcf | cut -f"$IDXS";
        tail -n+$((N_HEADER + 1)) uncompressed.vcf | cut -f"$IDXS" | fgrep -wf VIDs_to_keep.list;
      } \
        | /opt/sv-pipeline/scripts/vcf_qc/perSample_vcf_parsing_helper.R \
            /dev/stdin \
            ~{samples_list} \
            "~{prefix}_perSample_VIDs/"

      # Gzip all output lists
      for FILE in ~{prefix}_perSample_VIDs/*.VIDs_genotypes.txt; do
        gzip -f "$FILE"
        rm -f "$FILE"
      done

      # Check if one file per sample is present
      NUM_GENOTYPE_FILES=$(find "~{prefix}_perSample_VIDs/" -name "*.VIDs_genotypes.txt.gz" | wc -l)
      NUM_SAMPLES=$(sort ~{samples_list} | uniq | wc -l)
      if [ $NUM_GENOTYPE_FILES -lt $NUM_SAMPLES ]; then
        echo "ERROR IN TASK collect_VIDs_perSample! FEWER PER-SAMPLE GENOTYPE FILES LOCATED THAN NUMBER OF INPUT SAMPLES"
        exit 1
      fi
    fi
  >>>

  output {
    Array[File] vid_lists = glob("~{prefix}_perSample_VIDs/*.VIDs_genotypes.txt.gz")
  }
}
