version 1.0

import "Structs.wdl"

workflow MergeDepth {
  input {
    Array[String] samples
    Array[File] genotyped_segments_vcfs
    Array[File] contig_ploidy_calls
    File std_cnmops_del
    File std_cnmops_dup
    File large_cnmops_del
    File large_cnmops_dup
    String batch
    String sv_base_mini_docker
    String sv_pipeline_docker
    Int gcnv_qs_cutoff
    Float? defragment_max_dist
    RuntimeAttr? runtime_attr_merge_sample
    RuntimeAttr? runtime_attr_merge_set
    RuntimeAttr? runtime_attr_convert_gcnv
  }

  scatter (i in range(length(samples))) {
    call GcnvVcfToBed {
      input:
        sample_id = samples[i],
        vcf = genotyped_segments_vcfs[i],
        contig_ploidy_call_tar = contig_ploidy_calls[i],
        sv_pipeline_docker = sv_pipeline_docker,
        qs_cutoff = gcnv_qs_cutoff,
        runtime_attr_override = runtime_attr_convert_gcnv
    }
  }

  scatter (i in range(length(samples))) {
    call MergeSample as MergeSample_del {
      input:
        sample_id = samples[i],
        gcnv = GcnvVcfToBed.del_bed[i],
        cnmops = [std_cnmops_del, large_cnmops_del],
        max_dist = defragment_max_dist,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_merge_sample
    }
  }

  scatter (i in range(length(samples))) {
    call MergeSample as MergeSample_dup {
      input:
        sample_id = samples[i],
        gcnv = GcnvVcfToBed.dup_bed[i],
        cnmops = [std_cnmops_dup, large_cnmops_dup],
        max_dist = defragment_max_dist,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_merge_sample
    }
  }

  call MergeSet as MergeSet_del {
    input:
      beds = MergeSample_del.sample_bed,
      svtype = "DEL",
      batch = batch,
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_merge_set
  }

  call MergeSet as MergeSet_dup {
    input:
      beds = MergeSample_dup.sample_bed,
      svtype = "DUP",
      batch = batch,
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_merge_set
  }
  output{
    File del = MergeSet_del.out
    File del_index = MergeSet_del.out_idx
    File dup = MergeSet_dup.out
    File dup_index = MergeSet_dup.out_idx
  }
}

task GcnvVcfToBed {
  input {
    File vcf
    File contig_ploidy_call_tar
    String sample_id
    String sv_pipeline_docker
    Int qs_cutoff
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 3.75,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File del_bed = "~{sample_id}.del.bed"
    File dup_bed = "~{sample_id}.dup.bed"
  }
  command <<<

    set -e
    tar xzf ~{contig_ploidy_call_tar}
    tabix ~{vcf}
    python /opt/WGD/bin/convert_gcnv.py \
      --cutoff ~{qs_cutoff} \
      contig_ploidy.tsv \
      ~{vcf} \
      ~{sample_id}

  >>>
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

task MergeSample {
  input {
    File gcnv
    Array[File] cnmops
    Float? max_dist
    String sample_id
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 3.75,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File sample_bed = "~{sample_id}.merged.defrag.sorted.bed"
  }

  command <<<

    set -euo pipefail
    zcat ~{sep=" " cnmops} | awk -F "\t" -v OFS="\t" '{if ($5=="~{sample_id}") print}' > cnmops.cnv
    cat ~{gcnv} cnmops.cnv | sort -k1,1V -k2,2n > ~{sample_id}.bed
    bedtools merge -i ~{sample_id}.bed -d 0 -c 4,5,6,7 -o distinct > ~{sample_id}.merged.bed
    /opt/sv-pipeline/00_preprocessing/scripts/defragment_cnvs.py \
      --max-dist ~{if defined(max_dist) then max_dist else "0.25"} ~{sample_id}.merged.bed ~{sample_id}.merged.defrag.bed
    sort -k1,1V -k2,2n ~{sample_id}.merged.defrag.bed > ~{sample_id}.merged.defrag.sorted.bed
    
  >>>
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

task MergeSet {
  input {
    Array[File] beds
    String svtype
    String batch
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 3.75,
    disk_gb: 10,
    boot_disk_gb: 20,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File out = "~{batch}.~{svtype}.bed.gz"
    File out_idx = "~{batch}.~{svtype}.bed.gz.tbi"
  }
  command <<<

    set -euo pipefail

    # TODO: fail fast if localization failed. This is a Cromwell issue.
    while read file; do
      if [ ! -f $file ]; then
        echo "Localization failed: ${file}" >&2
        exit 1
      fi
    done < ~{write_lines(beds)};

    zcat -f ~{sep=' ' beds} \
      | sort -k1,1V -k2,2n \
      | awk -v OFS="\t" -v svtype=~{svtype} -v batch=~{batch} '{$4=batch"_"svtype"_"NR; print}' \
      | cat <(echo -e "#chr\\tstart\\tend\\tname\\tsample\\tsvtype\\tsources") - \
      | bgzip -c > ~{batch}.~{svtype}.bed.gz;
    tabix -p bed ~{batch}.~{svtype}.bed.gz
		
  >>>
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

