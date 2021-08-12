version 1.0

import "Structs.wdl"

workflow EvidenceMerging {
  input {
    Array[String] samples
    Array[File?]? BAF_files
    Array[File] PE_files
    Array[File] SR_files
    File inclusion_bed
    File genome_file
    String batch

    #Single file size estimates
    Int? BAF_size_mb
    Int? PE_size_mb
    Int? SR_size_mb

    Int? disk_overhead_gb   # Fixed extra disk
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_set_sample

    RuntimeAttr? runtime_attr_shard_baf
    RuntimeAttr? runtime_attr_merge_baf
    RuntimeAttr? runtime_attr_shard_pe
    RuntimeAttr? runtime_attr_merge_pe
    RuntimeAttr? runtime_attr_shard_sr
    RuntimeAttr? runtime_attr_merge_sr
  }

  Array[String] contigs = transpose(read_tsv(genome_file))[0]

  if (defined(BAF_files)) {
    Array[File?] BAF_files_value = select_first([BAF_files])
    scatter (i in range(length(samples))) {
      if (defined(BAF_files_value[i])) {
        call SetSampleId as SetSampleIdBAF {
           input:
             file = select_first([BAF_files_value[i]]),
             sample_id = samples[i],
             output_name = "BAF_GatherBatchEvidence.~{samples[i]}.txt.gz",
             sample_column_index = 4,
             sv_base_mini_docker = sv_base_mini_docker,
             runtime_attr_override = runtime_attr_set_sample
        }
      }
    }

    scatter (contig in contigs) {
      call MergeEvidenceFilesByContig as MergeBAFFilesByContig {
          input:
            files = select_all(SetSampleIdBAF.out),
            indexes = select_all(SetSampleIdBAF.out_index),
            contig = contig,
            batch = batch,
            evidence = "BAF",
            inclusion_bed = inclusion_bed,
            disk_size_factor = 5,
            sv_base_mini_docker = sv_base_mini_docker,
            runtime_attr_override = runtime_attr_shard_baf
      }
    }

    call MergeEvidenceShards as MergeBAFShards {
      input:
        files = MergeBAFFilesByContig.merged,
        batch = batch,
        evidence = "BAF",
        sv_base_mini_docker = sv_base_mini_docker,
        runtime_attr_override = runtime_attr_merge_baf
    }
  }

  scatter (i in range(length(samples))) {
         call SetSampleId as SetSampleIdSR {
             input:
             file = SR_files[i],
             sample_id = samples[i],
             output_name = "SR_GatherBatchEvidence.~{samples[i]}.txt.gz",
             sample_column_index = 5,
             sv_base_mini_docker = sv_base_mini_docker,
             runtime_attr_override = runtime_attr_set_sample
        }

         call SetSampleId as SetSampleIdPE {
             input:
             file = PE_files[i],
             sample_id = samples[i],
             output_name = "PE_GatherBatchEvidence.~{samples[i]}.txt.gz",
             sample_column_index = 7,
             sv_base_mini_docker = sv_base_mini_docker,
             runtime_attr_override = runtime_attr_set_sample
        }
  }

  scatter (contig in contigs) {
    call MergeEvidenceFilesByContig as MergeSRFilesByContig {
        input:
          files = SetSampleIdSR.out,
          indexes = SetSampleIdSR.out_index,
          contig = contig,
          batch = batch,
          evidence = "SR",
          inclusion_bed = inclusion_bed,
          disk_size_factor = 5,
          sv_base_mini_docker = sv_base_mini_docker,
          runtime_attr_override = runtime_attr_shard_sr
    }

    call MergeEvidenceFilesByContig as MergePEFilesByContig {
        input:
          files = SetSampleIdPE.out,
          indexes = SetSampleIdPE.out_index,
          contig = contig,
          batch = batch,
          evidence = "PE",
          inclusion_bed = inclusion_bed,
          disk_size_factor = 5,
          sv_base_mini_docker = sv_base_mini_docker,
          runtime_attr_override = runtime_attr_shard_pe
    }
  }

  call MergeEvidenceShards as MergeSRShards {
    input:
      files = MergeSRFilesByContig.merged,
      batch = batch,
      evidence = "SR",
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_merge_sr
  }

  call MergeEvidenceShards as MergePEShards {
    input:
      files = MergePEFilesByContig.merged,
      batch = batch,
      evidence = "PE",
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_merge_pe
  }

  output {
    File? merged_BAF = MergeBAFShards.out
    File? merged_BAF_idx = MergeBAFShards.out_idx
    File merged_SR = MergeSRShards.out
    File merged_SR_idx = MergeSRShards.out_idx
    File merged_PE = MergePEShards.out
    File merged_PE_idx = MergePEShards.out_idx
  }
}

task SetSampleId {
  input {
    File file
    String sample_id
    String output_name
    Int sample_column_index
    String sv_base_mini_docker
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

  output {
    File out = "~{output_name}"
    File out_index = "~{output_name}.tbi"
  }
  command <<<

    set -euo pipefail
    zcat ~{file} | awk -F "\t" -v OFS="\t" '$~{sample_column_index}="~{sample_id}"' | bgzip -c > ~{output_name}
    tabix -s1 -b2 -e2 ~{output_name}

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

task MergeEvidenceFilesByContig {
  input {
    Array[File] files
    Array[File] indexes
    String contig
    String batch
    String evidence
    File inclusion_bed
    Float disk_size_factor
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  Int disk_size_gb = 10 + ceil(size(files, "GB") * disk_size_factor)

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 3.75,
    disk_gb: disk_size_gb,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File merged = "~{batch}.~{evidence}.~{contig}.txt.gz"
  }
  command <<<

    set -euo pipefail

    mkdir data
    awk -F "\t" '{if ($1=="~{contig}") print}' ~{inclusion_bed} > regions.bed
    while read file; do
      filename=`basename $file`
      tabix -h -R regions.bed $file > data/$filename.txt
    done < ~{write_lines(files)}

    mkdir tmp
    sort -m -k1,1V -k2,2n -T tmp data/*.txt | bgzip -c > ~{batch}.~{evidence}.~{contig}.txt.gz
  
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

task MergeEvidenceShards {
  input {
    Array[File] files
    String batch
    String evidence
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  Int disk_size = 10 + ceil(size(files, "GiB") * 2)

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 3.75,
    disk_gb: disk_size,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File out = "~{batch}.~{evidence}.txt.gz"
    File out_idx = "~{batch}.~{evidence}.txt.gz.tbi"
  }
  command <<<

    set -euo pipefail
    cat ~{sep=" " files} > ~{batch}.~{evidence}.txt.gz
    tabix -f -s1 -b2 -e2 ~{batch}.~{evidence}.txt.gz

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