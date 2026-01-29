version 1.0

import "Structs.wdl"

workflow BatchEvidenceMerging {
  input {
    Array[String] samples
    Array[File?]? BAF_files
    Array[File] PE_files
    Array[File] SR_files
    Array[File]? SD_files
    File? sd_locs_vcf
    File reference_dict
    File primary_contigs_fai
    Boolean subset_primary_contigs  # If true, input PE/SR/BAF files will be subsetted to primary contigs only
    Boolean rename_samples  # If true, rename samples to IDs in the input array
    String batch
    String gatk_docker
    RuntimeAttr? runtime_attr_override
  }

  call MergeEvidence as MergeSREvidence {
    input:
      files = SR_files,
      batch = batch,
      evidence = "sr",
      samples = samples,
      reference_dict = reference_dict,
      primary_contigs_fai=primary_contigs_fai,
      subset_primary_contigs=subset_primary_contigs,
      rename_samples=rename_samples,
      gatk_docker = gatk_docker,
      runtime_attr_override = runtime_attr_override
  }

  call MergeEvidence as MergePEEvidence {
    input:
      files = PE_files,
      batch = batch,
      evidence = "pe",
      samples = samples,
      reference_dict = reference_dict,
      primary_contigs_fai=primary_contigs_fai,
      subset_primary_contigs=subset_primary_contigs,
      rename_samples=rename_samples,
      gatk_docker = gatk_docker,
      runtime_attr_override = runtime_attr_override
  }

  if (defined(BAF_files)) {
    call MergeEvidence as MergeBAFEvidence {
      input:
        files = select_all(select_first([BAF_files])),
        batch = batch,
        evidence = "baf",
        samples = samples,
        reference_dict = reference_dict,
        primary_contigs_fai=primary_contigs_fai,
        subset_primary_contigs=subset_primary_contigs,
        rename_samples=rename_samples,
        gatk_docker = gatk_docker,
        runtime_attr_override = runtime_attr_override
    }
  }
  if (!defined(BAF_files)) {
    call SDtoBAF {
      input:
        SD_files = select_first([SD_files]),
        sd_locs_vcf = select_first([sd_locs_vcf]),
        batch = batch,
        samples = samples,
        rename_samples=rename_samples,
        reference_dict = reference_dict,
        gatk_docker = gatk_docker,
        runtime_attr_override = runtime_attr_override
    }
  }

  output {
    File merged_BAF = select_first([MergeBAFEvidence.out, SDtoBAF.out])
    File merged_BAF_index = select_first([MergeBAFEvidence.out_index, SDtoBAF.out_index])
    File merged_SR = MergeSREvidence.out
    File merged_SR_index = MergeSREvidence.out_index
    File merged_PE = MergePEEvidence.out
    File merged_PE_index = MergePEEvidence.out_index
  }
}

task MergeEvidence {
  input {
    Array[File] files
    String batch
    String evidence
    Array[String] samples
    File primary_contigs_fai
    Boolean subset_primary_contigs
    Boolean rename_samples
    File reference_dict
    String gatk_docker
    RuntimeAttr? runtime_attr_override
  }

  Float file_size = size(files, "GiB")
  Float subset_disk = if (subset_primary_contigs || rename_samples) then file_size else 0
  Int disk_size = 10 + ceil(file_size * 3 + subset_disk)
  Int java_heap_size_mb = round(42.0 * length(files) + 1024.0)
  Float mem_size_gb = java_heap_size_mb / 1024.0 + 2.5

  RuntimeAttr default_attr = object {
    cpu_cores: 2,
    mem_gb: mem_size_gb,
    disk_gb: disk_size,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File out = "~{batch}.~{evidence}.txt.gz"
    File out_index = "~{batch}.~{evidence}.txt.gz.tbi"
  }
  command <<<

    set -euo pipefail
    ulimit -n 100000

    mv ~{write_lines(files)} evidence.list
    mv ~{write_lines(samples)} samples.list

    # For legacy evidence files that were not dictionary sorted, removing non-primary contigs fixes the GATK error
    # BAF records will be deduplicated by contig/coordinate
    if ~{subset_primary_contigs} || ~{rename_samples}; then
      mkdir evidence
      touch evidence.tmp
      cut -f1 ~{primary_contigs_fai} > contigs.list
      while read fil sample; do
        FILENAME=$(basename $fil)
        OUT="evidence/$FILENAME"
        if [[ "~{evidence}" == "pe" ]]; then
          zcat $fil \
            | awk -F'\t' -v OFS='\t' -v SAMPLE="$sample" '!second_file{chroms[$1]; next} {~{if subset_primary_contigs then "if ($1 in chroms && $4 in chroms)" else ""} print ~{if rename_samples then "$1,$2,$3,$4,$5,$6,SAMPLE" else "$0"} }' contigs.list second_file=1 - \
            | bgzip > $OUT
        elif [[ "~{evidence}" == "sr" ]]; then
          zcat $fil \
            | awk -F'\t' -v OFS='\t' -v SAMPLE="$sample" '!second_file{chroms[$1]; next} {~{if subset_primary_contigs then "if ($1 in chroms)" else ""} print ~{if rename_samples then "$1,$2,$3,$4,SAMPLE" else "$0"} }' contigs.list second_file=1 - \
            | bgzip > $OUT
        else
          # baf - also uniquify records from old files
          zcat $fil \
            | awk -F'\t' -v OFS='\t' -v SAMPLE="$sample" '!second_file{chroms[$1]; next} {~{if subset_primary_contigs then "if ($1 in chroms)" else ""} print ~{if rename_samples then "$1,$2,$3,SAMPLE" else "$0"} }' contigs.list second_file=1 - \
            | awk -F'\t' -v OFS='\t' '!_[$1"_"$2]++' - \
            | bgzip > $OUT
        fi
        echo "$OUT" >> evidence.tmp
      done < <(paste evidence.list samples.list)
      mv evidence.tmp evidence.list
    fi

    awk '/txt\.gz$/' evidence.list | while read fil; do
      tabix -f -0 -s1 -b2 -e2 $fil
    done

    /gatk/gatk --java-options "-Xmx~{java_heap_size_mb}m" PrintSVEvidence -F evidence.list --sample-names samples.list --sequence-dictionary ~{reference_dict} -O "~{batch}.~{evidence}.txt.gz"

    tabix -f -0 -s1 -b2 -e2 "~{batch}.~{evidence}.txt.gz"
  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: gatk_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    noAddress: true
  }
}

task SDtoBAF {
  input {
    Array[File] SD_files
    File sd_locs_vcf
    String batch
    Array[String] samples
    Boolean rename_samples
    File reference_dict
    Float min_het_probability = 0.05
    String gatk_docker
    RuntimeAttr? runtime_attr_override
  }

  Float file_size = size(SD_files, "GiB")
  Float subset_disk = if rename_samples then file_size else 0
  Int disk_size = ceil(10 + file_size * 3 + subset_disk)
  Int java_heap_size_mb = round(42.0 * length(SD_files) + 1024.0)
  Float mem_size_gb = java_heap_size_mb / 1024.0 + 2.5

  RuntimeAttr default_attr = object {
    cpu_cores: 2,
    mem_gb: mem_size_gb,
    disk_gb: disk_size,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File out = "~{batch}.baf.txt.gz"
    File out_index = "~{batch}.baf.txt.gz.tbi"
  }
  command <<<

    set -euo pipefail
    ulimit -n 100000

    mv ~{write_lines(SD_files)} inputs.list
    mv ~{write_lines(samples)} samples.list

    # Rename samples
    if ~{rename_samples}; then
      mkdir evidence
      touch evidence.tmp
      while read fil sample; do
        FILENAME=$(basename $fil)
        OUT="evidence/$FILENAME"
        zcat $fil \
          | awk -F'\t' -v OFS='\t' -v SAMPLE="$sample" '{print $1,$2,SAMPLE,$4,$5,$6,$7}' - \
          | bgzip > $OUT
        echo "$OUT" >> evidence.tmp
      done < <(paste inputs.list samples.list)
      mv evidence.tmp inputs.list
    fi

    awk '/txt\.gz$/' inputs.list | while read fil; do
      tabix -f -s1 -b2 -e2 $fil
    done

    /gatk/gatk --java-options "-Xmx~{java_heap_size_mb}m" SiteDepthtoBAF \
        -F inputs.list \
        --sample-names samples.list \
        --sequence-dictionary "~{reference_dict}" \
        --baf-sites-vcf "~{sd_locs_vcf}" \
        --min-het-probability "~{min_het_probability}" \
        -O "~{batch}.baf.txt.gz"

  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: gatk_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    noAddress: true
  }
}
