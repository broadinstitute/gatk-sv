##########################################################################################

## Base script:   https://portal.firecloud.org/#methods/Talkowski-SV/00_batch_evidence_merging/15/wdl

## Github commit: talkowski-lab/gatk-sv-v1:<ENTER HASH HERE IN FIRECLOUD>

##########################################################################################

version 1.0

import "Structs.wdl"

workflow BincovMatrix {
  input {
    Array[String] samples
    Array[File] count_files
    Array[String]? bincov_matrix_samples
    File? bincov_matrix
    String batch
    Int? binsize
    Int? disk_overhead_gb
    Int? bincov_size_mb
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  call MakeBincovMatrix {
    input:
      samples = samples,
      count_files = count_files,
      bincov_matrix_samples = bincov_matrix_samples,
      bincov_matrix = bincov_matrix,
      batch = batch,
      binsize = binsize,
      disk_overhead_gb = disk_overhead_gb,
      bincov_size_mb = bincov_size_mb,
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_override
  }

  output {
    File merged_bincov = MakeBincovMatrix.merged_bincov
    File merged_bincov_idx = MakeBincovMatrix.merged_bincov_idx
  }
}

task MakeBincovMatrix {
  input {
    Array[String] samples
    Array[File] count_files
    Array[String]? bincov_matrix_samples # Required if bincov_matrix is supplied and vice versa
    File? bincov_matrix
    String batch
    Int? binsize
    Int? disk_overhead_gb
    Int? bincov_size_mb
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  Array[String] all_samples = flatten(select_all([samples, bincov_matrix_samples]))
  Array[File] all_count_files = flatten([count_files, select_all([bincov_matrix])])

  Int overhead_gb = select_first([disk_overhead_gb, 10])
  Int bincov_matrix_size_gb = if defined(bincov_matrix) then ceil(15*size(bincov_matrix, "GB")) else 0
  Int count_files_size_gb = ceil((length(count_files) * select_first([bincov_size_mb, 750])) / 1000)
  Int disk_gb = overhead_gb + bincov_matrix_size_gb + count_files_size_gb

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 3.75, 
    disk_gb: disk_gb,
    boot_disk_gb: 20,
    preemptible_tries: 0,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File merged_bincov = "~{batch}.bincov.bed.gz"
    File merged_bincov_idx = "~{batch}.bincov.bed.gz.tbi"
  }
  command <<<

    set -eu

    # make the CollectReadCounts output consistent with the old bincov code
    # determine what format this is
    firstchar=$(gunzip -c ~{count_files[0]} | head -c 1)
    set -o pipefail
    if [ $firstchar == '@' ]; then
      shift=1  # GATK CollectReadCounts (to convert from 1-based closed intervals)
    else
      shift=0  # bincov sample or matrix
    fi

    # kill the dictionary | kill the header | adjust to bed format: 0-based half-open intervals
    zcat ~{all_count_files[0]} \
      | sed '/^@/d' \
      | sed '/^CONTIG	START	END	COUNT$/d' \
      | sed '/^#/d' \
      | awk -v x="${shift}" 'BEGIN{OFS="\t"}{$2=$2-x; print $1,$2,$3}' > locs

    # determine bin size, and drop all bins not exactly equal to this size
    if ~{!defined(binsize)}; then
      sed -n '1,1000p' locs | awk '{ print $3-$2 }' \
      | sort | uniq -c | sort -nrk1,1 \
      | sed -n '1p' | awk '{ print $2 }' \
      > most_common_binsize.txt
      binsize=$( cat most_common_binsize.txt )
    else
      binsize=~{binsize}
    fi
    awk -v FS="\t" -v b=$binsize '{ if ($3-$2==b) print $0 }' locs > locs2
    mv locs2 locs

    mkdir cargo
    fileNo=0
    for fil in ~{sep=' ' all_count_files}
    do
      set +o pipefail
      firstchar=$(gunzip -c $fil | head -c 1)
      set -o pipefail
      if [ $firstchar == '@' ]; then
        shift=1
      else
        shift=0
      fi
      zcat $fil \
        | sed '/^@/d' \
        | sed '/^CONTIG	START	END	COUNT$/d' \
        | sed '/^#/d' \
        | awk -v x="${shift}" -v b=$binsize \
          'BEGIN{OFS="\t"}{$2=$2-x; if ($3-$2==b) print $0}' > fil.bincov.bed
      if ! cut -f1-3 fil.bincov.bed | cmp locs; then
        echo $fil has different intervals than ~{all_count_files[0]}
        exit 1
      fi
      cut -f4- fil.bincov.bed > cargo/`printf "%08d" $fileNo`
      ((++fileNo))
    done

    echo "#Chr	Start	End	~{sep='	' all_samples}" > ~{batch}.bincov.bed
    paste locs cargo/* >> ~{batch}.bincov.bed
    bgzip ~{batch}.bincov.bed
    tabix ~{batch}.bincov.bed.gz

  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + runtime_attr.disk_gb + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_base_mini_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

