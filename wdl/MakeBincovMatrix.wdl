version 1.0

import "Structs.wdl"

workflow MakeBincovMatrix {
  input {
    Array[String] samples
    Array[File] count_files
    Array[String]? bincov_matrix_samples
    File? bincov_matrix
    String batch
    Int? binsize
    Int? disk_overhead_gb
    String sv_base_mini_docker
    String sv_base_docker
    RuntimeAttr? runtime_attr_override
  }

  call SetBins {
    input:
      count_file = all_count_files[0],
      binsize = binsize,
      bincov_matrix_samples = bincov_matrix_samples,
      sv_base_mini_docker = sv_base_mini_docker,
      disk_overhead_gb = disk_overhead_gb,
      runtime_attr_override = runtime_attr_override
  }
  if(defined(bincov_matrix_samples)) {
    String bincov_matrix_header = read_lines(SetBins.bincov_matrix_header_file)[0]
  }

  Array[String]+ all_samples = flatten([samples, select_all([bincov_matrix_header])])
  Array[File]+ all_count_files = flatten([count_files, select_all([bincov_matrix])])

  scatter(i in range(length(all_count_files))) {
    call MakeBincovMatrixColumns {
      input:
        count_file = all_count_files[i],
        sample = all_samples[i],
        binsize = SetBins.out_binsize,
        bin_locs = SetBins.bin_locs,
        disk_overhead_gb = disk_overhead_gb,
        sv_base_mini_docker = sv_base_mini_docker,
        runtime_attr_override = runtime_attr_override
    }
  }

  call ZPaste {
    input:
      column_files = flatten([[SetBins.bin_locs], MakeBincovMatrixColumns.bincov_bed]),
      matrix_file_name = "~{batch}.RD.txt.gz",
      disk_overhead_gb = disk_overhead_gb,
      sv_base_docker = sv_base_docker,
      runtime_attr_override = runtime_attr_override
  }

  output {
    File merged_bincov = ZPaste.matrix_file
    File merged_bincov_idx = ZPaste.matrix_file_idx
  }
}

task SetBins {
  input {
    File count_file
    Int? binsize
    Array[String]? bincov_matrix_samples
    Int? disk_overhead_gb
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  String bincov_header_file_name = "bincov_header_file.tsv"

  String bin_file_name = "locs.bed.gz"
  String binsize_output_file_name = "binsize.txt"

  Float disk_scale_factor = 10.0
  Int disk_gb = select_first([disk_overhead_gb, 10]) + ceil(disk_scale_factor * size(count_file, "GiB"))
  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 2.0,
    disk_gb: disk_gb,
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File bin_locs = bin_file_name
    Int out_binsize = read_int(binsize_output_file_name)
    File bincov_matrix_header_file = bincov_header_file_name
  }

  command <<<
    set -Eeu

    # make the CollectReadCounts output consistent with the old bincov code
    # determine what format this is
    firstchar=$(gunzip -c ~{count_file} | head -c 1)
    set -o pipefail
    if [ $firstchar == '@' ]; then
      shift=1  # GATK CollectReadCounts (to convert from 1-based closed intervals)
    else
      shift=0  # bincov sample or matrix
    fi

    # kill the dictionary | kill the header | adjust to bed format: 0-based half-open intervals
    zcat ~{count_file} \
      | sed '/^@/d' \
      | sed '/^CONTIG	START	END	COUNT$/d' \
      | sed '/^#/d' \
      | awk -v x="${shift}" 'BEGIN{OFS="\t"}{$2=$2-x; print $1,$2,$3}' > tmp_locs

    # determine bin size, and drop all bins not exactly equal to this size
    if ~{defined(binsize)}; then
      # use the provided bin size
      binsize=~{binsize}
    else
      # use the most common bin size from the bins
      binsize=$(
        sed -n '1,1000p' tmp_locs | awk '{ print $3-$2 }' \
        | sort | uniq -c | sort -nrk1,1 \
        | sed -n '1p' | awk '{ print $2 }'
      )
    fi
    # store binsize where cromwell can read it
    echo $binsize > ~{binsize_output_file_name}

    # write final bed file with header, and compress it
    awk -v FS="\t" -v b=$binsize 'BEGIN{ print "#Chr\tStart\tEnd" } { if ($3-$2==b) print $0 }' tmp_locs \
        | bgzip -c \
        > "~{bin_file_name}"

    # if bincov_matrix_samples was passed, convert to tab-separated string
    if ~{defined(bincov_matrix_samples)}; then
      mv "~{write_tsv([select_first([bincov_matrix_samples, ["dummy"]])])}" "~{bincov_header_file_name}"
    else
      touch "~{bincov_header_file_name}"
    fi
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

task MakeBincovMatrixColumns {
  input {
    File count_file
    String sample
    Int binsize
    File bin_locs
    Int? disk_overhead_gb = 10
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }
  Float disk_scale_factor = 10.0
  Int disk_gb = select_first([disk_overhead_gb, 10])+ ceil(disk_scale_factor * (size(count_file, "GiB") + size(bin_locs, "GiB")))
  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 2.0,
    disk_gb: disk_gb,
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  String bincov_file_name = "~{basename(count_file, '.tsv.gz')}.RD.txt.gz"

  output {
    File bincov_bed = bincov_file_name
  }

  command <<<
    set -Eeu
    firstchar=$(gunzip -c "~{count_file}" | head -c 1)
    set -o pipefail
    if [ $firstchar == '@' ]; then
      shift=1
    else
      shift=0
    fi

    TMP_BED="$(basename "~{count_file}").tmp.bed"
    printf "#Chr\tStart\tEnd\t%s\n" "~{sample}" > $TMP_BED
    zcat "~{count_file}" \
      | sed '/^@/d' \
      | sed '/^CONTIG	START	END	COUNT$/d' \
      | sed '/^#/d' \
      | awk -v x=$shift -v b=~{binsize} \
        'BEGIN{OFS="\t"}{$2=$2-x; if ($3-$2==b) print $0}' \
      >> "$TMP_BED"

    if ! cut -f1-3 "$TMP_BED" | cmp <(bgzip -cd "~{bin_locs}"); then
      echo "~{count_file} has different intervals than ~{bin_locs}"
      exit 1
    fi
    cut -f4- "$TMP_BED" | bgzip -c >> "~{bincov_file_name}"
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


task ZPaste {
  input {
    Array[File]+ column_files
    String matrix_file_name
    Int? disk_overhead_gb = 10
    Float mem_overhead_gb = 1.0
    String sv_base_docker
    RuntimeAttr? runtime_attr_override
  }

  # Only compressed files are stored (if localization_optional, then only output file is stored),
  # so this is a reasonably conservative estimate for disk:
  Int disk_gb = select_first([disk_overhead_gb, 10]) + ceil(3.0 * size(column_files, "GiB"))
  # Some memory is used up by the named pipes. Not a lot, but allocate in case the batch is huge:
  Float mem_gb = mem_overhead_gb + 0.003 * length(column_files)
  RuntimeAttr default_attr = object {
    cpu_cores: 4,
    mem_gb: mem_gb,
    disk_gb: disk_gb,
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File matrix_file = matrix_file_name
    File matrix_file_idx = "~{matrix_file_name}.tbi"
  }

  command <<<
    set -Eeu -o pipefail

    # Use loops rather than WDL sep feature in case there are enough samples to exceed bash line limits
    # Use named pipes to stream unzipped column files in memory
    mkdir -p column_file_fifos
    FILE_NUM=0
    while read -r COLUMN_FILE; do
      FIFO=$(printf "column_file_fifos/%08d" $FILE_NUM)
      mkfifo "$FIFO"
      bgzip -@$(nproc) -cd "$COLUMN_FILE" > "$FIFO" &
      ((++FILE_NUM))
    done < ~{write_lines(column_files)}

    # paste unzipped files and compress
    paste column_file_fifos/* | bgzip -@$(nproc) -c > "~{matrix_file_name}"
    tabix -p bed "~{matrix_file_name}"
  >>>

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_base_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}
