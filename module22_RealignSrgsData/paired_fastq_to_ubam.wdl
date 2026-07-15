version 1.0

workflow ConvertPairedFastQsToUnmappedBamWf {

  input {
    String        sample_name
    File          fastq_1
    File          fastq_2
    Array[String] library_names
    Array[String] platform_names
    Array[String] sequencing_centers
    Array[String]? run_dates

    File   ref_dict

    String samtools_docker = "staphb/samtools:1.17"
    String gatk_docker     = "broadinstitute/gatk:4.5.0.0"

    Int additional_disk_space_gb = 50
    Int machine_mem_gb           = 7
    Int preemptible_attempts     = 3
  }

  call SplitPairedFastqByReadGroupAndBuildMetadata {
    input:
      sample_name              = sample_name,
      fastq_1                  = fastq_1,
      fastq_2                  = fastq_2,
      library_names            = library_names,
      platform_names           = platform_names,
      sequencing_centers       = sequencing_centers,
      run_dates                = run_dates,
      ref_dict                 = ref_dict,
      docker                   = samtools_docker,
      additional_disk_space_gb = additional_disk_space_gb,
      machine_mem_gb           = machine_mem_gb,
      preemptible_attempts     = preemptible_attempts
  }

  scatter (i in range(length(SplitPairedFastqByReadGroupAndBuildMetadata.platform_units))) {
    call ImportReadGroupFastqPairToUbam {
      input:
        split_fastq_1             = SplitPairedFastqByReadGroupAndBuildMetadata.split_fastq_1s[i],
        split_fastq_2             = SplitPairedFastqByReadGroupAndBuildMetadata.split_fastq_2s[i],
        readgroup_name            = SplitPairedFastqByReadGroupAndBuildMetadata.readgroup_names[i],
        sample_name               = sample_name,
        library_name              = SplitPairedFastqByReadGroupAndBuildMetadata.effective_library_names[i],
        platform_unit             = SplitPairedFastqByReadGroupAndBuildMetadata.platform_units[i],
        run_date                  = SplitPairedFastqByReadGroupAndBuildMetadata.effective_run_dates[i],
        platform_name             = SplitPairedFastqByReadGroupAndBuildMetadata.effective_platform_names[i],
        sequencing_center         = SplitPairedFastqByReadGroupAndBuildMetadata.effective_sequencing_centers[i],
        readgroup_index           = i,
        docker                    = samtools_docker,
        additional_disk_space_gb  = additional_disk_space_gb,
        machine_mem_gb            = machine_mem_gb,
        preemptible_attempts      = preemptible_attempts
    }
  }

  call MergeReadGroupUbams {
    input:
      sample_name              = sample_name,
      readgroup_ubams          = ImportReadGroupFastqPairToUbam.output_unmapped_bam,
      header_sam               = SplitPairedFastqByReadGroupAndBuildMetadata.full_header_sam,
      docker                   = samtools_docker,
      additional_disk_space_gb = additional_disk_space_gb,
      machine_mem_gb           = machine_mem_gb,
      preemptible_attempts     = preemptible_attempts
  }

  call SortSam {
    input:
      input_bam                = MergeReadGroupUbams.output_unmapped_bam,
      sample_name              = sample_name,
      docker                   = gatk_docker,
      machine_mem_gb           = machine_mem_gb,
      additional_disk_space_gb = additional_disk_space_gb,
      preemptible_attempts     = preemptible_attempts
  }

  output {
    File output_unmapped_bam = SortSam.sorted_bam
    Array[String] readgroup_names = SplitPairedFastqByReadGroupAndBuildMetadata.readgroup_names
    Array[String] platform_units = SplitPairedFastqByReadGroupAndBuildMetadata.platform_units
    Array[String] effective_run_dates = SplitPairedFastqByReadGroupAndBuildMetadata.effective_run_dates
  }
}

task SplitPairedFastqByReadGroupAndBuildMetadata {
  input {
    String        sample_name
    File          fastq_1
    File          fastq_2
    Array[String] library_names
    Array[String] platform_names
    Array[String] sequencing_centers
    Array[String]? run_dates
    File          ref_dict

    Int    additional_disk_space_gb = 50
    Int    machine_mem_gb           = 7
    Int    preemptible_attempts     = 3
    String docker
  }

  Array[String] provided_run_dates = select_first([run_dates, []])
  Int disk_space_gb = ceil((size(fastq_1, "GB") + size(fastq_2, "GB")) * 6) + additional_disk_space_gb

  command <<<
    set -euo pipefail

    mkdir -p split_fastq

    python3 <<'PY'
import gzip
import re

fastq_1 = "~{fastq_1}"
fastq_2 = "~{fastq_2}"

def open_fastq(path, mode):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)

def parse_platform_unit(header):
    token = header.strip().split()[0]
    if token.startswith("@"):
        token = token[1:]
    fields = token.split(":")
    if len(fields) < 4:
        raise RuntimeError("FASTQ header does not match Illumina format: " + header.strip())
    return fields[2] + "." + fields[3]

writers = {}
counts = {}
paths = {}

try:
    with open_fastq(fastq_1, "rt") as r1, open_fastq(fastq_2, "rt") as r2:
        while True:
            rec1 = [r1.readline() for _ in range(4)]
            rec2 = [r2.readline() for _ in range(4)]

            if not rec1[0] and not rec2[0]:
                break
            if bool(rec1[0]) != bool(rec2[0]):
                raise RuntimeError("R1 and R2 have different number of records")
            if any(x == "" for x in rec1 + rec2):
                raise RuntimeError("Truncated FASTQ record encountered")

            pu1 = parse_platform_unit(rec1[0])
            pu2 = parse_platform_unit(rec2[0])
            if pu1 != pu2:
                raise RuntimeError("R1/R2 read-group mismatch: " + pu1 + " vs " + pu2)

            if pu1 not in writers:
                safe_key = re.sub(r"[^A-Za-z0-9._-]", "_", pu1)
                fq1_out = "split_fastq/" + safe_key + "_R1.fastq.gz"
                fq2_out = "split_fastq/" + safe_key + "_R2.fastq.gz"
                writers[pu1] = (gzip.open(fq1_out, "wt"), gzip.open(fq2_out, "wt"))
                paths[pu1] = (fq1_out, fq2_out)
                counts[pu1] = 0

            w1, w2 = writers[pu1]
            w1.writelines(rec1)
            w2.writelines(rec2)
            counts[pu1] += 1
finally:
    for w1, w2 in writers.values():
        w1.close()
        w2.close()

if not counts:
    raise RuntimeError("No reads found in FASTQ files")

pu_keys = sorted(counts.keys())
with open("platform_units.txt", "w") as f:
    for pu in pu_keys:
        f.write(pu + "\n")
with open("split_fastq_1.list", "w") as f1:
    for pu in pu_keys:
        f1.write(paths[pu][0] + "\n")
with open("split_fastq_2.list", "w") as f2:
    for pu in pu_keys:
        f2.write(paths[pu][1] + "\n")
PY

    mapfile -t PLATFORM_UNITS < platform_units.txt
    mapfile -t LIB_NAMES < ~{write_lines(library_names)}
    mapfile -t PLAT_NAMES < ~{write_lines(platform_names)}
    mapfile -t SEQ_CENTERS < ~{write_lines(sequencing_centers)}
    mapfile -t RUN_DATES < ~{write_lines(provided_run_dates)}

    n_rg=${#PLATFORM_UNITS[@]}
    if [ "$n_rg" -eq 0 ]; then
      echo "ERROR: No read groups detected from FASTQ headers." >&2
      exit 1
    fi

    check_array_length() {
      local arr_name="$1"
      declare -n arr_ref="$arr_name"
      local n="${#arr_ref[@]}"
      if [ "$n" -ne 1 ] && [ "$n" -ne "$n_rg" ]; then
        echo "ERROR: $arr_name length ($n) must be either 1 or number of read groups ($n_rg)." >&2
        exit 1
      fi
    }

    pick_value() {
      local arr_name="$1"
      local idx="$2"
      declare -n arr_ref="$arr_name"
      local n="${#arr_ref[@]}"
      if [ "$n" -eq 1 ]; then
        printf "%s" "${arr_ref[0]}"
      else
        printf "%s" "${arr_ref[$idx]}"
      fi
    }

    check_array_length LIB_NAMES
    check_array_length PLAT_NAMES
    check_array_length SEQ_CENTERS
    if [ "${#RUN_DATES[@]}" -ne 0 ]; then
      check_array_length RUN_DATES
    fi

    JOB_DATE=$(date -u +"%Y-%m-%dT%H:%M:%SZ")

    printf "@HD\tVN:1.6\tSO:queryname\n" > new_header.sam
    grep "^@SQ" ~{ref_dict} >> new_header.sam
    : > readgroup_names.txt
    : > effective_run_dates.txt
    : > effective_library_names.txt
    : > effective_platform_names.txt
    : > effective_sequencing_centers.txt

    for (( i=0; i<n_rg; i++ )); do
      RG_ID="~{sample_name}.${PLATFORM_UNITS[$i]}"
      LB=$(pick_value LIB_NAMES "$i")
      PL=$(pick_value PLAT_NAMES "$i")
      CN=$(pick_value SEQ_CENTERS "$i")
      if [ "${#RUN_DATES[@]}" -eq 0 ]; then
        DT="$JOB_DATE"
      else
        DT=$(pick_value RUN_DATES "$i")
      fi

      printf "@RG\tID:%s\tSM:%s\tLB:%s\tPU:%s\tDT:%s\tPL:%s\tCN:%s\n" \
        "$RG_ID" \
        "~{sample_name}" \
        "$LB" \
        "${PLATFORM_UNITS[$i]}" \
        "$DT" \
        "$PL" \
        "$CN" \
        >> new_header.sam

      echo "$RG_ID" >> readgroup_names.txt
      echo "$DT" >> effective_run_dates.txt
      echo "$LB" >> effective_library_names.txt
      echo "$PL" >> effective_platform_names.txt
      echo "$CN" >> effective_sequencing_centers.txt
    done
  >>>

  output {
    Array[String] readgroup_names = read_lines("readgroup_names.txt")
    Array[String] platform_units = read_lines("platform_units.txt")
    Array[String] effective_run_dates = read_lines("effective_run_dates.txt")
    Array[String] effective_library_names = read_lines("effective_library_names.txt")
    Array[String] effective_platform_names = read_lines("effective_platform_names.txt")
    Array[String] effective_sequencing_centers = read_lines("effective_sequencing_centers.txt")
    Array[File] split_fastq_1s = read_lines("split_fastq_1.list")
    Array[File] split_fastq_2s = read_lines("split_fastq_2.list")
    File full_header_sam = "new_header.sam"
  }

  runtime {
    docker:      docker
    memory:      machine_mem_gb + " GB"
    disks:       "local-disk " + disk_space_gb + " HDD"
    preemptible: preemptible_attempts
    maxRetries:  preemptible_attempts
  }
}

task ImportReadGroupFastqPairToUbam {
  input {
    File split_fastq_1
    File split_fastq_2
    String readgroup_name
    String sample_name
    String library_name
    String platform_unit
    String run_date
    String platform_name
    String sequencing_center
    Int readgroup_index

    Int    additional_disk_space_gb = 50
    Int    machine_mem_gb           = 7
    Int    preemptible_attempts     = 3
    String docker
  }

  Int disk_space_gb = ceil((size(split_fastq_1, "GB") + size(split_fastq_2, "GB")) * 4) + additional_disk_space_gb

  command <<<
    set -euo pipefail

    RG_FIELD=$(printf "ID:%s\tSM:%s\tLB:%s\tPU:%s\tDT:%s\tPL:%s\tCN:%s" \
      "~{readgroup_name}" \
      "~{sample_name}" \
      "~{library_name}" \
      "~{platform_unit}" \
      "~{run_date}" \
      "~{platform_name}" \
      "~{sequencing_center}")

    samtools import \
      -1 ~{split_fastq_1} \
      -2 ~{split_fastq_2} \
      -r "$RG_FIELD" \
      -O BAM \
      -o rg.~{readgroup_index}.unmapped.bam
  >>>

  output {
    File output_unmapped_bam = "rg.~{readgroup_index}.unmapped.bam"
  }

  runtime {
    docker:      docker
    memory:      machine_mem_gb + " GB"
    disks:       "local-disk " + disk_space_gb + " HDD"
    preemptible: preemptible_attempts
    maxRetries:  preemptible_attempts
  }
}

task MergeReadGroupUbams {
  input {
    String      sample_name
    Array[File] readgroup_ubams
    File        header_sam

    Int    additional_disk_space_gb = 50
    Int    machine_mem_gb           = 7
    Int    preemptible_attempts     = 3
    String docker
  }

  Int disk_space_gb = ceil(size(readgroup_ubams, "GB") * 4) + additional_disk_space_gb

  command <<<
    set -euo pipefail

    if [ "~{length(readgroup_ubams)}" -eq 1 ]; then
      cp ~{readgroup_ubams[0]} merged.raw.unmapped.bam
    else
      samtools merge -f -O BAM merged.raw.unmapped.bam ~{sep=' ' readgroup_ubams}
    fi

    samtools reheader ~{header_sam} merged.raw.unmapped.bam > ~{sample_name}.unmapped.bam
  >>>

  output {
    File output_unmapped_bam = "~{sample_name}.unmapped.bam"
  }

  runtime {
    docker:      docker
    memory:      machine_mem_gb + " GB"
    disks:       "local-disk " + disk_space_gb + " HDD"
    preemptible: preemptible_attempts
    maxRetries:  preemptible_attempts
  }
}

task SortSam {
  input {
    File   input_bam
    String sample_name
    String docker
    Int    machine_mem_gb           = 7
    Int    additional_disk_space_gb = 50
    Int    preemptible_attempts     = 3
  }

  Int command_mem_gb = machine_mem_gb - 1
  Int disk_space_gb  = ceil(size(input_bam, "GB") * 4) + additional_disk_space_gb

  command <<<
    set -euo pipefail

    mkdir -p tmp

    gatk --java-options "-Xmx~{command_mem_gb}g -Djava.io.tmpdir=./tmp" \
      SortSam \
      -I ~{input_bam} \
      -O ~{sample_name}.query_sorted.unmapped.bam \
      --SORT_ORDER queryname \
      --TMP_DIR ./tmp
  >>>

  output {
    File sorted_bam = "~{sample_name}.query_sorted.unmapped.bam"
  }

  runtime {
    docker:      docker
    memory:      machine_mem_gb + " GB"
    disks:       "local-disk " + disk_space_gb + " HDD"
    preemptible: preemptible_attempts
    maxRetries:  preemptible_attempts
  }
}
