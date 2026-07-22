version 1.0

workflow ConvertPairedFastQsToPerReadgroupUbamWf {
  input {
    String sample_name
    File fastq_1
    File fastq_2

    # If provided, skip ExtractUniqueReadGroups and use these as the readgroup units
    Array[String]? read_groups

    Array[String]? library_names
    Array[String]? platform_units
    Array[String]? run_dates
    Array[String]? platform_names
    Array[String]? sequencing_centers

    String samtools_docker = "staphb/samtools:1.17"

    Int additional_disk_space_gb = 50
    Int machine_mem_gb           = 16
    Int machine_cpu_cores        = 8
    Int preemptible_attempts     = 1
  }

  # Only run extraction when read_groups are not supplied
  if (!defined(read_groups)) {
    call ExtractUniqueReadGroups {
      input:
        sample_name              = sample_name,
        fastq_1                  = fastq_1,
        library_names            = library_names,
        platform_units           = platform_units,
        run_dates                = run_dates,
        platform_names           = platform_names,
        sequencing_centers       = sequencing_centers,
        docker                   = samtools_docker,
        additional_disk_space_gb = additional_disk_space_gb,
        machine_mem_gb           = machine_mem_gb,
        machine_cpu_cores        = machine_cpu_cores,
        preemptible_attempts     = preemptible_attempts
    }
  }

  # Resolve the readgroup units: provided input takes precedence over extracted
  Array[String] resolved_readgroup_units = select_first([read_groups,
                                                         ExtractUniqueReadGroups.detected_readgroup_units])

  call SplitPairedFastqByReadGroup {
    input:
      fastq_1                  = fastq_1,
      fastq_2                  = fastq_2,
      readgroup_units          = resolved_readgroup_units,
      docker                   = samtools_docker,
      additional_disk_space_gb = additional_disk_space_gb,
      machine_mem_gb           = machine_mem_gb,
      machine_cpu_cores        = machine_cpu_cores,
      preemptible_attempts     = preemptible_attempts
  }

  # Resolve per-RG metadata: from ExtractUniqueReadGroups when available, otherwise
  # fall back to the resolved_readgroup_units (used as RG names when extraction was skipped)
  Array[String] resolved_readgroup_names    = select_first([ExtractUniqueReadGroups.readgroup_names,
                                                            resolved_readgroup_units])
  Array[String] resolved_effective_lib      = select_first([ExtractUniqueReadGroups.effective_library_names,
                                                            resolved_readgroup_units])
  Array[String] resolved_effective_pu       = select_first([ExtractUniqueReadGroups.effective_platform_units,
                                                            resolved_readgroup_units])
  Array[String] resolved_effective_dt       = select_first([ExtractUniqueReadGroups.effective_run_dates,
                                                            resolved_readgroup_units])
  Array[String] resolved_effective_pl       = select_first([ExtractUniqueReadGroups.effective_platform_names,
                                                            resolved_readgroup_units])
  Array[String] resolved_effective_cn       = select_first([ExtractUniqueReadGroups.effective_sequencing_centers,
                                                            resolved_readgroup_units])

  scatter (i in range(length(SplitPairedFastqByReadGroup.split_fastq_1s))) {
    call ConvertReadGroupFastqToUbam {
      input:
        sample_name                   = sample_name,
        readgroup_name                = resolved_readgroup_names[i],
        effective_library_name        = resolved_effective_lib[i],
        effective_platform_unit       = resolved_effective_pu[i],
        effective_run_date            = resolved_effective_dt[i],
        effective_platform_name       = resolved_effective_pl[i],
        effective_sequencing_center   = resolved_effective_cn[i],
        split_fastq_1                 = SplitPairedFastqByReadGroup.split_fastq_1s[i],
        split_fastq_2                 = SplitPairedFastqByReadGroup.split_fastq_2s[i],
        readgroup_index               = i,
        docker                        = samtools_docker,
        additional_disk_space_gb      = additional_disk_space_gb,
        machine_mem_gb                = machine_mem_gb,
        machine_cpu_cores             = machine_cpu_cores,
        preemptible_attempts          = preemptible_attempts
    }

    call SortReadGroupUbam {
      input:
        input_ubam                    = ConvertReadGroupFastqToUbam.unmapped_ubam,
        readgroup_index               = i,
        docker                        = samtools_docker,
        additional_disk_space_gb      = additional_disk_space_gb,
        machine_mem_gb                = machine_mem_gb,
        machine_cpu_cores             = machine_cpu_cores,
        preemptible_attempts          = preemptible_attempts
    }
  }

  call MergeSortedReadGroupUbams {
    input:
      sample_name               = sample_name,
      input_ubams               = SortReadGroupUbam.sorted_ubam,
      docker                    = samtools_docker,
      additional_disk_space_gb  = additional_disk_space_gb,
      machine_mem_gb            = machine_mem_gb,
      machine_cpu_cores         = machine_cpu_cores,
      preemptible_attempts      = preemptible_attempts
  }

  output {
    Array[File] split_readgroup_fastq_1s   = SplitPairedFastqByReadGroup.split_fastq_1s
    Array[File] split_readgroup_fastq_2s   = SplitPairedFastqByReadGroup.split_fastq_2s
    Array[File] per_readgroup_sorted_ubams = SortReadGroupUbam.sorted_ubam
    File output_united_sorted_ubam         = MergeSortedReadGroupUbams.merged_ubam

    Array[String] resolved_readgroup_units_out = resolved_readgroup_units
    Array[String] readgroup_names              = resolved_readgroup_names
    Array[String] effective_platform_units     = resolved_effective_pu
    Array[String] effective_run_dates          = resolved_effective_dt
  }
}

task ExtractUniqueReadGroups {
  input {
    String sample_name
    File fastq_1

    Array[String]? library_names
    Array[String]? platform_units
    Array[String]? run_dates
    Array[String]? platform_names
    Array[String]? sequencing_centers

    Int additional_disk_space_gb = 50
    Int machine_mem_gb           = 16
    Int machine_cpu_cores        = 8
    Int preemptible_attempts     = 1
    String docker
  }

  Array[String] provided_library_names = select_first([library_names, []])
  Array[String] provided_platform_units = select_first([platform_units, []])
  Array[String] provided_run_dates = select_first([run_dates, []])
  Array[String] provided_platform_names = select_first([platform_names, []])
  Array[String] provided_sequencing_centers = select_first([sequencing_centers, []])
  Int disk_space_gb = ceil(size(fastq_1, "GB") * 4) + additional_disk_space_gb

  command <<<
    set -euo pipefail

    mapfile -t LIB_NAMES < ~{write_lines(provided_library_names)}
    mapfile -t PU_INPUTS < ~{write_lines(provided_platform_units)}
    mapfile -t RUN_DATES < ~{write_lines(provided_run_dates)}
    mapfile -t PLAT_NAMES < ~{write_lines(provided_platform_names)}
    mapfile -t SEQ_CENTERS < ~{write_lines(provided_sequencing_centers)}

    FASTQ1="~{fastq_1}"
    if [[ "$FASTQ1" == *.gz ]]; then
      gzip -cd "$FASTQ1"
    else
      cat "$FASTQ1"
    fi | awk 'NR % 4 == 1 {
      h=$1
      sub(/^@/, "", h)
      n=split(h, a, ":")
      if (n < 4) {
        print "ERROR: FASTQ header does not match Illumina format: " $0 > "/dev/stderr"
        exit 1
      }
      print a[3] "." a[4]
    }' | sort -u > detected_readgroup_units.txt

    n_rg=$(wc -l < detected_readgroup_units.txt)
    if [ "$n_rg" -eq 0 ]; then
      echo "ERROR: no read groups found in FASTQ header." >&2
      exit 1
    fi

    check_array_length() {
      local arr_name="$1"
      declare -n arr_ref="$arr_name"
      local n="${#arr_ref[@]}"
      if [ "$n" -ne 0 ] && [ "$n" -ne 1 ] && [ "$n" -ne "$n_rg" ]; then
        echo "ERROR: $arr_name length ($n) must be 0, 1, or number of read groups ($n_rg)." >&2
        exit 1
      fi
    }

    check_array_length LIB_NAMES
    check_array_length PU_INPUTS
    check_array_length RUN_DATES
    check_array_length PLAT_NAMES
    check_array_length SEQ_CENTERS

    pick_or_default() {
      local arr_name="$1"
      local idx="$2"
      local default_val="$3"
      declare -n arr_ref="$arr_name"
      local n="${#arr_ref[@]}"
      if [ "$n" -eq 0 ]; then
        printf "%s" "$default_val"
      elif [ "$n" -eq 1 ]; then
        printf "%s" "${arr_ref[0]}"
      else
        printf "%s" "${arr_ref[$idx]}"
      fi
    }

    JOB_DATE=$(date -u +"%Y-%m-%dT%H:%M:%SZ")
    : > readgroup_names.txt
    : > effective_library_names.txt
    : > effective_platform_units.txt
    : > effective_run_dates.txt
    : > effective_platform_names.txt
    : > effective_sequencing_centers.txt

    mapfile -t DETECTED_UNITS < detected_readgroup_units.txt
    for (( i=0; i<n_rg; i++ )); do
      DETECTED_UNIT="${DETECTED_UNITS[$i]}"
      echo "~{sample_name}.${DETECTED_UNIT}" >> readgroup_names.txt
      echo "$(pick_or_default LIB_NAMES "$i" "~{sample_name}")" >> effective_library_names.txt
      echo "$(pick_or_default PU_INPUTS "$i" "ILLUMINA")" >> effective_platform_units.txt
      echo "$(pick_or_default RUN_DATES "$i" "$JOB_DATE")" >> effective_run_dates.txt
      echo "$(pick_or_default PLAT_NAMES "$i" "ILLUMINA")" >> effective_platform_names.txt
      echo "$(pick_or_default SEQ_CENTERS "$i" "ILLUMINA")" >> effective_sequencing_centers.txt
    done
  >>>

  output {
    Array[String] detected_readgroup_units = read_lines("detected_readgroup_units.txt")
    Array[String] readgroup_names = read_lines("readgroup_names.txt")
    Array[String] effective_library_names = read_lines("effective_library_names.txt")
    Array[String] effective_platform_units = read_lines("effective_platform_units.txt")
    Array[String] effective_run_dates = read_lines("effective_run_dates.txt")
    Array[String] effective_platform_names = read_lines("effective_platform_names.txt")
    Array[String] effective_sequencing_centers = read_lines("effective_sequencing_centers.txt")
  }

  runtime {
    docker:      docker
    cpu:         machine_cpu_cores
    memory:      machine_mem_gb + " GB"
    disks:       "local-disk " + disk_space_gb + " HDD"
    preemptible: preemptible_attempts
    maxRetries:  preemptible_attempts
  }
}

task SplitPairedFastqByReadGroup {
  input {
    File fastq_1
    File fastq_2
    Array[String] readgroup_units

    Int additional_disk_space_gb = 50
    Int machine_mem_gb           = 16
    Int machine_cpu_cores        = 8
    Int preemptible_attempts     = 1
    String docker
  }

  Int disk_space_gb = ceil((size(fastq_1, "GB") + size(fastq_2, "GB")) * 6) + additional_disk_space_gb

  command <<<
    set -euo pipefail
    mkdir -p split_fastqs

    python3 <<'PY'
import gzip

fastq_1 = "~{fastq_1}"
fastq_2 = "~{fastq_2}"
rg_units_file = "~{write_lines(readgroup_units)}"

def open_fastq(path, mode):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)

def parse_readgroup(header):
    token = header.strip().split()[0]
    if token.startswith("@"):
        token = token[1:]
    fields = token.split(":")
    if len(fields) < 4:
        raise RuntimeError("FASTQ header does not match expected Illumina format: " + header.strip())
    return fields[2] + "." + fields[3]

with open(rg_units_file, "r") as f:
    rg_units = [line.strip() for line in f if line.strip()]

if not rg_units:
    raise RuntimeError("No read-group units were provided")

rg_to_index = {rg: i for i, rg in enumerate(rg_units)}

w1 = []
w2 = []
for i in range(len(rg_units)):
    w1.append(gzip.open(f"split_fastqs/rg_{i:06d}_R1.fastq.gz", "wt"))
    w2.append(gzip.open(f"split_fastqs/rg_{i:06d}_R2.fastq.gz", "wt"))

counts = [0] * len(rg_units)
total = 0

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

            rg1 = parse_readgroup(rec1[0])
            rg2 = parse_readgroup(rec2[0])
            if rg1 != rg2:
                raise RuntimeError("R1/R2 read-group mismatch: " + rg1 + " vs " + rg2)
            if rg1 not in rg_to_index:
                raise RuntimeError("Read group found in FASTQ but not in readgroup_units list: " + rg1)

            i = rg_to_index[rg1]
            w1[i].writelines(rec1)
            w2[i].writelines(rec2)
            counts[i] += 1
            total += 1
finally:
    for fh in w1:
        fh.close()
    for fh in w2:
        fh.close()

if total == 0:
    raise RuntimeError("No reads found in FASTQ files")

for i, c in enumerate(counts):
    if c == 0:
        raise RuntimeError(f"No reads assigned to readgroup index {i} ({rg_units[i]})")
PY
  >>>

  output {
    Array[File] split_fastq_1s = glob("split_fastqs/*_R1.fastq.gz")
    Array[File] split_fastq_2s = glob("split_fastqs/*_R2.fastq.gz")
  }

  runtime {
    docker:      docker
    cpu:         machine_cpu_cores
    memory:      machine_mem_gb + " GB"
    disks:       "local-disk " + disk_space_gb + " HDD"
    preemptible: preemptible_attempts
    maxRetries:  preemptible_attempts
  }
}

task ConvertReadGroupFastqToUbam {
  input {
    String sample_name
    String readgroup_name
    String effective_library_name
    String effective_platform_unit
    String effective_run_date
    String effective_platform_name
    String effective_sequencing_center
    File split_fastq_1
    File split_fastq_2
    Int readgroup_index

    Int additional_disk_space_gb = 50
    Int machine_mem_gb           = 16
    Int machine_cpu_cores        = 8
    Int preemptible_attempts     = 1
    String docker
  }

  Int disk_space_gb = ceil((size(split_fastq_1, "GB") + size(split_fastq_2, "GB")) * 5) + additional_disk_space_gb

  command <<<
    set -euo pipefail
    mkdir -p tmp

    RG_FIELD=$(printf "ID:%s\tSM:%s\tLB:%s\tPU:%s\tDT:%s\tPL:%s\tCN:%s" \
      "~{readgroup_name}" \
      "~{sample_name}" \
      "~{effective_library_name}" \
      "~{effective_platform_unit}" \
      "~{effective_run_date}" \
      "~{effective_platform_name}" \
      "~{effective_sequencing_center}")

    samtools import \
      -@ ~{machine_cpu_cores} \
      -1 ~{split_fastq_1} \
      -2 ~{split_fastq_2} \
      -r "$RG_FIELD" \
      -O BAM \
      -o rg.~{readgroup_index}.unmapped.bam
  >>>

  output {
    File unmapped_ubam = "rg.~{readgroup_index}.unmapped.bam"
  }

  runtime {
    docker:      docker
    cpu:         machine_cpu_cores
    memory:      machine_mem_gb + " GB"
    disks:       "local-disk " + disk_space_gb + " HDD"
    preemptible: preemptible_attempts
    maxRetries:  preemptible_attempts
  }
}

task SortReadGroupUbam {
  input {
    File input_ubam
    Int readgroup_index

    Int additional_disk_space_gb = 50
    Int machine_mem_gb           = 16
    Int machine_cpu_cores        = 8
    Int preemptible_attempts     = 1
    String docker
  }

  Int disk_space_gb = ceil(size(input_ubam, "GB") * 4) + additional_disk_space_gb

  command <<<
    set -euo pipefail
    mkdir -p tmp

    samtools sort \
      -n \
      -@ ~{machine_cpu_cores} \
      -m 1G \
      -T tmp/queryname_sort \
      -O BAM \
      -o rg.~{readgroup_index}.query_sorted.unmapped.bam \
      ~{input_ubam}
  >>>

  output {
    File sorted_ubam = "rg.~{readgroup_index}.query_sorted.unmapped.bam"
  }

  runtime {
    docker:      docker
    cpu:         machine_cpu_cores
    memory:      machine_mem_gb + " GB"
    disks:       "local-disk " + disk_space_gb + " HDD"
    preemptible: preemptible_attempts
    maxRetries:  preemptible_attempts
  }
}


task MergeSortedReadGroupUbams {
  input {
    String sample_name
    Array[File] input_ubams

    Int additional_disk_space_gb = 50
    Int machine_mem_gb           = 16
    Int machine_cpu_cores        = 8
    Int preemptible_attempts     = 1
    String docker
  }

  Int disk_space_gb = ceil(size(input_ubams, "GB") * 4) + additional_disk_space_gb

  command <<<
    set -euo pipefail

    if [ "~{length(input_ubams)}" -eq 1 ]; then
      cp ~{input_ubams[0]} ~{sample_name}.query_sorted.united.unmapped.bam
    else
      samtools merge         -n         -@ ~{machine_cpu_cores}         -f         -O BAM         ~{sample_name}.query_sorted.united.unmapped.bam         ~{sep=' ' input_ubams}
    fi
  >>>

  output {
    File merged_ubam = "~{sample_name}.query_sorted.united.unmapped.bam"
  }

  runtime {
    docker:      docker
    cpu:         machine_cpu_cores
    memory:      machine_mem_gb + " GB"
    disks:       "local-disk " + disk_space_gb + " HDD"
    preemptible: preemptible_attempts
    maxRetries:  preemptible_attempts
  }
}
