version 1.0

workflow ConvertPairedFastQsToUnmappedBamWf {
  input {
    String sample_name
    File fastq_1
    File fastq_2

    Array[String]? library_names
    Array[String]? platform_names
    Array[String]? sequencing_centers
    Array[String]? run_dates

    File ref_dict

    Int reads_per_chunk = 2000000
    String samtools_docker = "staphb/samtools:1.17"

    Int additional_disk_space_gb = 50
    Int machine_mem_gb           = 16
    Int machine_cpu_cores        = 8
    Int preemptible_attempts     = 1
  }

  call DiscoverReadGroupsAndMetadata {
    input:
      sample_name               = sample_name,
      fastq_1                   = fastq_1,
      library_names             = library_names,
      platform_names            = platform_names,
      sequencing_centers        = sequencing_centers,
      run_dates                 = run_dates,
      ref_dict                  = ref_dict,
      docker                    = samtools_docker,
      additional_disk_space_gb  = additional_disk_space_gb,
      machine_mem_gb            = machine_mem_gb,
      machine_cpu_cores         = machine_cpu_cores,
      preemptible_attempts      = preemptible_attempts
  }

  call SplitPairedFastqIntoChunks {
    input:
      fastq_1                   = fastq_1,
      fastq_2                   = fastq_2,
      reads_per_chunk           = reads_per_chunk,
      docker                    = samtools_docker,
      additional_disk_space_gb  = additional_disk_space_gb,
      machine_mem_gb            = machine_mem_gb,
      machine_cpu_cores         = machine_cpu_cores,
      preemptible_attempts      = preemptible_attempts
  }

  scatter (i in range(length(SplitPairedFastqIntoChunks.chunk_fastq_1s))) {
    call ConvertChunkToUbamByReadGroup {
      input:
        chunk_fastq_1              = SplitPairedFastqIntoChunks.chunk_fastq_1s[i],
        chunk_fastq_2              = SplitPairedFastqIntoChunks.chunk_fastq_2s[i],
        chunk_index                = i,
        sample_name                = sample_name,
        platform_units             = DiscoverReadGroupsAndMetadata.platform_units,
        readgroup_names            = DiscoverReadGroupsAndMetadata.readgroup_names,
        effective_library_names    = DiscoverReadGroupsAndMetadata.effective_library_names,
        effective_platform_names   = DiscoverReadGroupsAndMetadata.effective_platform_names,
        effective_sequencing_centers = DiscoverReadGroupsAndMetadata.effective_sequencing_centers,
        effective_run_dates        = DiscoverReadGroupsAndMetadata.effective_run_dates,
        docker                     = samtools_docker,
        additional_disk_space_gb   = additional_disk_space_gb,
        machine_mem_gb             = machine_mem_gb,
        machine_cpu_cores          = machine_cpu_cores,
        preemptible_attempts       = preemptible_attempts
    }
  }

  call MergeChunkUbams {
    input:
      sample_name               = sample_name,
      chunk_ubams               = ConvertChunkToUbamByReadGroup.output_chunk_ubam,
      header_sam                = DiscoverReadGroupsAndMetadata.full_header_sam,
      docker                    = samtools_docker,
      additional_disk_space_gb  = additional_disk_space_gb,
      machine_mem_gb            = machine_mem_gb,
      machine_cpu_cores         = machine_cpu_cores,
      preemptible_attempts      = preemptible_attempts
  }

  call SortSam {
    input:
      input_bam                 = MergeChunkUbams.output_unmapped_bam,
      sample_name               = sample_name,
      docker                    = samtools_docker,
      machine_mem_gb            = machine_mem_gb,
      machine_cpu_cores         = machine_cpu_cores,
      additional_disk_space_gb  = additional_disk_space_gb,
      preemptible_attempts      = preemptible_attempts
  }

  output {
    File output_unmapped_bam = SortSam.sorted_bam
    Array[String] readgroup_names = DiscoverReadGroupsAndMetadata.readgroup_names
    Array[String] platform_units = DiscoverReadGroupsAndMetadata.platform_units
    Array[String] effective_run_dates = DiscoverReadGroupsAndMetadata.effective_run_dates
    Array[File] split_chunk_fastq_1s = SplitPairedFastqIntoChunks.chunk_fastq_1s
    Array[File] split_chunk_fastq_2s = SplitPairedFastqIntoChunks.chunk_fastq_2s
  }
}

task DiscoverReadGroupsAndMetadata {
  input {
    String sample_name
    File fastq_1
    Array[String]? library_names
    Array[String]? platform_names
    Array[String]? sequencing_centers
    Array[String]? run_dates
    File ref_dict

    Int additional_disk_space_gb = 50
    Int machine_mem_gb           = 16
    Int machine_cpu_cores        = 8
    Int preemptible_attempts     = 1
    String docker
  }

  Array[String] provided_library_names = select_first([library_names, []])
  Array[String] provided_platform_names = select_first([platform_names, []])
  Array[String] provided_sequencing_centers = select_first([sequencing_centers, []])
  Array[String] provided_run_dates = select_first([run_dates, []])
  Int disk_space_gb = ceil(size(fastq_1, "GB") * 4) + additional_disk_space_gb

  command <<<
    set -euo pipefail

    mapfile -t LIB_NAMES < ~{write_lines(provided_library_names)}
    mapfile -t PLAT_NAMES < ~{write_lines(provided_platform_names)}
    mapfile -t SEQ_CENTERS < ~{write_lines(provided_sequencing_centers)}
    mapfile -t RUN_DATES < ~{write_lines(provided_run_dates)}

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
    }' | sort -u > platform_units.txt

    n_rg=$(wc -l < platform_units.txt)
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
    check_array_length PLAT_NAMES
    check_array_length SEQ_CENTERS
    check_array_length RUN_DATES

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
    : > effective_run_dates.txt
    : > effective_library_names.txt
    : > effective_platform_names.txt
    : > effective_sequencing_centers.txt
    printf "@HD\tVN:1.6\tSO:queryname\n" > new_header.sam
    grep "^@SQ" ~{ref_dict} >> new_header.sam

    mapfile -t PUS < platform_units.txt
    for (( i=0; i<n_rg; i++ )); do
      PU="${PUS[$i]}"
      RG_ID="~{sample_name}.${PU}"
      LB=$(pick_or_default LIB_NAMES "$i" "~{sample_name}")
      PL=$(pick_or_default PLAT_NAMES "$i" "ILLUMINA")
      CN=$(pick_or_default SEQ_CENTERS "$i" "ILLUMINA")
      DT=$(pick_or_default RUN_DATES "$i" "$JOB_DATE")

      printf "@RG\tID:%s\tSM:%s\tLB:%s\tPU:%s\tDT:%s\tPL:%s\tCN:%s\n" \
        "$RG_ID" "~{sample_name}" "$LB" "$PU" "$DT" "$PL" "$CN" >> new_header.sam
      echo "$RG_ID" >> readgroup_names.txt
      echo "$DT" >> effective_run_dates.txt
      echo "$LB" >> effective_library_names.txt
      echo "$PL" >> effective_platform_names.txt
      echo "$CN" >> effective_sequencing_centers.txt
    done
  >>>

  output {
    Array[String] platform_units = read_lines("platform_units.txt")
    Array[String] readgroup_names = read_lines("readgroup_names.txt")
    Array[String] effective_run_dates = read_lines("effective_run_dates.txt")
    Array[String] effective_library_names = read_lines("effective_library_names.txt")
    Array[String] effective_platform_names = read_lines("effective_platform_names.txt")
    Array[String] effective_sequencing_centers = read_lines("effective_sequencing_centers.txt")
    File full_header_sam = "new_header.sam"
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

task SplitPairedFastqIntoChunks {
  input {
    File fastq_1
    File fastq_2
    Int reads_per_chunk = 2000000

    Int additional_disk_space_gb = 50
    Int machine_mem_gb           = 16
    Int machine_cpu_cores        = 8
    Int preemptible_attempts     = 1
    String docker
  }

  Int disk_space_gb = ceil((size(fastq_1, "GB") + size(fastq_2, "GB")) * 6) + additional_disk_space_gb

  command <<<
    set -euo pipefail
    mkdir -p chunks

    python3 <<'PY'
import gzip

fastq_1 = "~{fastq_1}"
fastq_2 = "~{fastq_2}"
reads_per_chunk = int("~{reads_per_chunk}")

def open_fastq(path, mode):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)

if reads_per_chunk <= 0:
    raise RuntimeError("reads_per_chunk must be > 0")

chunk_idx = 0
reads_in_chunk = 0
w1 = None
w2 = None
chunk_paths_1 = []
chunk_paths_2 = []
total_reads = 0

def open_chunk(i):
    p1 = f"chunks/chunk_{i:06d}_R1.fastq.gz"
    p2 = f"chunks/chunk_{i:06d}_R2.fastq.gz"
    return gzip.open(p1, "wt"), gzip.open(p2, "wt"), p1, p2

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

        if w1 is None:
            w1, w2, p1, p2 = open_chunk(chunk_idx)
            chunk_paths_1.append(p1)
            chunk_paths_2.append(p2)
            chunk_idx += 1

        w1.writelines(rec1)
        w2.writelines(rec2)
        reads_in_chunk += 1
        total_reads += 1

        if reads_in_chunk >= reads_per_chunk:
            w1.close()
            w2.close()
            w1 = None
            w2 = None
            reads_in_chunk = 0

if w1 is not None:
    w1.close()
    w2.close()

if total_reads == 0:
    raise RuntimeError("No reads found in FASTQ files")

with open("chunk_fastq_1.list", "w") as f1:
    for p in chunk_paths_1:
        f1.write(p + "\n")
with open("chunk_fastq_2.list", "w") as f2:
    for p in chunk_paths_2:
        f2.write(p + "\n")
PY
  >>>

  output {
    Array[File] chunk_fastq_1s = glob("chunks/*_R1.fastq.gz")
    Array[File] chunk_fastq_2s = glob("chunks/*_R2.fastq.gz")
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

task ConvertChunkToUbamByReadGroup {
  input {
    File chunk_fastq_1
    File chunk_fastq_2
    Int chunk_index
    String sample_name

    Array[String] platform_units
    Array[String] readgroup_names
    Array[String] effective_library_names
    Array[String] effective_platform_names
    Array[String] effective_sequencing_centers
    Array[String] effective_run_dates

    Int additional_disk_space_gb = 50
    Int machine_mem_gb           = 16
    Int machine_cpu_cores        = 8
    Int preemptible_attempts     = 1
    String docker
  }

  Int disk_space_gb = ceil((size(chunk_fastq_1, "GB") + size(chunk_fastq_2, "GB")) * 5) + additional_disk_space_gb

  command <<<
    set -euo pipefail
    mkdir -p chunk_split per_rg

    mapfile -t PLATFORM_UNITS < ~{write_lines(platform_units)}
    mapfile -t RG_NAMES < ~{write_lines(readgroup_names)}
    mapfile -t LIB_NAMES < ~{write_lines(effective_library_names)}
    mapfile -t PLAT_NAMES < ~{write_lines(effective_platform_names)}
    mapfile -t SEQ_CENTERS < ~{write_lines(effective_sequencing_centers)}
    mapfile -t RUN_DATES < ~{write_lines(effective_run_dates)}

    n_rg=${#PLATFORM_UNITS[@]}
    for arr_name in RG_NAMES LIB_NAMES PLAT_NAMES SEQ_CENTERS RUN_DATES; do
      declare -n arr="$arr_name"
      if [ "${#arr[@]}" -ne "$n_rg" ]; then
        echo "ERROR: metadata array length mismatch for $arr_name" >&2
        exit 1
      fi
    done

    declare -A PU_TO_INDEX
    for ((i=0; i<n_rg; i++)); do
      PU_TO_INDEX["${PLATFORM_UNITS[$i]}"]="$i"
    done

    python3 <<'PY'
import gzip
import re

fastq_1 = "~{chunk_fastq_1}"
fastq_2 = "~{chunk_fastq_2}"

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
paths = {}
counts = {}

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
                safe = re.sub(r"[^A-Za-z0-9._-]", "_", pu1)
                out1 = "chunk_split/" + safe + "_R1.fastq.gz"
                out2 = "chunk_split/" + safe + "_R2.fastq.gz"
                writers[pu1] = (gzip.open(out1, "wt"), gzip.open(out2, "wt"))
                paths[pu1] = (out1, out2)
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
    raise RuntimeError("No reads found in chunk")

pus = sorted(counts.keys())
with open("local_platform_units.txt", "w") as f:
    for pu in pus:
        f.write(pu + "\n")
with open("local_fastq_1.list", "w") as f1:
    for pu in pus:
        f1.write(paths[pu][0] + "\n")
with open("local_fastq_2.list", "w") as f2:
    for pu in pus:
        f2.write(paths[pu][1] + "\n")
PY

    mapfile -t LOCAL_PUS < local_platform_units.txt
    mapfile -t LOCAL_FQ1 < local_fastq_1.list
    mapfile -t LOCAL_FQ2 < local_fastq_2.list

    n_local=${#LOCAL_PUS[@]}
    if [ "$n_local" -eq 0 ]; then
      echo "ERROR: no local read groups found in chunk." >&2
      exit 1
    fi

    for ((j=0; j<n_local; j++)); do
      pu="${LOCAL_PUS[$j]}"
      if [ -z "${PU_TO_INDEX[$pu]+x}" ]; then
        echo "ERROR: chunk read group $pu not found in global metadata" >&2
        exit 1
      fi
      idx="${PU_TO_INDEX[$pu]}"

      RG_FIELD=$(printf "ID:%s\tSM:%s\tLB:%s\tPU:%s\tDT:%s\tPL:%s\tCN:%s" \
        "${RG_NAMES[$idx]}" \
        "~{sample_name}" \
        "${LIB_NAMES[$idx]}" \
        "${PLATFORM_UNITS[$idx]}" \
        "${RUN_DATES[$idx]}" \
        "${PLAT_NAMES[$idx]}" \
        "${SEQ_CENTERS[$idx]}")

      samtools import \
        -@ ~{machine_cpu_cores} \
        -1 "${LOCAL_FQ1[$j]}" \
        -2 "${LOCAL_FQ2[$j]}" \
        -r "$RG_FIELD" \
        -O BAM \
        -o "per_rg/rg_${j}.bam"
    done

    if [ "$n_local" -eq 1 ]; then
      cp per_rg/rg_0.bam chunk.~{chunk_index}.unmapped.bam
    else
      samtools merge -@ ~{machine_cpu_cores} -f -O BAM chunk.~{chunk_index}.unmapped.bam per_rg/rg_*.bam
    fi
  >>>

  output {
    File output_chunk_ubam = "chunk.~{chunk_index}.unmapped.bam"
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

task MergeChunkUbams {
  input {
    String sample_name
    Array[File] chunk_ubams
    File header_sam

    Int additional_disk_space_gb = 50
    Int machine_mem_gb           = 16
    Int machine_cpu_cores        = 8
    Int preemptible_attempts     = 1
    String docker
  }

  Int disk_space_gb = ceil(size(chunk_ubams, "GB") * 4) + additional_disk_space_gb

  command <<<
    set -euo pipefail

    if [ "~{length(chunk_ubams)}" -eq 1 ]; then
      cp ~{chunk_ubams[0]} merged.raw.unmapped.bam
    else
      samtools merge -@ ~{machine_cpu_cores} -f -O BAM merged.raw.unmapped.bam ~{sep=' ' chunk_ubams}
    fi

    samtools reheader ~{header_sam} merged.raw.unmapped.bam > ~{sample_name}.unmapped.bam
  >>>

  output {
    File output_unmapped_bam = "~{sample_name}.unmapped.bam"
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

task SortSam {
  input {
    File input_bam
    String sample_name
    String docker
    Int machine_mem_gb           = 16
    Int machine_cpu_cores        = 8
    Int additional_disk_space_gb = 50
    Int preemptible_attempts     = 1
  }

  Int disk_space_gb  = ceil(size(input_bam, "GB") * 4) + additional_disk_space_gb

  command <<<
    set -euo pipefail
    mkdir -p tmp

    samtools sort \
      -n \
      -@ ~{machine_cpu_cores} \
      -m 1G \
      -T tmp/queryname_sort \
      -O BAM \
      -o ~{sample_name}.query_sorted.unmapped.bam \
      ~{input_bam}
  >>>

  output {
    File sorted_bam = "~{sample_name}.query_sorted.unmapped.bam"
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
