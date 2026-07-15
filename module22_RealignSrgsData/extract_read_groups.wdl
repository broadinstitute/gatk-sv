version 1.0

workflow ExtractReadGroupsWorkflow {
  input {
    File fastq_gz
    String sample_id
  }

  call ExtractReadGroups {
    input:
      fastq_gz = fastq_gz,
      sample_id = sample_id
  }

  output {
    Array[String] PU = ExtractReadGroups.PU
    Array[String] ID = ExtractReadGroups.ID
    String DT = ExtractReadGroups.DT
    String PL = ExtractReadGroups.PL
  }
}

task ExtractReadGroups {
  input {
    File fastq_gz
    String sample_id
  }

  command <<<
    set -euo pipefail

    # Extract flowcell + lane from every read header, then get unique combinations
    # Illumina header format: @<instrument>:<run>:<flowcell>:<lane>:<tile>:<x>:<y> <read>:<filtered>:<control>:<index>
    zcat ~{fastq_gz} | awk 'NR % 4 == 1' | awk -F: '{print $3"."$4}' | sort -u > lane_combos.txt

    # Build PU (flowcell.lane) and ID (sample.flowcell.lane) for each unique combo
    > pu_list.txt
    > id_list.txt
    while IFS= read -r combo; do
      echo "${combo}" >> pu_list.txt
      echo "~{sample_id}.${combo}" >> id_list.txt
    done < lane_combos.txt

    # DT: date the task is run, in ISO8601 format
    date -u +"%Y-%m-%dT%H:%M:%SZ" > dt.txt
  >>>

  output {
    Array[String] PU = read_lines("pu_list.txt")
    Array[String] ID = read_lines("id_list.txt")
    String DT = read_string("dt.txt")
    String PL = "ILLUMINA"
  }

  runtime {
    docker: "ubuntu:22.04"
    cpu: 1
    memory: "2 GB"
  }
}
