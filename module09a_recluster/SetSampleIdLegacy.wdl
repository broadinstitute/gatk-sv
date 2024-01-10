version 1.0

import "Structs.wdl"

workflow SetSampleIdLegacy {
  input {
    String sample_name
    File? BAF_file
    File? BAF_file_index
    File? PE_file
    File? PE_file_index
    File? SR_file
    File? SR_file_index
    File? SD_file
    File? SD_file_index
    File reference_dict
    String gatk_docker
    RuntimeAttr? runtime_attr_override
  }

  if (defined(BAF_file)) {
    call SetSampleId as SetSampleIdBAF {
      input:
        evidence_file = select_first([BAF_file]),
        evidence_file_index = BAF_file_index,
        file_type = "baf",
        sample_name = sample_name,
        sample_column = 4,
        reference_dict = reference_dict,
        gatk_docker = gatk_docker,
        runtime_attr_override = runtime_attr_override
    }
  }

  if (defined(PE_file)) {
    call SetSampleId as SetSampleIdPE {
      input:
        evidence_file = select_first([PE_file]),
        evidence_file_index = PE_file_index,
        file_type = "pe",
        sample_name = sample_name,
        sample_column = 7,
        reference_dict = reference_dict,
        gatk_docker = gatk_docker,
        runtime_attr_override = runtime_attr_override
    }
  }

  if (defined(SR_file)) {
    call SetSampleId as SetSampleIdSR {
      input:
        evidence_file = select_first([SR_file]),
        evidence_file_index = SR_file_index,
        file_type = "sr",
        sample_name = sample_name,
        sample_column = 5,
        reference_dict = reference_dict,
        gatk_docker = gatk_docker,
        runtime_attr_override = runtime_attr_override
    }
  }

  if (defined(SD_file)) {
    call SetSampleId as SetSampleIdSD {
      input:
        evidence_file = select_first([SD_file]),
        evidence_file_index = SD_file_index,
        file_type = "sd",
        sample_name = sample_name,
        sample_column = 3,
        reference_dict = reference_dict,
        gatk_docker = gatk_docker,
        runtime_attr_override = runtime_attr_override
    }
  }

  output {
    File? BAF_out = SetSampleIdBAF.out
    File? BAF_out_index = SetSampleIdBAF.out_index
    File? PE_out = SetSampleIdPE.out
    File? PE_out_index = SetSampleIdPE.out_index
    File? SR_out = SetSampleIdSR.out
    File? SR_out_index = SetSampleIdSR.out_index
    File? SD_out = SetSampleIdSD.out
    File? SD_out_index = SetSampleIdSD.out_index
  }
}

task SetSampleId {
  input {
    File evidence_file
    File? evidence_file_index
    String file_type
    String sample_name
    Int sample_column
    File reference_dict
    String gatk_docker
    RuntimeAttr? runtime_attr_override
  }

  Int disk_size = 10 + ceil(size(evidence_file, "GiB") * 2)

  RuntimeAttr default_attr = object {
    cpu_cores: 2,
    mem_gb: 3.75,
    disk_gb: disk_size,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File out = "~{sample_name}.~{file_type}.txt.gz"
    File out_index = "~{sample_name}.~{file_type}.txt.gz.tbi"
  }
  command <<<

    set -euo pipefail

    fifo_name="~{sample_name}.~{file_type}.txt"
    output_name="~{sample_name}.~{file_type}.txt.gz"

    if [ ! -f "~{evidence_file}.tbi" ]; then
        tabix -s1 -b2 -e2 ~{evidence_file}
    fi

    mkfifo $fifo_name
    /gatk/gatk --java-options "-Xmx2000m" PrintSVEvidence -F ~{evidence_file} --sequence-dictionary ~{reference_dict} -O $fifo_name &
    awk '{$~{sample_column}="~{sample_name}"}' < $fifo_name | bgzip -c > $output_name
    tabix -s1 -b2 -e2 $output_name

  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: gatk_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}
