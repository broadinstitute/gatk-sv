##########################################################################################

## Base script:   https://portal.firecloud.org/#methods/Talkowski-SV/pesr_collection/2/wdl

## Github commit: talkowski-lab/gatk-sv-v1:<ENTER HASH HERE IN FIRECLOUD>

##########################################################################################

version 1.0

import "Structs.wdl"

# Workflow to run PE/SR collection on a single sample
workflow PESRCollection {
  input {
    File cram
    File cram_index
    String sample_id
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  call RunPESRCollection {
    input:
      cram = cram,
      cram_index = cram_index,
      sample_id = sample_id,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_override
  }

  output {
    File disc_out = RunPESRCollection.disc_out
    File disc_out_index = RunPESRCollection.disc_out_index
    File split_out = RunPESRCollection.split_out
    File split_out_index = RunPESRCollection.split_out_index
  }
}

# Task to run collect-pesr on a single sample
task RunPESRCollection {
  input {
    File cram
    File cram_index
    String sample_id
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float cram_size = size(cram, "GiB")
  Int vm_disk_size = ceil(cram_size + 50)

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 3.75, 
    disk_gb: vm_disk_size,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File split_out = "${sample_id}.split.txt.gz"
    File split_out_index = "${sample_id}.split.txt.gz.tbi"
    File disc_out = "${sample_id}.disc.txt.gz"
    File disc_out_index = "${sample_id}.disc.txt.gz.tbi"
  }
  command <<<

    set -euo pipefail
    svtk collect-pesr \
      ~{cram} \
      ~{sample_id} \
      ~{sample_id}.split.unsorted.txt \
      ~{sample_id}.disc.unsorted.txt

    # As of 11/14/2018, svtk does not sort the split file
    tmpdir=$(mktemp -d)
    sort -k1,1V -k2,2n -T $tmpdir ~{sample_id}.split.unsorted.txt > ~{sample_id}.split.txt
    sort -k1,1V -k2,2n -T $tmpdir ~{sample_id}.disc.unsorted.txt > ~{sample_id}.disc.txt

    bgzip ~{sample_id}.disc.txt
    bgzip ~{sample_id}.split.txt
    
    tabix -f -s1 -b 2 -e 2 ~{sample_id}.disc.txt.gz
    tabix -f -s1 -b 2 -e 2 ~{sample_id}.split.txt.gz
  
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

