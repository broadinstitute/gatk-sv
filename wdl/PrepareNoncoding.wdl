version 1.0

import "Structs.wdl"

# Prepare input BED file list containing intervals of noncoding regions and output one BED file for use in annotation sub-module
workflow PrepareNoncoding {

  input {
    
    File noncoding_bed_list

    String sv_base_mini_docker

    RuntimeAttr? runtime_attr_clean_noncoding_bed
    RuntimeAttr? runtime_attr_make_noncoding_bed
  }

  Array[File] noncoding_beds = read_lines(noncoding_bed_list)

  scatter (bed in noncoding_beds) {
    
    call CleanNoncodingBed {
      input:

        bed = bed,

        sv_base_mini_docker = sv_base_mini_docker, 
        
        runtime_attr_override = runtime_attr_clean_noncoding_bed
    }
  }

  call MakeNoncodingBed {
    input:

      beds = CleanNoncodingBed.cleaned_bed,
      
      runtime_attr_override = runtime_attr_make_noncoding_bed
  }

  output {
    File noncoding_bed = MakeNoncodingBed.noncoding_bed
  }
}

task CleanNoncodingBed {

  input {

    File bed

    String sv_base_mini_docker
    
    RuntimeAttr? runtime_attr_override
  }
  
  String name = basename(bed, ".bed")
  
  output {
    File cleaned_bed = "${name}.bed"
  }

  command <<<

    set -euo pipefail
    cat ~{bed} \
      | cut -f -3 \
      | sort -k1,1V -k2,2n \
      | bedtools merge -i stdin \
      | awk -v OFS="\t" '{print $0, "~{name}"}' \
      > ~{name}.bed
  
  >>>
  
  #########################
  RuntimeAttr default_attr = object {
    cpu_cores:          1, 
    mem_gb:             3.75, 
    disk_gb:            10,
    boot_disk_gb:       10,
    preemptible_tries:  3,
    max_retries:        0
    # docker:             "gatksv/sv-base-mini:v0.1"
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  runtime {
    cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
    memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
    disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
    preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
    docker:                 sv_base_mini_docker # select_first([runtime_attr.docker,            default_attr.docker])
  }
}

task MakeNoncodingBed {
  
  input {
    
    Array[File] beds

    RuntimeAttr? runtime_attr_override
  }
  
  output {
    File noncoding_bed = "noncoding_elements.bed"
  }

  command <<<

    sort -k1,1V -k2,2n -m ~{sep=" "  beds} > noncoding_elements.bed
  
  >>>
  
  #########################
  # note here we did not make docker an input because it is really a small task (only sort), I actually suspect alpine is enough
  RuntimeAttr default_attr = object {
    cpu_cores:          1, 
    mem_gb:             3.75, 
    disk_gb:            10,
    boot_disk_gb:       10,
    preemptible_tries:  3,
    max_retries:        0
    # docker:             "ubuntu:18.10"
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  runtime {
    cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
    memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
    disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
    preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
    docker:                 "ubuntu:18.10" # select_first([runtime_attr.docker,            default_attr.docker])
  }
}
