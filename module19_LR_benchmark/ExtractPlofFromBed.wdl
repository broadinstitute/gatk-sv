version 1.0

workflow ExtractPlofFromBed {
  input {
    Array[File] bed_files   # each file ends in .bed.gz
    String sv_base_mini_docker
  }

  scatter (bed in bed_files) {
    call FilterBed { 
    input: 
      bed = bed,
      sv_base_mini_docker = sv_base_mini_docker
    }
  }

  call ConcatBeds {
    input:
      filtered_beds = FilterBed.filtered_bed,
      sv_base_mini_docker = sv_base_mini_docker
  }

  output {
    Array[File] per_file_outputs = FilterBed.filtered_bed
    File merged_bed = ConcatBeds.merged_bed
  }
}


task FilterBed {
  input {
    File bed   # .bed.gz
    String sv_base_mini_docker
  }

  String prefix = basename(bed, ".bed.gz")

  command <<<
    set -euo pipefail

    zcat ~{bed} \
      | cut -f1-6,41,594 \
      | awk '{if ($NF=="PASS") print}' \
      | awk '{if ($(NF-1)!="NA") print}' \
      | bgzip > ~{prefix}.pLoF.bed.gz
  >>>

  output {
    File filtered_bed = "~{prefix}.pLoF.bed.gz"
  }

  runtime { 
    docker: 
      sv_base_mini_docker }
}


task ConcatBeds {
  input {
    Array[File] filtered_beds
    String sv_base_mini_docker
  }

  command <<<
    set -euo pipefail

    # Concatenate all filtered beds
    zcat ~{sep=" " filtered_beds} | bgzip > merged.pLoF.bed.gz
  >>>

  output {
    File merged_bed = "merged.pLoF.bed.gz"
  }

  runtime { docker: sv_base_mini_docker}
}