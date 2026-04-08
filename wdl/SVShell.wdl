version 1.0

import "Structs.wdl"

workflow SVShell {
  input {
    File gcnv_model_tars_list
    File ref_pesr_split_files_list
    File ref_pesr_disc_files_list
    File ref_pesr_sd_files_list
    Array[File] genome_tracks
  }

  Array[File] gcnv_model_tars = read_lines(gcnv_model_tars_list)

  Array[File] ref_pesr_split_files = read_lines(ref_pesr_split_files_list)
  scatter (ref_pesr_split_file in ref_pesr_split_files) {
    File ref_pesr_split_file_index = ref_pesr_split_file + ".tbi"
  }

  Array[File] ref_pesr_disc_files = read_lines(ref_pesr_disc_files_list)
  scatter (ref_pesr_disc_file in ref_pesr_disc_files) {
    File ref_pesr_disc_file_index = ref_pesr_disc_file + ".tbi"
  }

  Array[File] ref_pesr_sd_files = read_lines(ref_pesr_sd_files_list)
  scatter (ref_pesr_sd_file in ref_pesr_sd_files) {
    File ref_pesr_sd_file_index = ref_pesr_sd_file + ".tbi"
  }

  scatter (genome_track in genome_tracks) {
    File genome_track_index = genome_track + ".tbi"
  }

  call RunSVShell {
    input:
      gcnv_model_tars = gcnv_model_tars,
      ref_pesr_split_files = ref_pesr_split_files,
      ref_pesr_split_files_indices = ref_pesr_split_file_index,
      ref_pesr_disc_files = ref_pesr_disc_files,
      ref_pesr_disc_files_indices = ref_pesr_disc_file_index,
      ref_pesr_sd_files = ref_pesr_sd_files,
      ref_pesr_sd_files_indices = ref_pesr_sd_file_index,
      genome_tracks = genome_tracks,
      genome_tracks_indices = genome_track_index,
  }


  output {
    File inputs_json = RunSVShell.inputs_json
  }
}

task RunSVShell {
  input {
    String batch
    Array[File] gcnv_model_tars
    Array[File] ref_pesr_split_files
    Array[File] ref_pesr_split_files_indices
    Array[File] ref_pesr_disc_files
    Array[File] ref_pesr_disc_files_indices
    Array[File] ref_pesr_sd_files
    Array[File] ref_pesr_sd_files_indices
    Array[File] genome_tracks
    Array[File] genome_tracks_indices

    String sv_shell_docker
    RuntimeAttr? runtime_attr_override
  }

  command <<<
    set -Exeuo pipefail

    export BASE_DIR="${PWD}"
    export SV_SHELL_BASE_DIR="${PWD}/wd"
    export TMPDIR="${PWD}/wd/tmp"
    mkdir -p "${PWD}/wd/tmp"

    touch single_sample_pipeline_inputs.json
  >>>

  output {
    File inputs_json = "single_sample_pipeline_inputs.json"
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 4,
    mem_gb: 16,
    disk_gb: 400,
    boot_disk_gb: 30,
    preemptible_tries: 0,
    max_retries: 0
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " SSD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_shell_docker
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}