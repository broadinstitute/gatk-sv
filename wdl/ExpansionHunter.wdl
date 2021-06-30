## Runs ExpansionHunter denovo (ehdn)

version 1.0

import "Structs.wdl"

struct FilenamePostfixes {
  String locus
  String motif
  String profile
  String merged_profile
  Int profile_len
}

workflow ComputeSTRProfiles {

  input {
    Array[File] case_reads_filenames
    Array[File] case_indexes_filenames
    Array[File] control_reads_filenames
    Array[File] control_indexes_filenames
    File reference_filename
    Int min_anchor_mapq
    Int max_irr_mapq
    String ehdn_docker
    RuntimeAttr? runtime_attr_str_profile
    RuntimeAttr? runtime_attr_merge
  }

  parameter_meta {
    case_reads_filenames: ""
    control_reads_filenames: ""
    manifest_filename: ""
    reference_filename: ""
    min_anchor_mapq: ""
    max_irr_mapq: ""
    ehdn_docker: ""
  }

  # The values of the variables are based on
  # ehdn's current latest hard-coded postfixes.
  # Do not change them unless they are changed
  # in ehdn.
  FilenamePostfixes postfixes = object {
    locus: ".locus.tsv",
    motif: ".motif.tsv",
    profile: ".str_profile.json",
    merged_profile: ".multisample_profile.json",

    # This the length of `profile` postfix without
    # the filename extension. It is used to remove
    # postfix in order to extract sample name from
    # ehdn generated output.
    # e.g., extract `sample1` from `sample1.str_profile.json`.
    profile_len: 12
  }

  String locus_filename_postfix = ".locus.tsv"
  String motif_filename_postfix = ".motif.tsv"
  String profile_filename_postfix = ".str_profile.json"

  RuntimeAttr runtime_attr_str_profile_default = object {
    cpu_cores: 1,
    mem_gb: 4,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }
  RuntimeAttr runtime_attr_str_profile = select_first([
    runtime_attr_str_profile,
    runtime_attr_str_profile_default])

  RuntimeAttr runtime_attr_merge_default = object {
    cpu_cores: 1,
    mem_gb: 4,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }
  RuntimeAttr runtime_attr_merge = select_first([
    runtime_attr_merge,
    runtime_attr_merge_default])

  scatter(pair in zip(case_reads_filenames, case_indexes_filenames)) {
    String case_sample_filename = basename(pair.left, ".bam")
    call ComputeSTRProfile as CasesProfiles {
      input:
        filename = case_sample_filename,
        reads_filename = pair.left,
        index_filename = pair.right,
        reference_filename = reference_filename,
        min_anchor_mapq = min_anchor_mapq,
        max_irr_mapq = max_irr_mapq,
        ehdn_docker = ehdn_docker,
        runtime_attr = runtime_attr_str_profile,
        postfixes = postfixes
    }
  }

  scatter(pair in zip(control_reads_filenames, control_indexes_filenames)) {
    String control_sample_filename = basename(pair.left, ".bam")
    call ComputeSTRProfile as ControlsProfiles {
      input:
        filename = control_sample_filename,
        reads_filename = pair.left,
        index_filename = pair.right,
        reference_filename = reference_filename,
        min_anchor_mapq = min_anchor_mapq,
        max_irr_mapq = max_irr_mapq,
        ehdn_docker = ehdn_docker,
        runtime_attr = runtime_attr_str_profile,
        postfixes = postfixes
    }
  }

  call Merge {
    input:
      cases = CasesProfiles.str_profile,
      controls = ControlsProfiles.str_profile,
      reference_filename = reference_filename,
      ehdn_docker = ehdn_docker,
      runtime_attr = runtime_attr_merge,
      postfixes = postfixes
  }

  output {
    Array[File] cases_locus = CasesProfiles.locus
    Array[File] cases_motif = CasesProfiles.motif
    Array[File] cases_str_profile = CasesProfiles.str_profile
    Array[File] controls_locus = ControlsProfiles.locus
    Array[File] controls_motif = ControlsProfiles.motif
    Array[File] controls_str_profile = ControlsProfiles.str_profile
    File multisample_profile = Merge.multisample_profile
  }
}

task ComputeSTRProfile {
  input {
    String filename
    File reads_filename
    File index_filename
    File reference_filename
    Int min_anchor_mapq
    Int max_irr_mapq
    String ehdn_docker
    RuntimeAttr runtime_attr
    FilenamePostfixes postfixes
  }

  output {
    File locus = "${filename}${postfixes.locus}"
    File motif = "${filename}${postfixes.motif}"
    File str_profile = "${filename}${postfixes.profile}"
  }

  command <<<

    ExpansionHunterDenovo profile \
    --reads ~{reads_filename} \
    --reference ~{reference_filename} \
    --output-prefix ~{filename} \
    --min-anchor-mapq ~{min_anchor_mapq} \
    --max-irr-mapq ~{max_irr_mapq}

  >>>

  runtime {
    docker: ehdn_docker
    cpu: runtime_attr.cpu_cores
    memory: runtime_attr.mem_gb + " GiB"
    disks: "local-disk " + runtime_attr.disk_gb + " HDD"
    bootDiskSizeGb: runtime_attr.boot_disk_gb
    preemptible: runtime_attr.preemptible_tries
    maxRetries: runtime_attr.max_retries
  }
}

task Merge {
  input {
    Array[File] cases
    Array[File] controls
    File reference_filename
    String ehdn_docker
    RuntimeAttr runtime_attr
    FilenamePostfixes postfixes
  }

  output {
    File multisample_profile = "${output_prefix}${postfixes.merged_profile}"
  }

  Int cases_length = length(cases)
  Int controls_length = length(controls)
  String output_prefix = "merged"

  # This shell script uses ehdn's `merge` command
  # to merge STR profiles on all the individual samples
  # into a single JSON file.
  #
  # The script is composed of two parts:
  # - Create a "manifest" file based on the given inputs.
  #   The manifest file is composed of three columns:
  #   (1) Sample name (infered from the filename);
  #   (2) case/control (infered from files in the inputs
  #       `cases` and `controls`);
  #   (3) file path.
  #
  # - Call ehdn's `merge` method using the input and
  #   the generaged manifest file.
  command <<<
    get_sample_name()
    {
      filename=$(basename -- "$1")
      filename="${filename%.*}"
      echo "${filename::-~{postfixes.profile_len}}"
    }
    manifest_filename="manifest.tsv"

    echo ~{cases_length}
    echo ~{controls_length}
    if [ ~{cases_length} -ne 0 ]; then
      for i in "~{sep=" " cases}"; do
        echo -e "$( get_sample_name "$i" )\tcase\t$i" >> $manifest_filename
      done
    fi
    if [ ~{controls_length} -ne 0 ]; then
      for i in "~{sep=" " controls}"; do
        echo -e "$( get_sample_name "$i" )\tcontrol\t$i" >> $manifest_filename
      done
    fi
    cat $manifest_filename

    ExpansionHunterDenovo merge \
      --reference ~{reference_filename} \
      --manifest $manifest_filename \
      --output-prefix ~{output_prefix}
  >>>

  runtime {
    docker: ehdn_docker
    cpu: runtime_attr.cpu_cores
    memory: runtime_attr.mem_gb + " GiB"
    disks: "local-disk " + runtime_attr.disk_gb + " HDD"
    bootDiskSizeGb: runtime_attr.boot_disk_gb
    preemptible: runtime_attr.preemptible_tries
    maxRetries: runtime_attr.max_retries
  }
}
