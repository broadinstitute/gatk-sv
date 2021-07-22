##  ExpansionHunter denovo (EHdn)
##
##  This WDL implements workflows for EHdn's two execution modes:
##  - case-control analysis (read EHdn docs on this:
##    https://github.com/Illumina/ExpansionHunterDenovo/blob/master/documentation/03_Case_control_quickstart.md);
##  - outlier analysis (read EHdn docs on this:
##    https://github.com/Illumina/ExpansionHunterDenovo/blob/master/documentation/04_Outlier_quickstart.md).

version 1.0

import "Structs.wdl"

struct FilenamePostfixes {
  String locus
  String motif
  String profile
  String merged_profile
  Int profile_len
}

workflow EHdnSTRAnalysis {

  input {
    String analysis_type
    String str_comparison_type
    Array[File] sample_bams_or_crams
    Array[File]? sample_bams_or_crams_indexes
    Array[String] samples_status
    File reference_fasta
    File? reference_fasta_index
    Int min_anchor_mapq
    Int max_irr_mapq
    String ehdn_docker
    RuntimeAttr? runtime_attr_str_profile
    RuntimeAttr? runtime_attr_merge
    RuntimeAttr? runtime_attr_analysis
  }

  parameter_meta {
    analysis_type: "Sets the analysis type; accepted values are: `casecontrol`, `outlier`, and `both`."
    str_comparison_type: "Set the STR comparison type; accepted values are: `locus`, `motif`, and `both`."
    case_indexes_filenames: "A list of index files (.bam.bai) for the `case` samples. Files should be in the same order as the case samples"
    manifest_filename: ""
    reference_filename: ""
    min_anchor_mapq: ""
    max_irr_mapq: ""
    ehdn_docker: ""
  }

  # The values of the variables are based on
  # EHdn's current latest hard-coded postfixes.
  # Do not change them unless they are changed
  # in EHdn.
  FilenamePostfixes postfixes = object {
    locus: ".locus.tsv",
    motif: ".motif.tsv",
    profile: ".str_profile.json",
    merged_profile: ".multisample_profile.json",

    # This the length of `profile` postfix without
    # the filename extension. It is used to remove
    # postfix in order to extract sample name from
    # EHdn generated output.
    # e.g., extract `sample1` from `sample1.str_profile.json`.
    profile_len: 12
  }

  scatter (i in range(length(sample_bams_or_crams))) {
    File sample_bam_or_cram_ = sample_bams_or_crams[i]
    File bam_or_cram_index_ =
      if defined(sample_bams_or_crams_indexes) then
        select_first([sample_bams_or_crams_indexes])[i]
      else
        sample_bam_or_cram_ +
        if  basename(sample_bam_or_cram_, ".bam") + ".bam" ==
            basename(sample_bam_or_cram_) then
          ".bai"
        else
          ".crai"

    String filename =
      if  basename(sample_bam_or_cram_, ".bam") + ".bam" ==
          basename(sample_bam_or_cram_) then
        basename(sample_bam_or_cram_, ".bam")
      else
        basename(sample_bam_or_cram_, ".cram")

    call ComputeSTRProfile {
      input:
        filename = filename,
        bam_or_cram = sample_bam_or_cram_,
        bam_or_cram_index = bam_or_cram_index_,
        reference_fasta = reference_fasta,
        reference_fasta_index = reference_fasta_index,
        min_anchor_mapq = min_anchor_mapq,
        max_irr_mapq = max_irr_mapq,
        ehdn_docker = ehdn_docker,
        runtime_attr_override = runtime_attr_str_profile,
        postfixes = postfixes
    }
  }

  scatter (i in range(length(samples_status))) {
    if (samples_status[i] == "case") {
      File case_str_profile_json = ComputeSTRProfile.str_profile[i]
    }
    if (samples_status[i] == "control") {
      File control_str_profile_json = ComputeSTRProfile.str_profile[i]
    }
  }

  Array[File] cases_str_profile_json = select_all(case_str_profile_json)
  Array[File] controls_str_profile_json = select_all(control_str_profile_json)

  call Merge {
    input:
      cases = cases_str_profile_json,
      controls = controls_str_profile_json,
      reference_fasta = reference_fasta,
      reference_fasta_index = reference_fasta_index,
      ehdn_docker = ehdn_docker,
      runtime_attr_override = runtime_attr_merge,
      postfixes = postfixes
  }

  call STRAnalyze {
    input:
      analysis_type = analysis_type,
      str_comparison_type = str_comparison_type,
      manifest = Merge.manifest,
      multisample_profile = Merge.multisample_profile,
      ehdn_docker = ehdn_docker,
      runtime_attr_override = runtime_attr_analysis
  }

  output {
    Array[File] cases_locus = ComputeSTRProfile.locus
    Array[File] cases_motif = ComputeSTRProfile.motif
    Array[File] cases_str_profile = ComputeSTRProfile.str_profile
    File multisample_profile = Merge.multisample_profile
    Array[File] analysis_results = STRAnalyze.results
  }
}

task ComputeSTRProfile {
  input {
    String filename
    File bam_or_cram
    File bam_or_cram_index
    File reference_fasta
    File? reference_fasta_index
    Int min_anchor_mapq
    Int max_irr_mapq
    String ehdn_docker
    RuntimeAttr? runtime_attr_override
    FilenamePostfixes postfixes
  }

  output {
    File locus = "${filename}${postfixes.locus}"
    File motif = "${filename}${postfixes.motif}"
    File str_profile = "${filename}${postfixes.profile}"
  }

  # Defining this varialbe is requires since it
  # will trigger localization of the file. The
  # localized file is used by the `profile` subcommand
  # of ExpansionHunterDenovo.
  File reference_fasta_index_ = select_first([
    reference_fasta_index,
    reference_fasta + ".fai"])

  command <<<
    set -euxo pipefail

    ExpansionHunterDenovo profile \
    --reads ~{bam_or_cram} \
    --reference ~{reference_fasta} \
    --output-prefix ~{filename} \
    --min-anchor-mapq ~{min_anchor_mapq} \
    --max-irr-mapq ~{max_irr_mapq}

  >>>

  RuntimeAttr runtime_attr_str_profile_default = object {
    cpu_cores: 1,
    mem_gb: 4,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([
    runtime_attr_override,
    runtime_attr_str_profile_default])

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
    File reference_fasta
    File? reference_fasta_index
    String ehdn_docker
    RuntimeAttr? runtime_attr_override
    FilenamePostfixes postfixes
  }

  output {
    File manifest = "manifest.tsv"
    File multisample_profile = "${output_prefix}${postfixes.merged_profile}"
  }

  Int cases_length = length(cases)
  Int controls_length = length(controls)
  String output_prefix = "merged"

  # Defining this varialbe is requires since it
  # will trigger localization of the file. The
  # localized file is used by the `profile` subcommand
  # of ExpansionHunterDenovo.
  File reference_fasta_index_ = select_first([
    reference_fasta_index,
    reference_fasta + ".fai"])

  # This shell script uses EHdn's `merge` command
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
  # - Call EHdn's `merge` method using the input and
  #   the generaged manifest file.
  command <<<
    set -euxo pipefail

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
      cases_arr=(~{sep=" " cases})
      for i in "${cases_arr[@]}"; do
        echo -e "$( get_sample_name "$i" )\tcase\t$i" >> $manifest_filename
      done
    fi
    if [ ~{controls_length} -ne 0 ]; then
      controls_arr=(~{sep=" " controls})
      for i in "${controls_arr[@]}"; do
        echo -e "$( get_sample_name "$i" )\tcontrol\t$i" >> $manifest_filename
      done
    fi
    cat $manifest_filename

    ExpansionHunterDenovo merge \
      --reference ~{reference_fasta} \
      --manifest $manifest_filename \
      --output-prefix ~{output_prefix}
  >>>

  RuntimeAttr runtime_attr_merge_default = object {
    cpu_cores: 1,
    mem_gb: 4,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([
    runtime_attr_override,
    runtime_attr_merge_default])

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

task STRAnalyze {
  input {
    String analysis_type
    String str_comparison_type
    File manifest
    File multisample_profile
    String ehdn_docker
    RuntimeAttr? runtime_attr_override
  }

  output {
    Array[File] results = glob("result_*.tsv")
  }

  command <<<
    set -euxo pipefail

    analysis_types=()
    comparison_types=()
    if [ ~{analysis_type} == "both" ]; then
      analysis_types+=("casecontrol")
      analysis_types+=("outlier")
    else
      analysis_types+=("~{analysis_type}")
    fi

    if [ ~{str_comparison_type} == "both" ]; then
      comparison_types+=("locus")
      comparison_types+=("motif")
    else
      comparison_types+=("~{str_comparison_type}")
    fi

    for analysis_type in "${analysis_types[@]}"; do
      for comparison_type in "${comparison_types[@]}"; do
        python ${SCRIPTS_DIR}/${analysis_type}.py ${comparison_type} \
          --manifest ~{manifest} \
          --multisample-profile ~{multisample_profile} \
          --output result_${analysis_type}_${comparison_type}.tsv
      done
    done

    for script in "${scripts[@]}"
    do
      echo "$script"
    done
  >>>

  RuntimeAttr runtime_attr_analysis_default = object {
    cpu_cores: 1,
    mem_gb: 4,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([
    runtime_attr_override,
    runtime_attr_analysis_default])

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
