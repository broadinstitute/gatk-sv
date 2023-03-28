version 1.0

import "Structs.wdl"

workflow MakeGqRecalibratorTrainingSetFromVapor {

  input {
    File vcf
    String? output_prefix
    Array[String] vapor_sample_ids
    Array[File] vapor_files

    # Optional: array intensity ratio files
    Array[File]? irs_sample_batches
    Array[File]? irs_test_reports
    File? irs_contigs_fai

    Float? vapor_min_precision = 0.99
    Int? vapor_max_cnv_size = 5000
    String vapor_strategy = "READS"
    Int? vapor_pos_read_threshold = 2
    Int? vapor_neg_read_threshold = 0
    Int? vapor_neg_cov_read_threshold = 5
    Int? irs_min_cnv_size = 50000
    Float? irs_good_pvalue_threshold = 0.001
    Float? irs_bad_pvalue_threshold = 0.6
    Int? irs_min_probes = 5

    # For debugging
    File? script

    String sv_utils_docker
    RuntimeAttr? runtime_attr
  }

  String output_prefix_ =
    if defined(output_prefix) then
      select_first([output_prefix])
    else
        basename(vcf, ".vcf.gz")

  call GetVariantListsFromVaporAndIRS {
    input:
      vcf=vcf,
      output_prefix=output_prefix_,
      vapor_sample_ids=vapor_sample_ids,
      vapor_files=vapor_files,
      irs_sample_batches=irs_sample_batches,
      irs_test_reports=irs_test_reports,
      irs_contigs_fai=irs_contigs_fai,
      vapor_max_cnv_size=vapor_max_cnv_size,
      vapor_min_precision=vapor_min_precision,
      vapor_strategy=vapor_strategy,
      vapor_pos_read_threshold=vapor_pos_read_threshold,
      vapor_neg_read_threshold=vapor_neg_read_threshold,
      vapor_neg_cov_read_threshold=vapor_neg_cov_read_threshold,
      irs_min_cnv_size=irs_min_cnv_size,
      irs_good_pvalue_threshold=irs_good_pvalue_threshold,
      irs_bad_pvalue_threshold=irs_bad_pvalue_threshold,
      irs_min_probes=irs_min_probes,
      script=script,
      sv_utils_docker=sv_utils_docker,
      runtime_attr_override=runtime_attr
  }

  output {
    File vapor_and_irs_output_json = GetVariantListsFromVaporAndIRS.output_json
  }
}

task GetVariantListsFromVaporAndIRS {
  input {
    File vcf
    String output_prefix
    Array[String] vapor_sample_ids
    Array[File] vapor_files

    Array[File]? irs_sample_batches
    Array[File]? irs_test_reports
    File? irs_contigs_fai

    Int? vapor_max_cnv_size
    Float? vapor_min_precision
    String vapor_strategy
    Int? vapor_pos_read_threshold
    Int? vapor_neg_read_threshold
    Int? vapor_neg_cov_read_threshold
    Int? irs_min_cnv_size
    Float? irs_good_pvalue_threshold
    Float? irs_bad_pvalue_threshold
    Int? irs_min_probes

    File? script

    String sv_utils_docker
    RuntimeAttr? runtime_attr_override
  }

  String vapor_json = "vapor_data.json"

  output {
    File output_json = "${output_prefix}.gq_recalibrator_labels.from_vapor_and_irs.json"
  }


  command <<<
    set -euxo pipefail

    # JSON prepprocessing taken from GetTruthOverlap.wdl
    # construct a vapor JSON file if vapor data was passed
    VAPOR_SAMPLE_IDS=~{write_lines(select_first([vapor_sample_ids]))}
    VAPOR_FILES=~{write_lines(select_first([vapor_files]))}
    # this is all horrible, but it's just taking text processing to turn the sample IDs and vapor files arrays
    # into a json file that contains a map with sample IDs as keys, and corresponding vapor files as values
    {
      echo '{'
      paste -d: $VAPOR_SAMPLE_IDS $VAPOR_FILES \
        | sed -e 's/^/\"/' -e 's/:/\":\"/' -e 's/$/\",/'
      echo '}'
    } \
      | tr -d '\n' \
      | sed 's/,}/}/' \
      | sed -e 's/\(,\|{\)/\1\n/g' -e 's/"}/"\n}\n/' \
      | sed 's/^"/  "/g' \
      > ~{vapor_json}
    printf "~{vapor_json}: "
    cat ~{vapor_json}

    python ~{default="/opt/sv_utils/src/sv_utils/get_confident_variant_lists_from_vapor_and_irs_test.py" script} \
      --vapor-json ~{vapor_json}  \
      --vcf ~{vcf} \
      ~{if defined(irs_sample_batches) then "--irs-sample-batch-lists " + write_lines(select_first([irs_sample_batches])) else ""} \
      ~{if defined(irs_test_reports) then "--irs-test-report-list " + write_lines(select_first([irs_test_reports])) else ""} \
      ~{"--irs-contigs-file " + irs_contigs_fai} \
      ~{"--vapor-max-cnv-size " + vapor_max_cnv_size} \
      ~{"--vapor-min-precision " + vapor_min_precision} \
      ~{"--irs-min-cnv-size " + irs_min_cnv_size} \
      ~{"--irs-min-probes " + irs_min_probes} \
      ~{"--irs-good-pvalue-threshold " + irs_good_pvalue_threshold} \
      ~{"--irs-bad-pvalue-threshold " + irs_bad_pvalue_threshold} \
      --vapor-strategy ~{vapor_strategy} \
      ~{"--vapor-read-support-pos-thresh " + vapor_pos_read_threshold} \
      ~{"--vapor-read-support-neg-thresh " + vapor_neg_read_threshold} \
      ~{"--vapor-read-support-neg-cov-thresh " + vapor_neg_cov_read_threshold} \
      -O ~{output_prefix}.gq_recalibrator_labels.from_vapor_and_irs.json
  >>>

  RuntimeAttr runtime_attr_str_profile_default = object {
    cpu_cores: 1,
    mem_gb: 4,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1,
    disk_gb: 25 + ceil(size([
      vcf], "GiB"))
  }
  RuntimeAttr runtime_attr = select_first([
    runtime_attr_override,
    runtime_attr_str_profile_default])

  runtime {
    docker: sv_utils_docker
    cpu: runtime_attr.cpu_cores
    memory: runtime_attr.mem_gb + " GiB"
    disks: "local-disk " + runtime_attr.disk_gb + " HDD"
    bootDiskSizeGb: runtime_attr.boot_disk_gb
    preemptible: runtime_attr.preemptible_tries
    maxRetries: runtime_attr.max_retries
  }
}
