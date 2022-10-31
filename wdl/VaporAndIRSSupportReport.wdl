version 1.0

import "Structs.wdl"

workflow VaporAndIRSSupportReport {

  input {
    File vcf
    Array[String] vapor_sample_ids
    Array[File] vapor_files
    Array[File] irs_sample_batches
    Array[File] irs_test_reports
    File? irs_contigs_fai
    Float? vapor_min_precision = 0.99
    Int? vapor_max_cnv_size = 5000
    Int? irs_min_cnv_size = 50000
    Float? irs_good_pvalue_threshold = 0.001
    Int? irs_min_probes = 5
    String? output_prefix
    String sv_utils_docker
    RuntimeAttr? runtime_attr
  }

  String output_prefix_ =
    if defined(output_prefix) then
      select_first([output_prefix])
    else
        basename(vcf, ".vcf.gz")

  call VaporAndIRSSupportReport {
    input:
      vcf=vcf,
      vapor_sample_ids=vapor_sample_ids,
      vapor_files=vapor_files,
      irs_sample_batches=irs_sample_batches,
      irs_test_reports=irs_test_reports,
      irs_contigs_fai=irs_contigs_fai,
      vapor_max_cnv_size=vapor_max_cnv_size,
      vapor_min_precision=vapor_min_precision,
      irs_min_cnv_size=irs_min_cnv_size,
      irs_good_pvalue_threshold=irs_good_pvalue_threshold,
      irs_min_probes=irs_min_probes,
      output_prefix=output_prefix_,
      sv_utils_docker=sv_utils_docker,
      runtime_attr_override=runtime_attr
  }

  output {
    File summary_report = VaporAndIRSSupportReport.summary
    File detail_report = VaporAndIRSSupportReport.detail
  }
}

task VaporAndIRSSupportReport {
  input {
    File vcf
    Array[String] vapor_sample_ids
    Array[File] vapor_files
    Array[File] irs_sample_batches
    Array[File] irs_test_reports
    File? irs_contigs_fai
    Int? vapor_max_cnv_size
    Float? vapor_min_precision
    Int? irs_min_cnv_size
    Float? irs_good_pvalue_threshold
    Int? irs_min_probes
    String output_prefix
    String sv_utils_docker
    RuntimeAttr? runtime_attr_override
  }

  String vapor_json = "vapor_data.json"

  output {
    File summary = "${output_prefix}.irs_vapor_support.summary.tsv"
    File detail = "${output_prefix}.irs_vapor_support.detail.tsv.gz"
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

    /opt/sv_utils/src/sv_utils/report_confident_irs_vapor_variants.py \
      --vapor-json ~{vapor_json}  \
      --vcf ~{vcf} \
      --irs-sample-batch-lists ~{write_lines(irs_sample_batches)} \
      --irs-test-report-list ~{write_lines(irs_test_reports)} \
      ~{"--irs-contigs-file " + irs_contigs_fai} \
      ~{"--vapor-max-cnv-size " + vapor_max_cnv_size} \
      ~{"--vapor-min-precision " + vapor_min_precision} \
      ~{"--irs-min-cnv-size " + irs_min_cnv_size} \
      ~{"--irs-min-probes " + irs_min_probes} \
      ~{"--irs-good-pvalue-threshold " + irs_good_pvalue_threshold} \
      --output-summary ~{output_prefix}.irs_vapor_support.summary.tsv \
      --output-detail ~{output_prefix}.irs_vapor_support.detail.tsv.gz
  >>>

  RuntimeAttr runtime_attr_report_default = object {
    cpu_cores: 1,
    mem_gb: 16,
    boot_disk_gb: 10,
    preemptible_tries: 1,
    max_retries: 1,
    disk_gb: 100 + ceil(size([
      vcf], "GiB")) * 2
  }
  RuntimeAttr runtime_attr = select_first([
    runtime_attr_override,
    runtime_attr_report_default])

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
