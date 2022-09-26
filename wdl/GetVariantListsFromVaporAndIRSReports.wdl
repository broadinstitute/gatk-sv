version 1.0

import "Structs.wdl"

workflow GetVariantListsFromVaporAndIRSReports {

  input {
    File vcf
    Array[String] vapor_sample_ids
    Array[File] vapor_files
    File irs_test_report
    Float? vapor_min_precision = 0.01
    Int? vapor_max_cnv_size = 5000
    Int? irs_min_cnv_size = 50000
    Float? irs_min_pvalue = 0.001
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

  call GetVariantListsFromVaporAndIRSReports {
    input:
      vcf=vcf,
      vapor_sample_ids=vapor_sample_ids,
      vapor_files=vapor_files,
      irs_test_report=irs_test_report,
      vapor_max_cnv_size=vapor_max_cnv_size,
      vapor_min_precision=vapor_min_precision,
      irs_min_cnv_size=irs_min_cnv_size,
      irs_min_pvalue=irs_min_pvalue,
      irs_min_probes=irs_min_probes,
      output_prefix=output_prefix_,
      sv_utils_docker=sv_utils_docker,
      runtime_attr_override=runtime_attr
  }

  output {
    File output_json = GetVariantListsFromVaporAndIRSReports.output_json
  }
}

task GetVariantListsFromVaporAndIRSReports {
  input {
    File vcf
    Array[String] vapor_sample_ids
    Array[File] vapor_files
    File irs_test_report
    Int? vapor_max_cnv_size
    Float? vapor_min_precision
    Int? irs_min_cnv_size
    Float? irs_min_pvalue
    Int? irs_min_probes
    String output_prefix
    String sv_utils_docker
    RuntimeAttr? runtime_attr_override
  }

  String vapor_json = "vapor_data.json"

  output {
    File output_json = "${output_prefix}.json"
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

    /opt/sv_utils/src/sv_utils/get_confident_variant_lists_from_vapor_and_irs_test.py \
      --vapor-json ~{vapor_json}  \
      --vcf ~{vcf} \
      --irs-test-report ~{irs_test_report} \
      ~{"--vapor-max-cnv-size " + vapor_max_cnv_size} \
      ~{"--vapor-min-precision " + vapor_min_precision} \
      ~{"--irs-min-cnv-size " + irs_min_cnv_size} \
      ~{"--irs-min-probes " + irs_min_probes} \
      ~{"--irs-pvalue-threshold " + irs_min_pvalue} \
      -O ~{output_prefix}.json
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
