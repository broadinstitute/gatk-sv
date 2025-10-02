version 1.0

import "Structs.wdl"
import "SVConcordancePacBioSample.wdl" as concordance
import "Utils.wdl" as utils
import "TasksMakeCohortVcf.wdl" as mini_tasks

workflow MakeGqRecalibratorTrainingSetFromPacBio {

  input {
    # Cleaned GATK-formatted vcf
    # SVConcordance should be run first if the training set is a proper subset of the cohort
    # This can be either a single whole-genome vcf or multiple vcf shards.
    # Assumes all vcfs have indexes, i.e. at {VCF_PATH}.tbi
    Array[File] vcfs

    Array[File]? vapor_files
    Array[File]? pbsv_vcfs
    Array[File]? pav_vcfs
    Array[File]? sniffles_vcfs
    Array[File]? hifi_cnv_vcfs
    Array[File]? hapdiff_vcfs
    
    File training_sample_ids  # Sample IDs with PacBio or array data
    Array[String] pacbio_sample_ids  # Corresponding to files below (must be a subset of training_sample_ids)
    
    String? output_prefix
    File ploidy_table

    # Optional: array intensity ratio files
    Array[File]? irs_sample_batches
    Array[File]? irs_test_reports
    File? irs_contigs_fai

    Boolean write_detail_report = true
    Int? vapor_max_cnv_size = 5000
    Int? vapor_max_non_cnv_size = 5000
    Float? vapor_min_precision = 0.999
    Int? vapor_pos_read_threshold = 2
    Int? irs_min_cnv_size = 10000
    Float? irs_good_pvalue_threshold = 0.000001
    Int? irs_min_probes = 5

    Float? pesr_interval_overlap_strict = 0.1
    Float? pesr_size_similarity_strict = 0.5
    Int? pesr_breakend_window_strict = 5000

    Float? depth_interval_overlap_strict = 0.8
    Float? depth_size_similarity_strict = 0.0
    Int? depth_breakend_window_strict = 10000000

    Float? mixed_interval_overlap_strict = 0.8
    Float? mixed_size_similarity_strict = 0.0
    Int? mixed_breakend_window_strict = 10000000

    File? clustering_config_strict
    File? stratification_config_strict
    Array[String]? track_names_strict
    Array[File]? track_intervals_strict

    Float? pesr_interval_overlap_loose = 0
    Float? pesr_size_similarity_loose = 0
    Int? pesr_breakend_window_loose = 5000

    Float? depth_interval_overlap_loose = 0.5
    Float? depth_size_similarity_loose = 0.0
    Int? depth_breakend_window_loose = 10000000

    Float? mixed_interval_overlap_loose = 0.5
    Float? mixed_size_similarity_loose = 0.0
    Int? mixed_breakend_window_loose = 10000000

    File? clustering_config_loose
    File? stratification_config_loose
    Array[String]? track_names_loose
    Array[File]? track_intervals_loose

    File reference_dict

    String sv_utils_docker
    String gatk_docker
    String sv_base_mini_docker
    String sv_pipeline_docker
    String linux_docker
  }

  Array[String] tool_names = flatten([
    if defined(pbsv_vcfs) then ["pbsv"] else [],
    if defined(pav_vcfs) then ["pav"] else [],
    if defined(sniffles_vcfs) then ["sniffles"] else [],
    if defined(hifi_cnv_vcfs) then ["hifi_cnv"] else [],
    if defined(hapdiff_vcfs) then ["hapdiff"] else []
  ])

  Array[Array[File]] tool_vcfs_list = select_all([pbsv_vcfs, pav_vcfs, sniffles_vcfs, hifi_cnv_vcfs, hapdiff_vcfs])

  Array[Array[File]] pacbio_vcfs = transpose(tool_vcfs_list)

  String output_prefix_ = if defined(output_prefix) then select_first([output_prefix]) else basename(vcfs[0], ".vcf.gz")

  scatter (i in range(length(vcfs))) {
    call utils.SubsetVcfBySamplesList as SubsetTrainingSamples {
      input:
        vcf = vcfs[i],
        vcf_idx = vcfs[i] + ".tbi",
        list_of_samples = training_sample_ids,
        outfile_name = "~{output_prefix_}.training_samples.shard_~{i}.vcf.gz",
        remove_samples = false,
        remove_private_sites = true,
        keep_af = true,
        sv_base_mini_docker = sv_base_mini_docker
    }
  }

  call mini_tasks.ConcatVcfs as ConcatTrainingSampleVcfs {
    input:
      vcfs=SubsetTrainingSamples.vcf_subset,
      vcfs_idx=SubsetTrainingSamples.vcf_subset_index,
      naive=true,
      outfile_prefix="~{output_prefix_}.concat_training_sample_vcfs",
      sv_base_mini_docker=sv_base_mini_docker
  }

  if (defined(vapor_files)) {
    call GetVariantListsFromVaporAndIRS {
      input:
        vcf=ConcatTrainingSampleVcfs.concat_vcf,
        output_prefix=output_prefix_,
        vapor_sample_ids=pacbio_sample_ids,
        vapor_files=select_first([vapor_files, []]),
        irs_sample_batches=select_first([irs_sample_batches, []]),
        irs_test_reports=select_first([irs_test_reports, []]),
        irs_contigs_fai=irs_contigs_fai,
        vapor_max_cnv_size=vapor_max_cnv_size,
        vapor_min_precision=vapor_min_precision,
        vapor_pos_read_threshold=vapor_pos_read_threshold,
        irs_min_cnv_size=irs_min_cnv_size,
        irs_good_pvalue_threshold=irs_good_pvalue_threshold,
        irs_min_probes=irs_min_probes,
        sv_utils_docker=sv_utils_docker
    }

    call VaporAndIRSSupportReport {
      input:
        vcf=ConcatTrainingSampleVcfs.concat_vcf,
        output_prefix=output_prefix_,
        vapor_sample_ids=pacbio_sample_ids,
        vapor_files=select_first([vapor_files, []]),
        write_detail_report=write_detail_report,
        irs_sample_batches=select_first([irs_sample_batches, []]),
        irs_test_reports=select_first([irs_test_reports, []]),
        irs_contigs_fai=irs_contigs_fai,
        vapor_max_cnv_size=vapor_max_cnv_size,
        vapor_min_precision=vapor_min_precision,
        vapor_pos_read_threshold=vapor_pos_read_threshold,
        irs_min_cnv_size=irs_min_cnv_size,
        irs_good_pvalue_threshold=irs_good_pvalue_threshold,
        irs_min_probes=irs_min_probes,
        sv_utils_docker=sv_utils_docker
    }
  }

  call utils.WriteLines as WritePacBioSampleIds {
    input:
      lines = pacbio_sample_ids,
      output_filename = "~{output_prefix_}.pacbio_sample_ids.list",
      linux_docker = linux_docker
  }

  scatter (i in range(length(vcfs))) {
    call utils.SubsetVcfBySamplesList as SubsetPacBioSamples {
      input:
        vcf = vcfs[i],
        vcf_idx = vcfs[i] + ".tbi",
        list_of_samples = WritePacBioSampleIds.out,
        outfile_name = "~{output_prefix_}.pacbio_samples.shard_~{i}.vcf.gz",
        remove_samples = false,
        remove_private_sites = true,
        keep_af = true,
        sv_base_mini_docker = sv_base_mini_docker
    }
  }

  call mini_tasks.ConcatVcfs as ConcatPacbioSampleVcfs {
    input:
      vcfs=SubsetPacBioSamples.vcf_subset,
      vcfs_idx=SubsetPacBioSamples.vcf_subset_index,
      naive=true,
      outfile_prefix="~{output_prefix_}.concat_pacbio_sample_vcfs",
      sv_base_mini_docker=sv_base_mini_docker
  }

  scatter (i in range(length(pacbio_sample_ids))) {
    call PrepSampleVcf {
      input:
        sample_id=pacbio_sample_ids[i],
        vcf=ConcatPacbioSampleVcfs.concat_vcf,
        vcf_index=ConcatPacbioSampleVcfs.concat_vcf_idx,
        max_variant_size=vapor_max_cnv_size,
        output_prefix=output_prefix_,
        sv_pipeline_docker=sv_pipeline_docker
    }

    call concordance.SVConcordancePacBioSample as SVConcordanceLoose {
      input:
        sample_id=pacbio_sample_ids[i],
        sample_vcf=PrepSampleVcf.out,
        pacbio_sample_vcfs=pacbio_vcfs[i],
        tool_names=tool_names,
        prefix="~{output_prefix_}.loose",
        ploidy_table=ploidy_table,
        reference_dict=reference_dict,
        pesr_interval_overlap=pesr_interval_overlap_loose,
        pesr_size_similarity=pesr_size_similarity_loose,
        pesr_breakend_window=pesr_breakend_window_loose,
        depth_interval_overlap=depth_interval_overlap_loose,
        depth_size_similarity=depth_size_similarity_loose,
        depth_breakend_window=depth_breakend_window_loose,
        mixed_interval_overlap=mixed_interval_overlap_loose,
        mixed_size_similarity=mixed_size_similarity_loose,
        mixed_breakend_window=mixed_breakend_window_loose,
        clustering_config=clustering_config_loose,
        stratification_config=stratification_config_loose,
        track_names=track_names_loose,
        track_intervals=track_intervals_loose,
        sv_base_mini_docker=sv_base_mini_docker,
        sv_pipeline_docker=sv_pipeline_docker,
        gatk_docker=gatk_docker,
        linux_docker=linux_docker
    }

    call concordance.SVConcordancePacBioSample as SVConcordanceStrict {
      input:
        sample_id=pacbio_sample_ids[i],
        sample_vcf=PrepSampleVcf.out,
        pacbio_sample_vcfs=pacbio_vcfs[i],
        tool_names=tool_names,
        prefix="~{output_prefix_}.strict",
        ploidy_table=ploidy_table,
        reference_dict=reference_dict,
        pesr_interval_overlap=pesr_interval_overlap_strict,
        pesr_size_similarity=pesr_size_similarity_strict,
        pesr_breakend_window=pesr_breakend_window_strict,
        depth_interval_overlap=depth_interval_overlap_strict,
        depth_size_similarity=depth_size_similarity_strict,
        depth_breakend_window=depth_breakend_window_strict,
        mixed_interval_overlap=mixed_interval_overlap_strict,
        mixed_size_similarity=mixed_size_similarity_strict,
        mixed_breakend_window=mixed_breakend_window_strict,
        clustering_config=clustering_config_strict,
        stratification_config=stratification_config_strict,
        track_names=track_names_strict,
        track_intervals=track_intervals_strict,
        sv_base_mini_docker=sv_base_mini_docker,
        sv_pipeline_docker=sv_pipeline_docker,
        gatk_docker=gatk_docker,
        linux_docker=linux_docker
    }

    if (defined(vapor_files)) {
      call RefineSampleLabels {
        input:
          sample_id=pacbio_sample_ids[i],
          main_vcf=PrepSampleVcf.out,
          vapor_json=select_first([GetVariantListsFromVaporAndIRS.vapor_json]),
          tool_names=tool_names,
          loose_concordance_vcfs=SVConcordanceLoose.pacbio_concordance_vcfs,
          strict_concordance_vcfs=SVConcordanceStrict.pacbio_concordance_vcfs,
          vapor_max_cnv_size=vapor_max_cnv_size,
          vapor_max_non_cnv_size=vapor_max_non_cnv_size,
          output_prefix="~{output_prefix_}.gq_training_labels.~{pacbio_sample_ids[i]}",
          sv_pipeline_docker=sv_pipeline_docker
      }
    }
  }

  if (defined(vapor_files)) {
    call MergeJsons {
      input:
        jsons=select_all(flatten([[GetVariantListsFromVaporAndIRS.irs_json], RefineSampleLabels.out_json])),
        output_prefix="~{output_prefix_}.gq_training_labels",
        sv_pipeline_docker=sv_pipeline_docker
    }

    call mini_tasks.ConcatHeaderedTextFiles {
      input:
        text_files=select_all(RefineSampleLabels.out_table),
        gzipped=false,
        output_filename="~{output_prefix_}.gq_training_labels.tsv",
        linux_docker=linux_docker
    }
  }

  output {
    File? gq_recalibrator_training_json = MergeJsons.out
    File? pacbio_support_summary_table = ConcatHeaderedTextFiles.out

    File training_sample_vcf = ConcatTrainingSampleVcfs.concat_vcf
    File training_sample_vcf_index = ConcatTrainingSampleVcfs.concat_vcf_idx
    File pacbio_sample_vcf = ConcatPacbioSampleVcfs.concat_vcf
    File pacbio_sample_vcf_index = ConcatPacbioSampleVcfs.concat_vcf_idx

    File? vapor_output_json = GetVariantListsFromVaporAndIRS.vapor_json
    File? irs_output_json = GetVariantListsFromVaporAndIRS.irs_json
    File? vapor_and_irs_summary_report = VaporAndIRSSupportReport.summary
    File? vapor_and_irs_detail_report = VaporAndIRSSupportReport.detail
    Array[File] loose_pacbio_concordance_vcf_tars = SVConcordanceLoose.pacbio_concordance_vcfs_tar
    Array[File] strict_pacbio_concordance_vcf_tars = SVConcordanceStrict.pacbio_concordance_vcfs_tar
  }
}

task GetVariantListsFromVaporAndIRS {
  input {
    File vcf
    String output_prefix
    Array[String] vapor_sample_ids
    Array[File] vapor_files

    Array[File] irs_sample_batches
    Array[File] irs_test_reports
    File? irs_contigs_fai

    Int? vapor_max_cnv_size
    Float? vapor_min_precision
    String vapor_strategy = "READS"
    Int? vapor_pos_read_threshold
    Int? vapor_neg_read_threshold
    Int? vapor_neg_cov_read_threshold
    Int? irs_min_cnv_size
    Float? irs_good_pvalue_threshold
    Float? irs_bad_pvalue_threshold
    Int? irs_min_probes

    File? script
    String? additional_args

    String sv_utils_docker
    RuntimeAttr? runtime_attr_override
  }

  output {
    File vapor_json = "~{output_prefix}.gq_recalibrator_labels.vapor.json"
    File irs_json = "~{output_prefix}.gq_recalibrator_labels.irs.json"
  }

  command <<<
    set -euxo pipefail

    # JSON preprocessing taken from GetTruthOverlap.wdl
    # construct a vapor JSON file if vapor data was passed
    VAPOR_SAMPLE_IDS=~{write_lines(vapor_sample_ids)}
    VAPOR_FILES=~{write_lines(vapor_files)}
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
      > vapor_data.json
    printf "vapor_data.json: "
    cat vapor_data.json

    IRS_SAMPLE_BATCHES=~{write_lines(irs_sample_batches)}
    IRS_TEST_REPORTS=~{write_lines(irs_test_reports)}
    python ~{default="/opt/sv_utils/src/sv_utils/get_confident_variant_lists_from_vapor_and_irs_test.py" script} \
      --vapor-json vapor_data.json  \
      --vcf ~{vcf} \
      ~{if length(irs_sample_batches) > 0 then "--irs-sample-batch-lists $IRS_SAMPLE_BATCHES" else ""} \
      ~{if length(irs_test_reports) > 0 then "--irs-test-report-list $IRS_TEST_REPORTS" else ""} \
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
      --irs-output ~{output_prefix}.gq_recalibrator_labels.irs.json \
      --vapor-output ~{output_prefix}.gq_recalibrator_labels.vapor.json \
      ~{additional_args}
  >>>

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 4,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1,
    disk_gb: 25 + ceil(size([
      vcf], "GiB"))
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_utils_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

task PrepSampleVcf {
  input {
    String sample_id
    File vcf
    File vcf_index
    Int? max_variant_size
    File? script
    String output_prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 3.75,
                               disk_gb: ceil(50 + size(vcf, "GB")),
                               boot_disk_gb: 10,
                               preemptible_tries: 1,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File out = "~{output_prefix}.~{sample_id}.prepped_for_pb.vcf.gz"
    File out_index = "~{output_prefix}.~{sample_id}.prepped_for_pb.vcf.gz.tbi"
  }
  command <<<
    set -euxo pipefail

    # Subset to DEL/DUP/INS and convert DUP to INS
    bcftools view -s "~{sample_id}" ~{vcf} | bcftools view --no-update --min-ac 1 -Oz -o tmp1.vcf.gz
    python ~{default="/opt/sv-pipeline/scripts/preprocess_gatk_for_pacbio_eval.py" script} \
      tmp1.vcf.gz \
      --max-variant-size ~{max_variant_size} \
      | bcftools sort -Oz -o ~{output_prefix}.~{sample_id}.prepped_for_pb.vcf.gz
    tabix ~{output_prefix}.~{sample_id}.prepped_for_pb.vcf.gz
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


task VaporAndIRSSupportReport {
  input {
    File vcf
    String output_prefix
    Array[String] vapor_sample_ids
    Array[File] vapor_files
    Array[File] irs_sample_batches
    Array[File] irs_test_reports
    Boolean write_detail_report
    File? irs_contigs_fai
    Int? vapor_max_cnv_size
    Float? vapor_min_precision
    Int? vapor_pos_read_threshold
    Int? irs_min_cnv_size
    Float? irs_good_pvalue_threshold
    Int? irs_min_probes

    File? script

    String sv_utils_docker
    RuntimeAttr? runtime_attr_override
  }

  output {
    File summary = "${output_prefix}.irs_vapor_support.summary.tsv"
    File detail = "${output_prefix}.irs_vapor_support.detail.tsv.gz"
  }

  command <<<
    set -euxo pipefail

    # JSON preprocessing taken from GetTruthOverlap.wdl
    # construct a vapor JSON file if vapor data was passed
    VAPOR_SAMPLE_IDS=~{write_lines(vapor_sample_ids)}
    VAPOR_FILES=~{write_lines(vapor_files)}
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
      > vapor_data.json
    printf "vapor_data.json: "
    cat vapor_data.json

    touch ~{output_prefix}.irs_vapor_support.detail.tsv.gz
    IRS_SAMPLE_BATCHES=~{write_lines(irs_sample_batches)}
    IRS_TEST_REPORTS=~{write_lines(irs_test_reports)}
    python ~{default="/opt/sv_utils/src/sv_utils/report_confident_irs_vapor_variants.py" script} \
      --vapor-json vapor_data.json  \
      --vcf ~{vcf} \
      ~{if length(irs_sample_batches) > 0 then "--irs-sample-batch-lists $IRS_SAMPLE_BATCHES" else ""} \
      ~{if length(irs_test_reports) > 0 then "--irs-test-report-list $IRS_TEST_REPORTS" else ""} \
      ~{"--irs-contigs-file " + irs_contigs_fai} \
      ~{"--vapor-max-cnv-size " + vapor_max_cnv_size} \
      ~{"--vapor-min-precision " + vapor_min_precision} \
      ~{"--irs-min-cnv-size " + irs_min_cnv_size} \
      ~{"--irs-min-probes " + irs_min_probes} \
      ~{"--irs-good-pvalue-threshold " + irs_good_pvalue_threshold} \
      --output-summary ~{output_prefix}.irs_vapor_support.summary.tsv \
      ~{"--vapor-read-support-pos-thresh " + vapor_pos_read_threshold} \
      ~{if write_detail_report then "--output-detail ~{output_prefix}.irs_vapor_support.detail.tsv.gz" else ""}
  >>>

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 16,
    boot_disk_gb: 10,
    preemptible_tries: 1,
    max_retries: 1,
    disk_gb: 100 + ceil(size([
      vcf], "GiB")) * 2
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_utils_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}


task RefineSampleLabels {
  input {
    String sample_id
    File main_vcf
    File vapor_json
    Array[String] tool_names
    Array[File] loose_concordance_vcfs
    Array[File] strict_concordance_vcfs
    Int? vapor_max_cnv_size
    Int? vapor_max_non_cnv_size
    String? additional_args
    File? script
    String output_prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 3.75,
                               disk_gb: 10,
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  String additional_args_ = if defined(additional_args) then select_first([additional_args]) else ""

  output {
    File out_json = "~{output_prefix}.json"
    File out_table = "~{output_prefix}.tsv"
  }
  command <<<
    set -euo pipefail
    python ~{default="/opt/sv-pipeline/scripts/refine_training_set.py" script} \
      --loose-concordance-vcfs ~{sep=" " loose_concordance_vcfs} \
      --strict-concordance-vcfs ~{sep=" " strict_concordance_vcfs} \
      --main-vcf ~{main_vcf} \
      --vapor-json ~{vapor_json} \
      --sample-id ~{sample_id} \
      --cnv-size-cutoff ~{vapor_max_cnv_size} \
      --non-cnv-size-cutoff ~{vapor_max_non_cnv_size} \
      --json-out ~{output_prefix}.json \
      --table-out ~{output_prefix}.tsv \
      --truth-algorithms ~{sep="," tool_names} \
      ~{additional_args_}
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

task MergeJsons {
  input {
    Array[File] jsons
    String output_prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 7.5,
                               disk_gb: ceil(10 + size(jsons, "GB") * 2),
                               boot_disk_gb: 10,
                               preemptible_tries: 1,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File out = "~{output_prefix}.json"
  }
  command <<<
    set -euo pipefail
    python3 <<CODE
    import json
    with open('~{write_lines(jsons)}') as f:
        paths = [line.strip() for line in f]
    data = {}
    for p in paths:
        with open(p) as f:
          for key1, value1 in json.load(f).items():
              if key1 not in data:
                  data[key1] = value1
              else:
                  for key2, value2 in value1.items():
                    data[key1][key2] = data[key1][key2] + value2
    with open('~{output_prefix}.json', 'w') as f:
        f.write(json.dumps(data))
    CODE
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
