version 1.0

import "Structs.wdl"
import "TasksClusterBatch.wdl" as tasks_cluster
import "TasksMakeCohortVcf.wdl" as tasks_cohort

workflow RefinePacbioTruthset {
  input {

    String cohort
    File cleaned_vcf

    # Sample ids as they appear in the cleaned vcf
    Array[String] sample_ids
    Array[File] pbsv_vcfs
    Array[File] pav_vcfs
    Array[File] sniffles_vcfs
    File gq_recalibrator_training_json

    Float? pesr_interval_overlap
    Float? pesr_size_similarity
    Int? pesr_breakend_window

    File ploidy_table

    File reference_fasta
    File reference_fasta_fai
    File reference_dict

    File? main_vcf_preprocess_script
    File? truth_vcf_preprocess_script
    File? label_refinement_script

    String gatk_docker
    String sv_base_mini_docker
    String sv_pipeline_docker

    Float? java_mem_fraction

    RuntimeAttr? runtime_attr_combine_truth
    RuntimeAttr? runtime_attr_svcluster
    RuntimeAttr? runtime_attr_refine_labels
    RuntimeAttr? runtime_attr_merge_jsons
  }

  Array[String] tool_names = ["pbsv", "pav", "sniffles"]

  scatter (i in range(length(sample_ids))) {
    call PrepSampleVcfs {
      input:
        sample_id=sample_ids[i],
        main_vcf=cleaned_vcf,
        main_vcf_index=cleaned_vcf + ".tbi",
        truth_vcfs=[pbsv_vcfs[i], pav_vcfs[i], sniffles_vcfs[i]],
        tool_names=tool_names,
        ploidy_table=ploidy_table,
        output_prefix="~{cohort}.prep_sample_vcfs.~{sample_ids[i]}",
        main_script=main_vcf_preprocess_script,
        truth_script=truth_vcf_preprocess_script,
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_attr_combine_truth
    }
    call tasks_cluster.SVCluster {
      input:
        vcfs=[PrepSampleVcfs.truth_out, PrepSampleVcfs.main_out],
        ploidy_table=ploidy_table,
        output_prefix="~{cohort}.svcluster.~{sample_ids[i]}",
        fast_mode=true,
        algorithm="SINGLE_LINKAGE",
        pesr_sample_overlap=0,
        pesr_interval_overlap=select_first([pesr_interval_overlap, 0]),
        pesr_size_similarity=select_first([pesr_size_similarity, 0]),
        pesr_breakend_window=select_first([pesr_breakend_window, 1000]),
        reference_fasta=reference_fasta,
        reference_fasta_fai=reference_fasta_fai,
        reference_dict=reference_dict,
        java_mem_fraction=java_mem_fraction,
        gatk_docker=gatk_docker,
        runtime_attr_override=runtime_attr_svcluster
    }
    call RefineLabels {
      input:
        sample_id=sample_ids[i],
        main_vcf=cleaned_vcf,
        clustered_vcf=SVCluster.out,
        training_json=gq_recalibrator_training_json,
        tool_names=tool_names,
        output_prefix="~{cohort}.refine_labels.~{sample_ids[i]}",
        script=label_refinement_script,
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_attr_refine_labels
    }
  }

  call MergeJsons {
    input:
      jsons=RefineLabels.out,
      output_prefix="~{cohort}.refined_labels",
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_attr_merge_jsons
  }

  output {
    File refined_labels = MergeJsons.out
  }
}

task PrepSampleVcfs {
  input {
    String sample_id
    File main_vcf
    File main_vcf_index
    Array[File] truth_vcfs
    Array[String] tool_names
    File ploidy_table
    File? main_script
    File? truth_script
    String output_prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 3.75,
                               disk_gb: ceil(50 + size(truth_vcfs, "GB") * 3 + size(main_vcf, "GB")),
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File truth_out = "~{output_prefix}.truth.vcf.gz"
    File truth_out_index = "~{output_prefix}.truth.vcf.gz.tbi"
    File main_out = "{output_prefix}.main.vcf.gz"
    File main_out_index = "{output_prefix}.main.vcf.gz.tbi"
  }
  command <<<
    set -euo pipefail

    # Subset to DEL/DUP/INS and convert DUP to INS
    bcftools view -s "~{sample_id}" ~{main_vcf} -Oz -o tmp.vcf.gz
    python ~{default="/opt/sv-pipeline/scripts/preprocess_gatk_for_pacbio_eval.py" main_script} tmp.vcf.gz \
      | bgzip \
      > ~{output_prefix}.main.vcf.gz
    tabix ~{output_prefix}.main.vcf.gz

    # Convert format
    mkdir tmp/
    while read path algorithm; do
      bcftools reheader -s "~{sample_id}" $path | bgzip > tmp1.vcf.gz
      python ~{default="/opt/sv-pipeline/scripts/format_pb_for_gatk.py" truth_script} \
        --vcf tmp1.vcf.gz \
        --algorithm $algorithm \
        --out tmp2.vcf.gz \
        --ploidy-table ~{ploidy_table}
      bcftools sort tmp2.vcf.gz -Oz -o tmp/$algorithm.vcf.gz
      tabix tmp/$algorithm.vcf.gz
    done < <(paste ~{write_lines(truth_vcfs)} ~{write_lines(tool_names)})

    bcftools concat --allow-overlaps tmp/*.vcf.gz -Oz -o ~{output_prefix}.truth.vcf.gz
    tabix ~{output_prefix}.truth.vcf.gz
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

task RefineLabels {
  input {
    String sample_id
    File main_vcf
    File clustered_vcf
    File training_json
    Array[String] tool_names
    String? additional_args
    File? script
    String output_prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 3.75,
                               disk_gb: ceil(50 + size(clustered_vcf, "GB") + size(main_vcf, "GB")),
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  String additional_args_ = if defined(additional_args) then additional_args else ""

  output {
    File out = "~{output_prefix}.json"
  }
  command <<<
    set -euo pipefail
    python ~{default="/opt/sv-pipeline/scripts/refine_training_set.py" script}
      --clustered-vcf ~{clustered_vcf} \
      --main-vcf ~{main_vcf} \
      --truth-json ~{training_json} \
      --sample-id ~{sample_id} \
      --out ~{output_prefix}.json \
      --truth-algorithms ~{sep="," tool_names}
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
                               mem_gb: 3.75,
                               disk_gb: ceil(50 + size(jsons, "GB") * 2),
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
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
            data.update(json.load(f))
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