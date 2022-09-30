version 1.0

import "Structs.wdl"
import "SVConcordance.wdl" as conc
import "TasksMakeCohortVcf.wdl" as tasks_cohort

workflow SVConcordancePacBio {
  input {
    String cohort
    File vcf

    # Sample ids as they appear in the vcfs. Pacbio vcf samples will be renamed
    Array[String] sample_ids
    Array[File] pacbio_sample_vcfs
    String pacbio_tool_name
    Int? pacbio_min_size

    File ploidy_table

    File reference_fasta
    File reference_fasta_fai
    File reference_dict

    File? svtk_to_gatk_formatter_script
    File? main_vcf_preprocess_script
    File? pacbio_vcf_preprocess_script

    String gatk_docker
    String sv_pipeline_docker

    Float? java_mem_fraction

    RuntimeAttr? runtime_attr_combine_truth
    RuntimeAttr? runtime_attr_sv_concordance
  }

  scatter (i in range(length(sample_ids))) {
    call PrepSampleVcfs {
      input:
        sample_id=sample_ids[i],
        main_vcf=vcf,
        main_vcf_index=vcf + ".tbi",
        truth_vcf=pacbio_sample_vcfs[i],
        min_size=select_first([pacbio_min_size, 25]),
        truth_tool_name=pacbio_tool_name,
        ploidy_table=ploidy_table,
        output_prefix="~{cohort}.prep_sample_vcfs.~{sample_ids[i]}",
        svtk_to_gatk_formatter_script=svtk_to_gatk_formatter_script,
        main_script=main_vcf_preprocess_script,
        truth_script=pacbio_vcf_preprocess_script,
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_attr_combine_truth
    }
    call conc.SVConcordanceTask {
      input:
        eval_vcf=PrepSampleVcfs.main_out,
        truth_vcf=PrepSampleVcfs.truth_out,
        output_prefix="~{cohort}.concordance.~{sample_ids[i]}",
        reference_dict=reference_dict,
        java_mem_fraction=java_mem_fraction,
        gatk_docker=gatk_docker,
        runtime_attr_override=runtime_attr_sv_concordance
    }
  }

  output {
    Array[File] pacbio_sample_concordance_vcfs = SVConcordanceTask.out
    Array[File] pacbio_sample_concordance_vcf_indexes = SVConcordanceTask.out_index
  }
}

task PrepSampleVcfs {
  input {
    String sample_id
    File main_vcf
    File main_vcf_index
    File truth_vcf
    String truth_tool_name
    Int min_size
    File ploidy_table
    File? svtk_to_gatk_formatter_script
    File? main_script
    File? truth_script
    String output_prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 3.75,
                               disk_gb: ceil(50 + size(truth_vcf, "GB") * 3 + size(main_vcf, "GB")),
                               boot_disk_gb: 10,
                               preemptible_tries: 1,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File truth_out = "~{output_prefix}.~{truth_tool_name}.vcf.gz"
    File truth_out_index = "~{output_prefix}.~{truth_tool_name}.vcf.gz.tbi"
    File main_out = "~{output_prefix}.main.vcf.gz"
    File main_out_index = "~{output_prefix}.main.vcf.gz.tbi"
  }
  command <<<
    set -euxo pipefail

    # Subset to DEL/DUP/INS and convert DUP to INS
    bcftools view -s "~{sample_id}" ~{main_vcf} | bcftools view --min-ac 1 -Oz -o tmp1.vcf.gz
    python ~{default="/opt/sv-pipeline/scripts/format_svtk_vcf_for_gatk.py" svtk_to_gatk_formatter_script} \
      --vcf tmp1.vcf.gz \
      --out tmp2.vcf.gz \
      --ploidy-table ~{ploidy_table}
    python ~{default="/opt/sv-pipeline/scripts/preprocess_gatk_for_pacbio_eval.py" main_script} tmp2.vcf.gz \
      | bcftools sort -Oz -o ~{output_prefix}.main.vcf.gz
    tabix ~{output_prefix}.main.vcf.gz

    # Convert format
    bcftools reheader -s <(echo "~{sample_id}") ~{truth_vcf} > tmp1.vcf.gz
    python ~{default="/opt/sv-pipeline/scripts/format_pb_for_gatk.py" truth_script} \
      --vcf tmp1.vcf.gz \
      --algorithm ~{truth_tool_name} \
      --min-size ~{min_size} \
      --out tmp2.vcf.gz \
      --ploidy-table ~{ploidy_table}
    bcftools sort tmp2.vcf.gz \
      | bcftools annotate --set-id '~{truth_tool_name}_%CHROM\_%POS\_%END\_%SVTYPE\_%SVLEN' -Oz -o ~{output_prefix}.~{truth_tool_name}.vcf.gz
    tabix ~{output_prefix}.~{truth_tool_name}.vcf.gz
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
