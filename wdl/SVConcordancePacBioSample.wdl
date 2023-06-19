version 1.0

import "Structs.wdl"
import "SVConcordance.wdl" as conc
import "TasksMakeCohortVcf.wdl" as tasks_cohort
import "Utils.wdl" as utils

workflow SVConcordancePacBioSample {
  input {
    String sample_id  # Sample IDs will be renamed to this in the PacBio vcfs
    File sample_vcf  # Single-sample vcf, usually derived from the cleaned vcf. Sample ID should match sample_id.
    Array[File] pacbio_sample_vcfs  # Raw single-sample vcfs from each tool. Sample ID will be changed to sample_id.
    Array[String] tool_names  # Names of PacBio tools in same order as pacbio_sample_vcfs
    String prefix
    File ploidy_table

    Int? pacbio_min_size

    Float? pesr_interval_overlap
    Float? pesr_size_similarity
    Int? pesr_breakend_window

    File reference_dict

    # For debugging
    File? preprocess_script

    String gatk_docker
    String sv_pipeline_docker
    String linux_docker

    Float? java_mem_fraction

    RuntimeAttr? runtime_attr_combine_truth
    RuntimeAttr? runtime_attr_sv_concordance
    RuntimeAttr? runtime_attr_tar
  }

  scatter (i in range(length(tool_names))) {
    call PrepPacBioVcf {
      input:
        sample_id=sample_id,
        vcf=pacbio_sample_vcfs[i],
        min_size=pacbio_min_size,
        truth_tool_name=tool_names[i],
        ploidy_table=ploidy_table,
        output_prefix="~{prefix}.prep_~{tool_names[i]}_vcf.~{sample_id}",
        script=preprocess_script,
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_attr_combine_truth
    }
    call conc.SVConcordanceTask {
      input:
        truth_vcf=PrepPacBioVcf.out,
        eval_vcf=sample_vcf,
        output_prefix="~{prefix}.concordance.~{tool_names[i]}.~{sample_id}",
        additional_args="--pesr-interval-overlap ~{pesr_interval_overlap} --pesr-size-similarity ~{pesr_size_similarity} --pesr-breakend-window ~{pesr_breakend_window}",
        reference_dict=reference_dict,
        java_mem_fraction=java_mem_fraction,
        gatk_docker=gatk_docker,
        runtime_attr_override=runtime_attr_sv_concordance
    }
  }

  call utils.TarFiles {
    input:
      files=flatten([SVConcordanceTask.out, SVConcordanceTask.out_index]),
      prefix="~{prefix}.pacbio_concordance_vcfs.~{sample_id}",
      linux_docker=linux_docker,
      runtime_attr_override=runtime_attr_tar
  }

  output {
    Array[File] pacbio_concordance_vcfs = SVConcordanceTask.out
    File pacbio_concordance_vcfs_tar = TarFiles.out
  }
}

task PrepPacBioVcf {
  input {
    String sample_id
    File vcf
    String truth_tool_name
    Int min_size = 25
    File ploidy_table
    File? script
    String output_prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 1.5,
                               disk_gb: ceil(10 + 3 * size(vcf, "GB")),
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File out = "~{output_prefix}.vcf.gz"
    File out_index = "~{output_prefix}.vcf.gz.tbi"
  }
  command <<<
    set -euxo pipefail

    # Convert format
    bcftools reheader -s <(echo "~{sample_id}") ~{vcf} > tmp1.vcf.gz
    python ~{default="/opt/sv-pipeline/scripts/format_pb_for_gatk.py" script} \
      --vcf tmp1.vcf.gz \
      --algorithm ~{truth_tool_name} \
      --min-size ~{min_size} \
      --out tmp2.vcf.gz \
      --ploidy-table ~{ploidy_table}
    bcftools sort tmp2.vcf.gz \
      | bcftools annotate --set-id '~{truth_tool_name}_%CHROM\_%POS\_%END\_%SVTYPE\_%SVLEN' -Oz -o ~{output_prefix}.vcf.gz
    tabix ~{output_prefix}.vcf.gz
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

