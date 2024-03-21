version 1.0

import "Module07XfBatchEffect.wdl" as batch_effect

workflow BatchEffectAcrossContigs {
  input {
    Array[File] vcfs
    File sample_batch_assignments
    File batches_list
    File excludesamples_list #empty file if need be
    File sample_pop_assignments
    File famfile
    File contiglist
    File? par_bed
    Int? onevsall_cutoff=2
    String prefix

    String sv_pipeline_docker
    String sv_base_mini_docker
  }

  Array[String] contigs = read_lines(contiglist)

  Array[String] batches = read_lines(batches_list)

  scatter ( vcf in vcfs)  {
    call IndexVcf {
      input:
        vcf=vcf,
        sv_base_mini_docker=sv_base_mini_docker
    }
  }

  scatter ( batch in batches ) {
    # Get list of samples to include & exclude per batch
    call batch_effect.GetBatchSamplesList {
      input:
        vcf=IndexVcf.vcf_out[0],
        vcf_idx=IndexVcf.vcf_out_index[0],
        batch=batch,
        sample_batch_assignments=sample_batch_assignments,
        probands_list=excludesamples_list,
        sv_pipeline_docker=sv_pipeline_docker
    }
  }

  scatter ( i in range(length(vcfs)) ) {

    call batch_effect.XfBatchEffect {
      input:
        vcf=IndexVcf.vcf_out[i],
        vcf_idx=IndexVcf.vcf_out_index[i],
        sample_batch_assignments=sample_batch_assignments,
        batches_list=batches_list,
        batch_include_lists = GetBatchSamplesList.include_samples_list,
        sample_pop_assignments=sample_pop_assignments,
        famfile=famfile,
        par_bed=par_bed,
        onevsall_cutoff=onevsall_cutoff,
        prefix="~{prefix}.~{contigs[i]}",
        sv_base_mini_docker=sv_base_mini_docker,
        sv_pipeline_docker=sv_pipeline_docker
    }
  }

  output {
    Array[File] batch_effects_labeled_vcfs = XfBatchEffect.labeled_vcf
    Array[File] batch_effects_labeled_vcf_indexes = XfBatchEffect.labeled_vcf_idx
    Array[File] batch_effect_reclassification_tables = XfBatchEffect.reclassification_table
  }
}

task IndexVcf {
  input {
    File vcf
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  String outfile_name = basename(vcf)

  # Disk must be scaled proportionally to the size of the VCF
  Float input_size = size(vcf, "GiB")
  RuntimeAttr default_attr = object {
    mem_gb: 3.75,
    disk_gb: ceil(10.0 + (input_size * 2)),
    cpu_cores: 1,
    preemptible_tries: 3,
    max_retries: 1,
    boot_disk_gb: 10
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<

    set -euo pipefail

    mv ~{vcf} .
    tabix -p vcf ~{outfile_name}

  >>>

  output {
    File vcf_out = outfile_name
    File vcf_out_index = "~{outfile_name}.tbi"
  }

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_base_mini_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}
