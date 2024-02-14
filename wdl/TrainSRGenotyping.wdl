version 1.0

import "Structs.wdl"
import "TasksGenotypeBatch.wdl" as tasksgenotypebatch

workflow TrainSRGenotyping {
  input {
    File batch_vcf
    File splitfile
    File? splitfile_index
    Int n_per_split
    File medianfile
    Array[String] samples
    String batch_ID
    File RF_cutoffs
    File RD_melted_genotypes
    File PE_train
    File PE_genotypes
    File ref_dict

    String sv_base_mini_docker
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_split_vcf
    RuntimeAttr? runtime_attr_count_sr
    RuntimeAttr? runtime_attr_merge_counts
    RuntimeAttr? runtime_attr_genotype
  }

  call tasksgenotypebatch.SplitVcf as SplitVcf {
    input:
      vcf = batch_vcf,
      n_per_split = n_per_split,
      evidence_type = "sr",
      bgzip = true,
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_split_vcf
  }

  scatter (vcf in SplitVcf.vcfs) {
    call tasksgenotypebatch.CountSR as CountSR {
      input:
        vcf = vcf,
        splitfile = splitfile,
        splitfile_index = splitfile_index,
        medianfile = medianfile,
        samples = samples,
        ref_dict = ref_dict,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_count_sr
    }
  }

  call tasksgenotypebatch.MergePESRCounts as MergeSRCounts {
    input:
      count_list = CountSR.sr_counts,
      sum_list = CountSR.sr_sum,
      evidence_type = "sr",
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_merge_counts
  }

  call GenotypeSRPart1 {
    input:
      vcf = batch_vcf,
      SR_counts = MergeSRCounts.counts,
      SR_sum = MergeSRCounts.sum,
      RD_melted_genotypes = RD_melted_genotypes,
      RF_cutoffs = RF_cutoffs,
      samples = samples,
      PE_train = PE_train,
      PE_genotypes = PE_genotypes,
      batch = batch_ID,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_genotype
  }

  output {
    File SR_metrics = GenotypeSRPart1.SR_metrics
  }
}

task GenotypeSRPart1 {
  input {
    File vcf
    File SR_counts
    File SR_sum
    File RD_melted_genotypes
    File RF_cutoffs
    Array[String] samples
    File PE_train
    File PE_genotypes
    String batch
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 6,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File SR_metrics = "~{batch}.sr_metric_file.txt"
  }
  command <<<

    /opt/sv-pipeline/04_variant_resolution/scripts/SR_genotype.opt_part1.sh \
      ~{vcf} \
      ~{SR_counts} \
      ~{SR_sum} \
      ~{RD_melted_genotypes} \
      ~{RF_cutoffs} \
      ~{write_lines(samples)} \
      ~{PE_train} \
      ~{PE_genotypes} \
      ~{batch}
  
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
