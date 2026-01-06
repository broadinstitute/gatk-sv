version 1.0

import "Structs.wdl"
import "TasksGenotypeBatch.wdl" as tasksgenotypebatch

workflow TrainPEGenotyping {
  input {
    File batch_vcf    # variants from just the batch in question
    Array[String] samples
    File discfile
    File? discfile_index
    File medianfile
    Int n_per_split
    String batch_ID
    File RF_cutoffs
    File RD_genotypes
    File RD_melted_genotypes
    File exclude_list
    File ref_dict

    String sv_base_mini_docker
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_split_vcf
    RuntimeAttr? runtime_attr_make_batch_bed
    RuntimeAttr? runtime_attr_count_pe
    RuntimeAttr? runtime_attr_merge_counts
    RuntimeAttr? runtime_attr_genotype
  }

  call tasksgenotypebatch.SplitVcf as SplitVcf {
    input:
      vcf = batch_vcf,
      n_per_split = n_per_split,
      evidence_type = "pe",
      bgzip = false,
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_split_vcf
  }

  call VcfToBed {
    input:
      vcf = batch_vcf,
      prefix = batch_ID,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_make_batch_bed
  }

  scatter (vcf in SplitVcf.vcfs) {
    call tasksgenotypebatch.CountPE as CountPE {
      input:
        vcf = vcf,
        discfile = discfile,
        discfile_index = discfile_index,
        medianfile = medianfile,
        samples = samples,
        ref_dict = ref_dict,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_count_pe
    }
  }

  call tasksgenotypebatch.MergePESRCounts as MergePECounts {
    input:
      count_list = CountPE.pe_counts,
      sum_list = [],
      evidence_type = "pe",
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_merge_counts
  }

  call GenotypePEPart1 {
    input:
      bed = VcfToBed.bed,
      RF_cutoffs = RF_cutoffs,
      PE_counts = MergePECounts.counts,
      RD_genotypes = RD_genotypes,
      RD_melted_genotypes = RD_melted_genotypes,
      exclude_list = exclude_list,
      batch = batch_ID,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_genotype
  }

  output {
    File PE_varGQ = GenotypePEPart1.varGQ
    File PE_train = GenotypePEPart1.PE_train
    File PE_genotypes = GenotypePEPart1.genotypes
    File PE_metrics = GenotypePEPart1.PE_metrics
  }
}

task VcfToBed {
  input {
    File vcf
    String prefix
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

  output {
    File bed = "${prefix}.bed"
  }
  command <<<

    svtk vcf2bed ~{vcf} -i ALGORITHMS ~{prefix}.bed
  
  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    noAddress: true
  }
}

task GenotypePEPart1 {
  input {
    File bed
    File RF_cutoffs
    File PE_counts
    File RD_genotypes
    File RD_melted_genotypes
    File exclude_list
    String batch
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

  output {
    File PE_train = "~{batch}.pe.train.include.txt"
    File PE_metrics = "~{batch}.pe_metric_file.txt"
    File genotypes = "~{batch}.pe.geno.withquality.txt.gz"
    File varGQ = "~{batch}.pe.variant.quality.final.txt.gz"
  }
  command <<<

    /opt/sv-pipeline/04_variant_resolution/scripts/PE_genotype.sh \
      ~{bed} \
      ~{PE_counts} \
      ~{RD_genotypes} \
      ~{RD_melted_genotypes} \
      ~{RF_cutoffs} \
      ~{exclude_list} \
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
    noAddress: true
  }
}
