version 1.0

import "Structs.wdl"
import "TasksGenotypeBatch.wdl" as tasksgenotypebatch
import "TasksMakeCohortVcf.wdl" as tasksmakecohortvcf

workflow GenotypePESRPart2 {
  input {
    File bin_exclude
    File cohort_vcf
    File RD_pesr_sepcutoff
    File RD_depth_sepcutoff
    File PE_metrics
    File SR_metrics
    Int n_per_split
    Int n_RdTest_bins
    String batch
    File ref_dict

    File medianfile
    Array[String] samples

    File coveragefile
    File? coveragefile_index
    File discfile
    File? discfile_index
    File splitfile
    File? splitfile_index

    # SR genotyping parameters
    Int? sr_median_hom_ins
    Float? sr_hom_cutoff_multiplier

    String sv_pipeline_docker
    String sv_base_mini_docker
    String linux_docker
    String sv_pipeline_rdtest_docker
    RuntimeAttr? runtime_attr_split_variants
    RuntimeAttr? runtime_attr_make_subset_vcf
    RuntimeAttr? runtime_attr_count_pe
    RuntimeAttr? runtime_attr_genotype_pe
    RuntimeAttr? runtime_attr_count_sr
    RuntimeAttr? runtime_attr_genotype_sr
    RuntimeAttr? runtime_attr_rdtest_genotype
    RuntimeAttr? runtime_attr_integrate_gq
    RuntimeAttr? runtime_attr_integrate_pesr_gq
    RuntimeAttr? runtime_attr_add_genotypes
    RuntimeAttr? runtime_attr_triple_stream_cat
    RuntimeAttr? runtime_attr_concat_vcfs
  }

  File bin_exclude_idx = bin_exclude + ".tbi"

  call tasksgenotypebatch.SplitVariants as SplitVariants {
    input:
      vcf = cohort_vcf,
      n_per_split = n_per_split,
      generate_bca = true,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_split_variants
  }

  scatter (lt5kb_bed in SplitVariants.lt5kb_beds) {
    call tasksgenotypebatch.MakeSubsetVcf as MakeSubsetVcfUnder5kb {
      input:
        vcf = cohort_vcf,
        bed = lt5kb_bed,
        sv_base_mini_docker = sv_base_mini_docker,
        runtime_attr_override = runtime_attr_make_subset_vcf
    }

    call tasksgenotypebatch.CountPE as CountPEUnder5kb {
      input:
        vcf = MakeSubsetVcfUnder5kb.subset_vcf,
        discfile = discfile,
        discfile_index = discfile_index,
        medianfile = medianfile,
        samples = samples,
        ref_dict = ref_dict,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_count_pe
    }

    call GenotypePEPart2 as GenotypePEPart2Under5kb {
      input:
        PE_counts = CountPEUnder5kb.pe_counts,
        PE_metrics = PE_metrics,
        sv_pipeline_rdtest_docker = sv_pipeline_rdtest_docker,
        runtime_attr_override = runtime_attr_genotype_pe
    }

    call tasksgenotypebatch.CountSR as CountSRUnder5kb {
      input:
        vcf = MakeSubsetVcfUnder5kb.subset_vcf,
        splitfile = splitfile,
        splitfile_index = splitfile_index,
        medianfile = medianfile,
        samples = samples,
        ref_dict = ref_dict,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_count_sr
    }

    call GenotypeSRPart2 as GenotypeSRPart2Under5kb {
      input:
        vcf = MakeSubsetVcfUnder5kb.subset_vcf,
        SR_counts = CountSRUnder5kb.sr_counts,
        SR_sum = CountSRUnder5kb.sr_sum,
        SR_metrics = SR_metrics,
        svtype = "CNV_LT_5KBP",
        hom_cutoff_multiplier = sr_hom_cutoff_multiplier,
        sv_pipeline_rdtest_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_genotype_sr
    }

    call tasksgenotypebatch.RDTestGenotype as RDTestGenotypeUnder5kb {
      input:
        bin_exclude=bin_exclude,
        bin_exclude_idx=bin_exclude_idx,
        bed = lt5kb_bed,
        coveragefile = coveragefile,
        coveragefile_index = coveragefile_index,
        medianfile = medianfile,
        samples = samples,
        gt_cutoffs = RD_pesr_sepcutoff,
        n_bins = n_RdTest_bins,
        prefix = basename(lt5kb_bed, ".bed"),
        generate_melted_genotypes = true,
        ref_dict = ref_dict,
        sv_pipeline_rdtest_docker = sv_pipeline_rdtest_docker,
        runtime_attr_override = runtime_attr_rdtest_genotype
    }

    call IntegrateGQ as IntegrateGQUnder5kb {
      input:
        vcf = MakeSubsetVcfUnder5kb.subset_vcf,
        RD_melted_genotypes = RDTestGenotypeUnder5kb.melted_genotypes,
        RD_vargq = RDTestGenotypeUnder5kb.varGQ,
        PE_genotypes = GenotypePEPart2Under5kb.genotypes,
        PE_vargq = GenotypePEPart2Under5kb.varGQ,
        SR_genotypes = GenotypeSRPart2Under5kb.genotypes,
        SR_vargq = GenotypeSRPart2Under5kb.varGQ,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_integrate_gq
    }

    call tasksgenotypebatch.AddGenotypes as AddGenotypesUnder5kb {
      input:
        vcf = MakeSubsetVcfUnder5kb.subset_vcf,
        genotypes = IntegrateGQUnder5kb.genotypes,
        varGQ = IntegrateGQUnder5kb.varGQ,
        prefix = basename(lt5kb_bed, ".bed"),
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_add_genotypes
    }
  }

  scatter (gt5kb_bed in SplitVariants.gt5kb_beds) {
    call tasksgenotypebatch.MakeSubsetVcf as MakeSubsetVcfOver5kb {
      input:
        vcf = cohort_vcf,
        bed = gt5kb_bed,
        sv_base_mini_docker = sv_base_mini_docker,
        runtime_attr_override = runtime_attr_make_subset_vcf
    }

    call tasksgenotypebatch.CountPE as CountPEOver5kb {
      input:
        vcf = MakeSubsetVcfOver5kb.subset_vcf,
        discfile = discfile,
        discfile_index = discfile_index,
        medianfile = medianfile,
        samples = samples,
        ref_dict = ref_dict,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_count_pe
    }

    call GenotypePEPart2 as GenotypePEPart2Over5kb {
      input:
        PE_counts = CountPEOver5kb.pe_counts,
        PE_metrics = PE_metrics,
        sv_pipeline_rdtest_docker = sv_pipeline_rdtest_docker,
        runtime_attr_override = runtime_attr_genotype_pe
    }

    call tasksgenotypebatch.CountSR as CountSROver5kb {
      input:
        vcf = MakeSubsetVcfOver5kb.subset_vcf,
        splitfile = splitfile,
        splitfile_index = splitfile_index,
        medianfile = medianfile,
        samples = samples,
        ref_dict = ref_dict,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_count_sr
    }

    call GenotypeSRPart2 as GenotypeSRPart2Over5kb {
      input:
        vcf = MakeSubsetVcfOver5kb.subset_vcf,
        SR_counts = CountSROver5kb.sr_counts,
        SR_sum = CountSROver5kb.sr_sum,
        SR_metrics = SR_metrics,
        svtype = "CNV_GT_5KBP",
        hom_cutoff_multiplier = sr_hom_cutoff_multiplier,
        sv_pipeline_rdtest_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_genotype_sr
    }

    call tasksgenotypebatch.RDTestGenotype as RDTestGenotypeOver5kb {
      input:
        bin_exclude=bin_exclude,
        bin_exclude_idx=bin_exclude_idx,
        bed = gt5kb_bed,
        coveragefile = coveragefile,
        coveragefile_index = coveragefile_index,
        medianfile = medianfile,
        samples = samples,
        gt_cutoffs = RD_depth_sepcutoff,
        n_bins = n_RdTest_bins,
        prefix = basename(gt5kb_bed),
        generate_melted_genotypes = true,
        ref_dict = ref_dict,
        sv_pipeline_rdtest_docker = sv_pipeline_rdtest_docker,
        runtime_attr_override = runtime_attr_rdtest_genotype
    }

    call IntegrateGQ as IntegrateGQOver5kb {
      input:
        vcf = MakeSubsetVcfOver5kb.subset_vcf,
        RD_melted_genotypes = RDTestGenotypeOver5kb.melted_genotypes,
        RD_vargq = RDTestGenotypeOver5kb.varGQ,
        PE_genotypes = GenotypePEPart2Over5kb.genotypes,
        PE_vargq = GenotypePEPart2Over5kb.varGQ,
        SR_genotypes = GenotypeSRPart2Over5kb.genotypes,
        SR_vargq = GenotypeSRPart2Over5kb.varGQ,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_integrate_gq
    }

    call tasksgenotypebatch.AddGenotypes as AddGenotypesOver5kb {
      input:
        vcf = MakeSubsetVcfOver5kb.subset_vcf,
        genotypes = IntegrateGQOver5kb.genotypes,
        varGQ = IntegrateGQOver5kb.varGQ,
        prefix = basename(gt5kb_bed, ".bed"),
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_add_genotypes
    }
  }

  scatter (bca_bed in SplitVariants.bca_beds) {
    call tasksgenotypebatch.MakeSubsetVcf as MakeSubsetVcfBca {
      input:
        vcf = cohort_vcf,
        bed = bca_bed,
        sv_base_mini_docker = sv_base_mini_docker,
        runtime_attr_override = runtime_attr_make_subset_vcf
    }

    call tasksgenotypebatch.CountPE as CountPEBca {
      input:
        vcf = MakeSubsetVcfBca.subset_vcf,
        discfile = discfile,
        discfile_index = discfile_index,
        medianfile = medianfile,
        samples = samples,
        ref_dict = ref_dict,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_count_pe
    }

    call GenotypePEPart2 as GenotypePEPart2Bca {
      input:
        PE_counts = CountPEBca.pe_counts,
        PE_metrics = PE_metrics,
        sv_pipeline_rdtest_docker = sv_pipeline_rdtest_docker,
        runtime_attr_override = runtime_attr_genotype_pe
    }

    call tasksgenotypebatch.CountSR as CountSRBca {
      input:
        vcf = MakeSubsetVcfBca.subset_vcf,
        splitfile = splitfile,
        splitfile_index = splitfile_index,
        medianfile = medianfile,
        samples = samples,
        ref_dict = ref_dict,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_count_sr
    }

    call GenotypeSRPart2 as GenotypeSRPart2Bca {
      input:
        vcf = MakeSubsetVcfBca.subset_vcf,
        SR_counts = CountSRBca.sr_counts,
        SR_sum = CountSRBca.sr_sum,
        SR_metrics = SR_metrics,
        svtype = "BCA",
        hom_cutoff_multiplier = sr_hom_cutoff_multiplier,
        sv_pipeline_rdtest_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_genotype_sr
    }

    call IntegratePesrGQ as IntegratePesrGQBca {
      input:
        vcf = MakeSubsetVcfBca.subset_vcf,
        PE_genotypes = GenotypePEPart2Bca.genotypes,
        PE_vargq = GenotypePEPart2Bca.varGQ,
        SR_genotypes = GenotypeSRPart2Bca.genotypes,
        SR_vargq = GenotypeSRPart2Bca.varGQ,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_integrate_pesr_gq
    }

    call tasksgenotypebatch.AddGenotypes as AddGenotypesBca {
      input:
        vcf = MakeSubsetVcfBca.subset_vcf,
        genotypes = IntegratePesrGQBca.genotypes,
        varGQ = IntegratePesrGQBca.varGQ,
        prefix = basename(bca_bed, ".bed"),
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_add_genotypes
    }
  }

  scatter (ins_bed in SplitVariants.ins_beds) {
    call tasksgenotypebatch.MakeSubsetVcf as MakeSubsetVcfIns {
      input:
        vcf = cohort_vcf,
        bed = ins_bed,
        sv_base_mini_docker = sv_base_mini_docker,
        runtime_attr_override = runtime_attr_make_subset_vcf
    }

    call tasksgenotypebatch.CountPE as CountPEIns {
      input:
        vcf = MakeSubsetVcfIns.subset_vcf,
        discfile = discfile,
        discfile_index = discfile_index,
        medianfile = medianfile,
        samples = samples,
        ref_dict = ref_dict,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_count_pe
    }

    call GenotypePEPart2 as GenotypePEPart2Ins {
      input:
        PE_counts = CountPEIns.pe_counts,
        PE_metrics = PE_metrics,
        sv_pipeline_rdtest_docker = sv_pipeline_rdtest_docker,
        runtime_attr_override = runtime_attr_genotype_pe
    }

    call tasksgenotypebatch.CountSR as CountSRIns {
      input:
        vcf = MakeSubsetVcfIns.subset_vcf,
        splitfile = splitfile,
        splitfile_index = splitfile_index,
        medianfile = medianfile,
        samples = samples,
        ref_dict = ref_dict,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_count_sr
    }

    call GenotypeSRPart2 as GenotypeSRPart2Ins {
      input:
        vcf = MakeSubsetVcfIns.subset_vcf,
        SR_counts = CountSRIns.sr_counts,
        SR_sum = CountSRIns.sr_sum,
        SR_metrics = SR_metrics,
        svtype = "INS",
        hom_cutoff_multiplier = sr_hom_cutoff_multiplier,
        median_hom_ins = sr_median_hom_ins,
        sv_pipeline_rdtest_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_genotype_sr
    }

    call IntegratePesrGQ as IntegratePesrGQIns {
      input:
        vcf = MakeSubsetVcfIns.subset_vcf,
        PE_genotypes = GenotypePEPart2Ins.genotypes,
        PE_vargq = GenotypePEPart2Ins.varGQ,
        SR_genotypes = GenotypeSRPart2Ins.genotypes,
        SR_vargq = GenotypeSRPart2Ins.varGQ,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_integrate_pesr_gq
    }

    call tasksgenotypebatch.AddGenotypes as AddGenotypesIns {
      input:
        vcf = MakeSubsetVcfIns.subset_vcf,
        genotypes = IntegratePesrGQIns.genotypes,
        varGQ = IntegratePesrGQIns.varGQ,
        prefix = basename(ins_bed, ".bed"),
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_add_genotypes
    }
  }

  call CatFiles as CatFilesFail {
    input:
      files = flatten([GenotypeSRPart2Under5kb.background_fail, GenotypeSRPart2Over5kb.background_fail, GenotypeSRPart2Bca.background_fail, GenotypeSRPart2Ins.background_fail]),
      outfile = "~{batch}.genotype_SR_part2_background_fail.txt",
      linux_docker = linux_docker,
      runtime_attr_override = runtime_attr_triple_stream_cat
  }

  call CatFiles as CatFilesPass {
    input:
      files = flatten([GenotypeSRPart2Under5kb.bothside_pass, GenotypeSRPart2Over5kb.bothside_pass, GenotypeSRPart2Bca.bothside_pass, GenotypeSRPart2Ins.bothside_pass]),
      outfile = "~{batch}.genotype_SR_part2_bothside_pass.txt",
      linux_docker = linux_docker,
      runtime_attr_override = runtime_attr_triple_stream_cat
  }

  call tasksmakecohortvcf.ConcatVcfs {
    input:
      vcfs=flatten([AddGenotypesUnder5kb.genotyped_vcf, AddGenotypesOver5kb.genotyped_vcf, AddGenotypesBca.genotyped_vcf, AddGenotypesIns.genotyped_vcf]),
      vcfs_idx=flatten([AddGenotypesUnder5kb.genotyped_vcf_index, AddGenotypesOver5kb.genotyped_vcf_index, AddGenotypesBca.genotyped_vcf_index, AddGenotypesIns.genotyped_vcf_index]),
      allow_overlaps=true,
      outfile_prefix="~{batch}.genotyped_pesr",
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_attr_concat_vcfs
  }

  output {
    File bothside_pass = CatFilesPass.merged_file
    File background_fail = CatFilesFail.merged_file
    File genotyped_vcf = ConcatVcfs.concat_vcf
    File genotyped_vcf_index = ConcatVcfs.concat_vcf_idx
  }
}

task IntegrateGQ {
  input {
    File vcf
    File RD_melted_genotypes
    File RD_vargq
    File PE_genotypes
    File PE_vargq
    File SR_genotypes
    File SR_vargq
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
    File genotypes = "genotype.indiv.txt.gz"
    File varGQ = "genotype.variant.txt.gz"
  }
  command <<<

    /opt/sv-pipeline/04_variant_resolution/scripts/IntegrateGQ.sh \
      ~{vcf} \
      ~{RD_melted_genotypes} \
      ~{RD_vargq} \
      ~{PE_genotypes} \
      ~{PE_vargq} \
      ~{SR_genotypes} \
      ~{SR_vargq}
  
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

task IntegratePesrGQ {
  input {
    File vcf
    File PE_genotypes
    File PE_vargq
    File SR_genotypes
    File SR_vargq
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
    File genotypes = "genotype.indiv.txt.gz"
    File varGQ = "genotype.variant.txt.gz"
  }
  command <<<

    /opt/sv-pipeline/04_variant_resolution/scripts/IntegrateGQ_PESR.sh \
      ~{vcf} \
      ~{PE_genotypes} \
      ~{PE_vargq} \
      ~{SR_genotypes} \
      ~{SR_vargq};
  
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

task CatFiles {
  input {
    Array[File] files
    String outfile
    String linux_docker
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
    File merged_file = "${outfile}"
  }
  command <<<

    set -euo pipefail
    cat ~{sep=" " files} \
      | sort -Vk1,1 \
      | uniq > ~{outfile}
  
  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: linux_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

task GenotypePEPart2 {
  input {
    File PE_counts
    File PE_metrics
    String sv_pipeline_rdtest_docker
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
    File genotypes = "pe.geno.withquality.txt.gz"
    File varGQ = "pe.variant.quality.final.txt.gz"
  }
  command <<<

    /opt/sv-pipeline/04_variant_resolution/scripts/PE_genotype.opt_part2.sh \
      ~{PE_counts} \
      ~{PE_metrics}
  
  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_rdtest_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

task GenotypeSRPart2 {
  input {
    File vcf
    File SR_counts
    File SR_sum
    File SR_metrics
    File? script
    Float hom_cutoff_multiplier = 1.6
    Int median_hom_ins = 105
    String svtype
    String sv_pipeline_rdtest_docker
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
    File genotypes = "sr.geno.withquality.txt.gz"
    File varGQ = "sr.variant.quality.final.txt.gz"
    File background_fail = "background.variant.fail.txt"
    File bothside_pass = "both.pass.txt"
  }
  command <<<

    bash ~{default="/opt/sv-pipeline/04_variant_resolution/scripts/SR_genotype.opt_part2.sh" script} \
      ~{vcf} \
      ~{SR_counts} \
      ~{SR_sum} \
      ~{SR_metrics} \
      ~{svtype} \
      ~{hom_cutoff_multiplier} \
      ~{median_hom_ins}
  
  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_rdtest_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}
