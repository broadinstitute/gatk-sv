##########################################################################################

## Base script:   https://portal.firecloud.org/#methods/Talkowski-SV/04_v2_genotype_depth_part2/14/wdl

## Github commit: talkowski-lab/gatk-sv-v1:<ENTER HASH HERE IN FIRECLOUD>

##########################################################################################

version 1.0

import "Structs.wdl"
import "Tasks04.wdl" as tasks04

workflow GenotypeDepthPart2 {
  input {
    File bin_exclude
    File cohort_vcf
    File RD_pesr_sepcutoff
    File RD_depth_sepcutoff
    Int n_per_split
    Int n_RdTest_bins
    String batch
    File ref_dict
    File medianfile
    File famfile
    Array[String] samples

    File coveragefile

    String sv_pipeline_docker
    String sv_base_mini_docker
    String sv_pipeline_rdtest_docker
    RuntimeAttr? runtime_attr_split_variants
    RuntimeAttr? runtime_attr_rdtest_genotype
    RuntimeAttr? runtime_attr_make_subset_vcf
    RuntimeAttr? runtime_attr_integrate_depth_gq
    RuntimeAttr? runtime_attr_add_genotypes
    RuntimeAttr? runtime_attr_concat_vcfs
  }

  File bin_exclude_idx = bin_exclude + ".tbi"

  call tasks04.SplitVariants as SplitVariants {
    input:
      vcf = cohort_vcf,
      n_per_split = n_per_split,
      generate_bca = false,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_split_variants
  }

  scatter (gt5kb_bed in SplitVariants.gt5kb_beds) {

    call tasks04.MakeSubsetVcf as MakeSubsetVcfOver5kb {
      input:
        vcf = cohort_vcf,
        bed = gt5kb_bed,
        sv_base_mini_docker = sv_base_mini_docker,
        runtime_attr_override = runtime_attr_make_subset_vcf
    }

    call tasks04.RDTestGenotype as RDTestGenotypeOver5kb {
      input:
        bin_exclude=bin_exclude,
        bin_exclude_idx=bin_exclude_idx,
        bed = gt5kb_bed,
        coveragefile = coveragefile,
        medianfile = medianfile,
        famfile = famfile,
        samples = samples,
        gt_cutoffs = RD_depth_sepcutoff,
        n_bins = n_RdTest_bins,
        prefix = basename(gt5kb_bed),
        generate_melted_genotypes = true,
        ref_dict = ref_dict,
        sv_pipeline_rdtest_docker = sv_pipeline_rdtest_docker,
        runtime_attr_override = runtime_attr_rdtest_genotype
    }

    call tasks04.IntegrateDepthGq as IntegrateDepthGqOver5kb {
      input:
        vcf = MakeSubsetVcfOver5kb.subset_vcf,
        RD_melted_genotypes = RDTestGenotypeOver5kb.melted_genotypes,
        RD_vargq = RDTestGenotypeOver5kb.varGQ,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_integrate_depth_gq
    }

    call tasks04.AddGenotypes as AddGenotypesOver5kb {
      input:
        vcf = MakeSubsetVcfOver5kb.subset_vcf,
        genotypes = IntegrateDepthGqOver5kb.genotypes,
        varGQ = IntegrateDepthGqOver5kb.varGQ,
        prefix = basename(gt5kb_bed, ".bed"),
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_add_genotypes
    }
  }

  scatter (lt5kb_bed in SplitVariants.lt5kb_beds) {

    call tasks04.RDTestGenotype as RDTestGenotypeUnder5kb {
      input:
        bin_exclude=bin_exclude,
        bin_exclude_idx=bin_exclude_idx,
        bed = lt5kb_bed,
        coveragefile = coveragefile,
        medianfile = medianfile,
        famfile = famfile,
        samples = samples,
        gt_cutoffs = RD_pesr_sepcutoff,
        n_bins = n_RdTest_bins,
        prefix = basename(lt5kb_bed, ".bed"),
        generate_melted_genotypes = true,
        ref_dict = ref_dict,
        sv_pipeline_rdtest_docker = sv_pipeline_rdtest_docker,
        runtime_attr_override = runtime_attr_rdtest_genotype
    }

    call tasks04.MakeSubsetVcf as MakeSubsetVcfUnder5kb {
      input:
        vcf = cohort_vcf,
        bed = lt5kb_bed,
        sv_base_mini_docker = sv_base_mini_docker,
        runtime_attr_override = runtime_attr_make_subset_vcf
    }

    call tasks04.IntegrateDepthGq as IntegrateDepthGqUnder5kb {
      input:
        vcf = MakeSubsetVcfUnder5kb.subset_vcf,
        RD_melted_genotypes = RDTestGenotypeUnder5kb.melted_genotypes,
        RD_vargq = RDTestGenotypeUnder5kb.varGQ,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_integrate_depth_gq
    }
    call tasks04.AddGenotypes as AddGenotypesUnder5kb {
      input:
        vcf = MakeSubsetVcfUnder5kb.subset_vcf,
        genotypes = IntegrateDepthGqUnder5kb.genotypes,
        varGQ = IntegrateDepthGqUnder5kb.varGQ,
        prefix = basename(lt5kb_bed, ".bed"),
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_add_genotypes
    }
  }

  call tasks04.ConcatGenotypedVcfs as ConcatGenotypedVcfs {
    input:
      lt5kb_vcfs = AddGenotypesUnder5kb.genotyped_vcf,
      gt5kb_vcfs = AddGenotypesOver5kb.genotyped_vcf,
      bca_vcfs = [],
      batch = batch,
      evidence_type = "depth",
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_concat_vcfs
  }
  output {
    File genotyped_vcf = ConcatGenotypedVcfs.genotyped_vcf
    File genotyped_vcf_index = ConcatGenotypedVcfs.genotyped_vcf_index
    Array[File] regeno_coverage_medians = RDTestGenotypeOver5kb.copy_states
  }
}

task IntegrateDepthGq {
  input {
    File vcf
    File RD_melted_genotypes
    File RD_vargq
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
    File genotypes = "genotype.indiv.depth.txt.gz"
    File varGQ = "genotype.variant.depth.txt.gz"
  }
  command <<<

    /opt/sv-pipeline/04_variant_resolution/scripts/IntegrateGQ_depthonly.sh \
      ~{vcf} \
      ~{RD_melted_genotypes} \
      ~{RD_vargq}
  
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
