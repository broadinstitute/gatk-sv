version 1.0

import "TasksGenerateBatchMetrics.wdl" as tasksbatchmetrics

workflow PETestChromosome {
  input {
    String chrom
    String algorithm
    File vcf
    File discfile
    File medianfile
    File ref_dict
    String batch
    Int split_size
    File ped_file
    Int? suffix_len
    File male_samples
    File female_samples
    File male_only_variant_ids
    File samples
    Boolean allosome
    Int common_cnv_size_cutoff
    File? outlier_sample_ids

    String sv_base_mini_docker
    String linux_docker
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_split_vcf
    RuntimeAttr? runtime_attr_petest
    RuntimeAttr? runtime_attr_merge_allo
    RuntimeAttr? runtime_attr_merge_stats
  }

  File discfile_idx = discfile + ".tbi"

  call tasksbatchmetrics.SplitVCF as SplitVCF {
    input:
      vcf = vcf,
      batch = batch,
      algorithm = algorithm,
      chrom = chrom,
      split_size = split_size,
      suffix_len = select_first([suffix_len, 4]),
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_split_vcf
  }

  scatter (split in SplitVCF.split_vcfs) {
    if (allosome) {
      call PETest as PETestMale {
        input:
          vcf = split,
          discfile = discfile,
          medianfile = medianfile,
          discfile_idx = discfile_idx,
          include_list = male_samples,
          prefix = basename(split),
          ref_dict = ref_dict,
          common_model = false,
          outlier_sample_ids = outlier_sample_ids,
          sv_pipeline_docker = sv_pipeline_docker,
          runtime_attr_override = runtime_attr_petest
      }

      call PETest as PETestFemale {
        input:
          vcf = split,
          discfile = discfile,
          medianfile = medianfile,
          discfile_idx = discfile_idx,
          include_list = female_samples,
          prefix = basename(split),
          ref_dict = ref_dict,
          common_model = false,
          outlier_sample_ids = outlier_sample_ids,
          sv_pipeline_docker = sv_pipeline_docker,
          runtime_attr_override = runtime_attr_petest
      }

      call tasksbatchmetrics.MergeAllosomes as MergeAllosomes {
        input:
          male_test = PETestMale.stats,
          female_test = PETestFemale.stats,
          male_only_ids_list = male_only_variant_ids,
          chrom = chrom,
          sv_pipeline_docker = sv_pipeline_docker,
          runtime_attr_override = runtime_attr_merge_allo,
          male_only_expr = "females.name.isin(male_only_ids)"
      }
    }
    if (!allosome) {
      call PETest as PETestAutosome {
        input:
          vcf = split,
          discfile = discfile,
          medianfile = medianfile,
          discfile_idx = discfile_idx,
          include_list = samples,
          prefix = basename(split),
          ref_dict = ref_dict,
          common_model = false,
          outlier_sample_ids = outlier_sample_ids,
          sv_pipeline_docker = sv_pipeline_docker,
          runtime_attr_override = runtime_attr_petest
      }
    }
  }

  if (!allosome) {
    call tasksbatchmetrics.GetCommonVCF {
      input:
        vcf = vcf,
        cnv_size_cutoff = common_cnv_size_cutoff,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_split_vcf
    }

    call tasksbatchmetrics.SplitVCF as SplitCommonVCF {
      input:
        vcf = GetCommonVCF.common_vcf,
        batch = batch,
        algorithm = algorithm,
        chrom = chrom,
        split_size = split_size,
        suffix_len = select_first([suffix_len, 4]),
        sv_base_mini_docker = sv_base_mini_docker,
        runtime_attr_override = runtime_attr_split_vcf
    }

    scatter (split in SplitCommonVCF.split_vcfs) {
      call PETest as PETestAutosomeCommon {
        input:
          vcf = split,
          discfile = discfile,
          medianfile = medianfile,
          discfile_idx = discfile_idx,
          include_list = samples,
          ref_dict = ref_dict,
          prefix = basename(split),
          common_model = true,
          outlier_sample_ids = outlier_sample_ids,
          sv_pipeline_docker = sv_pipeline_docker,
          runtime_attr_override = runtime_attr_petest
      }
    }
  }

  Array[File] unmerged_stats = if allosome then select_all(MergeAllosomes.merged_test) else select_all(PETestAutosome.stats)
  Array[File] unmerged_stats_common = select_first([PETestAutosomeCommon.stats, []])

  call tasksbatchmetrics.MergeStats as MergeStats {
    input:
      stats = unmerged_stats,
      prefix = "${batch}.${algorithm}.${chrom}",
      linux_docker = linux_docker,
      runtime_attr_override = runtime_attr_merge_stats
  }

  if (!allosome) {
    call tasksbatchmetrics.MergeStats as MergeStatsCommon {
      input:
        stats = unmerged_stats_common,
        prefix = "${batch}.${algorithm}.${chrom}.common",
        linux_docker = linux_docker,
        runtime_attr_override = runtime_attr_merge_stats
    }
  }

  output {
    File stats = MergeStats.merged_stats
    File? stats_common = MergeStatsCommon.merged_stats
  }
}

task PETest {
  input {
    File vcf
    File discfile
    File medianfile
    File discfile_idx
    File include_list
    Boolean common_model
    File ref_dict
    String prefix
    File? outlier_sample_ids
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Int window = 1000
  String common_arg = if common_model then "--common" else ""

  parameter_meta {
    discfile: {
      localization_optional: true
    }
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

  Float mem_gb = select_first([runtime_attr.mem_gb, default_attr.mem_gb])
  Int java_mem_mb = ceil(mem_gb * 1000 * 0.8)

  output {
    File stats = "${prefix}.stats"
  }
  command <<<

    set -euo pipefail
    svtk vcf2bed --split-bnd --no-header ~{vcf} test.bed
    awk -v OFS="\t" '{if ($2-~{window}>0){print $1,$2-~{window},$2+~{window}}else{print $1,0,$2+~{window}}}' test.bed  >> region.bed
    awk -v OFS="\t" '{if ($3-~{window}>0){print $1,$3-~{window},$3+~{window}}else{print $1,0,$3+~{window}}}' test.bed  >> region.bed
    sort -k1,1 -k2,2n region.bed > region.sorted.bed
    bedtools merge -d 16384 -i region.sorted.bed > region.merged.bed

    if [ -s region.merged.bed ]; then
      java -Xmx~{java_mem_mb}M -jar ${GATK_JAR} PrintSVEvidence \
        --sequence-dictionary ~{ref_dict} \
        --evidence-file ~{discfile} \
        -L region.merged.bed \
        -O local.PE.txt.gz

      tabix -f -0 -s1 -b2 -e2 local.PE.txt.gz
    else
      touch local.PE.txt
      bgzip local.PE.txt
      tabix -0 -s1 -b2 -e2 local.PE.txt.gz
    fi

    svtk pe-test -o ~{window} ~{common_arg} --medianfile ~{medianfile} --samples ~{include_list} ~{vcf} local.PE.txt.gz ~{prefix}.stats ~{if defined(outlier_sample_ids) then "--outlier-sample-ids ~{outlier_sample_ids}" else ""}
  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: mem_gb + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

