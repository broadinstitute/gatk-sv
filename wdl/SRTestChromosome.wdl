version 1.0

import "TasksGenerateBatchMetrics.wdl" as tasksbatchmetrics

workflow SRTestChromosome {
  input {
    File vcf
    String chrom
    String batch
    Int split_size
    String algorithm
    File medianfile
    File splitfile
    File ped_file
    Int? suffix_len
    File male_samples
    File female_samples
    File male_only_variant_ids
    File samples
    File ref_dict
    Boolean allosome
    Boolean run_common
    Int? common_cnv_size_cutoff
    File? outlier_sample_ids

    String sv_pipeline_docker
    String linux_docker
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_split_vcf
    RuntimeAttr? runtime_attr_srtest
    RuntimeAttr? runtime_attr_merge_allo
    RuntimeAttr? runtime_attr_merge_stats
  }
  
  File splitfile_idx = splitfile + ".tbi"

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
      call SRTest as SRTestFemale {
        input:
          vcf = split,
          splitfile = splitfile,
          medianfile = medianfile,
          splitfile_idx = splitfile_idx,
          include_list = female_samples,
          common_model = false,
          prefix = basename(split),
          ref_dict = ref_dict,
          outlier_sample_ids = outlier_sample_ids,
          sv_pipeline_docker = sv_pipeline_docker,
          runtime_attr_override = runtime_attr_srtest
      }

      call SRTest as SRTestMale {
        input:
          vcf = split,
          splitfile = splitfile,
          medianfile = medianfile,
          splitfile_idx = splitfile_idx,
          include_list = male_samples,
          common_model = false,
          prefix = basename(split),
          ref_dict = ref_dict,
          outlier_sample_ids = outlier_sample_ids,
          sv_pipeline_docker = sv_pipeline_docker,
          runtime_attr_override = runtime_attr_srtest
      }

      call tasksbatchmetrics.MergeAllosomes as MergeAllosomes {
        input:
          male_test = SRTestMale.stats,
          female_test = SRTestFemale.stats,
          male_only_ids_list = male_only_variant_ids,
          chrom = chrom,
          runtime_attr_override = runtime_attr_merge_allo,
          sv_pipeline_docker = sv_pipeline_docker,
          male_only_expr = "females.name.isin(male_only_ids)"
      }
    } 

    if (!allosome) {
      call SRTest as SRTestAutosome {
        input:
          vcf = split,
          splitfile = splitfile,
          medianfile = medianfile,
          splitfile_idx = splitfile_idx,
          include_list = samples,
          common_model = false,
          prefix = basename(split),
          ref_dict = ref_dict,
          outlier_sample_ids = outlier_sample_ids,
          sv_pipeline_docker = sv_pipeline_docker,
          runtime_attr_override = runtime_attr_srtest
      }
    }
  }

  if (run_common && !allosome) {
    call tasksbatchmetrics.GetCommonVCF {
      input:
        vcf = vcf,
        cnv_size_cutoff = select_first([common_cnv_size_cutoff]),
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
      call SRTest as SRTestAutosomeCommon {
        input:
          vcf = split,
          splitfile = splitfile,
          medianfile = medianfile,
          splitfile_idx = splitfile_idx,
          include_list = samples,
          common_model = true,
          prefix = basename(split),
          ref_dict = ref_dict,
          outlier_sample_ids = outlier_sample_ids,
          sv_pipeline_docker = sv_pipeline_docker,
          runtime_attr_override = runtime_attr_srtest
      }
    }
  }
  
  Array[File] unmerged_stats = if allosome then select_all(MergeAllosomes.merged_test) else select_all(SRTestAutosome.stats)
  Array[File] unmerged_stats_common = select_first([SRTestAutosomeCommon.stats, []])

  call tasksbatchmetrics.MergeStats as MergeStats {
    input:
      stats = unmerged_stats,
      prefix = "${batch}.${algorithm}.${chrom}",
      linux_docker = linux_docker,
      runtime_attr_override = runtime_attr_merge_stats
  }

  if (run_common && !allosome) {
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

task SRTest {
  input {
    File vcf
    File splitfile
    File medianfile
    File splitfile_idx
    File include_list
    Boolean common_model
    String prefix
    File ref_dict
    File? outlier_sample_ids
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  String common_arg = if common_model then "--common" else ""

  parameter_meta {
    splitfile: {
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
    awk -v OFS="\t" '{if ($2-250>0){print $1,$2-250,$2+250}else{print $1,0,$2+250}}' test.bed  >> region.bed
    awk -v OFS="\t" '{if ($3-250>0){print $1,$3-250,$3+250}else{print $1,0,$3+250}}' test.bed  >> region.bed
    sort -k1,1 -k2,2n region.bed > region.sorted.bed
    bedtools merge -d 16384 -i region.sorted.bed > region.merged.bed

    if [ -s region.merged.bed ]; then
      java -Xmx~{java_mem_mb}M -jar ${GATK_JAR} PrintSVEvidence \
        --sequence-dictionary ~{ref_dict} \
        --evidence-file ~{splitfile} \
        -L region.merged.bed \
        -O local.SR.txt.gz
    else
      touch local.SR.txt
      bgzip local.SR.txt
      tabix -0 -s1 -b2 -e2 local.SR.txt.gz
    fi
    
    svtk sr-test -w 50 --log ~{common_arg} --medianfile ~{medianfile} --samples ~{include_list} ~{vcf} local.SR.txt.gz ~{prefix}.stats ~{if defined(outlier_sample_ids) then "--outlier-sample-ids ~{outlier_sample_ids}" else ""}
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


