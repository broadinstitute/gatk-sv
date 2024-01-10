version 1.0

import "Structs.wdl"
import "TasksMakeCohortVcf.wdl" as MiniTasks
import "FormatVcfForGatk.wdl" as fvcf
import "CleanVcf1b.wdl" as c1b
import "CleanVcf5.wdl" as c5
import "HailMerge.wdl" as HailMerge

workflow CleanVcfChromosome {
  input {
    File vcf
    String contig
    File background_list
    File ped_file
    File allosome_fai
    String prefix
    Int max_shards_per_chrom_step1
    File bothsides_pass_list
    Int min_records_per_shard_step1
    Int samples_per_step2_shard
    Int clean_vcf1b_records_per_shard
    Int clean_vcf5_records_per_shard
    Int? clean_vcf5_threads_per_task
    File? outlier_samples_list
    Int? max_samples_per_shard_step3

    File ploidy_table
    String chr_x
    String chr_y

    File? svtk_to_gatk_script  # For debugging

    Boolean use_hail
    String? gcs_project

    String linux_docker
    String sv_base_mini_docker
    String sv_pipeline_docker
    String sv_pipeline_hail_docker
    String sv_pipeline_updates_docker

    # overrides for local tasks
    RuntimeAttr? runtime_override_clean_vcf_1a
    RuntimeAttr? runtime_override_clean_vcf_2
    RuntimeAttr? runtime_override_clean_vcf_3
    RuntimeAttr? runtime_override_clean_vcf_4
    RuntimeAttr? runtime_override_clean_vcf_5_scatter
    RuntimeAttr? runtime_override_clean_vcf_5_make_cleangq
    RuntimeAttr? runtime_override_clean_vcf_5_find_redundant_multiallelics
    RuntimeAttr? runtime_override_clean_vcf_5_polish
    RuntimeAttr? runtime_override_stitch_fragmented_cnvs
    RuntimeAttr? runtime_override_final_cleanup

    # Clean vcf 1b
    RuntimeAttr? runtime_attr_override_subset_large_cnvs_1b
    RuntimeAttr? runtime_attr_override_sort_bed_1b
    RuntimeAttr? runtime_attr_override_intersect_bed_1b
    RuntimeAttr? runtime_attr_override_build_dict_1b
    RuntimeAttr? runtime_attr_override_scatter_1b
    RuntimeAttr? runtime_attr_override_filter_vcf_1b
    RuntimeAttr? runtime_override_concat_vcfs_1b
    RuntimeAttr? runtime_override_cat_multi_cnvs_1b

    RuntimeAttr? runtime_override_preconcat_step1
    RuntimeAttr? runtime_override_hail_merge_step1
    RuntimeAttr? runtime_override_fix_header_step1

    RuntimeAttr? runtime_override_preconcat_drc
    RuntimeAttr? runtime_override_hail_merge_drc
    RuntimeAttr? runtime_override_fix_header_drc

    # overrides for MiniTasks
    RuntimeAttr? runtime_override_split_vcf_to_clean
    RuntimeAttr? runtime_override_combine_step_1_sex_chr_revisions
    RuntimeAttr? runtime_override_split_include_list
    RuntimeAttr? runtime_override_combine_clean_vcf_2
    RuntimeAttr? runtime_override_combine_revised_4
    RuntimeAttr? runtime_override_combine_multi_ids_4
    RuntimeAttr? runtime_override_drop_redundant_cnvs
    RuntimeAttr? runtime_override_combine_step_1_vcfs
    RuntimeAttr? runtime_override_sort_drop_redundant_cnvs
    RuntimeAttr? runtime_attr_format

  }

  call MiniTasks.SplitVcf as SplitVcfToClean {
    input:
      vcf=vcf,
      contig=contig,
      prefix="~{prefix}.shard_",
      n_shards=max_shards_per_chrom_step1,
      min_vars_per_shard=min_records_per_shard_step1,
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_override_split_vcf_to_clean
  }

  scatter ( i in range(length(SplitVcfToClean.vcf_shards)) ) {
    call CleanVcf1a {
      input:
        vcf=SplitVcfToClean.vcf_shards[i],
        prefix="~{prefix}.clean_vcf_1.shard_~{i}",
        background_fail_list=background_list,
        bothsides_pass_list=bothsides_pass_list,
        ped_file=ped_file,
        allosome_fai=allosome_fai,
        chr_x=chr_x,
        chr_y=chr_y,
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_override_clean_vcf_1a
    }
  }

  if (use_hail) {
    call HailMerge.HailMerge as CombineStep1VcfsHail {
      input:
        vcfs=CleanVcf1a.intermediate_vcf,
        prefix="~{prefix}.combine_step_1_vcfs",
        gcs_project=gcs_project,
        sv_base_mini_docker=sv_base_mini_docker,
        sv_pipeline_docker=sv_pipeline_docker,
        sv_pipeline_hail_docker=sv_pipeline_hail_docker,
        runtime_override_preconcat=runtime_override_preconcat_step1,
        runtime_override_hail_merge=runtime_override_hail_merge_step1,
        runtime_override_fix_header=runtime_override_fix_header_step1
    }
  }
  if (!use_hail) {
    call MiniTasks.ConcatVcfs as CombineStep1Vcfs {
      input:
        vcfs=CleanVcf1a.intermediate_vcf,
        vcfs_idx=CleanVcf1a.intermediate_vcf_idx,
        naive=true,
        generate_index=false,
        outfile_prefix="~{prefix}.combine_step_1_vcfs",
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_attr_override=runtime_override_combine_step_1_vcfs
    }
  }

  call MiniTasks.CatUncompressedFiles as CombineStep1SexChrRevisions {
    input:
      shards=CleanVcf1a.sex,
      outfile_name="~{prefix}.combine_step_1_sex_chr_revisions.txt",
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_override_combine_step_1_sex_chr_revisions
  }

  call c1b.CleanVcf1b {
    input:
      intermediate_vcf=select_first([CombineStep1Vcfs.concat_vcf, CombineStep1VcfsHail.merged_vcf]),
      prefix="~{prefix}.clean_vcf_1b",
      records_per_shard=clean_vcf1b_records_per_shard,
      sv_pipeline_docker=sv_pipeline_docker,
      sv_pipeline_updates_docker=sv_pipeline_updates_docker,
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override_subset_large_cnvs=runtime_attr_override_subset_large_cnvs_1b,
      runtime_attr_override_sort_bed=runtime_attr_override_sort_bed_1b,
      runtime_attr_override_intersect_bed=runtime_attr_override_intersect_bed_1b,
      runtime_attr_override_build_dict=runtime_attr_override_build_dict_1b,
      runtime_attr_override_scatter=runtime_attr_override_scatter_1b,
      runtime_attr_override_filter_vcf=runtime_attr_override_filter_vcf_1b,
      runtime_override_concat_vcfs=runtime_override_concat_vcfs_1b,
      runtime_override_cat_multi_cnvs=runtime_override_cat_multi_cnvs_1b
  }

  call MiniTasks.SplitUncompressed as SplitIncludeList {
    input:
      whole_file=CleanVcf1a.include_list[0],
      lines_per_shard=samples_per_step2_shard,
      shard_prefix="~{prefix}.split_include_list.",
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_override_split_include_list
  }

  scatter ( i in range(length(SplitIncludeList.shards)) ){
    call CleanVcf2 {
      input:
        normal_revise_vcf=CleanVcf1b.normal,
        prefix="~{prefix}.clean_vcf_2.shard_~{i}",
        include_list=SplitIncludeList.shards[i],
        multi_cnvs=CleanVcf1b.multi,
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_override_clean_vcf_2
      }
  }

  call MiniTasks.CatUncompressedFiles as CombineCleanVcf2 {
    input:
      shards=CleanVcf2.out,
      outfile_name="~{prefix}.combine_clean_vcf_2.txt",
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_override_combine_clean_vcf_2
  }

  call CleanVcf3 {
    input:
      rd_cn_revise=CombineCleanVcf2.outfile,
      max_samples_shard = max_samples_per_shard_step3,
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_override_clean_vcf_3
  }

  scatter ( i in range(length(CleanVcf3.shards)) ){
    call CleanVcf4 {
      input:
        rd_cn_revise=CleanVcf3.shards[i],
        normal_revise_vcf=CleanVcf1b.normal,
        prefix="~{prefix}.clean_vcf_4.shard_~{i}",
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_override_clean_vcf_4
    }
  }

  call MiniTasks.CatUncompressedFiles as CombineRevised4 {
    input:
      shards=CleanVcf4.out,
      outfile_name="~{prefix}.combine_revised_4.txt.gz",
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_override_combine_revised_4
  }

  call MiniTasks.CatUncompressedFiles as CombineMultiIds4 {
    input:
      shards=CleanVcf4.multi_ids,
      outfile_name="~{prefix}.combine_multi_ids_4.txt.gz",
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_override_combine_multi_ids_4
  }

  call c5.CleanVcf5 {
    input:
      revise_vcf_lines=CombineRevised4.outfile,
      normal_revise_vcf=CleanVcf1b.normal,
      ped_file=ped_file,
      sex_chr_revise=CombineStep1SexChrRevisions.outfile,
      multi_ids=CombineMultiIds4.outfile,
      outlier_samples_list=outlier_samples_list,
      contig=contig,
      prefix="~{prefix}.clean_vcf_5",
      records_per_shard=clean_vcf5_records_per_shard,
      threads_per_task=clean_vcf5_threads_per_task,
      sv_pipeline_docker=sv_pipeline_updates_docker,
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override_scatter=runtime_override_clean_vcf_5_scatter,
      runtime_attr_override_make_cleangq=runtime_override_clean_vcf_5_make_cleangq,
      runtime_attr_override_find_redundant_multiallelics=runtime_override_clean_vcf_5_find_redundant_multiallelics,
      runtime_attr_override_polish=runtime_override_clean_vcf_5_polish
  }

  call DropRedundantCnvs {
    input:
      vcf=CleanVcf5.polished,
      prefix="~{prefix}.drop_redundant_cnvs",
      contig=contig,
      sv_pipeline_docker=sv_pipeline_updates_docker,
      runtime_attr_override=runtime_override_drop_redundant_cnvs
  }

  if (use_hail) {
    call HailMerge.HailMerge as SortDropRedundantCnvsHail {
      input:
        vcfs=[DropRedundantCnvs.out],
        prefix="~{prefix}.drop_redundant_cnvs.sorted",
        gcs_project=gcs_project,
        reset_cnv_gts=true,
        sv_base_mini_docker=sv_base_mini_docker,
        sv_pipeline_docker=sv_pipeline_docker,
        sv_pipeline_hail_docker=sv_pipeline_hail_docker,
        runtime_override_preconcat=runtime_override_preconcat_drc,
        runtime_override_hail_merge=runtime_override_hail_merge_drc,
        runtime_override_fix_header=runtime_override_fix_header_drc
    }
  }
  if (!use_hail) {
    call MiniTasks.SortVcf as SortDropRedundantCnvs {
      input:
        vcf=DropRedundantCnvs.out,
        outfile_prefix="~{prefix}.drop_redundant_cnvs.sorted",
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_attr_override=runtime_override_sort_drop_redundant_cnvs
    }
  }

  call StitchFragmentedCnvs {
    input:
      vcf=select_first([SortDropRedundantCnvs.out, SortDropRedundantCnvsHail.merged_vcf]),
      prefix="~{prefix}.stitch_fragmented_cnvs",
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_override_stitch_fragmented_cnvs
  }

  call FinalCleanup {
    input:
      vcf=StitchFragmentedCnvs.stitched_vcf_shard,
      contig=contig,
      prefix="~{prefix}.final_cleanup",
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_override_final_cleanup
  }

  call fvcf.FormatVcf {
    input:
      vcf=FinalCleanup.final_cleaned_shard,
      ploidy_table=ploidy_table,
      args="--scale-down-gq",
      output_prefix="~{prefix}.final_format",
      script=svtk_to_gatk_script,
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_attr_format
  }
  
  output {
    File out = FormatVcf.out
    File out_idx = FormatVcf.out_index
  }
}


task CleanVcf1a {
  input {
    File vcf
    String prefix
    File background_fail_list
    File bothsides_pass_list
    File ped_file
    File allosome_fai
    String chr_x
    String chr_y
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size([vcf, background_fail_list, bothsides_pass_list], "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 3.75,
                                  disk_gb: ceil(10.0 + input_size * 2),
                                  cpu_cores: 1,
                                  preemptible_tries: 3,
                                  max_retries: 1,
                                  boot_disk_gb: 10
                                }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euo pipefail

    touch ~{prefix}.includelist.txt
    touch ~{prefix}.sexchr.revise.txt

    # outputs
    # includelist.txt: the names of all the samples in the input vcf
    # sexchr.revise.txt: the names of the events where genotypes got tweaked on allosomes
    # stdout: a revised vcf
    java -jar $CLEAN_VCF_PART_1_JAR \
      ~{vcf} \
      ~{ped_file} \
      ~{chr_x} \
      ~{chr_y} \
      ~{background_fail_list} \
      ~{bothsides_pass_list} \
      ~{prefix}.includelist.txt \
      ~{prefix}.sexchr.revise.txt \
      | bgzip \
      > ~{prefix}.vcf.gz
    tabix ~{prefix}.vcf.gz
  >>>

  output {
    File include_list="~{prefix}.includelist.txt"
    File sex="~{prefix}.sexchr.revise.txt"
    File intermediate_vcf="~{prefix}.vcf.gz"
    File intermediate_vcf_idx="~{prefix}.vcf.gz.tbi"
  }
}

task CleanVcf2 {
  input {
    File normal_revise_vcf
    String prefix
    File include_list
    File multi_cnvs
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  # generally assume working disk size is ~2 * inputs, and outputs are ~2 *inputs, and inputs are not removed
  # generally assume working memory is ~3 * inputs
  Float input_size = size([normal_revise_vcf, include_list, multi_cnvs], "GB")
  Float base_disk_gb = 10.0
  Float input_disk_scale = 3.0
  RuntimeAttr runtime_default = object {
    mem_gb: 2.0,
    disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
    cpu_cores: 1,
    preemptible_tries: 3,
    max_retries: 1,
    boot_disk_gb: 10
  }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -eu -o pipefail

    bcftools index ~{normal_revise_vcf}
    /opt/sv-pipeline/04_variant_resolution/scripts/clean_vcf_part2.sh \
      ~{normal_revise_vcf} \
      ~{include_list} \
      ~{multi_cnvs} \
      "~{prefix}.txt"
  >>>

  output {
    File out="~{prefix}.txt"
  }
}


task CleanVcf3 {
  input {
    File rd_cn_revise
    Int? max_samples_shard
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }
  Int max_samples_shard_ = select_first([max_samples_shard, 7000])
  # generally assume working disk size is ~2 * inputs, and outputs are ~2 *inputs, and inputs are not removed
  # generally assume working memory is ~3 * inputs
  Float input_size = size(rd_cn_revise, "GB")
  RuntimeAttr runtime_default = object {
    mem_gb: 3.75,
    disk_gb: ceil(10.0 + input_size * 2.0),
    cpu_cores: 1,
    preemptible_tries: 3,
    max_retries: 1,
    boot_disk_gb: 10
  }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euo pipefail
    python /opt/sv-pipeline/04_variant_resolution/scripts/clean_vcf_part3.py ~{rd_cn_revise} -s ~{max_samples_shard_}
    # Ensure there is at least one shard
    touch shards/out.0_0.txt
  >>>

  output {
     Array[File] shards = glob("shards/*")
  }
}


task CleanVcf4 {
  input {
    File rd_cn_revise
    File normal_revise_vcf
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size([rd_cn_revise, normal_revise_vcf], "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 2.0,
                                  disk_gb: 50,
                                  cpu_cores: 1,
                                  preemptible_tries: 3,
                                  max_retries: 1,
                                  boot_disk_gb: 10
                                }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euo pipefail
    python3 <<CODE
    import pysam
    import os

    # Inputs
    REGENO_FILE="~{rd_cn_revise}"
    VCF_FILE="~{normal_revise_vcf}"

    # Build map of variants to regenotype
    with open(REGENO_FILE) as f:
      vid_sample_cn_map = {}
      for line in f:
        tokens = line.strip().split('\t')
        vid = tokens[0]
        if vid not in vid_sample_cn_map:
          vid_sample_cn_map[vid] = []
        vid_sample_cn_map[vid].append(tuple(tokens[1:]))

    # Traverse VCF and replace genotypes
    with open("~{prefix}.revise_vcf_lines.txt", "w") as f:
      vcf = pysam.VariantFile(VCF_FILE)
      num_vcf_records = 0
      for record in vcf:
        num_vcf_records += 1
        if record.id not in vid_sample_cn_map:
          continue
        for entry in vid_sample_cn_map[record.id]:
          s = record.samples[entry[0]]
          s['GT'] = (0, 1)
          s['RD_CN'] = int(entry[1])
        f.write(str(record))
      vcf.close()

    # Get batch size
    regeno_file_name_tokens = os.path.basename(REGENO_FILE).split('.')[1].split('_')
    batch_num = max(int(regeno_file_name_tokens[0]), 1)
    total_batch = max(int(regeno_file_name_tokens[1]), 1)
    segments = num_vcf_records / float(total_batch)
    print("{} {} {}".format(batch_num, total_batch, segments))

    vcf = pysam.VariantFile(VCF_FILE)
    # Max sample count with PE or SR GT over 3
    max_vf = max(len(vcf.header.samples) * 0.01, 2)
    record_start = (batch_num - 1) * segments
    record_end = batch_num * segments
    record_idx = 0
    print("{} {} {}".format(max_vf, record_start, record_end))
    multi_geno_ids = set([])
    for record in vcf:
      record_idx += 1
      if record_idx < record_start:
        continue
      elif record_idx > record_end:
        break
      num_gt_over_2 = 0
      for sid in record.samples:
        s = record.samples[sid]
        # Pick best GT
        if s['PE_GT'] is None:
          continue
        elif s['SR_GT'] is None:
          gt = s['PE_GT']
        elif s['PE_GT'] > 0 and s['SR_GT'] == 0:
          gt = s['PE_GT']
        elif s['PE_GT'] == 0:
          gt = s['SR_GT']
        elif s['PE_GQ'] >= s['SR_GQ']:
          gt = s['PE_GT']
        else:
          gt = s['SR_GT']
        if gt > 2:
          num_gt_over_2 += 1
      if num_gt_over_2 > max_vf:
        multi_geno_ids.add(record.id)
    vcf.close()

    multi_geno_ids = sorted(list(multi_geno_ids))
    with open("~{prefix}.multi_geno_ids.txt", "w") as f:
      for vid in multi_geno_ids:
        f.write(vid + "\n")
    CODE

    bgzip ~{prefix}.revise_vcf_lines.txt
    gzip ~{prefix}.multi_geno_ids.txt
  >>>

  output {
    File out="~{prefix}.revise_vcf_lines.txt.gz"
    File multi_ids="~{prefix}.multi_geno_ids.txt.gz"
  }
}


# Remove CNVs that are redundant with CPX events or other CNVs
task DropRedundantCnvs {
  input {
    File vcf
    String prefix
    String contig
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size(vcf, "GiB")
  # disk is cheap, read/write speed is proportional to disk size, and disk IO is a significant time factor:
  # in tests on large VCFs, memory usage is ~1.0 * input VCF size
  # the biggest disk usage is at the end of the task, with input + output VCF on disk
  Int cpu_cores = 2 # speed up compression / decompression of VCFs
  RuntimeAttr runtime_default = object {
    mem_gb: 3.75 + input_size * 1.5,
    disk_gb: ceil(100.0 + input_size * 2.0),
    cpu_cores: cpu_cores,
    preemptible_tries: 3,
    max_retries: 1,
    boot_disk_gb: 10
  }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euo pipefail
    /opt/sv-pipeline/04_variant_resolution/scripts/resolve_cpx_cnv_redundancies.py \
      ~{vcf} ~{prefix}.vcf.gz --temp-dir ./tmp
  >>>

  output {
    File out = "~{prefix}.vcf.gz"
  }
}


# Stitch fragmented RD-only calls found in 100% of the same samples
task StitchFragmentedCnvs {
  input {
    File vcf
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size(vcf, "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 7.5,
                                  disk_gb: ceil(10.0 + input_size * 2),
                                  cpu_cores: 1,
                                  preemptible_tries: 3,
                                  max_retries: 1,
                                  boot_disk_gb: 10
                                }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  Float mem_gb = select_first([runtime_override.mem_gb, runtime_default.mem_gb])
  Int java_mem_mb = ceil(mem_gb * 1000 * 0.8)

  runtime {
    memory: "~{mem_gb} GB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euo pipefail
    echo "First pass..."
    java -Xmx~{java_mem_mb}M -jar ${STITCH_JAR} 0.2 200000 0.2 ~{vcf} \
      | bgzip \
      > tmp.vcf.gz
    rm ~{vcf}
    echo "Second pass..."
    java -Xmx~{java_mem_mb}M -jar ${STITCH_JAR} 0.2 200000 0.2 tmp.vcf.gz \
      | bgzip \
      > ~{prefix}.vcf.gz
  >>>

  output {
    File stitched_vcf_shard = "~{prefix}.vcf.gz"
  }
}


# Final VCF cleanup
task FinalCleanup {
  input {
    File vcf
    String contig
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  # generally assume working disk size is ~2 * inputs, and outputs are ~2 *inputs, and inputs are not removed
  # generally assume working memory is ~3 * inputs
  Float input_size = size(vcf, "GB")
  Float base_disk_gb = 10.0
  Float base_mem_gb = 2.0
  Float input_mem_scale = 3.0
  Float input_disk_scale = 5.0
  RuntimeAttr runtime_default = object {
    mem_gb: base_mem_gb + input_size * input_mem_scale,
    disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
    cpu_cores: 1,
    preemptible_tries: 3,
    max_retries: 1,
    boot_disk_gb: 10
  }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -eu -o pipefail
    
    /opt/sv-pipeline/04_variant_resolution/scripts/rename_after_vcfcluster.py \
      --chrom ~{contig} \
      --prefix ~{prefix} \
      ~{vcf} stdout \
      | bcftools annotate --no-version -e 'SVTYPE=="CNV" && SVLEN<5000' -x INFO/MEMBERS -Oz -o ~{prefix}.vcf.gz
    tabix ~{prefix}.vcf.gz
  >>>

  output {
    File final_cleaned_shard = "~{prefix}.vcf.gz"
    File final_cleaned_shard_idx = "~{prefix}.vcf.gz.tbi"
  }
}