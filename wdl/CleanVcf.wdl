version 1.0

import "Structs.wdl"
import "Tasks0506.wdl" as MiniTasks

workflow CleanVcf {
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
    File? outlier_samples_list

    String sv_base_mini_docker
    String sv_pipeline_docker

    # overrides for local tasks
    RuntimeAttr? runtime_override_clean_vcf_1a
    RuntimeAttr? runtime_override_clean_vcf_1b
    RuntimeAttr? runtime_override_clean_vcf_2
    RuntimeAttr? runtime_override_clean_vcf_3
    RuntimeAttr? runtime_override_clean_vcf_4
    RuntimeAttr? runtime_override_clean_vcf_5
    RuntimeAttr? runtime_override_drop_redundant_cnvs
    RuntimeAttr? runtime_override_stitch_fragmented_cnvs
    RuntimeAttr? runtime_override_final_cleanup

    # overrides for MiniTasks
    RuntimeAttr? runtime_override_split_vcf_to_clean
    RuntimeAttr? runtime_override_combine_step_1_vcfs
    RuntimeAttr? runtime_override_combine_step_1_sex_chr_revisions
    RuntimeAttr? runtime_override_split_include_list
    RuntimeAttr? runtime_override_combine_clean_vcf_2
    RuntimeAttr? runtime_override_combine_revised_4
    RuntimeAttr? runtime_override_combine_multi_ids_4

  }

  call MiniTasks.SplitVcf as SplitVcfToClean {
    input:
      vcf=vcf,
      contig=contig,
      prefix="~{prefix}.~{contig}.shard_",
      n_shards=max_shards_per_chrom_step1,
      min_vars_per_shard=min_records_per_shard_step1,
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_override_split_vcf_to_clean
  }

  scatter ( vcf_shard in SplitVcfToClean.vcf_shards ) {
    call CleanVcf1a {
      input:
        vcf=vcf_shard,
        background_list=background_list,
        ped_file=ped_file,
        sv_pipeline_docker=sv_pipeline_docker,
        bothsides_pass_list=bothsides_pass_list,
        allosome_fai=allosome_fai,
        runtime_attr_override=runtime_override_clean_vcf_1a
    }
  }

  call MiniTasks.ConcatVcfs as CombineStep1Vcfs {
    input:
      vcfs=CleanVcf1a.intermediate_vcf,
      vcfs_idx=CleanVcf1a.intermediate_vcf_idx,
      outfile_prefix=prefix + ".cleanVCF_step1.intermediate_vcf.merged",
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_override_combine_step_1_vcfs
  }

  call MiniTasks.CatUncompressedFiles as CombineStep1SexChrRevisions {
    input:
      shards=CleanVcf1a.sex,
      outfile_name=prefix + ".cleanVCF_step1.sexchr_revise.merged.txt",
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_override_combine_step_1_sex_chr_revisions
  }

  call CleanVcf1b {
    input:
      intermediate_vcf=CombineStep1Vcfs.concat_vcf,
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_override_clean_vcf_1b
  }

  call MiniTasks.SplitUncompressed as SplitIncludeList {
    input:
      whole_file=CleanVcf1a.include_list[0],
      lines_per_shard=samples_per_step2_shard,
      shard_prefix="includeexclude.",
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_override_split_include_list
  }

  scatter ( included_interval in SplitIncludeList.shards ){
    call CleanVcf2{
      input:
        normal_revise_vcf=CleanVcf1b.normal,
        include_list=included_interval,
        multi_cnvs=CleanVcf1b.multi,
        vcftools_idx=CleanVcf1b.vcftools_idx,
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_override_clean_vcf_2
      }
  }

  call MiniTasks.CatUncompressedFiles as CombineCleanVcf2 {
    input:
      shards=CleanVcf2.out,
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_override_combine_clean_vcf_2
  }

  call CleanVcf3 {
    input:
      rd_cn_revise=CombineCleanVcf2.outfile,
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_override_clean_vcf_3
  }

  scatter ( rd_cn_revise in CleanVcf3.shards ){
    call CleanVcf4 {
      input:
        rd_cn_revise=rd_cn_revise,
        normal_revise_vcf=CleanVcf1b.normal,
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_override_clean_vcf_4
    }
  }

  call MiniTasks.CatUncompressedFiles as CombineRevised4 {
    input:
      shards=CleanVcf4.out,
      outfile_name="revise.vcf.lines.txt.gz",
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_override_combine_revised_4
  }

  call MiniTasks.CatUncompressedFiles as CombineMultiIds4 {
    input:
      shards=CleanVcf4.multi_ids,
      outfile_name="multi.geno.ids.txt.gz",
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_override_combine_multi_ids_4
  }

  call CleanVcf5 {
    input:
      revise_vcf_lines=CombineRevised4.outfile,
      normal_revise_vcf=CleanVcf1b.normal,
      ped_file=ped_file,
      sex_chr_revise=CombineStep1SexChrRevisions.outfile,
      multi_ids=CombineMultiIds4.outfile,
      outlier_samples_list=outlier_samples_list,
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_override_clean_vcf_5
  }

  call DropRedundantCnvs {
    input:
      vcf=CleanVcf5.polished,
      contig=contig,
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_override_drop_redundant_cnvs
  }

  call StitchFragmentedCnvs {
    input:
      vcf=DropRedundantCnvs.cleaned_vcf_shard,
      contig=contig,
      prefix=prefix,
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_override_stitch_fragmented_cnvs
  }

  call FinalCleanup {
    input:
      vcf=StitchFragmentedCnvs.stitched_vcf_shard,
      contig=contig,
      prefix=prefix,
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_override_final_cleanup

  }
  
  output {
    File out=FinalCleanup.final_cleaned_shard
    File out_idx=FinalCleanup.final_cleaned_shard_idx
  }
}


#CleanVCF 1a is sharded
task CleanVcf1a {
  input {
    File vcf
    File background_list
    File ped_file
    String sv_pipeline_docker
    File bothsides_pass_list
    File allosome_fai
    RuntimeAttr? runtime_attr_override
  }

  # generally assume working disk size is ~2 * inputs, and outputs are ~2 *inputs, and inputs are not removed
  # generally assume working memory is ~3 * inputs
  Float shard_size = size([vcf, background_list, ped_file], "GB")
  Float base_disk_gb = 10.0
  Float base_mem_gb = 2.0
  Float input_mem_scale = 3.0
  Float input_disk_scale = 5.0
  RuntimeAttr runtime_default = object {
    mem_gb: base_mem_gb + shard_size * input_mem_scale,
    disk_gb: ceil(base_disk_gb + shard_size * input_disk_scale),
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
    
    /opt/sv-pipeline/04_variant_resolution/scripts/clean_vcf_part1.sh ~{vcf} ~{background_list} ~{ped_file} ~{allosome_fai}
    /opt/sv-pipeline/04_variant_resolution/scripts/add_bothsides_support_filter.py \
      --bgzip \
      --outfile int.w_bothsides.vcf.gz \
      int.vcf.gz \
      ~{bothsides_pass_list}
    tabix int.w_bothsides.vcf.gz
  >>>

  output {
    File include_list="includelist.txt"
    File sex="sexchr.revise.txt"
    File intermediate_vcf="int.w_bothsides.vcf.gz"
    File intermediate_vcf_idx="int.w_bothsides.vcf.gz.tbi"
  }
}


task CleanVcf1b {
  input {
    File intermediate_vcf
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  # generally assume working disk size is ~2 * inputs, and outputs are ~2 *inputs, and inputs are not removed
  # generally assume working memory is ~3 * inputs
  Float input_size = size(intermediate_vcf, "GB")
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
    set -euo pipefail
    
    /opt/sv-pipeline/04_variant_resolution/scripts/clean_vcf_part1b.sh ~{intermediate_vcf}
  >>>

  output {
    File multi="multi.cnvs.txt"
    File normal="normal.revise.vcf.gz"
    File vcftools_idx = "normal.revise.vcf.gz.csi"
  }
}


task CleanVcf2 {
  input {
    File normal_revise_vcf
    File include_list
    File multi_cnvs
    File vcftools_idx
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  # generally assume working disk size is ~2 * inputs, and outputs are ~2 *inputs, and inputs are not removed
  # generally assume working memory is ~3 * inputs
  Float input_size = size([normal_revise_vcf, include_list, multi_cnvs, vcftools_idx], "GB")
  Float base_disk_gb = 10.0
  Float base_mem_gb = 4.0
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

    /opt/sv-pipeline/04_variant_resolution/scripts/clean_vcf_part2.sh \
      ~{normal_revise_vcf} \
      ~{include_list} \
      ~{multi_cnvs} \
      "output.txt"
  >>>

  output {
    File out="output.txt"
  }
}


task CleanVcf3{
  input {
    File rd_cn_revise
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  # generally assume working disk size is ~2 * inputs, and outputs are ~2 *inputs, and inputs are not removed
  # generally assume working memory is ~3 * inputs
  Float input_size = size(rd_cn_revise, "GB")
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
    
    /opt/sv-pipeline/04_variant_resolution/scripts/clean_vcf_part3.sh ~{rd_cn_revise}

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
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  # generally assume working disk size is ~2 * inputs, and outputs are ~2 *inputs, and inputs are not removed
  # generally assume working memory is ~3 * inputs
  Float input_size = size([rd_cn_revise, normal_revise_vcf], "GB")
  RuntimeAttr runtime_default = object {
    mem_gb: 2.0 + input_size * 3.0,
    disk_gb: ceil(10 + input_size * 10),
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

    /opt/sv-pipeline/04_variant_resolution/scripts/clean_vcf_part4.sh ~{rd_cn_revise} ~{normal_revise_vcf}
  >>>

  output {
    File out="revise.vcf.lines.txt.gz"
    File multi_ids="multi.geno.ids.txt.gz"
  }
}


task CleanVcf5 {
  input {
    File revise_vcf_lines
    File normal_revise_vcf
    File ped_file
    File sex_chr_revise
    File multi_ids
    File? outlier_samples_list
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  # generally assume working disk size is ~2 * inputs, and outputs are ~2 *inputs, and inputs are not removed
  # generally assume working memory is ~3 * inputs
  Float input_size = size(
    select_all([revise_vcf_lines, normal_revise_vcf, ped_file, sex_chr_revise, multi_ids, outlier_samples_list]),
    "GB"
  )
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

    ~{if defined(outlier_samples_list) then "ln ~{outlier_samples_list} outliers.txt" else "touch outliers.txt"}

    /opt/sv-pipeline/04_variant_resolution/scripts/clean_vcf_part5.sh \
      ~{revise_vcf_lines} \
      ~{normal_revise_vcf} \
      ~{ped_file} \
      ~{sex_chr_revise} \
      ~{multi_ids} \
      outliers.txt
  >>>

  output {
    File polished="polished.vcf.gz"
  }
}


task DropRedundantCnvs {
  input {
    File vcf
    String contig
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  String outfile_name = contig + ".shard.no_CNV_redundancies.vcf.gz"

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
    
    /opt/sv-pipeline/04_variant_resolution/scripts/resolve_CPX_CNV_redundancies.sh \
      ~{vcf} \
      ~{outfile_name}
  >>>

  output {
    File cleaned_vcf_shard = outfile_name
  }
}


# Stitch fragmented RD-only calls found in 100% of the same samples
task StitchFragmentedCnvs {
  input {
    File vcf
    String contig
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }
  
  String stitched_vcf_name = contig + ".shard.fragmented_CNVs_stitched.vcf.gz"

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
    
    /opt/sv-pipeline/04_variant_resolution/scripts/stitch_fragmented_CNVs.sh \
      ~{vcf} \
      "tmp_~{stitched_vcf_name}"
    
    /opt/sv-pipeline/04_variant_resolution/scripts/stitch_fragmented_CNVs.sh \
      "tmp_~{stitched_vcf_name}" \
      "~{stitched_vcf_name}"
  >>>

  output {
    File stitched_vcf_shard = stitched_vcf_name
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
  
  String cleaned_shard_name = prefix + "." + contig + ".final_cleanup.vcf.gz"

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
      | fgrep -v "##INFO=<ID=HIGH_SR_BACKGROUND" \
      | /opt/sv-pipeline/04_variant_resolution/scripts/sanitize_filter_field.py stdin stdout \
      | fgrep -v "##INFO=<ID=MEMBERS,Number=.,Type=String," \
      | bgzip -c \
      > "~{cleaned_shard_name}"
    tabix ~{cleaned_shard_name}
  >>>

  output {
    File final_cleaned_shard = cleaned_shard_name
    File final_cleaned_shard_idx = cleaned_shard_name + ".tbi"
  }
}
