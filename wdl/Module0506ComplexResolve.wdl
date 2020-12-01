version 1.0

import "ResolveCpxSv.wdl" as ResolveComplexContig
import "Tasks0506.wdl" as MiniTasks

workflow Module0506ComplexResolve {
  input {
    String cohort_name

    Boolean merge_vcfs = false

    Array[File] cluster_vcfs
    Array[File] cluster_bothside_pass_lists
    Array[File] cluster_background_fail_lists

    Array[File] disc_files
    Array[File] rf_cutoff_files

    File contig_list
    Int max_shards_per_chrom
    Int min_variants_per_shard
    File cytobands
    File mei_bed
    File pe_exclude_list
    File ref_dict

    String sv_base_mini_docker
    String sv_pipeline_docker

    # overrides for local tasks
    RuntimeAttr? runtime_override_update_sr_list
    RuntimeAttr? runtime_override_breakpoint_overlap_filter
    RuntimeAttr? runtime_override_integrate_resolved_vcfs
    RuntimeAttr? runtime_override_rename_variants

    # overrides for mini tasks
    RuntimeAttr? runtime_override_subset_inversions
    RuntimeAttr? runtime_override_merge_fam_file_list
    RuntimeAttr? runtime_override_concat

    # overrides for ResolveComplexContig
    RuntimeAttr? runtime_override_get_se_cutoff
    RuntimeAttr? runtime_override_shard_vcf_cpx
    RuntimeAttr? runtime_override_resolve_prep
    RuntimeAttr? runtime_override_resolve_cpx_per_shard
    RuntimeAttr? runtime_override_restore_unresolved_cnv_per_shard
    RuntimeAttr? runtime_override_concat_resolved_per_shard
  }

  #Scatter per chromosome
  Array[String] contigs = transpose(read_tsv(contig_list))[0]
  scatter ( i in range(length(contigs)) ) {

    String contig = contigs[i]

    #Subset inversions from PESR+RD VCF
    call MiniTasks.FilterVcf as SubsetInversions {
      input:
        vcf=cluster_vcfs[i],
        outfile_prefix="~{cohort_name}.~{contig}.inversions_only",
        records_filter="fgrep SVTYPE=INV",
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_attr_override=runtime_override_subset_inversions
    }

    #Resolve inversion-only VCF
    call ResolveComplexContig.ResolveComplexSv as ResolveCpxInv {
      input:
        vcf=SubsetInversions.filtered_vcf,
        prefix="~{cohort_name}.inv_only",
        contig=contig,
        max_shards_per_chrom=max_shards_per_chrom,
        min_variants_per_shard=100,
        cytobands=cytobands,
        disc_files=disc_files,
        mei_bed=mei_bed,
        pe_exclude_list=pe_exclude_list,
        rf_cutoff_files=rf_cutoff_files,
        inv_only=true,
        ref_dict=ref_dict,
        sv_pipeline_docker=sv_pipeline_docker,
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_override_get_se_cutoff=runtime_override_get_se_cutoff,
        runtime_override_shard_vcf_cpx=runtime_override_shard_vcf_cpx,
        runtime_override_resolve_prep=runtime_override_resolve_prep,
        runtime_override_resolve_cpx_per_shard=runtime_override_resolve_cpx_per_shard,
        runtime_override_restore_unresolved_cnv_per_shard=runtime_override_restore_unresolved_cnv_per_shard,
        runtime_override_concat_resolved_per_shard=runtime_override_concat_resolved_per_shard
    }

    #Run same-bp overlap filter on full vcf
    call BreakpointOverlapFilter {
      input:
        vcf=cluster_vcfs[i],
        prefix="~{cohort_name}.~{contig}",
        bothside_pass=cluster_bothside_pass_lists[i],
        background_fail=cluster_background_fail_lists[i],
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_override_breakpoint_overlap_filter
    }

    #Resolve all-variants VCF after same-bp overlap filter
    call ResolveComplexContig.ResolveComplexSv as ResolveCpxAll {
      input:
        vcf=BreakpointOverlapFilter.bp_filtered_vcf,
        prefix="~{cohort_name}.all_variants",
        contig=contig,
        max_shards_per_chrom=max_shards_per_chrom,
        min_variants_per_shard=100,
        cytobands=cytobands,
        disc_files=disc_files,
        mei_bed=mei_bed,
        pe_exclude_list=pe_exclude_list,
        rf_cutoff_files=rf_cutoff_files,
        inv_only=false,
        ref_dict=ref_dict,
        sv_pipeline_docker=sv_pipeline_docker,
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_override_get_se_cutoff=runtime_override_get_se_cutoff,
        runtime_override_shard_vcf_cpx=runtime_override_shard_vcf_cpx,
        runtime_override_resolve_prep=runtime_override_resolve_prep,
        runtime_override_resolve_cpx_per_shard=runtime_override_resolve_cpx_per_shard,
        runtime_override_restore_unresolved_cnv_per_shard=runtime_override_restore_unresolved_cnv_per_shard,
        runtime_override_concat_resolved_per_shard=runtime_override_concat_resolved_per_shard
    }

    #Integrate inv-only and all-variants resolved VCFs
    call IntegrateResolvedVcfs {
      input:
        inv_res_vcf=ResolveCpxInv.resolved_vcf_merged,
        all_res_vcf=ResolveCpxAll.resolved_vcf_merged,
        prefix="~{cohort_name}.resolved.~{contig}",
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_override_integrate_resolved_vcfs
    }

    #Apply consistent variant naming scheme to integrated VCF
    call RenameVariants {
      input:
        vcf=IntegrateResolvedVcfs.integrated_vcf,
        prefix="~{cohort_name}.~{contig}",
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_override_rename_variants
    }

    #Update SR background fail & bothside pass files
    call MiniTasks.UpdateSrList as UpdateBothsidePass {
      input:
        vcf=RenameVariants.renamed_vcf,
        original_list=cluster_bothside_pass_lists[i],
        outfile="sr_bothside_pass.~{contig}.updated3.txt",
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_override_update_sr_list
    }

    #Update SR background fail & bothside pass files
    call MiniTasks.UpdateSrList as UpdateBackgroundFail {
      input:
        vcf=RenameVariants.renamed_vcf,
        original_list=cluster_background_fail_lists[i],
        outfile="sr_background_fail.~{contig}.updated3.txt",
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_override_update_sr_list
    }
  }

  #Merge resolved vcfs for QC
  if (merge_vcfs) {
    call MiniTasks.ConcatVcfs {
      input:
        vcfs=RenameVariants.renamed_vcf,
        vcfs_idx=RenameVariants.renamed_vcf_index,
        merge_sort=true,
        outfile_prefix="~{cohort_name}.0506_complex",
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_attr_override=runtime_override_concat
    }
  }

  #Final outputs
  output {
    Array[File] complex_resolve_vcfs = RenameVariants.renamed_vcf
    Array[File] complex_resolve_vcf_indexes = RenameVariants.renamed_vcf_index
    Array[File] complex_resolve_bothside_pass_lists = UpdateBothsidePass.updated_list
    Array[File] complex_resolve_background_fail_lists = UpdateBackgroundFail.updated_list
    File? merged_vcf = ConcatVcfs.concat_vcf
    File? merged_vcf_index = ConcatVcfs.concat_vcf_idx
  }
}

#Run Harrison's overlapping breakpoint filter prior to complex resolution
task BreakpointOverlapFilter {
  input {
    File vcf
    String prefix
    File bothside_pass
    File background_fail
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  String temp_output_file = "non_redundant.vcf.gz"
  String output_file = prefix + "." + temp_output_file

  Float input_size = size([vcf, bothside_pass, background_fail], "GiB")
  Float base_mem_gb = 2.0
  Float base_disk_gb = 5.0
  RuntimeAttr runtime_default = object {
                                  mem_gb: base_mem_gb,
                                  disk_gb: ceil(base_disk_gb + input_size * 3.0),
                                  cpu_cores: 1,
                                  preemptible_tries: 3,
                                  max_retries: 1,
                                  boot_disk_gb: 10
                                }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GiB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -eu -o pipefail

    /opt/sv-pipeline/04_variant_resolution/scripts/overlapbpchange.sh \
    ~{vcf} \
    ~{background_fail} \
    ~{bothside_pass}

    mv "~{temp_output_file}" "~{output_file}"
    tabix -p vcf -f "~{output_file}"
  >>>

  output {
    File bp_filtered_vcf = output_file
    File bp_filtered_vcf_idx = output_file + ".tbi"
  }
}


#Merge inversion-only and all-variant cpx-resolved outputs
task IntegrateResolvedVcfs {
  input {
    File inv_res_vcf
    File all_res_vcf
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size([inv_res_vcf, all_res_vcf], "GiB")
  Float base_mem_gb = 2.0
  Float base_disk_gb = 5.0
  RuntimeAttr runtime_default = object {
                                  mem_gb: base_mem_gb,
                                  disk_gb: ceil(base_disk_gb + input_size * 3.0),
                                  cpu_cores: 1,
                                  preemptible_tries: 3,
                                  max_retries: 1,
                                  boot_disk_gb: 10
                                }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GiB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -eu -o pipefail

    /opt/sv-pipeline/04_variant_resolution/scripts/Complex_Inversion_Integration.sh \
    ~{inv_res_vcf} \
    ~{all_res_vcf} \
    ~{prefix}.integrated_resolved.vcf.gz

    tabix -p vcf -f "~{prefix}.integrated_resolved.vcf.gz"
  >>>

  output {
    File integrated_vcf = "~{prefix}.integrated_resolved.vcf.gz"
    File integrated_vcf_idx = "~{prefix}.integrated_resolved.vcf.gz.tbi"
  }
}


# Rename variants in VCF
task RenameVariants {
  input {
    File vcf
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size(vcf, "GiB")
  Float base_mem_gb = 2.0
  Float base_disk_gb = 5.0
  RuntimeAttr runtime_default = object {
                                  mem_gb: base_mem_gb,
                                  disk_gb: ceil(base_disk_gb + input_size * 2.0),
                                  cpu_cores: 1,
                                  preemptible_tries: 3,
                                  max_retries: 1,
                                  boot_disk_gb: 10
                                }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GiB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -eu -o pipefail

    /opt/sv-pipeline/04_variant_resolution/scripts/rename.py \
    --prefix ~{prefix} ~{vcf} - \
    | bgzip -c > "~{prefix}.04_renamed.vcf.gz"
    tabix ~{prefix}.04_renamed.vcf.gz
  >>>

  output {
    File renamed_vcf = "~{prefix}.04_renamed.vcf.gz"
    File renamed_vcf_index = "~{prefix}.04_renamed.vcf.gz.tbi"
  }
}
