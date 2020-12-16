version 1.0

import "ResolveCpxSv.wdl" as ResolveComplexContig
import "TasksMakeCohortVcf.wdl" as MiniTasks

workflow ResolveComplexVariants {
  input {
    String cohort_name

    Boolean merge_vcfs = false

    Array[File] cluster_vcfs
    Array[File] cluster_bothside_pass_lists
    Array[File] cluster_background_fail_lists

    Array[File] disc_files
    Array[File] rf_cutoff_files

    File contig_list
    Int max_shard_size
    File cytobands
    File mei_bed
    File pe_exclude_list
    File ref_dict

    File hail_script
    String project

    String sv_base_mini_docker
    String sv_pipeline_docker

    # overrides for local tasks
    RuntimeAttr? runtime_override_update_sr_list_pass
    RuntimeAttr? runtime_override_update_sr_list_fail
    RuntimeAttr? runtime_override_integrate_resolved_vcfs
    RuntimeAttr? runtime_override_rename_variants

    RuntimeAttr? runtime_override_breakpoint_overlap_filter

    # overrides for mini tasks
    RuntimeAttr? runtime_override_subset_inversions
    RuntimeAttr? runtime_override_concat

    # overrides for ResolveComplexContig
    RuntimeAttr? runtime_override_get_se_cutoff
    RuntimeAttr? runtime_override_shard_vcf_cpx
    RuntimeAttr? runtime_override_shard_vids
    RuntimeAttr? runtime_override_resolve_prep
    RuntimeAttr? runtime_override_resolve_cpx_per_shard
    RuntimeAttr? runtime_override_restore_unresolved_cnv_per_shard
    RuntimeAttr? runtime_override_concat_resolved_per_shard
    RuntimeAttr? runtime_override_merge_resolve_inner
    RuntimeAttr? runtime_override_pull_vcf_shard
    RuntimeAttr? runtime_override_preconcat
    RuntimeAttr? runtime_override_hail_merge
    RuntimeAttr? runtime_override_fix_header

    RuntimeAttr? runtime_override_get_se_cutoff_inv
    RuntimeAttr? runtime_override_shard_vcf_cpx_inv
    RuntimeAttr? runtime_override_shard_vids_inv
    RuntimeAttr? runtime_override_resolve_prep_inv
    RuntimeAttr? runtime_override_resolve_cpx_per_shard_inv
    RuntimeAttr? runtime_override_restore_unresolved_cnv_per_shard_inv
    RuntimeAttr? runtime_override_concat_resolved_per_shard_inv
    RuntimeAttr? runtime_override_merge_resolve_inner_inv
    RuntimeAttr? runtime_override_pull_vcf_shard_inv
    RuntimeAttr? runtime_override_preconcat_inv
    RuntimeAttr? runtime_override_hail_merge_inv
    RuntimeAttr? runtime_override_fix_header_inv
  }

  #Scatter per chromosome
  Array[String] contigs = transpose(read_tsv(contig_list))[0]
  scatter ( i in range(length(contigs)) ) {

    String contig = contigs[i]

    #Subset inversions from PESR+RD VCF
    call MiniTasks.FilterVcf as SubsetInversions {
      input:
        vcf=cluster_vcfs[i],
        vcf_index=cluster_vcfs[i] + ".tbi",
        outfile_prefix="~{cohort_name}.~{contig}.inversions_only",
        records_filter='INFO/SVTYPE="INV"',
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_attr_override=runtime_override_subset_inversions
    }

    #Resolve inversion-only VCF
    call ResolveComplexContig.ResolveComplexSv as ResolveCpxInv {
      input:
        vcf=SubsetInversions.filtered_vcf,
        prefix="~{cohort_name}.~{contig}.inv_only",
        contig=contig,
        max_shard_size=max_shard_size,
        cytobands=cytobands,
        disc_files=disc_files,
        mei_bed=mei_bed,
        pe_exclude_list=pe_exclude_list,
        rf_cutoff_files=rf_cutoff_files,
        ref_dict=ref_dict,
        precluster_distance=50000000,
        precluster_overlap_frac=0.1,
        hail_script=hail_script,
        project=project,
        sv_pipeline_docker=sv_pipeline_docker,
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_override_get_se_cutoff=runtime_override_get_se_cutoff_inv,
        runtime_override_shard_vcf_cpx=runtime_override_shard_vcf_cpx_inv,
        runtime_override_shard_vids=runtime_override_shard_vids_inv,
        runtime_override_resolve_prep=runtime_override_resolve_prep_inv,
        runtime_override_resolve_cpx_per_shard=runtime_override_resolve_cpx_per_shard_inv,
        runtime_override_restore_unresolved_cnv_per_shard=runtime_override_restore_unresolved_cnv_per_shard_inv,
        runtime_override_concat_resolved_per_shard=runtime_override_concat_resolved_per_shard_inv,
        runtime_override_merge_resolve_inner=runtime_override_merge_resolve_inner_inv,
        runtime_override_pull_vcf_shard=runtime_override_pull_vcf_shard_inv,
        runtime_override_preconcat=runtime_override_preconcat_inv,
        runtime_override_hail_merge=runtime_override_hail_merge_inv,
        runtime_override_fix_header=runtime_override_fix_header_inv
    }

    #Run same-bp overlap filter on full vcf
    call BreakpointOverlap {
      input:
        vcf=cluster_vcfs[i],
        vcf_index=cluster_vcfs[i] + ".tbi",
        prefix="~{cohort_name}.~{contig}.breakpoint_overlap",
        bothside_pass_list=cluster_bothside_pass_lists[i],
        background_fail_list=cluster_background_fail_lists[i],
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_override_breakpoint_overlap_filter
    }

    #Resolve all-variants VCF after same-bp overlap filter
    call ResolveComplexContig.ResolveComplexSv as ResolveCpxAll {
      input:
        vcf=BreakpointOverlap.out,
        prefix="~{cohort_name}.~{contig}.all",
        contig=contig,
        max_shard_size=max_shard_size,
        cytobands=cytobands,
        disc_files=disc_files,
        mei_bed=mei_bed,
        pe_exclude_list=pe_exclude_list,
        rf_cutoff_files=rf_cutoff_files,
        ref_dict=ref_dict,
        precluster_distance=2000,
        precluster_overlap_frac=0.000000001,
        hail_script=hail_script,
        project=project,
        sv_pipeline_docker=sv_pipeline_docker,
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_override_get_se_cutoff=runtime_override_get_se_cutoff,
        runtime_override_shard_vcf_cpx=runtime_override_shard_vcf_cpx,
        runtime_override_shard_vids=runtime_override_shard_vids,
        runtime_override_resolve_prep=runtime_override_resolve_prep,
        runtime_override_resolve_cpx_per_shard=runtime_override_resolve_cpx_per_shard,
        runtime_override_restore_unresolved_cnv_per_shard=runtime_override_restore_unresolved_cnv_per_shard,
        runtime_override_concat_resolved_per_shard=runtime_override_concat_resolved_per_shard,
        runtime_override_merge_resolve_inner=runtime_override_merge_resolve_inner,
        runtime_override_pull_vcf_shard=runtime_override_pull_vcf_shard,
        runtime_override_preconcat=runtime_override_preconcat,
        runtime_override_hail_merge=runtime_override_hail_merge,
        runtime_override_fix_header=runtime_override_fix_header
    }

    #Integrate inv-only and all-variants resolved VCFs
    call IntegrateResolvedVcfs {
      input:
        inv_res_vcf=ResolveCpxInv.resolved_vcf_merged,
        all_res_vcf=ResolveCpxAll.resolved_vcf_merged,
        prefix="~{cohort_name}.~{contig}.resolved",
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_override_integrate_resolved_vcfs
    }

    #Apply consistent variant naming scheme to integrated VCF
    call RenameVariants {
      input:
        vcf=IntegrateResolvedVcfs.integrated_vcf,
        prefix="~{cohort_name}.~{contig}.renamed",
        chrom=contig,
        vid_prefix=cohort_name,
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_override_rename_variants
    }

    #Update SR background fail & bothside pass files
    call MiniTasks.UpdateSrList as UpdateBothsidePass {
      input:
        vcf=RenameVariants.renamed_vcf,
        original_list=cluster_bothside_pass_lists[i],
        outfile="~{cohort_name}.~{contig}.sr_bothside_pass.updated3.txt",
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_override_update_sr_list_pass
    }

    #Update SR background fail & bothside pass files
    call MiniTasks.UpdateSrList as UpdateBackgroundFail {
      input:
        vcf=RenameVariants.renamed_vcf,
        original_list=cluster_background_fail_lists[i],
        outfile="~{cohort_name}.~{contig}.sr_background_fail.updated3.txt",
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_override_update_sr_list_fail
    }
  }

  #Merge resolved vcfs for QC
  if (merge_vcfs) {
    call MiniTasks.ConcatVcfs {
      input:
        vcfs=RenameVariants.renamed_vcf,
        vcfs_idx=RenameVariants.renamed_vcf_index,
        allow_overlaps=true,
        outfile_prefix="~{cohort_name}.complex_resolve",
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
    Array[File] breakpoint_overlap_dropped_record_vcfs = BreakpointOverlap.dropped_record_vcf
    Array[File] breakpoint_overlap_dropped_record_vcf_indexes = BreakpointOverlap.dropped_record_vcf_index
    File? merged_vcf = ConcatVcfs.concat_vcf
    File? merged_vcf_index = ConcatVcfs.concat_vcf_idx
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
  RuntimeAttr runtime_default = object {
                                  mem_gb: 3.75,
                                  disk_gb: ceil(10 + input_size * 20),
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
    set -euxo pipefail

    ##make bed of the inversion resolved vcf##
    zcat ~{inv_res_vcf} \
      |awk '{if ($8!~"UNRESOLVED" || $1~"#") print}' \
      |svtk vcf2bed stdin inv.resolve.bed --no-samples -i MEMBERS

    ##make beds of the fully resolved vcf##
    zcat ~{all_res_vcf} \
      |awk '{if ($8~"UNRESOLVED" || $1~"#") print}' \
      |svtk vcf2bed stdin all.unresolved.inv.bed --no-samples -i MEMBERS

    zcat ~{all_res_vcf} \
      |awk '{if ($8!~"UNRESOLVED" || $1~"#") print}' \
      |svtk vcf2bed stdin all.resolved.inv.bed --no-samples -i MEMBERS

    ##get unresolved variants from full vcf that are resolved in inversion resolved vcf###
    zcat ~{inv_res_vcf} \
      |fgrep -v "#" \
      |awk '{if ($8!~"UNRESOLVED") print}' \
      |awk -F'\t' -v OFS='\t' 'ARGIND==1{inFileA[$1]; next} {if (!($3 in inFileA)) print }' \
        <(awk '{if ($NF!="MEMBERS") print $NF}' all.resolved.inv.bed | tr ',' '\n') - \
      >add.vcf.lines.txt

    ##get unresolved variants id from full vcf to strip since they are resolved in inversion resolved vcf###
    ##inversions that cluster were other variants (rare) are kept as unresolved though they will also be part of a resolved variant in add.vcf.lines.txt##
    awk '{if ($NF!="MEMBERS") print $NF}' inv.resolve.bed \
      |tr ',' '\n'\
      |awk -F'\t' -v OFS='\t' 'ARGIND==1{inFileA[$4]; next} {if ($4 in inFileA) print }' all.resolved.inv.bed - \
      |awk '{if ($NF!~",")print $4}' \
      >remove.unresolved.vcf.ids.txt

    mkdir temp
    zcat ~{all_res_vcf} \
      |awk -F'\t' -v OFS='\t' 'ARGIND==1{inFileA[$1]; next} {if (!($3 in inFileA)) print }' remove.unresolved.vcf.ids.txt - \
      |cat - add.vcf.lines.txt \
      |bcftools sort - -O z -T temp \
      > ~{prefix}.vcf.gz

    tabix ~{prefix}.vcf.gz
  >>>

  output {
    File integrated_vcf = "~{prefix}.vcf.gz"
    File integrated_vcf_idx = "~{prefix}.vcf.gz.tbi"
  }
}


task BreakpointOverlap {
  input {
    File vcf
    File vcf_index
    File bothside_pass_list
    File background_fail_list
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
    set -euo pipefail
    python /opt/sv-pipeline/04_variant_resolution/scripts/overlap_breakpoint_filter.py \
      ~{vcf} \
      ~{bothside_pass_list} \
      ~{background_fail_list} \
      ~{prefix}.dropped_records.vcf.gz \
      | bgzip \
      > ~{prefix}.vcf.gz
    tabix ~{prefix}.vcf.gz
    tabix ~{prefix}.dropped_records.vcf.gz
  >>>

  output {
    File out = "~{prefix}.vcf.gz"
    File out_index = "~{prefix}.vcf.gz.tbi"
    File dropped_record_vcf = "~{prefix}.dropped_records.vcf.gz"
    File dropped_record_vcf_index = "~{prefix}.dropped_records.vcf.gz.tbi"
  }
}

# Rename variants in VCF
task RenameVariants {
  input {
    File vcf
    String vid_prefix
    String prefix
    String chrom
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
    set -euo pipefail
    /opt/sv-pipeline/04_variant_resolution/scripts/rename.py --chrom ~{chrom} --prefix ~{vid_prefix} ~{vcf} - \
      | bgzip \
      > ~{prefix}.vcf.gz
    tabix ~{prefix}.vcf.gz
  >>>

  output {
    File renamed_vcf = "~{prefix}.vcf.gz"
    File renamed_vcf_index = "~{prefix}.vcf.gz.tbi"
  }
}
