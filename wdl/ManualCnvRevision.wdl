version 1.0

import "Structs.wdl"
import "HailMerge.wdl" as HailMerge
import "TasksMakeCohortVcf.wdl" as tasks
import "Utils.wdl" as utils

# Applies some post-hoc revisions to CTX and large CNV/CPX/INV. Intended to be run following genotype filtering.

workflow ManualCnvRevision {
  input {
    File vcf
    File ped_file
    String output_prefix

    # See src/sv-pipeline/scripts/manual_review.py for file descriptions and formats
    File? new_cnv_table
    File? new_ctx_table
    File? remove_vids_list
    File? multiallelic_vids_list
    File? add_call_table
    File? remove_call_table
    File? filter_call_table
    File? coords_table
    File? gd_table
    File? spanned_del_table

    File cytobands

    Int records_per_shard

    # For concatentation
    Boolean use_hail = false
    String? gcs_project

    File? apply_manual_review_script

    String sv_base_mini_docker
    String sv_pipeline_docker

    # Do not use
    File? NONE_FILE_

    RuntimeAttr? runtime_override_hail_preconcat
    RuntimeAttr? runtime_override_hail_merge
    RuntimeAttr? runtime_override_hail_fix_header
    RuntimeAttr? runtime_override_concat

    RuntimeAttr? runtime_attr_override_scatter
    RuntimeAttr? runtime_attr_ids_from_vcf
    RuntimeAttr? runtime_attr_override_spanned_cpx
    RuntimeAttr? runtime_attr_apply
    RuntimeAttr? runtime_attr_apply_cpx
  }

  call tasks.ScatterVcf {
    input:
      vcf=vcf,
      records_per_shard = records_per_shard,
      prefix = "~{output_prefix}.scatter_vcf",
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_attr_override_scatter
  }

  scatter ( i in range(length(ScatterVcf.shards)) ) {
    call ApplyManualReviewUpdates {
      input:
      vcf=ScatterVcf.shards[i],
      ped_file=ped_file,
      prefix="~{output_prefix}.manual_review.shard_~{i}",
      cytobands=cytobands,
      cytobands_index = "~{cytobands}.tbi",
      new_cnv_table=if i == 0 then new_cnv_table else NONE_FILE_,
      new_ctx_table=if i == 0 then new_ctx_table else NONE_FILE_,
      remove_vids_list=remove_vids_list,
      multiallelic_vids_list=multiallelic_vids_list,
      add_call_table=add_call_table,
      remove_call_table=remove_call_table,
      filter_call_table=filter_call_table,
      coords_table=coords_table,
      gd_table=gd_table,
      spanned_del_table=spanned_del_table,
      script=apply_manual_review_script,
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_attr_apply
    }
  }

  Array[File] vcf_shards = ApplyManualReviewUpdates.out
  Array[File] vcf_shard_indexes = ApplyManualReviewUpdates.out_index

  if (use_hail) {
    call HailMerge.HailMerge {
      input:
        vcfs=ApplyManualReviewUpdates.out,
        prefix="~{output_prefix}.manual_cnv_revision",
        gcs_project=gcs_project,
        reset_cnv_gts=true,
        sv_base_mini_docker=sv_base_mini_docker,
        sv_pipeline_docker=sv_pipeline_docker,
        sv_pipeline_hail_docker=sv_pipeline_docker,
        runtime_override_preconcat=runtime_override_hail_preconcat,
        runtime_override_hail_merge=runtime_override_hail_merge,
        runtime_override_fix_header=runtime_override_hail_fix_header
    }
  }
  if (!use_hail) {
    call tasks.ConcatVcfs {
      input:
        vcfs=vcf_shards,
        vcfs_idx=vcf_shard_indexes,
        allow_overlaps=true,
        outfile_prefix="~{output_prefix}.manual_cnv_revision",
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_attr_override=runtime_override_concat
    }
  }

  output {
    File manual_cnv_revision_vcf = select_first([ConcatVcfs.concat_vcf, HailMerge.merged_vcf])
    File manual_cnv_revision_vcf_index = select_first([ConcatVcfs.concat_vcf_idx, HailMerge.merged_vcf_index])
  }
}

task GetSpannedDeletionsFromComplexResolve {
  input {
    File vcf
    File vcf_index
    File vids_list
    File sample_list
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr runtime_default = object {
                                  mem_gb: 3.75,
                                  disk_gb: ceil(10.0 + size(vcf, "GB") * 2.0),
                                  cpu_cores: 1,
                                  preemptible_tries: 1,
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

    # Subset samples, subset variants, and remove EV fields which are still in Integer format
    bcftools view --no-version -S ~{sample_list} -i 'ID=@~{vids_list}' ~{vcf} \
      | bcftools annotate -x FORMAT/EV -Oz -o ~{prefix}.vcf.gz
    tabix ~{prefix}.vcf.gz
  >>>

  output {
    File out="~{prefix}.vcf.gz"
    File out_index="~{prefix}.vcf.gz.tbi"
  }
}

task ApplyManualReviewUpdates {
  input {
    File vcf
    File ped_file
    String prefix

    File cytobands
    File cytobands_index

    File? new_cnv_table
    File? new_ctx_table
    File? remove_vids_list
    File? multiallelic_vids_list
    File? add_call_table
    File? remove_call_table
    File? filter_call_table
    File? coords_table
    File? gd_table
    File? spanned_del_table
    File? spanned_del_vids_list

    # For debugging
    File? script

    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr runtime_default = object {
                                  mem_gb: 3.75,
                                  disk_gb: ceil(50.0 + size(vcf, "GB") * 5.0),
                                  cpu_cores: 1,
                                  preemptible_tries: 1,
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
    python ~{default="/opt/sv-pipeline/scripts/manual_review.py" script} \
      --vcf ~{vcf} \
      --out ~{prefix}.unsorted.vcf.gz \
      --ped-file ~{ped_file} \
      --cytobands ~{cytobands} \
      ~{"--new-cnv-table " + new_cnv_table} \
      ~{"--new-ctx-table " + new_ctx_table} \
      ~{"--remove-vids-list " + remove_vids_list} \
      ~{"--multiallelic-vids-list " + multiallelic_vids_list} \
      ~{"--add-call-table " + add_call_table} \
      ~{"--remove-call-table " + remove_call_table} \
      ~{"--filter-call-table " + filter_call_table} \
      ~{"--coords-table " + coords_table} \
      ~{"--gd-table " + gd_table} \
      ~{"--spanned-del-table " + spanned_del_table} \
      ~{"--spanned-del-vids-list " + spanned_del_vids_list}
    mkdir tmp
    bcftools sort -T ./tmp ~{prefix}.unsorted.vcf.gz -Oz -o ~{prefix}.vcf.gz
    tabix ~{prefix}.vcf.gz
  >>>

  output {
    File out="~{prefix}.vcf.gz"
    File out_index="~{prefix}.vcf.gz.tbi"
  }
}
