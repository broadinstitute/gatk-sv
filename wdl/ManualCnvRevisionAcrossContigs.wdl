version 1.0

import "ManualCnvRevision.wdl" as revise

workflow ManualCnvRevisionAcrossContigs {
  input {
    Array[File] vcfs
    File ped_file
    String output_prefix

    File primary_contigs_list

    # See src/sv-pipeline/scripts/manual_review.py for file descriptions and formats
    Array[File]? new_cnv_table
    Array[File]? new_ctx_table
    Array[File]? remove_vids_list
    Array[File]? multiallelic_vids_list
    Array[File]? add_call_table
    Array[File]? remove_call_table
    Array[File]? filter_call_table
    Array[File]? coords_table
    Array[File]? gd_table
    Array[File]? spanned_del_table

    File cytobands

    Int records_per_shard

    # For concatentation
    Boolean use_hail = false
    String? gcs_project

    File? apply_manual_review_script

    File? NONE_FILE_  # DO NOT USE

    String sv_base_mini_docker
    String sv_pipeline_docker

  }

  Array[String] contigs = read_lines(primary_contigs_list)

  scatter (i in range(length(vcfs))) {
    call revise.ManualCnvRevision {
      input:
        vcf = vcfs[i],
        ped_file=ped_file,
        output_prefix="~{output_prefix}.~{contigs[i]}",
        new_cnv_table=if defined(new_cnv_table) then select_first([new_cnv_table])[i] else NONE_FILE_,
        new_ctx_table=if defined(new_ctx_table) then select_first([new_ctx_table])[i] else NONE_FILE_,
        remove_vids_list=if defined(remove_vids_list) then select_first([remove_vids_list])[i] else NONE_FILE_,
        multiallelic_vids_list=if defined(multiallelic_vids_list) then select_first([multiallelic_vids_list])[i] else NONE_FILE_,
        add_call_table=if defined(add_call_table) then select_first([add_call_table])[i] else NONE_FILE_,
        remove_call_table=if defined(remove_call_table) then select_first([remove_call_table])[i] else NONE_FILE_,
        filter_call_table=if defined(filter_call_table) then select_first([filter_call_table])[i] else NONE_FILE_,
        coords_table=if defined(coords_table) then select_first([coords_table])[i] else NONE_FILE_,
        gd_table=if defined(gd_table) then select_first([gd_table])[i] else NONE_FILE_,
        spanned_del_table=if defined(spanned_del_table) then select_first([spanned_del_table])[i] else NONE_FILE_,
        cytobands=cytobands,
        records_per_shard=records_per_shard,
        use_hail=use_hail,
        gcs_project=gcs_project,
        apply_manual_review_script=apply_manual_review_script,
        sv_base_mini_docker=sv_base_mini_docker,
        sv_pipeline_docker=sv_pipeline_docker
    }
  }

  output {
    Array[File] manual_revision_vcfs = ManualCnvRevision.manual_cnv_revision_vcf
    Array[File] manual_revision_vcf_idxs = ManualCnvRevision.manual_cnv_revision_vcf_index
  }
}
