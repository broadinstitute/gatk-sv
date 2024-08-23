version 1.0

import "CombineSRBothsidePass.wdl" as CombineSRBothsidePass
import "FormatVcfForGatk.wdl" as GatkFormatting
import "TasksClusterBatch.wdl" as ClusterTasks
import "TasksGenotypeBatch.wdl" as GenotypeTasks
import "TasksMakeCohortVcf.wdl" as MiniTasks

workflow CombineBatches {
  input {
    String cohort_name
    Array[String] batches
    File ped_file

    Boolean merge_vcfs = false

    Array[File] pesr_vcfs
    Array[File] depth_vcfs
    # Set to true if using vcfs generated with a prior version, i.e. not ending in "_reformatted.vcf.gz"
    Boolean legacy_vcfs = false

    Array[File] raw_sr_bothside_pass_files
    Array[File] raw_sr_background_fail_files

    File contig_list
    Int localize_shard_size = 100000
    File pe_exclude_list
    File depth_exclude_list
    Float min_sr_background_fail_batches

    File reference_fasta
    File reference_fasta_fai
    File reference_dict
    String? chr_x
    String? chr_y

    File empty_file

    Boolean use_hail = false
    String? gcs_project

    Float? java_mem_fraction

    String gatk_docker
    String sv_base_mini_docker
    String sv_pipeline_docker

    # overrides for local tasks
    RuntimeAttr? runtime_override_update_sr_list
    RuntimeAttr? runtime_override_merge_pesr_depth
    RuntimeAttr? runtime_override_reheader
    RuntimeAttr? runtime_override_pull_header
    RuntimeAttr? runtime_attr_create_ploidy
    RuntimeAttr? runtime_attr_reformat_1
    RuntimeAttr? runtime_attr_reformat_2
    RuntimeAttr? runtime_attr_svcluster

    # overrides for mini tasks
    RuntimeAttr? runtime_attr_get_non_ref_vids
    RuntimeAttr? runtime_attr_calculate_support_frac
    RuntimeAttr? runtime_override_clean_background_fail
    RuntimeAttr? runtime_override_concat
    RuntimeAttr? runtime_override_concat_pesr_depth
    RuntimeAttr? runtime_override_update_fix_pesr_header
    RuntimeAttr? runtime_override_count_samples
    RuntimeAttr? runtime_override_preconcat_pesr_depth
    RuntimeAttr? runtime_override_hail_merge_pesr_depth
    RuntimeAttr? runtime_override_fix_header_pesr_depth
    RuntimeAttr? runtime_override_concat_large_pesr_depth

    # overrides for VcfClusterContig
    RuntimeAttr? runtime_override_localize_vcfs
    RuntimeAttr? runtime_override_join_vcfs
    RuntimeAttr? runtime_override_fix_multiallelic
    RuntimeAttr? runtime_override_fix_ev_tags
    RuntimeAttr? runtime_override_subset_bothside_pass
    RuntimeAttr? runtime_override_subset_background_fail
    RuntimeAttr? runtime_override_subset_sv_type
    RuntimeAttr? runtime_override_shard_clusters
    RuntimeAttr? runtime_override_shard_vids
    RuntimeAttr? runtime_override_pull_vcf_shard
    RuntimeAttr? runtime_override_svtk_vcf_cluster
    RuntimeAttr? runtime_override_get_vcf_header_with_members_info_line
    RuntimeAttr? runtime_override_concat_vcf_cluster
    RuntimeAttr? runtime_override_concat_svtypes
    RuntimeAttr? runtime_override_concat_sharded_cluster
    RuntimeAttr? runtime_override_make_sites_only
    RuntimeAttr? runtime_override_sort_merged_vcf_cluster
    RuntimeAttr? runtime_override_preconcat_sharded_cluster
    RuntimeAttr? runtime_override_hail_merge_sharded_cluster
    RuntimeAttr? runtime_override_fix_header_sharded_cluster

    # overerides for merge pesr depth
    RuntimeAttr? runtime_override_shard_clusters_mpd
    RuntimeAttr? runtime_override_shard_vids_mpd
    RuntimeAttr? runtime_override_pull_vcf_shard_mpd
    RuntimeAttr? runtime_override_merge_pesr_depth_mpd

    RuntimeAttr? runtime_override_sort_merged_vcf_mpd
    RuntimeAttr? runtime_override_subset_small_mpd
    RuntimeAttr? runtime_override_subset_large_mpd
    RuntimeAttr? runtime_override_make_sites_only_mpd
    RuntimeAttr? runtime_override_concat_large_pesr_depth_mpd
    RuntimeAttr? runtime_override_concat_shards_mpd

    RuntimeAttr? runtime_override_preconcat_large_pesr_depth_mpd
    RuntimeAttr? runtime_override_hail_merge_large_pesr_depth_mpd
    RuntimeAttr? runtime_override_fix_header_large_pesr_depth_mpd

    RuntimeAttr? runtime_override_preconcat_pesr_depth_shards_mpd
    RuntimeAttr? runtime_override_hail_merge_pesr_depth_shards_mpd
    RuntimeAttr? runtime_override_fix_header_pesr_depth_shards_mpd

  }

  # Preprocess some inputs
  call CombineSRBothsidePass.CombineSRBothsidePass {
    input:
      pesr_vcfs=pesr_vcfs,
      raw_sr_bothside_pass_files=raw_sr_bothside_pass_files,
      prefix="~{cohort_name}.sr_bothside_pass",
      sv_base_mini_docker=sv_base_mini_docker,
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_get_non_ref_vids=runtime_attr_get_non_ref_vids,
      runtime_attr_calculate_support_frac=runtime_attr_calculate_support_frac
  }

  Float min_background_fail_first_col = min_sr_background_fail_batches * length(raw_sr_background_fail_files)
  call MiniTasks.CatUncompressedFiles as CombineBackgroundFail {
    input:
      shards=raw_sr_background_fail_files,
      filter_command="sort | uniq -c | awk -v OFS='\\t' '{if($1 >= ~{min_background_fail_first_col}) print $2}'",
      outfile_name="~{cohort_name}.background_fail.txt",
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_override_clean_background_fail
  }

  call ClusterTasks.CreatePloidyTableFromPed {
    input:
      ped_file=ped_file,
      contig_list=contig_list,
      retain_female_chr_y=true,
      chr_x=chr_x,
      chr_y=chr_y,
      output_prefix="~{cohort_name}.ploidy",
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_attr_create_ploidy
  }

  Array[File] all_vcfs = flatten([pesr_vcfs, depth_vcfs])

  scatter (vcf in all_vcfs) {
    if (legacy_vcfs) {
      call GenotypeTasks.ReformatGenotypedVcf {
        input:
          vcf = vcf,
          output_prefix = basename(vcf, ".vcf.gz") + ".reformat_1",
          sv_pipeline_docker = sv_pipeline_docker,
          runtime_attr_override = runtime_attr_reformat_1
      }
    }
    File reformatted_vcf = select_first([ReformatGenotypedVcf.out, vcf])
    call GatkFormatting.FormatVcf {
      input:
        vcf=reformatted_vcf,
        ploidy_table=CreatePloidyTableFromPed.out,
        output_prefix=basename(vcf, ".vcf.gz") + ".reformat_2",
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_attr_reformat_2
    }
  }

  #Scatter per chromosome
  Array[String] contigs = transpose(read_tsv(contig_list))[0]
  scatter ( contig in contigs ) {

    # TODO: disable sorting
    call ClusterTasks.SVCluster {
      input:
        vcfs=FormatVcf.out,
        ploidy_table=CreatePloidyTableFromPed.out,
        output_prefix="~{cohort_name}.combine_batches.~{contig}.svcluster",
        contig=contig,
        fast_mode=false,
        pesr_sample_overlap=0.5,
        pesr_interval_overlap=0.1,
        pesr_breakend_window=300,
        depth_sample_overlap=0.5,
        depth_interval_overlap=0.5,
        depth_breakend_window=500000,
        mixed_sample_overlap=0.5,
        mixed_interval_overlap=0.5,
        mixed_breakend_window=1000000,
        reference_fasta=reference_fasta,
        reference_fasta_fai=reference_fasta_fai,
        reference_dict=reference_dict,
        java_mem_fraction=java_mem_fraction,
        variant_prefix="~{cohort_name}_~{contig}_",
        gatk_docker=gatk_docker,
        runtime_attr_override=runtime_attr_svcluster
    }

    #Subset bothside_pass & background_fail to chromosome of interest
    call SubsetVariantList as SubsetBothsidePass {
      input:
        vid_list=CombineSRBothsidePass.out,
        vid_col=2,
        vcf=SVCluster.out,
        outfile_name="~{cohort_name}.combine_batches.sr_bothside_pass.~{contig}.subset.list",
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_attr_override=runtime_override_subset_bothside_pass
    }
    call SubsetVariantList as SubsetBackgroundFail {
      input:
        vid_list=CombineBackgroundFail.outfile,
        vid_col=1,
        vcf=SVCluster.out,
        outfile_name="~{cohort_name}.combine_batches.sr_background_fail.~{contig}.subset.list",
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_attr_override=runtime_override_subset_background_fail
    }

    #Update SR background fail & bothside pass files (1)
    call MiniTasks.UpdateSrList as UpdateBackgroundFail {
      input:
        vcf=SVCluster.out,
        original_list=SubsetBothsidePass.filtered_vid_list,
        outfile="~{cohort_name}.combine_batches.sr_bothside_pass.~{contig}.list",
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_override_update_sr_list
    }
    call MiniTasks.UpdateSrList as UpdateBothsidePass {
      input:
        vcf=SVCluster.out,
        original_list=SubsetBackgroundFail.filtered_vid_list,
        outfile="~{cohort_name}.combine_batches.sr_background_fail.~{contig}.list",
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_override_update_sr_list
    }
  }

  # Merge resolved vcfs for QC
  if (merge_vcfs) {
    call MiniTasks.ConcatVcfs {
      input:
        vcfs=SVCluster.out,
        vcfs_idx=SVCluster.out_index,
        naive=true,
        outfile_prefix="~{cohort_name}.combine_batches",
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_attr_override=runtime_override_concat
    }
  }

  #Final outputs
  output {
    Array[File] combined_vcfs = SVCluster.out
    Array[File] combined_vcf_indexes = SVCluster.out_index
    Array[File] cluster_bothside_pass_lists = UpdateBothsidePass.updated_list
    Array[File] cluster_background_fail_lists = UpdateBackgroundFail.updated_list
    File? combine_batches_merged_vcf = ConcatVcfs.concat_vcf
    File? combine_batches_merged_vcf_index = ConcatVcfs.concat_vcf_idx
  }
}


# Find intersection of Variant IDs from vid_list with those present in vcf, return as filtered_vid_list
task SubsetVariantList {
  input {
    File vid_list
    Int vid_col
    File vcf
    String outfile_name
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  # when filtering/sorting/etc, memory usage will likely go up (much of the data will have to
  # be held in memory or disk while working, potentially in a form that takes up more space)
  RuntimeAttr runtime_default = object {
                                  mem_gb: 3.75,
                                  disk_gb: ceil(10.0 + size(vid_list, "GB") * 2.0 + size(vcf, "GB")),
                                  cpu_cores: 1,
                                  preemptible_tries: 3,
                                  max_retries: 1,
                                  boot_disk_gb: 10
                                }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: select_first([runtime_override.mem_gb, runtime_default.mem_gb]) + " GB"
    disks: "local-disk " + select_first([runtime_override.disk_gb, runtime_default.disk_gb]) + " HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_base_mini_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euo pipefail
    bcftools query -f '[%MEMBERS\n]' ~{vcf} | tr ',' '\n' | sort -u > valid_vids.list
    awk -F'\t' -v OFS='\t' 'ARGIND==1{inFileA[$1]; next} {if ($~{vid_col} in inFileA) print }' valid_vids.list ~{vid_list} \
      > ~{outfile_name}
  >>>

  output {
    File filtered_vid_list = outfile_name
  }
}
