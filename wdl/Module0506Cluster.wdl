version 1.0

import "VcfClusterSingleChromsome.wdl" as VcfClusterContig
import "Tasks0506.wdl" as MiniTasks

workflow Module0506Cluster {
  input {
    String cohort_name
    Array[String] batches

    Boolean merge_vcfs = false

    Array[File] pesr_vcfs
    Array[File] depth_vcfs

    Array[File] raw_sr_bothside_pass_files
    Array[File] raw_sr_background_fail_files

    File contig_list
    Int max_shards_per_chrom
    Int min_variants_per_shard
    File pe_exclude_list
    File depth_exclude_list
    Float min_sr_background_fail_batches

    File empty_file

    String sv_base_mini_docker
    String sv_pipeline_docker

    # overrides for local tasks
    RuntimeAttr? runtime_override_update_sr_list
    RuntimeAttr? runtime_override_merge_pesr_depth

    # overrides for mini tasks
    RuntimeAttr? runtime_override_clean_bothside_pass
    RuntimeAttr? runtime_override_clean_background_fail
    RuntimeAttr? runtime_override_concat
    RuntimeAttr? runtime_override_sort_pesr_depth_merged_vcf
    RuntimeAttr? runtime_override_concat_pesr_depth

    # overrides for VcfClusterContig
    RuntimeAttr? runtime_override_localize_vcfs
    RuntimeAttr? runtime_override_join_vcfs
    RuntimeAttr? runtime_override_fix_multiallelic
    RuntimeAttr? runtime_override_fix_ev_tags
    RuntimeAttr? runtime_override_subset_bothside_pass
    RuntimeAttr? runtime_override_subset_background_fail
    RuntimeAttr? runtime_override_subset_sv_type
    RuntimeAttr? runtime_override_concat_sv_types
    RuntimeAttr? runtime_override_shard_vcf_precluster
    RuntimeAttr? runtime_override_svtk_vcf_cluster
    RuntimeAttr? runtime_override_get_vcf_header_with_members_info_line
    RuntimeAttr? runtime_override_concat_shards
  }

  # Preprocess some inputs
  Int num_pass_lines=length(raw_sr_bothside_pass_files)
  call MiniTasks.CatUncompressedFiles as CleanBothsidePass {
    input:
      shards=raw_sr_bothside_pass_files,
      filter_command="sort -m | uniq -c | awk -v OFS='\\t' '{print $1/~{num_pass_lines}, $2}'",
      outfile_name="cohort_sr_genotyping_bothside_pass_list.txt",
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_override_clean_bothside_pass
  }
  File sr_bothend_pass = CleanBothsidePass.outfile

  Float min_background_fail_first_col = min_sr_background_fail_batches * length(raw_sr_background_fail_files)
  call MiniTasks.CatUncompressedFiles as CleanBackgroundFail {
    input:
      shards=raw_sr_background_fail_files,
      filter_command="sort -m | uniq -c | awk -v OFS='\\t' '{if($1 >= ~{min_background_fail_first_col}) print $2}'",
      outfile_name="cohort_sr_genotyping_background_fail_list.txt",
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_override_clean_background_fail
  }
  File sr_background_fail = CleanBackgroundFail.outfile

  #Scatter per chromosome
  Array[String] contigs = transpose(read_tsv(contig_list))[0]
  scatter ( contig in contigs ) {

    #Subset PESR VCFs to single chromosome & cluster
    #Note: also subsets bothside_pass and background_fail files to variants
    #present on chromosome of interest
    call VcfClusterContig.VcfClusterSingleChrom as ClusterPesr {
      input:
        vcfs=pesr_vcfs,
        batches=batches,
        prefix="~{cohort_name}.pesr",
        dist=300,
        frac=0.1,
        sample_overlap=0.5,
        exclude_list=pe_exclude_list,
        sv_size=50,
        sv_types=["DEL","DUP","INV","BND","INS"],
        contig=contig,
        max_shards_per_chrom_svtype=100,
        min_variants_per_shard_per_chrom_svtype=100,
        subset_sr_lists=true,
        bothside_pass=sr_bothend_pass,
        background_fail=sr_background_fail,
        empty_file=empty_file,
        sv_pipeline_docker=sv_pipeline_docker,
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_override_localize_vcfs = runtime_override_localize_vcfs,
        runtime_override_join_vcfs = runtime_override_join_vcfs,
        runtime_override_fix_multiallelic = runtime_override_fix_multiallelic,
        runtime_override_fix_ev_tags = runtime_override_fix_ev_tags,
        runtime_override_subset_bothside_pass=runtime_override_subset_bothside_pass,
        runtime_override_subset_background_fail=runtime_override_subset_background_fail,
        runtime_override_subset_sv_type=runtime_override_subset_sv_type,
        runtime_override_concat_sv_types=runtime_override_concat_sv_types,
        runtime_override_shard_vcf_precluster=runtime_override_shard_vcf_precluster,
        runtime_override_svtk_vcf_cluster=runtime_override_svtk_vcf_cluster,
        runtime_override_get_vcf_header_with_members_info_line=runtime_override_get_vcf_header_with_members_info_line,
        runtime_override_concat_shards=runtime_override_concat_shards
    }

    #Subset RD VCFs to single chromosome & cluster
    call VcfClusterContig.VcfClusterSingleChrom as ClusterDepth {
      input:
        vcfs=depth_vcfs,
        batches=batches,
        prefix="~{cohort_name}.depth",
        dist=500000,
        frac=0.5,
        sample_overlap=0.5,
        exclude_list=depth_exclude_list,
        sv_size=5000,
        sv_types=["DEL","DUP"],
        contig=contig,
        max_shards_per_chrom_svtype=100,
        min_variants_per_shard_per_chrom_svtype=100,
        subset_sr_lists=false,
        bothside_pass=sr_bothend_pass,
        background_fail=sr_background_fail,
        empty_file=empty_file,
        sv_pipeline_docker=sv_pipeline_docker,
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_override_localize_vcfs = runtime_override_localize_vcfs,
        runtime_override_join_vcfs = runtime_override_join_vcfs,
        runtime_override_fix_multiallelic = runtime_override_fix_multiallelic,
        runtime_override_fix_ev_tags = runtime_override_fix_ev_tags,
        runtime_override_subset_bothside_pass=runtime_override_subset_bothside_pass,
        runtime_override_subset_background_fail=runtime_override_subset_background_fail,
        runtime_override_subset_sv_type=runtime_override_subset_sv_type,
        runtime_override_concat_sv_types=runtime_override_concat_sv_types,
        runtime_override_shard_vcf_precluster=runtime_override_shard_vcf_precluster,
        runtime_override_svtk_vcf_cluster=runtime_override_svtk_vcf_cluster,
        runtime_override_get_vcf_header_with_members_info_line=runtime_override_get_vcf_header_with_members_info_line,
        runtime_override_concat_shards=runtime_override_concat_shards
    }

    #Update SR background fail & bothside pass files (1)
    call MiniTasks.UpdateSrList as UpdateBackgroundFailFirst {
      input:
        vcf=ClusterPesr.clustered_vcf,
        original_list=ClusterPesr.filtered_background_fail,
        outfile="sr_background_fail.~{contig}.updated.txt",
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_override_update_sr_list
    }
    call MiniTasks.UpdateSrList as UpdateBothsidePassFirst {
      input:
        vcf=ClusterPesr.clustered_vcf,
        original_list=ClusterPesr.filtered_bothside_pass,
        outfile="sr_bothside_pass.~{contig}.updated.txt",
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_override_update_sr_list
    }

    #Merge PESR & RD VCFs
    call MiniTasks.ConcatVcfs as ConcatPesrDepth {
      input:
        vcfs=[ClusterPesr.clustered_vcf, ClusterDepth.clustered_vcf],
        vcfs_idx=[ClusterPesr.clustered_vcf_idx, ClusterDepth.clustered_vcf_idx],
        merge_sort=true,
        outfile_prefix="all_batches.pesr_depth.~{contig}.unmerged",
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_attr_override=runtime_override_concat_pesr_depth
    }
    call MergePesrDepth {
      input:
        vcf=ConcatPesrDepth.concat_vcf,
        vcf_index=ConcatPesrDepth.concat_vcf_idx,
        contig=contig,
        prefix="all_batches.pesr_depth.~{contig}.unsorted",
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_override_merge_pesr_depth
    }
    call MiniTasks.SortVcf as SortMergePesrDepth {
      input:
        vcf = MergePesrDepth.merged_vcf,
        outfile_prefix = "all_batches.pesr_depth.~{contig}",
        sv_base_mini_docker = sv_base_mini_docker,
        runtime_attr_override = runtime_override_sort_pesr_depth_merged_vcf
    }

    #Update SR background fail & bothside pass files (2)
    call MiniTasks.UpdateSrList as UpdateBackgroundFailSecond {
      input:
        vcf=SortMergePesrDepth.out,
        original_list=UpdateBackgroundFailFirst.updated_list,
        outfile="sr_background_fail.~{contig}.updated2.txt",
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_override_update_sr_list
    }
    call MiniTasks.UpdateSrList as UpdateBothsidePassSecond {
      input:
        vcf=SortMergePesrDepth.out,
        original_list=UpdateBothsidePassFirst.updated_list,
        outfile="sr_bothside_pass.~{contig}.updated2.txt",
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_override_update_sr_list
    }
  }

  #Merge resolved vcfs for QC
  if (merge_vcfs) {
    call MiniTasks.ConcatVcfs {
      input:
        vcfs=SortMergePesrDepth.out,
        vcfs_idx=SortMergePesrDepth.out_index,
        merge_sort=true,
        outfile_prefix="~{cohort_name}.0506_clustered",
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_attr_override=runtime_override_concat
    }
  }

  #Final outputs
  output {
    Array[File] vcfs = SortMergePesrDepth.out
    Array[File] vcf_indexes = SortMergePesrDepth.out_index
    Array[File] cluster_bothside_pass_lists = UpdateBothsidePassSecond.updated_list
    Array[File] cluster_background_fail_lists = UpdateBackgroundFailSecond.updated_list
    File? merged_vcf = ConcatVcfs.concat_vcf
    File? merged_vcf_index = ConcatVcfs.concat_vcf_idx
  }
}


#Merge PESR + RD VCFs
task MergePesrDepth {
  input {
    File vcf
    File vcf_index
    String prefix
    String contig
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  String output_file = prefix + ".vcf.gz"

  # when filtering/sorting/etc, memory usage will likely go up (much of the data will have to
  # be held in memory or disk while working, potentially in a form that takes up more space)
  Float input_size = size(vcf, "GiB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 2.0 + 0.3 * input_size,
                                  disk_gb: ceil(10.0 + 0.5 * input_size),
                                  cpu_cores: 1,
                                  preemptible_tries: 0,
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
    /opt/sv-pipeline/04_variant_resolution/scripts/merge_pesr_depth.py \
      --prefix pesr_depth_merged_~{contig} \
      ~{vcf} \
      ~{output_file}
  >>>

  output {
    File merged_vcf = output_file
  }
}
