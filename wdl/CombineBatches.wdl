version 1.0

import "CombineSRBothsidePass.wdl" as CombineSRBothsidePassWorkflow
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

    # Reclustering parameters
    File clustering_config_part1
    File stratification_config_part1
    File clustering_config_part2
    File stratification_config_part2
    # These arrays give the names and intervals for reference contexts for stratification (same lengths)
    # Names must correspond to those in the stratification config files
    Array[String] context_names
    Array[File] context_bed_files

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

    RuntimeAttr? runtime_override_update_sr_list
    RuntimeAttr? runtime_attr_create_ploidy
    RuntimeAttr? runtime_attr_reformat_1
    RuntimeAttr? runtime_attr_reformat_2
    RuntimeAttr? runtime_attr_join_vcfs
    RuntimeAttr? runtime_attr_cluster_sites
    RuntimeAttr? runtime_attr_recluster_part1
    RuntimeAttr? runtime_attr_recluster_part2
    RuntimeAttr? runtime_override_subset_bothside_pass
    RuntimeAttr? runtime_override_subset_background_fail
    RuntimeAttr? runtime_attr_get_non_ref_vids
    RuntimeAttr? runtime_attr_calculate_support_frac
    RuntimeAttr? runtime_override_clean_background_fail
    RuntimeAttr? runtime_override_concat
  }

  # Preprocess some inputs
  call CombineSRBothsidePassWorkflow.CombineSRBothsidePass {
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
          output_prefix = basename(vcf, ".vcf.gz") + ".reformat_legacy",
          sv_pipeline_docker = sv_pipeline_docker,
          runtime_attr_override = runtime_attr_reformat_1
      }
    }
    File reformatted_vcf = select_first([ReformatGenotypedVcf.out, vcf])
    call GatkFormatting.FormatVcf {
      input:
        vcf=reformatted_vcf,
        ploidy_table=CreatePloidyTableFromPed.out,
        args="--fix-end",
        output_prefix=basename(vcf, ".vcf.gz") + ".reformat_gatk",
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_attr_reformat_2
    }
  }

  #Scatter per chromosome
  Array[String] contigs = transpose(read_tsv(contig_list))[0]
  scatter ( contig in contigs ) {

    # First round of clustering
    call ClusterTasks.SVCluster as JoinVcfs {
      input:
        vcfs=FormatVcf.out,
        ploidy_table=CreatePloidyTableFromPed.out,
        output_prefix="~{cohort_name}.combine_batches.~{contig}.join_vcfs",
        contig=contig,
        fast_mode=false,
        pesr_sample_overlap=0,
        pesr_interval_overlap=1,
        pesr_breakend_window=0,
        depth_sample_overlap=0,
        depth_interval_overlap=1,
        depth_breakend_window=0,
        mixed_sample_overlap=0,
        mixed_interval_overlap=1,
        mixed_breakend_window=0,
        reference_fasta=reference_fasta,
        reference_fasta_fai=reference_fasta_fai,
        reference_dict=reference_dict,
        java_mem_fraction=java_mem_fraction,
        gatk_docker=gatk_docker,
        runtime_attr_override=runtime_attr_join_vcfs
    }

    # First round of clustering
    call ClusterTasks.SVCluster as ClusterSites {
      input:
        vcfs=[JoinVcfs.out],
        ploidy_table=CreatePloidyTableFromPed.out,
        output_prefix="~{cohort_name}.combine_batches.~{contig}.cluster_sites",
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
        runtime_attr_override=runtime_attr_cluster_sites
    }

    #Subset bothside_pass & background_fail to chromosome of interest
    call SubsetVariantList as SubsetBothsidePass {
      input:
        vid_list=CombineSRBothsidePass.out,
        vid_col=2,
        vcf=ClusterSites.out,
        outfile_name="~{cohort_name}.combine_batches.~{contig}.sr_bothside_pass.subset.list",
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_attr_override=runtime_override_subset_bothside_pass
    }
    call SubsetVariantList as SubsetBackgroundFail {
      input:
        vid_list=CombineBackgroundFail.outfile,
        vid_col=1,
        vcf=ClusterSites.out,
        outfile_name="~{cohort_name}.combine_batches.~{contig}.sr_background_fail.subset.list",
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_attr_override=runtime_override_subset_background_fail
    }

    #Update SR background fail & bothside pass files (1)
    call MiniTasks.UpdateSrList as UpdateBackgroundFail1 {
      input:
        vcf=ClusterSites.out,
        original_list=SubsetBothsidePass.filtered_vid_list,
        outfile="~{cohort_name}.combine_batches.~{contig}.sr_bothside_pass1.list",
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_override_update_sr_list
    }
    call MiniTasks.UpdateSrList as UpdateBothsidePass1 {
      input:
        vcf=ClusterSites.out,
        original_list=SubsetBackgroundFail.filtered_vid_list,
        outfile="~{cohort_name}.combine_batches.~{contig}.sr_background_fail1.list",
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_override_update_sr_list
    }

    # Second round of clustering
    call GroupedSVClusterTask as GroupedSVClusterPart1 {
      input:
        vcf=ClusterSites.out,
        ploidy_table=CreatePloidyTableFromPed.out,
        output_prefix="~{cohort_name}.combine_batches.~{contig}.recluster_part_1",
        reference_fasta=reference_fasta,
        reference_fasta_fai=reference_fasta_fai,
        reference_dict=reference_dict,
        clustering_config=clustering_config_part1,
        stratification_config=stratification_config_part1,
        context_bed_files=context_bed_files,
        context_names=context_names,
        java_mem_fraction=java_mem_fraction,
        gatk_docker=gatk_docker,
        runtime_attr_override=runtime_attr_recluster_part1
    }

    #Update SR background fail & bothside pass files (2)
    call MiniTasks.UpdateSrList as UpdateBackgroundFail2 {
      input:
        vcf=GroupedSVClusterPart1.out,
        original_list=UpdateBackgroundFail1.updated_list,
        outfile="~{cohort_name}.combine_batches.~{contig}.sr_bothside_pass2.list",
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_override_update_sr_list
    }
    call MiniTasks.UpdateSrList as UpdateBothsidePass2 {
      input:
        vcf=GroupedSVClusterPart1.out,
        original_list=UpdateBothsidePass1.updated_list,
        outfile="~{cohort_name}.combine_batches.~{contig}.sr_background_fail2.list",
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_override_update_sr_list
    }

    # Third round of clustering
    call GroupedSVClusterTask as GroupedSVClusterPart2 {
      input:
        vcf=GroupedSVClusterPart1.out,
        ploidy_table=CreatePloidyTableFromPed.out,
        output_prefix="~{cohort_name}.combine_batches.~{contig}.recluster_part_2",
        reference_fasta=reference_fasta,
        reference_fasta_fai=reference_fasta_fai,
        reference_dict=reference_dict,
        clustering_config=clustering_config_part2,
        stratification_config=stratification_config_part2,
        context_bed_files=context_bed_files,
        context_names=context_names,
        java_mem_fraction=java_mem_fraction,
        gatk_docker=gatk_docker,
        runtime_attr_override=runtime_attr_recluster_part2
    }

    #Update SR background fail & bothside pass files (3)
    call MiniTasks.UpdateSrList as UpdateBackgroundFail3 {
      input:
        vcf=GroupedSVClusterPart2.out,
        original_list=UpdateBackgroundFail2.updated_list,
        outfile="~{cohort_name}.combine_batches.~{contig}.sr_bothside_pass3.list",
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_override_update_sr_list
    }
    call MiniTasks.UpdateSrList as UpdateBothsidePass3 {
      input:
        vcf=GroupedSVClusterPart2.out,
        original_list=UpdateBothsidePass2.updated_list,
        outfile="~{cohort_name}.combine_batches.~{contig}.sr_background_fail3.list",
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_override_update_sr_list
    }
  }

  # Merge resolved vcfs for QC
  if (merge_vcfs) {
    call MiniTasks.ConcatVcfs {
      input:
        vcfs=GroupedSVClusterPart2.out,
        vcfs_idx=GroupedSVClusterPart2.out_index,
        naive=true,
        outfile_prefix="~{cohort_name}.combine_batches.concat_all_contigs",
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_attr_override=runtime_override_concat
    }
  }

  #Final outputs
  output {
    Array[File] combined_vcfs = GroupedSVClusterPart2.out
    Array[File] combined_vcf_indexes = GroupedSVClusterPart2.out_index
    Array[File] cluster_background_fail_lists = UpdateBackgroundFail3.updated_list
    Array[File] cluster_bothside_pass_lists = UpdateBothsidePass3.updated_list
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

task GroupedSVClusterTask {
  input {
    File vcf
    File ploidy_table
    String output_prefix

    File reference_fasta
    File reference_fasta_fai
    File reference_dict

    File clustering_config
    File stratification_config
    Array[File] context_bed_files
    Array[String] context_names

    Float context_overlap_frac = 0
    Int context_num_breakpoint_overlaps = 1
    Int context_interchrom_num_breakpoint_overlaps = 1

    String? contig
    String? additional_args

    Float? java_mem_fraction
    String gatk_docker
    RuntimeAttr? runtime_attr_override
  }

  parameter_meta {
    vcf: {
           localization_optional: true
         }
  }

  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 3.75,
                               disk_gb: ceil(10 + size(vcf, "GB") * 2),
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File out = "~{output_prefix}.vcf.gz"
    File out_index = "~{output_prefix}.vcf.gz.tbi"
  }
  command <<<
    set -euo pipefail

    function getJavaMem() {
    # get JVM memory in MiB by getting total memory from /proc/meminfo
    # and multiplying by java_mem_fraction
    cat /proc/meminfo \
      | awk -v MEM_FIELD="$1" '{
        f[substr($1, 1, length($1)-1)] = $2
      } END {
        printf "%dM", f[MEM_FIELD] * ~{default="0.85" java_mem_fraction} / 1024
      }'
    }
    JVM_MAX_MEM=$(getJavaMem MemTotal)
    echo "JVM memory: $JVM_MAX_MEM"

    gatk --java-options "-Xmx${JVM_MAX_MEM}" GroupedSVCluster \
      ~{"-L " + contig} \
      --reference ~{reference_fasta} \
      --ploidy-table ~{ploidy_table} \
      -V ~{vcf} \
      -O ~{output_prefix}.vcf.gz \
      --clustering-config ~{clustering_config} \
      --stratify-config ~{stratification_config} \
      --context-intervals ~{sep=" --context-intervals " context_bed_files} \
      --context-name ~{sep=" --context-name " context_names} \
      --stratify-overlap-fraction ~{context_overlap_frac} \
      --stratify-num-breakpoint-overlaps ~{context_num_breakpoint_overlaps} \
      --stratify-num-breakpoint-overlaps-interchromosomal ~{context_interchrom_num_breakpoint_overlaps} \
      ~{additional_args}
  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: gatk_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}