version 1.0

import "Structs.wdl"
import "CalcAF.wdl" as calcAF
import "MinGQRocOpt.wdl" as roc_opt_sub
import "TasksMakeCohortVcf.wdl" as MiniTasks
import "Module07MinGQTasks.wdl" as minGQTasks
import "ReviseSVtypeINStoMEI.wdl" as ReviseSVtype

workflow Module07FilterGTsPart2 {
  input {
    Array[File] PCRMINUS_trio_tarballs
    File PCRMINUS_cleaned_trios_famfile

    Array[File]? PCRPLUS_trio_tarballs
    File? PCRPLUS_cleaned_trios_famfile

    Int roc_shards
    String optimize_minSizes
    String optimize_maxSizes
    String optimize_minFreqs
    String optimize_maxFreqs
    String optimize_includeSVTYPEs
    String optimize_includeFILTERs
    String optimize_excludeFILTERs
    String optimize_includeEV
    String optimize_excludeEV
    Int optimize_maxTrainingTrios = 1000
    Int optimize_maxSVperTrio = 1000
    Int min_sv_per_proband_per_condition = 10
    Float roc_max_fdr_PCRMINUS
    Float roc_max_fdr_PCRPLUS
    String optimize_metric = "GQ"
    Int roc_min_metric
    Int roc_max_metric
    Float roc_step_metric
    String prefix

    String sv_pipeline_docker
    String sv_pipeline_base_docker
    String? sv_pipeline_base_docker_buildTree
    String sv_base_mini_docker

    RuntimeAttr? runtime_attr_EnumerateConditions
    RuntimeAttr? runtime_attr_ConcatTarball
    RuntimeAttr? runtime_attr_SubsetTrioTarball
    RuntimeAttr? runtime_attr_roc_single
  }

  # Get table of all conditions to evaluate
  call minGQTasks.EnumerateConditions {
    input:
      prefix=prefix,
      condition_shards=roc_shards,
      optimize_minSizes=optimize_minSizes,
      optimize_maxSizes=optimize_maxSizes,
      optimize_minFreqs=optimize_minFreqs,
      optimize_maxFreqs=optimize_maxFreqs,
      optimize_includeSVTYPEs=optimize_includeSVTYPEs,
      optimize_includeFILTERs=optimize_includeFILTERs,
      optimize_excludeFILTERs=optimize_excludeFILTERs,
      optimize_includeEV=optimize_includeEV,
      optimize_excludeEV=optimize_excludeEV,
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_attr_EnumerateConditions
  }

  # Concatenate training data across multiple input shards (e.g., chromosomes), if necessary
  if (length(PCRMINUS_trio_tarballs) > 1){
    call minGQTasks.ConcatTarball as ConcatTarballPCRMINUS {
      input: 
        tarballs = PCRMINUS_trio_tarballs,
        sv_base_mini_docker = sv_base_mini_docker,
        runtime_attr_override = runtime_attr_ConcatTarball
    }
  }
  File PCRMINUS_merged_tarball = select_first([ConcatTarballPCRMINUS.tarball, 
                                               PCRMINUS_trio_tarballs[0]])
  if (defined(PCRPLUS_cleaned_trios_famfile)){
    call minGQTasks.ConcatTarball as ConcatTarballPCRPLUS {
      input: 
        tarballs = select_first([PCRPLUS_trio_tarballs, []]),
        sv_base_mini_docker = sv_base_mini_docker,
        runtime_attr_override = runtime_attr_ConcatTarball
    }
    File PCRPLUS_merged_tarball = select_first([ConcatTarballPCRPLUS.tarball, 
                                                select_first([PCRPLUS_trio_tarballs])[0]])
  }

  # Subset trios for training, if optioned
  if(defined(optimize_maxTrainingTrios)) { 
    call minGQTasks.SubsetTrioTarball as SubsetTriosPCRMINUS {
      input:
        tarball_in=PCRMINUS_merged_tarball,
        max_trios=select_first([optimize_maxTrainingTrios, 1000000]),
        prefix=basename(PCRMINUS_merged_tarball, ".tar.gz") + "subsetted",
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_attr_override=runtime_attr_SubsetTrioTarball
    }
    if (defined(PCRPLUS_cleaned_trios_famfile)){
      call minGQTasks.SubsetTrioTarball as SubsetTriosPCRPLUS {
        input:
          tarball_in=select_first([PCRPLUS_merged_tarball]),
          max_trios=select_first([optimize_maxTrainingTrios, 1000000]),
          prefix=basename(select_first([PCRPLUS_merged_tarball]), ".tar.gz") + "subsetted",
          sv_base_mini_docker=sv_base_mini_docker,
          runtime_attr_override=runtime_attr_SubsetTrioTarball
      }
    }
  }
  File PCRMINUS_training_tarball = select_first([SubsetTriosPCRMINUS.subsetted_tarball, 
                                                 PCRMINUS_merged_tarball])
  if (defined(PCRPLUS_cleaned_trios_famfile)) {
    File PCRPLUS_training_tarball = select_first([SubsetTriosPCRPLUS.subsetted_tarball, 
                                                   PCRPLUS_merged_tarball])
  }

  # Train PCR- filtering model
  scatter ( shard in EnumerateConditions.conditions_table_noHeader_shards ) {
    call roc_opt_sub.MinGQRocOpt as roc_opt_PCRMINUS {
      input:
        trio_tarball=PCRMINUS_training_tarball,
        prefix="~{prefix}.PCRMINUS",
        trios_list=PCRMINUS_cleaned_trios_famfile,
        conditions_table=shard,
        maxSVperTrio=optimize_maxSVperTrio,
        roc_max_fdr=roc_max_fdr_PCRMINUS,
        optimize_metric=optimize_metric,
        roc_min_metric=roc_min_metric,
        roc_max_metric=roc_max_metric,
        roc_step_metric=roc_step_metric,
        min_sv_per_proband_per_condition=min_sv_per_proband_per_condition,
        sv_base_mini_docker=sv_base_mini_docker,
        sv_pipeline_base_docker=sv_pipeline_base_docker,
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_roc_single=runtime_attr_roc_single
    }
  }
  call minGQTasks.CombineRocOptResults as combine_roc_optimal_PCRMINUS {
    input:
      shards=roc_opt_PCRMINUS.roc_optimal_merged,
      outfile="~{prefix}.PCRMINUS.min~{optimize_metric}_condition_opts.txt",
      sv_base_mini_docker=sv_base_mini_docker
  }
  call minGQTasks.CombineRocOptResults as combine_roc_stats_PCRMINUS {
    input:
      shards=roc_opt_PCRMINUS.distrib_stats_merged,
      outfile="~{prefix}.min~{optimize_metric}_condition_distrib_stats.txt",
      sv_base_mini_docker=sv_base_mini_docker
  }
  call minGQTasks.BuildFilterTree as build_tree_PCRMINUS {
    input:
      conditions_table=EnumerateConditions.conditions_table,
      condition_optimizations=combine_roc_optimal_PCRMINUS.merged_file,
      condition_distrib_stats=combine_roc_stats_PCRMINUS.merged_file,
      prefix="~{prefix}.PCRMINUS",
      optimize_metric=optimize_metric,
      sv_pipeline_base_docker=select_first([sv_pipeline_base_docker_buildTree, sv_pipeline_base_docker])
  }

  # Train PCR+ filtering model, if optioned
  if (defined(PCRPLUS_cleaned_trios_famfile)) {
    scatter ( shard in EnumerateConditions.conditions_table_noHeader_shards ) {
      call roc_opt_sub.MinGQRocOpt as roc_opt_PCRPLUS {
        input:
          trio_tarball=PCRPLUS_training_tarball,
          prefix="~{prefix}.PCRPLUS",
          trios_list=PCRPLUS_cleaned_trios_famfile,
          conditions_table=shard,
          maxSVperTrio=optimize_maxSVperTrio,
          roc_max_fdr=roc_max_fdr_PCRPLUS,
          optimize_metric=optimize_metric,
          roc_min_metric=roc_min_metric,
          roc_max_metric=roc_max_metric,
          roc_step_metric=roc_step_metric,
          min_sv_per_proband_per_condition=min_sv_per_proband_per_condition,
          sv_base_mini_docker=sv_base_mini_docker,
          sv_pipeline_base_docker=sv_pipeline_base_docker,
          sv_pipeline_docker=sv_pipeline_docker,
          runtime_attr_roc_single=runtime_attr_roc_single
      }
    }
    call minGQTasks.CombineRocOptResults as combine_roc_optimal_PCRPLUS {
      input:
        shards=roc_opt_PCRPLUS.roc_optimal_merged,
        outfile="~{prefix}.PCRPLUS.min~{optimize_metric}_condition_opts.txt",
        sv_base_mini_docker=sv_base_mini_docker
    }
    call minGQTasks.CombineRocOptResults as combine_roc_stats_PCRPLUS {
      input:
        shards=roc_opt_PCRPLUS.distrib_stats_merged,
        outfile="~{prefix}.min~{optimize_metric}_condition_distrib_stats.txt",
        sv_base_mini_docker=sv_base_mini_docker
    }
    call minGQTasks.BuildFilterTree as build_tree_PCRPLUS {
      input:
        conditions_table=EnumerateConditions.conditions_table,
        condition_optimizations=combine_roc_optimal_PCRPLUS.merged_file,
        condition_distrib_stats=combine_roc_stats_PCRPLUS.merged_file,
        prefix="~{prefix}.PCRPLUS",
        optimize_metric=optimize_metric,
        sv_pipeline_base_docker=select_first([sv_pipeline_base_docker_buildTree, sv_pipeline_base_docker])
    }
  }


  output {
    File PCRMINUS_lookup_table = build_tree_PCRMINUS.filter_lookup_table
    File? PCRPLUS_lookup_table = build_tree_PCRPLUS.filter_lookup_table
  }

}

