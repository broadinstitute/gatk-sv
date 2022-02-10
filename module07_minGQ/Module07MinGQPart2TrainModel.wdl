version 1.0

import "Structs.wdl"
import "CalcAF.wdl" as calcAF
import "MinGQRocOpt.wdl" as roc_opt_sub
import "TasksMakeCohortVcf.wdl" as MiniTasks
import "Module07MinGQTasks.wdl" as minGQTasks
import "ReviseSVtypeINStoMEI.wdl" as ReviseSVtype



workflow Module07MinGQPart2 {
    input {
        String prefix
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
        Int optimize_maxSVperTrio
        Float roc_max_fdr_PCRMINUS
        Float roc_max_fdr_PCRPLUS
        Int roc_min_gq
        Int roc_max_gq
        Int roc_step_gq
        Int min_sv_per_proband_per_condition

        Array[File] PCRMINUS_vcf_lists
        Array[File] PCRMINUS_vcf_idx_lists
        Array[File] PCRPLUS_vcf_lists
        Array[File] PCRPLUS_vcf_idx_lists

        File? pcrplus_samples_list

        Array[File] PCRMINUS_trio_tarball
        Array[File] PCRPLUS_trio_tarball
        File PCRMINUS_cleaned_trios_famfile
        File? PCRPLUS_cleaned_trios_famfile

        String sv_pipeline_docker
        String sv_base_mini_docker

        RuntimeAttr? runtime_attr_EnumerateConditions
        RuntimeAttr? runtime_attr_ConcatTarball
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
        sv_pipeline_docker=sv_pipeline_docker
    }

    ### PCRMINUS
    call minGQTasks.ConcatTarball as ConcatTarballPCRMINUS{
        input: 
            tarballs = PCRMINUS_trio_tarball,
            sv_base_mini_docker = sv_base_mini_docker,
            runtime_attr_override = runtime_attr_ConcatTarball
    }

    scatter ( shard in EnumerateConditions.minGQ_conditions_table_noHeader_shards ) {
        ### PCRMINUS
        call roc_opt_sub.MinGQRocOpt as roc_opt_PCRMINUS {
            input:
                trio_tarball=ConcatTarballPCRMINUS.tarball,
                prefix="~{prefix}.PCRMINUS",
                trios_list=PCRMINUS_cleaned_trios_famfile,
                conditions_table=shard,
                maxSVperTrio=optimize_maxSVperTrio,
                roc_max_fdr=roc_max_fdr_PCRMINUS,
                roc_min_gq=roc_min_gq,
                roc_max_gq=roc_max_gq,
                roc_step_gq=roc_step_gq,
                min_sv_per_proband_per_condition=min_sv_per_proband_per_condition,
                sv_base_mini_docker=sv_base_mini_docker,
                sv_pipeline_docker=sv_pipeline_docker
      }
    }

    # Merge ROC results to build minGQ filtering lookup tree
    call minGQTasks.CombineRocOptResults as combine_roc_optimal_PCRMINUS {
      input:
        shards=roc_opt_PCRMINUS.roc_optimal_merged,
        outfile="~{prefix}.PCRMINUS.minGQ_condition_opts.txt",
        sv_base_mini_docker=sv_base_mini_docker
    }

    call minGQTasks.CombineRocOptResults as combine_roc_stats_PCRMINUS {
      input:
        shards=roc_opt_PCRMINUS.distrib_stats_merged,
        outfile="~{prefix}.minGQ_condition_distrib_stats.txt",
        sv_base_mini_docker=sv_base_mini_docker
    }

    # Create final minGQ filtering tree
    call minGQTasks.BuildFilterTree as build_tree_PCRMINUS {
      input:
        conditions_table=EnumerateConditions.minGQ_conditions_table,
        condition_optimizations=combine_roc_optimal_PCRMINUS.merged_file,
        condition_distrib_stats=combine_roc_stats_PCRMINUS.merged_file,
        prefix="~{prefix}.PCRMINUS",
        sv_pipeline_docker=sv_pipeline_docker
    }


    ### PCRPLUS
    if (defined(pcrplus_samples_list)){
        call minGQTasks.ConcatTarball as ConcatTarballPCRPLUS{
            input: 
                tarballs = PCRPLUS_trio_tarball,
                sv_base_mini_docker = sv_base_mini_docker,
                runtime_attr_override = runtime_attr_ConcatTarball
        }

        scatter ( shard in EnumerateConditions.minGQ_conditions_table_noHeader_shards ) {
            call roc_opt_sub.MinGQRocOpt as roc_opt_PCRPLUS {
                input:
                    trio_tarball=ConcatTarballPCRPLUS.tarball,
                    prefix="~{prefix}.PCRPLUS",
                    trios_list=PCRPLUS_cleaned_trios_famfile,
                    conditions_table=shard,
                    maxSVperTrio=optimize_maxSVperTrio,
                    roc_max_fdr=roc_max_fdr_PCRPLUS,
                    roc_min_gq=roc_min_gq,
                    roc_max_gq=roc_max_gq,
                    roc_step_gq=roc_step_gq,
                    min_sv_per_proband_per_condition=min_sv_per_proband_per_condition,
                    sv_base_mini_docker=sv_base_mini_docker,
                    sv_pipeline_docker=sv_pipeline_docker
          }
        }

        # Merge ROC results to build minGQ filtering lookup tree
        call minGQTasks.CombineRocOptResults as combine_roc_optimal_PCRPLUS {
          input:
            shards=roc_opt_PCRPLUS.roc_optimal_merged,
            outfile="~{prefix}.PCRPLUS.minGQ_condition_opts.txt",
            sv_base_mini_docker=sv_base_mini_docker
        }

        call minGQTasks.CombineRocOptResults as combine_roc_stats_PCRPLUS {
          input:
            shards=roc_opt_PCRPLUS.distrib_stats_merged,
            outfile="~{prefix}.minGQ_condition_distrib_stats.txt",
            sv_base_mini_docker=sv_base_mini_docker
        }
     
        # Create final minGQ filtering tree
        call minGQTasks.BuildFilterTree as build_tree_PCRPLUS {
          input:
            conditions_table=EnumerateConditions.minGQ_conditions_table,
            condition_optimizations=combine_roc_optimal_PCRPLUS.merged_file,
            condition_distrib_stats=combine_roc_stats_PCRPLUS.merged_file,
            prefix="~{prefix}.PCRPLUS",
            sv_pipeline_docker=sv_pipeline_docker
        }
    }


    output {
        File PCRMINUS_lookup_table = build_tree_PCRMINUS.filter_lookup_table
        File? PCRPLUS_lookup_table = build_tree_PCRPLUS.filter_lookup_table
    }

}

