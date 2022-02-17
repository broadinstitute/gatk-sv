version 1.0

import "Structs.wdl"
import "CalcAF.wdl" as calcAF
import "MinGQRocOpt.wdl" as roc_opt_sub
import "TasksMakeCohortVcf.wdl" as MiniTasks
import "Module07MinGQTasks.wdl" as minGQTasks
import "ReviseSVtypeINStoMEI.wdl" as ReviseSVtype


workflow Module07MinGQPart3 {
    input {
        String prefix
        File PCRMINUS_lookup_table
        File? PCRPLUS_lookup_table
        Array[File] PCRMINUS_vcf_lists
        Array[File] PCRMINUS_vcf_idx_lists
        Array[File]? PCRPLUS_vcf_lists
        Array[File]? PCRPLUS_vcf_idx_lists

        Float max_noCallRate
        Int global_minGQ

        String sv_pipeline_docker
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_CombineVcfs
    }


    # Apply filter per chromosome
    scatter (i in range(length(PCRMINUS_vcf_lists))){
        ### PCRMINUS
        call minGQTasks.ApplyMinGQFilter as apply_filter_PCRMINUS {
          input:
            vcf=PCRMINUS_vcf_lists[i],
            minGQ_lookup_table=PCRMINUS_lookup_table,
            prefix="~{prefix}.PCRMINUS",
            PCR_status="PCRMINUS",
            maxNCR=max_noCallRate,
            global_minGQ=global_minGQ,
            sv_pipeline_docker=sv_pipeline_docker
        }

        if (defined(PCRPLUS_vcf_lists)){
          ### PCRPLUS
          call minGQTasks.ApplyMinGQFilter as apply_filter_PCRPLUS {
            input:
              vcf=select_first([PCRPLUS_vcf_lists])[i],
              minGQ_lookup_table= PCRPLUS_lookup_table,
              prefix="~{prefix}.PCRPLUS",
              PCR_status="PCRPLUS",
              maxNCR=max_noCallRate,
              global_minGQ=global_minGQ,
              sv_pipeline_docker=sv_pipeline_docker
          }

          call minGQTasks.MergePcrVCFs as merge_PCR_VCFs {
            input:
              PCRPLUS_vcf=apply_filter_PCRPLUS.filtered_vcf,
              PCRMINUS_vcf=apply_filter_PCRMINUS.filtered_vcf,
              prefix=prefix, 
              sv_pipeline_docker = sv_pipeline_docker
          }
        }

        File filtered_vcf_shards = select_first([apply_filter_PCRMINUS.filtered_vcf, merge_PCR_VCFs.merged_vcf ])
      }

    call MiniTasks.ConcatVcfs as CombineVcfs {
        input:
          vcfs=filtered_vcf_shards,
          naive=true,
          outfile_prefix=prefix,
          sv_base_mini_docker=sv_base_mini_docker,
          runtime_attr_override=runtime_attr_CombineVcfs
    }

  output {
    File filtered_vcf = CombineVcfs.concat_vcf
    File filtered_vcf_idx = CombineVcfs.concat_vcf_idx
  }
 }