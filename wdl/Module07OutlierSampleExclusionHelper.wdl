##########################
## EXPERIMENTAL WORKFLOW
##########################

# Based on : https://portal.firecloud.org/#methods/Talkowski-SV/filter_svcount_outlier_samples/17/wdl

version 1.0 

import "Structs.wdl"
import "TasksMakeCohortVcf.wdl" as Utils
import "Module07FilterOutlierSamples.wdl" as ParentWorkflow

# Helper workflow for Module07FilterOutlierSamplesMultiVcf.wdl 
# to scatter outlier sample exclusion step across an array of sharded VCFs
# prior to concatenating the shards

workflow ScatteredSampleExclusion {
  input {
    Array[File] vcfs
    Array[File] vcf_idxs
    File? plus_outliers_list
    File minus_outliers_list
    String prefix

    String sv_base_mini_docker

    RuntimeAttr? runtime_overide_concat_shards
  }
  Array[Pair[File, File]] vcf_info_pairs = zip(vcfs, vcf_idxs)

  # Exclude outliers from vcf
  scatter ( vcf_info in vcf_info_pairs ) {
    call ParentWorkflow.ExcludeOutliers as ShardExclude {
      input:
        vcf=vcf_info.left,
        vcf_idx=vcf_info.right,
        plus_outliers_list=plus_outliers_list,
        minus_outliers_list=minus_outliers_list,
        outfile=basename(vcf_info.left, "vcf.gz") + ".outliers_removed.vcf.gz",
        prefix=prefix,
        sv_base_mini_docker=sv_base_mini_docker
    }
  }

  # Combine filtered shards
  call Utils.ConcatVcfs as ConcatShards {
    input:
      vcfs=ShardExclude.vcf_no_outliers,
      vcfs_idx=ShardExclude.vcf_no_outliers_idx,
      allow_overlaps=true,
      outfile_prefix="~{prefix}.outliers_removed",
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_overide_concat_shards
  }

  output {
    File filtered_vcf = ConcatShards.concat_vcf
    File filtered_vcf_idx = ConcatShards.concat_vcf_idx
    File excluded_samples = ShardExclude.merged_outliers_list[0]
  }
}
