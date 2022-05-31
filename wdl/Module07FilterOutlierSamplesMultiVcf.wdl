##########################
## EXPERIMENTAL WORKFLOW
##########################

version 1.0 

import "Structs.wdl"
import "Module07FilterOutlierSamples.wdl" as FilterSingle

# This is an analysis WDL that wraps Module07FilterOutlierSamples.wdl and 
# applies it in parallel across multiple input VCFs.


workflow FilterOutlierSamplesPostHocMultiVcf {
  input {
    Array[File] vcfs
    Array[File] vcf_idxs
    File? pcrplus_samples_list
    Int? n_iqr_cutoff_pcrplus
    Int n_iqr_cutoff_pcrminus
    Int records_per_shard
    String prefix
    File autosomes_fai

    String sv_pipeline_docker
    String sv_pipeline_docker_CombineCounts #note: can remove this input once debugged
    String sv_base_mini_docker

    RuntimeAttr? runtime_overide_shard_vcf
  }

  Array[Pair[File, File]] vcf_pairs = zip(vcfs, vcf_idxs)
  Boolean PCRPLUS = defined(pcrplus_samples_list)

  # Get count of biallelic autosomal variants per sample for each VCF
  scatter ( vcf_info in vcf_pairs ) {
    call FilterSingle.FilterOutlierSamplesPostHoc as CollectData {
      input:
        vcf=vcf_info.left,
        vcf_idx=vcf_info.right,
        pcrplus_samples_list=pcrplus_samples_list,
        records_per_shard=records_per_shard,
        prefix=prefix,
        autosomes_fai=autosomes_fai,
        collect_data_only=true,
        sv_pipeline_docker=sv_pipeline_docker,
        sv_pipeline_docker_CombineCounts=sv_pipeline_docker_CombineCounts,
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_overide_shard_vcf=runtime_overide_shard_vcf
    }
  }
  call FilterSingle.CombineCounts as Combine {
    input:
      svcounts=CollectData.svcounts_per_sample_data,
      prefix=prefix,
      sv_pipeline_docker=sv_pipeline_docker
  }
  File all_samples_list = CollectData.all_samples_list[0]
  File plus_samples_list = CollectData.plus_samples_list[0]
  File minus_samples_list = CollectData.minus_samples_list[0]

  # Get outliers
  if (PCRPLUS) {
    call FilterSingle.IdentifyOutliers as IdentifyPlusOutliers {
      input:
        svcounts=Combine.summed_svcounts,
        n_iqr_cutoff=select_first([n_iqr_cutoff_pcrplus]),
        samples_list=plus_samples_list,
        prefix="~{prefix}.PCRPLUS",
        sv_pipeline_docker=sv_pipeline_docker
    }
  }
  call FilterSingle.IdentifyOutliers as IdentifyMinusOutliers {
    input:
      svcounts=Combine.summed_svcounts,
      n_iqr_cutoff=n_iqr_cutoff_pcrminus,
      samples_list=minus_samples_list,
      prefix="~{prefix}.PCRMINUS",
      sv_pipeline_docker=sv_pipeline_docker
  }

  # Exclude outliers from vcfs
  scatter (vcf_info in vcf_pairs) {
    call FilterSingle.ExcludeOutliers as Exclude {
        input:
          vcf=vcf_info.left,
          vcf_idx=vcf_info.right,
          plus_outliers_list=IdentifyPlusOutliers.outliers_list,
          minus_outliers_list=IdentifyMinusOutliers.outliers_list,
          outfile=basename(vcf_info.left, ".vcf.gz") + ".outliers_removed.vcf.gz",
          prefix=prefix,
          sv_pipeline_docker=sv_pipeline_docker
      }
  }

  # Write new list of samples without outliers
  call FilterSingle.FilterSampleList as FilterList {
    input:
      original_samples_list=all_samples_list,
      outlier_samples=Exclude.merged_outliers_list[0],
      prefix=prefix,
      sv_pipeline_docker=sv_pipeline_docker
  }

  # Final outputs
  output {
    Array[File] vcfs_noOutliers = Exclude.vcf_no_outliers
    Array[File] vcf_noOutliers_idx = Exclude.vcf_no_outliers_idx
    File nooutliers_samples_list = FilterList.filtered_samples_list
    File excluded_samples_list = Exclude.merged_outliers_list[0]
    File svcounts_per_sample_data = Combine.summed_svcounts
    File? svcounts_per_sample_plots_PCRPLUS = IdentifyPlusOutliers.svcount_distrib_plots
    File svcounts_per_sample_plots_PCRMINUS = IdentifyMinusOutliers.svcount_distrib_plots
  }
}
