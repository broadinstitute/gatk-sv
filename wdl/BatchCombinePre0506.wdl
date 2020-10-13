version 1.0

import "SelectBatchesAndColumns.wdl" as select

workflow BatchCombinePre0506 {
  input {
    File sample_sets_table_tsv
    Array[String]? batches
    String? linux_docker
    RuntimeAttr? runtime_override_batchcombine
  }
  Array[String] columns = ["trained_genotype_depth_depth_sepcutoff", "regenotyped_depth_vcfs", "merged_bincov", "merged_PE","median_cov", "ped_file_postOutlierExclusion", "genotyped_pesr_vcf", "sr_background_fail", "sr_bothside_pass", "cutoffs"]

  call select.SelectBatchesAndColumns {
    input:
      sample_sets_table_tsv = sample_sets_table_tsv,
      columns = columns,
      batches_opt = batches,
      linux_docker = linux_docker,
      runtime_attr_override = runtime_override_batchcombine
  }

  output {
    Array[String]? batches_combined_pre_0506 = batches
    Array[File] trained_genotype_depth_depth_sepcutoff = SelectBatchesAndColumns.columns_as_rows[0]
    Array[File] regenotyped_depth_vcfs = SelectBatchesAndColumns.columns_as_rows[1]
    Array[File] merged_bincov = SelectBatchesAndColumns.columns_as_rows[2]
    Array[File] merged_PE = SelectBatchesAndColumns.columns_as_rows[3]
    Array[File] median_cov = SelectBatchesAndColumns.columns_as_rows[4]
    Array[File] ped_file_postOutlierExclusion = SelectBatchesAndColumns.columns_as_rows[5]
    Array[File] genotyped_pesr_vcf = SelectBatchesAndColumns.columns_as_rows[6]
    Array[File] sr_background_fail = SelectBatchesAndColumns.columns_as_rows[7]
    Array[File] sr_bothside_pass = SelectBatchesAndColumns.columns_as_rows[8]
    Array[File] cutoffs = SelectBatchesAndColumns.columns_as_rows[9]
  }
  
}
