version 1.0

import "SelectBatchesAndColumns.wdl" as select

workflow BatchCombinePreMerge {
  input {
    File sample_sets_table_tsv
    Array[String]? batches
    String? python_docker
    RuntimeAttr? runtime_override_batchcombine
  }
  Array[String] columns = ["filtered_depth_vcf", "filtered_pesr_vcf"]
  
  call select.SelectBatchesAndColumns {
    input:
      sample_sets_table_tsv = sample_sets_table_tsv,
      columns = columns,
      batches_opt = batches,
      python_docker = python_docker,
      runtime_attr_override = runtime_override_batchcombine
  }

  output {
    Array[File] filtered_depth_vcf = SelectBatchesAndColumns.columns_as_rows[0]
    Array[File] filtered_pesr_vcf = SelectBatchesAndColumns.columns_as_rows[1]
    Array[String]? batches_combined_pre_merge = batches
  }
  
}


