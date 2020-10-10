version 1.0

import "SelectBatchesAndColumns.wdl" as select

workflow BatchCombinePreMerge {
  input {
    File sample_sets_table_tsv
    Array[String]? batches
    String? linux_docker
    RuntimeAttr? runtime_override_batchcombine
  }
  Array[String] columns = ["filtered_depth_vcf", "filtered_pesr_vcf"]
  
  scatter (column in columns) {
    call select.SelectBatchesAndColumns {
      input:
        sample_sets_table_tsv = sample_sets_table_tsv,
        column = column,
        batches_opt = batches,
        linux_docker = linux_docker,
        runtime_attr_override = runtime_override_batchcombine
    }
  }
  

  output {
    Array[File] filtered_depth_vcf = SelectBatchesAndColumns.column_as_row[0]
    Array[File] filtered_pesr_vcf = SelectBatchesAndColumns.column_as_row[1]
    Array[String]? batches_combined_pre_merge = batches
  }
  
}


