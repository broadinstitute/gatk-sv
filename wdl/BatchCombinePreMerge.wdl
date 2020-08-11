version 1.0

import "SelectBatchesAndColumns.wdl" as select

workflow BatchCombinePreMerge {
	input {
		File sample_sets_table_tsv
		Array[String]? batches
		String? docker
	}
	Array[String] columns=["filtered_depth_vcf","filtered_pesr_vcf"]
	
	scatter (column in columns) {
		call select.SelectBatchesAndColumns as SelectBatchesAndColumns {
			input:
				sample_sets_table_tsv=sample_sets_table_tsv,
				column=column,
				batches_opt=batches,
				docker=docker
		}
	}
	

	output {
		Array[File] filtered_depth_vcf=SelectBatchesAndColumns.column_as_row[0]
		Array[File] filtered_pesr_vcf=SelectBatchesAndColumns.column_as_row[1]
		Array[String]? batches_combined_pre_merge = batches
	}
	
}


