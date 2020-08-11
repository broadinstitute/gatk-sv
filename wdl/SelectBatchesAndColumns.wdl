version 1.0

task SelectBatchesAndColumns {
	input {
		File sample_sets_table_tsv
		Array[String]? batches_opt
		String column
		String? docker
	}
	
	Boolean select_batches=defined(batches_opt)
	Array[String] batches=select_first([batches_opt, ["all"]])

	String docker_set = select_first([docker, "ubuntu:18.04"])

	runtime {
		docker: docker_set
	}

	command <<<
		set -euo pipefail
		1>&2 echo "Selecting columns from sample sets table"
		FILE=~{sample_sets_table_tsv}

		COLNAME=~{column}
		1>&2 echo "COLUMN NAME=$COLNAME"

		COLNUM=$(head -1 $FILE | awk -v RS='\t' -v field="$COLNAME" '$0~field{print NR; exit}')
		1>&2 echo "COLUMN NUMBER=$COLNUM"

		# if batches specified, select rows containing given batches in 1st column then select specified column
		# otherwise, select specified column from all rows
		# delete quotes and brackets from any bracketed, comma-separated arrays and split by newlines
		for BATCH in '~{sep="' '" batches}'; do
			1>&2 echo "BATCH=$BATCH"
			~{if select_batches then "1>&2 echo 'Selecting batches from sample sets table'; awk -v FS='\t' -v bat=$BATCH '$1 == bat' $FILE" else "cat $FILE"} \
				| cut -f $COLNUM | tr -d []'"' | tr , '\n' >> ~{column}_selected.tsv
		done
		
	>>>

	output {
		# return transposed tsv such that each row is all the files for a given column of the sample set table
		Array[File] column_as_row = transpose(read_tsv("~{column}_selected.tsv"))[0]
	}
}
