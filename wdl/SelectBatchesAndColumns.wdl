version 1.0

task SelectBatchesAndColumns {
  input {
    File sample_sets_table_tsv
    Array[String]? batches_opt
    String column
    String? linux_docker
    RuntimeAttr? runtime_attr_override
  }
  
  Boolean select_batches = defined(batches_opt)
  Array[String] batches = select_first([batches_opt, ["all"]])

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 0.9,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: select_first([linux_docker, "ubuntu:18.04"])
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
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
