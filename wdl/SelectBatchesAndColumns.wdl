version 1.0

import "Structs.wdl"

task SelectBatchesAndColumns {
  input {
    File sample_sets_table_tsv
    Array[String]? batches_opt
    Array[String] columns
    String? python_docker
    RuntimeAttr? runtime_attr_override
  }
  
  Boolean select_batches = defined(batches_opt)
  Array[String] batches = select_first([batches_opt, ["EMPTY_LIST_PLACEHOLDER"]])

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
    docker: select_first([python_docker, "python:3.7-slim"])
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

  command <<<
    set -euo pipefail
    python3 <<CODE
    import sys
    columns = ["~{sep='\",\"' columns}"]
    batches_lst = ["~{sep='\",\"' batches}"]
    all_batches = False
    if len(batches_lst) == 1 and batches_lst[0] == "EMPTY_LIST_PLACEHOLDER":
      all_batches = True
    batches = {batch: 0 for batch in batches_lst}

    with open("~{sample_sets_table_tsv}", 'r') as table:
      # output format: one tab-separated line of files (one per selected column, in order) per selected batch
      # this can be transposed to the format WDL can interpret as an Array[Array[File]]
      # note: assumes one file per column per batch. Array[File] outputs from Module04b should only contain one item because 04b is run per batch in Terra
      with open("selected.tsv", 'w') as out:
        header = None
        for entry in table:
          if header is None:
            # check input file format
            if not entry.startswith("entity:sample_set_id\t"):
              sys.exit("Malformed file: sample_sets_table_tsv should be a tab-separated sample sets table downloaded from Terra, starting with entity:sample_set_id.")
            # get columns dictionary from header line
            header = {colname:colnum for colnum,colname in enumerate(entry.strip().split('\t'))}
          else:
            line = entry.lstrip().split('\t')
            # select desired batches
            batch = line[header["entity:sample_set_id"]]
            if all_batches or batch in batches:
              if not all_batches:
                batches[batch] += 1
              output_line = ""
              # select desired columns
              for column in columns:
                if output_line != "":
                  output_line += '\t'
                # if column does not exist, enter DNE string - because some files are optional inputs/outputs so not ok to exit
                # but need to fill the spot so indices of outputs will not be disrupted 
                # (this is why simple list comprehension + tab-join was not used to generate output_line)
                if column not in header:
                  output_line += "COLUMN_DNE"
                  continue
                # get entry from this batch for each desired column, remove brackets if an array of files, and save
                file = line[header[column]].lstrip('["').rstrip('"]')
                output_line += file
              out.write(output_line + '\n')

    # raise exception if batch not found
    for batch in batches_lst:
      if not all_batches and batches[batch] == 0:
        sys.exit("Batch " + batch + " not found in sample sets table.")
    CODE
  >>>

  output {
    # return transposed tsv such that each row is all the files for a given column of the sample set table
    Array[Array[File]] columns_as_rows = transpose(read_tsv("selected.tsv"))
  }
}
