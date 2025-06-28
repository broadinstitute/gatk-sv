version 1.0
import "Structs.wdl"

workflow ExtractFileByIndex {
  input {
    Array[Array[File]] nested_files  # Input: list of list of files
    Int index
  }

  # Extract the first file from each inner list
  scatter (file_list in nested_files) {
    File first_file = file_list[index]
  }

  # You can pass the extracted first files to a task
  call ProcessFiles {
    input:
      input_files = first_file
  }

  output {
    Array[File] extracted_first_files = first_file
    Array[File] processed_files = ProcessFiles.output_files
  }
}

task ProcessFiles {
  input {
    Array[File] input_files
  }

  command <<<
    mkdir processed
    for file in ~{sep=' ' input_files}; do
      cp "$file" processed/
    done
  >>>

  output {
    Array[File] output_files = glob("processed/*")
  }

  runtime {
    docker: "ubuntu:20.04"
    memory: "2G"
    cpu: 1
  }
}
