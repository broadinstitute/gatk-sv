version 1.0

import "TasksMakeCohortVcf.wdl" as tasks

workflow ConcatTextFiles {

  input {
    Array[File] text_files
    String output_prefix
    String output_suffix = "concat.txt"

    Boolean gzipped = false
    Boolean headered = false

    String linux_docker
    String sv_base_mini_docker
  }

  if (!headered) {
    # Disable filter command since input might be compressed
    call tasks.CatUncompressedFiles {
      input:
        shards=text_files,
        outfile_name="~{output_prefix}.~{output_suffix}",
        filter_command="",
        sv_base_mini_docker=sv_base_mini_docker
    }
  }

  if (headered) {
    call tasks.ConcatHeaderedTextFiles {
      input:
        text_files=text_files,
        gzipped=gzipped,
        output_filename="~{output_prefix}.~{output_suffix}",
        linux_docker=linux_docker
    }
  }

  output {
    File concatenated_files = select_first([ConcatHeaderedTextFiles.out, CatUncompressedFiles.outfile])
  }
}
