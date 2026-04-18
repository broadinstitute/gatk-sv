version 1.0

import "Structs.wdl"

workflow ExtractHprcSvBedAcrossContigs {
  input {
    Array[File] input_vcfs
    File extractor_script

    String source_value = "HPRC_SV_Integration"
    String output_prefix = "hprc_sv_integration"
    String python_docker = "python:3.11-slim"

    RuntimeAttr? runtime_attr_process
    RuntimeAttr? runtime_attr_concat
  }

  scatter (input_vcf in input_vcfs) {
    call ProcessContigVcf {
      input:
        input_vcf              = input_vcf,
        extractor_script       = extractor_script,
        source_value           = source_value,
        python_docker          = python_docker,
        runtime_attr_override  = runtime_attr_process
    }
  }

  call ConcatBedGz as concat_beds {
    input:
      bed_files              = ProcessContigVcf.output_bed,
      output_file            = output_prefix + ".all_contigs.bed.gz",
      python_docker          = python_docker,
      runtime_attr_override  = runtime_attr_concat
  }

  output {
    Array[File] filtered_vcfs = ProcessContigVcf.output_vcf
    Array[File] bed_files = ProcessContigVcf.output_bed
    File merged_bed = concat_beds.merged_bed
  }
}


task ProcessContigVcf {
  input {
    File input_vcf
    File extractor_script
    String source_value
    String python_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 4,
    disk_gb: 100,
    boot_disk_gb: 10,
    preemptible_tries: 1,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  String base0 = basename(input_vcf, ".vcf.gz")
  String base1 = if base0 == basename(input_vcf) then basename(input_vcf, ".vcf") else base0
  String base2 = if base1 == basename(input_vcf) then basename(input_vcf, ".bcf") else base1
  String output_vcf_name = base2 + ".hprc_sv_integration.only.vcf.gz"
  String output_bed_name = base2 + ".hprc_sv_integration.only.bed.gz"

  command <<<
    set -euo pipefail

    python3 ~{extractor_script} \
      --input-vcf ~{input_vcf} \
      --output-vcf ~{output_vcf_name} \
      --output-bed ~{output_bed_name} \
      --source-value ~{source_value}
  >>>

  output {
    File output_vcf = output_vcf_name
    File output_bed = output_bed_name
  }

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: python_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}


task ConcatBedGz {
  input {
    Array[File] bed_files
    String output_file
    String python_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 4,
    disk_gb: 50,
    boot_disk_gb: 10,
    preemptible_tries: 1,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail

    python3 << 'PY'
import gzip

inputs = "~{sep=' ' bed_files}".split()
out_path = "~{output_file}"

if not inputs:
    raise RuntimeError("No BED files provided for concatenation")

header_written = False
rows = 0

with gzip.open(out_path, 'wt') as out:
    for path in inputs:
        with gzip.open(path, 'rt') as fin:
            header = fin.readline()
            if not header:
                continue
            if not header_written:
                out.write(header)
                header_written = True
            for line in fin:
                if line.strip():
                    out.write(line)
                    rows += 1

print({"inputs": len(inputs), "rows": rows, "output": out_path})
PY
  >>>

  output {
    File merged_bed = output_file
  }

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: python_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}
