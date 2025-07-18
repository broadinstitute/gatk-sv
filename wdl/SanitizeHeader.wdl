version 1.0

import "Structs.wdl"

workflow SanitizeHeader {
  input {
    Array[File] vcfs
    String prefix
    String drop_fields
    File sample_id_rename_map
    File primary_contigs_list
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_sanitize_header
  }

  Array[String] contigs = read_lines(primary_contigs_list)
  scatter (i in range(length(vcfs))) {
    call SanitizeHeaderTask {
      input:
        vcf=vcfs[i],
        vcf_index="~{vcfs[i]}.tbi",
        prefix="~{prefix}.~{contigs[i]}",
        drop_fields=drop_fields,
        sample_id_rename_map=sample_id_rename_map,
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_attr_sanitize_header
    }
  }

  output {
    Array[File] vcf_header_sanitized = SanitizeHeaderTask.out
    Array[File] vcf_header_sanitized_index = SanitizeHeaderTask.out_index

  }
}

task SanitizeHeaderTask {
  input {
    File vcf
    File vcf_index
    String prefix
    File sample_id_rename_map
    String drop_fields
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr runtime_default = object {
                                  mem_gb: 3.75,
                                  disk_gb: ceil(10.0 + size(vcf, "GiB") * 2.0),
                                  cpu_cores: 1,
                                  preemptible_tries: 3,
                                  max_retries: 1,
                                  boot_disk_gb: 10
                                }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GiB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euxo pipefail

    bcftools query -l ~{vcf} | sort > samples.list
    python <<CODE
with open("~{sample_id_rename_map}", 'r') as rename, open("samples.list", 'r') as header:
  all_to_rename = set()
  for line in rename:
    all_to_rename.add(line.strip("\n").split("\t")[0])  # create set of sample IDs to rename
  for line in header:
    sample = line.strip("\n")
    if sample not in all_to_rename:
      raise ValueError(f"Sample {sample} is in the VCF header but not in the renaming map")
CODE

    bcftools view --no-version -h ~{vcf} > header.vcf

    grep -v -e ^"##bcftools" header.vcf \
      -e ^"##source=depth" \
      -e ^"##source=cleanvcf" \
      -e ^"##ALT=<ID=UNR" \
      -e "assembly=38" \
      | sed 's/Split read genotype quality/Split-read genotype quality/g' \
      | sed 's/##ALT=<ID=BND,Description="Translocation">/##ALT=<ID=BND,Description="Breakend">/g' > newheader.vcf

    bcftools reheader -h newheader.vcf ~{vcf} --samples ~{sample_id_rename_map} \
      | bcftools annotate -x ~{drop_fields} \
        --no-version \
        -O z \
        -o ~{prefix}.vcf.gz

    tabix ~{prefix}.vcf.gz
  >>>

  output {
    File out = "~{prefix}.vcf.gz"
    File out_index = "~{prefix}.vcf.gz.tbi"
  }
}
