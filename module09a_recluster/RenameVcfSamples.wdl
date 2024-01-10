version 1.0

import "Structs.wdl"

workflow RenameVcfSamples {
  input {
    File vcf
    File? vcf_index
    Array[String]? current_sample_ids
    Array[String]? new_sample_ids
    File? sample_id_rename_map  # TSV with old sample ID in first column, new sample ID in second column. Either provide this file or both of the string arrays
    String prefix
    Boolean? check_rename_all_samples  # if True, check that all sample IDs in header are in renaming map / current sample IDs list
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  call RenameVcfSamplesTask {
    input:
      vcf=vcf,
      vcf_idx=select_first([vcf_index, "~{vcf}.tbi"]),
      sample_id_rename_map = if (defined(current_sample_ids) && defined(new_sample_ids)) then
            write_tsv(transpose([select_first([current_sample_ids]), select_first([new_sample_ids])])) else
            select_first([sample_id_rename_map]),
      prefix="~{prefix}.renamed",
      check_rename_all_samples = check_rename_all_samples,
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_attr_override
  }

  output {
    File vcf_out = RenameVcfSamplesTask.out
    File vcf_out_index = RenameVcfSamplesTask.out_index
  }
}

task RenameVcfSamplesTask {
  input {
    File vcf
    File vcf_idx
    File sample_id_rename_map
    String prefix
    Boolean? check_rename_all_samples = true
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr runtime_default = object {
    mem_gb: 1.0,
    disk_gb: ceil(10 + size(vcf, "GiB") * 2.0),
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
    set -euo pipefail
    if ~{check_rename_all_samples}; then
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
    fi
    bcftools reheader --samples ~{sample_id_rename_map} -o ~{prefix}.vcf.gz ~{vcf}
    tabix ~{prefix}.vcf.gz
  >>>

  output {
    File out = "~{prefix}.vcf.gz"
    File out_index = "~{prefix}.vcf.gz.tbi"
  }
}
