version 1.0

import "Structs.wdl"
import "TasksMakeCohortVcf.wdl" as tasks

workflow ExtractVariantsInGeneList {
  input {
    Array[File] vcfs
    String prefix
    File gene_list  # one gene per line
    String svtypes  # comma-separated list of SVTYPE values to include in output, ie. "DEL,DUP"
    String consequences  # comma-separated list of functional consequence annotations to search for gene overlaps, ie. "PREDICTED_LOF,PREDICTED_COPY_GAIN"

    String sv_pipeline_docker
    String sv_base_mini_docker
  }

  scatter (vcf in vcfs) {
    call ExtractVariantsInGeneListTask {
      input:
        vcf=vcf,
        vcf_idx=vcf + ".tbi",
        gene_list=gene_list,
        svtypes=svtypes,
        consequences=consequences,
        sv_pipeline_docker=sv_pipeline_docker
    }
  }

  call tasks.ConcatVcfs {
    input:
      vcfs = ExtractVariantsInGeneListTask.extracted_variants_vcf,
      vcfs_idx = ExtractVariantsInGeneListTask.extracted_variants_vcf_idx,
      outfile_prefix = "~{prefix}.gene_list_extracted",
      sv_base_mini_docker = sv_base_mini_docker
  }

  call tasks.CatUncompressedFiles {
    input:
      shards = ExtractVariantsInGeneListTask.extracted_variants_table,
      outfile_name = "~{prefix}.vids_genes_table.extracted.tsv",
      sv_base_mini_docker = sv_base_mini_docker
  }

  output {
    File extracted_variants_vcf = ConcatVcfs.concat_vcf
    File extracted_variants_vcf_idxs = ConcatVcfs.concat_vcf_idx
    File extracted_variants_tables = CatUncompressedFiles.outfile
  }
}


task ExtractVariantsInGeneListTask {
  input {
    File vcf
    File vcf_idx
    File gene_list
    String svtypes
    String consequences
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  String prefix = basename(vcf, ".vcf.gz") + ".gene_list_extracted"

  RuntimeAttr runtime_default = object {
                                  mem_gb: 3.75,
                                  disk_gb: ceil(10.0 + size(vcf, "GB") * 2.0),
                                  cpu_cores: 1,
                                  preemptible_tries: 1,
                                  max_retries: 1,
                                  boot_disk_gb: 10
                                }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euo pipefail
    python <<CODE
import pysam
vcf = pysam.VariantFile("~{vcf}", 'r')
out = pysam.VariantFile("~{prefix}.vcf.gz", 'w', header=vcf.header)
svtypes = "~{svtypes}".split(",")
consequences = "~{consequences}".split(",")
with open("~{gene_list}", 'r') as gene_file:
    gene_set = set([line.strip('\n') for line in gene_file])
with open("~{prefix}.vids_genes_table.tsv", 'w') as tbl:
    for record in vcf:
        present = [x for x in consequences if x in record.info]
        genes = []
        if record.info['SVTYPE'] in svtypes and len(present) > 0:
            for consequence in present:
                for gene in record.info[consequence]:
                    if gene in gene_set:
                        genes.append(gene)
            if len(genes) > 0:
                out.write(record)
                tbl.write(f"{record.id}\t{record.info['AF'][0]}\t" + ",".join(genes) + "\n")
vcf.close()
out.close()
CODE
    tabix ~{prefix}.vcf.gz
  >>>

  output {
    File extracted_variants_vcf="~{prefix}.vcf.gz"
    File extracted_variants_vcf_idx="~{prefix}.vcf.gz.tbi"
    File extracted_variants_table="~{prefix}.vids_genes_table.tsv"
  }
}
