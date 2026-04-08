version 1.0

import "Structs.wdl"

workflow SparsifySD {
  input {
    Array[File] sd_files
    Array[File] sd_indices
    Array[String] sample_ids
    File sparse_sites_vcf
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_sparsify_sd
  }

  scatter (i in range(length(sample_ids))) {
    call SparsifySDTask {
      input:
        sd_file = sd_files[i],
        sd_index = sd_indices[i],
        sample_id = sample_ids[i],
        sparse_sites_vcf = sparse_sites_vcf,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_sparsify_sd
    }
  }

  output {
    Array[File] sparse_sd_files = SparsifySDTask.sparse_sd_file
    Array[File] sparse_sd_indices = SparsifySDTask.sparse_sd_index
  }
}

task SparsifySDTask {
  input {
    File sd_file
    File sd_index
    String sample_id
    File sparse_sites_vcf
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size([sd_file, sparse_sites_vcf], "GiB")
  Int vm_disk_size = ceil(input_size * 2 + 10)

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 3.75,
    disk_gb: vm_disk_size,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File sparse_sd_file = "~{sample_id}.sparse.sd.txt.gz"
    File sparse_sd_index = "~{sample_id}.sparse.sd.txt.gz.tbi"
  }

  command <<<
    set -euo pipefail

    python3 <<CODE
    import gzip
    import pysam
    import sys

    # Extract 0-based positions from sparse sites VCF, keyed by contig
    sites = dict()
    with pysam.VariantFile("~{sparse_sites_vcf}") as vcf:
        for record in vcf:
            # VCF POS is 1-based; SD file uses 0-based coordinates
            pos_0based = record.pos - 1
            sites.setdefault(record.chrom, set()).add(pos_0based)

    n_in = 0
    n_out = 0
    with gzip.open("~{sd_file}", "rt") as fin, \
         gzip.open("~{sample_id}.sparse.sd.txt.gz", "wt") as fout:
        for line in fin:
            n_in += 1
            fields = line.rstrip("\n").split("\t", 3)
            contig = fields[0]
            pos = int(fields[1])
            if contig in sites and pos in sites[contig]:
                fout.write(line)
                n_out += 1

    print(f"Sample ~{sample_id}: kept {n_out}/{n_in} records", file=sys.stderr)
    CODE

    tabix -f -0 -s1 -b2 -e2 "~{sample_id}.sparse.sd.txt.gz"
  >>>

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}
