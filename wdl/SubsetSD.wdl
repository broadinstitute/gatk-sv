version 1.0

import "Structs.wdl"

workflow SubsetSD {
  input {
    File sd_file
    File sites_vcf
    String sv_pipeline_docker
    String prefix = basename(sd_file, ".sd.txt.gz")
    Int stride = 1
    RuntimeAttr? runtime_attr_subset_sd
  }

  call SubsetSDTask {
    input:
      sd_file = sd_file,
      sites_vcf = sites_vcf,
      prefix = prefix,
      stride = stride,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_subset_sd
  }

  output {
    File subset_sd_file = SubsetSDTask.subset_sd_file
    File subset_sd_index = SubsetSDTask.subset_sd_index
  }
}

task SubsetSDTask {
  input {
    File sd_file
    File sites_vcf
    String prefix
    Int stride = 1
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  String subset_sd_txt = "~{prefix}.subset.sd.txt"
  String subset_sd_filename = subset_sd_txt + ".gz"

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 3.75,
    disk_gb: 100,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail

    python3 <<CODE
    from collections import defaultdict
    import gzip
    import pysam
    import sys

    stride = ~{stride}
    if stride < 1:
        raise ValueError(f"Stride must be >= 1; got {stride}")

    sites_by_contig = defaultdict(set)
    with pysam.VariantFile("~{sites_vcf}") as vcf:
        for record in vcf:
            # VCF POS is 1-based; SD files use 0-based coordinates.
            sites_by_contig[record.chrom].add(record.pos - 1)

    matched_records = 0
    kept_records = 0
    with gzip.open("~{sd_file}", "rt") as fin, open("~{subset_sd_txt}", "wt") as fout:
        for line in fin:
            fields = line.split("\t", 2)
            if len(fields) < 2:
                raise ValueError(f"Malformed SD line: {line.rstrip()}")

            contig = fields[0]
            pos = int(fields[1])
            contig_sites = sites_by_contig.get(contig)
            if contig_sites is None or pos not in contig_sites:
                continue

            if matched_records % stride == 0:
                fout.write(line)
                kept_records += 1
            matched_records += 1

    print(
        f"Matched {matched_records} SD records; kept {kept_records} with stride {stride}",
        file=sys.stderr,
    )
    CODE

    bgzip -f "~{subset_sd_txt}"
    tabix -f -0 -s1 -b2 -e2 "~{subset_sd_filename}"
  >>>

  output {
    File subset_sd_file = subset_sd_filename
    File subset_sd_index = subset_sd_filename + ".tbi"
  }

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    noAddress: true
  }
}