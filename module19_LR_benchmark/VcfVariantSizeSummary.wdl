version 1.0

import "Structs.wdl"

workflow VcfVariantSizeSummary {
  input {
    String prefix

    File vcf_gz
    File vcf_tbi

    Array[String] contig_list
    String sv_pipeline_base_docker

    RuntimeAttr? runtime_attr_detect_contigs
    RuntimeAttr? runtime_attr_count_contig
    RuntimeAttr? runtime_attr_merge
  }

  scatter (c in contig_list) {
    call SplitVcfByContig {
      input:
        vcf_gz = vcf_gz,
        vcf_tbi = vcf_tbi,
        contig = c,
        docker_image = sv_pipeline_base_docker,
        runtime_attr_override = runtime_attr_count_contig
      }

    call CountVariantsPerContig {
      input:
        vcf_gz = SplitVcfByContig.contig_vcf,
        vcf_tbi = SplitVcfByContig.contig_vcf_tbi,
        contig = c,
        docker_image = sv_pipeline_base_docker,
        runtime_attr_override = runtime_attr_count_contig
    }
  }

  call ConcatVariantSummaries {
    input:
      prefix = prefix, 
      contig_summary_tables = CountVariantsPerContig.count_table,
      docker_image = sv_pipeline_base_docker,
      runtime_attr_override = runtime_attr_merge
  }

  call SumVariantStatsAcrossContigs {
    input :
      prefix = prefix, 
      contig_summary_table = ConcatVariantSummaries.merged_summary,
      docker_image = sv_pipeline_base_docker,
      runtime_attr_override = runtime_attr_merge
     }

  output {
    File  summary_table = SumVariantStatsAcrossContigs.genome_wide_summary
  }
}

task SplitVcfByContig {
  input {
    File vcf_gz
    File vcf_tbi
    String contig
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 4,
    disk_gb: ceil(10 + size(vcf_gz, "GB") * 2),
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail

    bcftools view \
      -r ~{contig} \
      -f "PASS,SMALL_SV,." \
      -Oz \
      -o ~{contig}.vcf.gz \
      ~{vcf_gz}

    tabix -p vcf ~{contig}.vcf.gz
  >>>

  output {
    File contig_vcf = "~{contig}.vcf.gz"
    File contig_vcf_tbi = "~{contig}.vcf.gz.tbi"
  }

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: docker_image
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

task DetectContigs {
  input {
    File vcf_gz
    File vcf_tbi
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 4,
    disk_gb: ceil(10 + size(vcf_gz, "GB") * 2),
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail
    bcftools index -n ~{vcf_gz} > contigs.txt
  >>>

  output {
    Array[String] contigs = read_lines("contigs.txt")
  }

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: docker_image
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

task CountVariantsPerContig {
  input {
    File vcf_gz
    File vcf_tbi
    String contig
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 4,
    disk_gb: ceil(10 + size(vcf_gz, "GB") * 2),
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail

    python3 << 'EOF'
    import sys
    import gzip
    from collections import defaultdict

    vcf_file = "~{vcf_gz}"
    contig = "~{contig}"

    # counters
    counts = defaultdict(int)
    positions = defaultdict(set)

    def open_vcf(path):
        if path.endswith(".gz"):
            return gzip.open(path, "rt")
        return open(path, "r")

    with open_vcf(vcf_file) as f:
        for line in f:
            if line.startswith("#"):
                continue

            fields = line.rstrip().split("\t")
            chrom = fields[0]
            pos = fields[1]
            ref = fields[3]
            alt = fields[4]

            alts = alt.split(",")

            for a in alts:
                rlen = len(ref)
                alen = len(a)

                if rlen == 1 and alen == 1:
                    cat = "SNV"
                elif rlen > alen:
                    d = rlen - alen
                    cat = "DEL_1_49" if d <= 49 else "DEL_GT49"
                elif alen > rlen:
                    d = alen - rlen
                    cat = "INS_1_49" if d <= 49 else "INS_GT49"
                else:
                    continue

                counts[cat] += 1
                positions[cat].add(f"{chrom}:{pos}")

    categories = ["SNV", "DEL_1_49", "INS_1_49", "DEL_GT49", "INS_GT49"]

    with open("variant_summary.tsv", "w") as out:
        out.write("CONTIG\tCATEGORY\tVARIANT_COUNT\tUNIQUE_POS_COUNT\n")
        for c in categories:
            out.write(
                f"{contig}\t{c}\t{counts[c]}\t{len(positions[c])}\n"
            )
    EOF
  >>>

  output {
    File count_table = "variant_summary.tsv"
  }

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: docker_image
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

task ConcatVariantSummaries {
  input {
    Array[File] contig_summary_tables
    String prefix
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 2,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail

    # write header from first file
    head -n 1 ~{contig_summary_tables[0]} > "~{prefix}.variant_summary.all_contigs.tsv"

    # append data rows from all files (skip headers)
    for f in ~{sep=' ' contig_summary_tables}; do
      tail -n +2 "$f" >> "~{prefix}.variant_summary.all_contigs.tsv"
    done
  >>>

  output {
    File merged_summary = "~{prefix}.variant_summary.all_contigs.tsv"
  }

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: docker_image
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

task SumVariantStatsAcrossContigs {
  input {
    String prefix
    File contig_summary_table
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 2,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail

    python3 <<CODE
    import pandas as pd

    # Desired category order
    category_order = [
        "SNV",
        "DEL_1_49",
        "INS_1_49",
        "DEL_GT49",
        "INS_GT49"
    ]

    # Read input
    df = pd.read_csv("~{contig_summary_table}", sep="\t")

    # Pivot VARIANT_COUNT
    vc = df.pivot(index="CONTIG",
                  columns="CATEGORY",
                  values="VARIANT_COUNT")

    # Reindex to enforce order
    vc = vc.reindex(columns=category_order, fill_value=0)
    vc.columns = [f"{c}_all" for c in vc.columns]

    # Pivot UNIQUE_POS_COUNT
    up = df.pivot(index="CONTIG",
                  columns="CATEGORY",
                  values="UNIQUE_POS_COUNT")

    # Reindex to enforce order
    up = up.reindex(columns=category_order, fill_value=0)
    up.columns = [f"{c}_unique" for c in up.columns]

    # Merge wide tables
    out = vc.join(up)

    # ---- add summary row ----
    summary = out.sum(axis=0)
    summary.name = "all"
    out = pd.concat([out, summary.to_frame().T])

    # Reset index for output
    out = out.reset_index().rename(columns={"index": "CONTIG"})

    # Write output
    out.to_csv("~{prefix}.variant_summary.genome_wide.tsv", sep="\t", index=False)

    CODE

  >>>

  output {
    File genome_wide_summary = "~{prefix}.variant_summary.genome_wide.tsv"
  }

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: docker_image
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}