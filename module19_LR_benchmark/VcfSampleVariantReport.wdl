version 1.0

import "Structs.wdl"
import "VcfVariantSizeSummary.wdl" as VcfVariantSizeSummary

workflow SampleVariantReport {
  input {
    File vcf_gz
    File vcf_idx
    Array[String]? sample_list
    String docker_image

    RuntimeAttr? runtime_attr_split
    RuntimeAttr? runtime_attr_count
    RuntimeAttr? runtime_attr_annot
    RuntimeAttr? runtime_attr_merge
    RuntimeAttr? runtime_attr_count_variant_per_sample
  }


  if (!defined(sample_list)){
      call ExtractSamplesFromVcf {
        input:
          vcf = vcf_gz,
          vcf_idx = vcf_idx,
          docker_image = docker_image
      }

  }

  Array[String] sample_list_new = select_first([ExtractSamplesFromVcf.samples, sample_list])

  scatter (sample in sample_list_new) {

    call SplitSampleNonRef {
      input:
        vcf_gz = vcf_gz,
        vcf_idx = vcf_idx,
        sample = sample,
        docker_image = docker_image,
        runtime_attr_override = runtime_attr_split
    }

    call CountVariantsPerSample {
      input:
        vcf_gz = SplitSampleNonRef.sample_vcf,
        vcf_tbi = SplitSampleNonRef.sample_vcf_tbi,
        sample = sample,
        docker_image = docker_image,
        runtime_attr_override = runtime_attr_count_variant_per_sample
    }

    call CountVariantClasses {
      input:
        sample_vcf = SplitSampleNonRef.sample_vcf,
        sample = sample,
        docker_image = docker_image,
        runtime_attr_override = runtime_attr_count
    }

    call SummarizeAnnotations {
      input:
        sample_vcf = SplitSampleNonRef.sample_vcf,
        sample = sample,
        docker_image = docker_image,
        runtime_attr_override = runtime_attr_annot
    }

    call MergeSampleReport {
      input:
        class_counts = CountVariantClasses.class_table,
        annotation_summary = SummarizeAnnotations.annotation_table,
        sample = sample,
        docker_image = docker_image,
        runtime_attr_override = runtime_attr_merge
    }
  }

  output {
    Array[File] reports = MergeSampleReport.report
  }
}

# Task 1: extract sample list
task ExtractSamplesFromVcf {
  input {
    File vcf
    File vcf_idx
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 1,
    disk_gb: 5,
    boot_disk_gb: ceil(10 + size(vcf, "GB") * 2),
    preemptible_tries: 0,
    max_retries: 1
  }

  RuntimeAttr runtime_attr =
    select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail
    bcftools query -l ~{vcf} > samples.txt
  >>>

  output {
    Array[String] samples = read_lines("samples.txt")
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

# Task 1: Split per-sample non-ref VCF
task SplitSampleNonRef {
  input {
    File vcf_gz
    File vcf_idx
    String sample
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
    bcftools view -s ~{sample} -c 1 -Oz \
      -o ~{sample}.nonref.vcf.gz ~{vcf_gz}
    bcftools index ~{sample}.nonref.vcf.gz
  >>>

  output {
    File sample_vcf = "~{sample}.nonref.vcf.gz"
    File sample_vcf_tbi = "~{sample}.nonref.vcf.gz.tbi"
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

# Task 2: Count variant size classes + PREDICTED_LOF
task CountVariantClasses {
  input {
    File sample_vcf
    String sample
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 4,
    disk_gb: ceil(10 + size(sample_vcf, "GB") * 2),
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  String prefix = basename(sample_vcf , ".vcf.gz")

  command <<<
    set -euo pipefail

    python3 <<CODE

    import gzip
    import re
    from collections import defaultdict

    vcf_file = "~{sample_vcf}"
    out_file = "~{prefix}.SV_lof_summary.tsv"

    categories = [
        "SNV",
        "DEL_1_49",
        "INS_1_49",
        "DEL_GT49",
        "INS_GT49"
    ]

    total = defaultdict(int)
    lof_count = defaultdict(int)
    lof_genes = defaultdict(set)

    lof_re = re.compile(r"PREDICTED_LOF=([^;]+)")

    with gzip.open(vcf_file, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            ref = fields[3]
            alts = fields[4].split(",")
            info = fields[7]
            for alt in alts:
                rlen = len(ref)
                alen = len(alt)
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
                total[cat] += 1
                m = lof_re.search(info)
                if m:
                    lof_count[cat] += 1
                    genes = m.group(1).split(",")
                    for g in genes:
                        if g != "." and g != "":
                            lof_genes[cat].add(g)

    with open(out_file, "w") as out:
        out.write("CATEGORY\tVARIANT_COUNT\tLOF_VARIANT_COUNT\tLOF_GENES\n")
        for cat in categories:
            genes = ",".join(sorted(lof_genes[cat]))
            out.write(
                f"{cat}\t{total[cat]}\t{lof_count[cat]}\t{genes}\n"
            )

    CODE

  >>>

  output {
    File class_table = "~{prefix}.SV_lof_summary.tsv"
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

# Task 3: VEP consequence + coding-disruptive genes
task SummarizeAnnotations {
  input {
    File sample_vcf
    String sample
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 4,
    disk_gb: ceil(10 + size(sample_vcf, "GB") * 2),
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  String prefix = basename(sample_vcf , ".vcf.gz")

  command <<<
    set -euo pipefail

    python3 <<CODE

    import sys
    import gzip
    from collections import defaultdict

    vcf = "~{sample_vcf}"

    def open_vcf(path):
        if path.endswith(".gz"):
            return gzip.open(path, "rt")
        return open(path)

    def variant_class(ref, alt):
        r = len(ref)
        a = len(alt)
        if r == 1 and a == 1:
            return "SNV"
        if r > a:
            d = r - a
            return "DEL_1_49" if d <= 49 else "DEL_GT49"
        if a > r:
            d = a - r
            return "INS_1_49" if d <= 49 else "INS_GT49"
        return None

    # class -> consequence -> transcript -> set(variant_ids)
    counts = defaultdict(lambda: defaultdict(lambda: defaultdict(set)))
    # class -> consequence -> transcript -> set(genes)
    genes = defaultdict(lambda: defaultdict(lambda: defaultdict(set)))
    # transcript -> canonical flag
    canonical_flag = {}

    with open_vcf(vcf) as f:
        for line in f:
            if line.startswith("#"):
                continue

            fields = line.rstrip().split("\t")
            chrom, pos, vid, ref, alts, _, _, info = fields[:8]
            alts = alts.split(",")

            info_dict = {}
            for x in info.split(";"):
                if "=" in x:
                    k, v = x.split("=", 1)
                    info_dict[k] = v

            if "vep" not in info_dict:
                continue

            csq_entries = info_dict["vep"].split(",")

            for alt in alts:
                vclass = variant_class(ref, alt)
                if vclass is None:
                    continue

                var_id = f"{chrom}:{pos}:{ref}:{alt}"

                for csq in csq_entries:
                    parts = csq.split("|")
                    if len(parts) < 25:
                        continue

                    consequence_field = parts[1]
                    gene = parts[3]
                    transcript = parts[6]        # Feature (transcript ID)
                    canonical = parts[26]         # CANONICAL field

                    if transcript:
                        canonical_flag[transcript] = canonical if canonical else "NO"

                    for consequence in consequence_field.split("&"):
                        counts[vclass][consequence][transcript].add(var_id)
                        if gene:
                            genes[vclass][consequence][transcript].add(gene)

    with open("~{prefix}.vep_consequence_summary.tsv", "w") as out:
        out.write(
            "VARIANT_CLASS\tCONSEQUENCE\tTRANSCRIPT_ID\tCANONICAL\tVARIANT_COUNT\tGENES\n"
        )

        for vclass in sorted(counts):
            for cons in sorted(counts[vclass]):
                for transcript in sorted(counts[vclass][cons]):
                    n = len(counts[vclass][cons][transcript])
                    gene_list = ",".join(sorted(genes[vclass][cons][transcript]))
                    canonical = canonical_flag.get(transcript, "NO")

                    out.write(
                        f"{vclass}\t{cons}\t{transcript}\t{canonical}\t{n}\t{gene_list}\n"
                    )
    CODE
  >>>

  output {
    File annotation_table = "~{prefix}.vep_consequence_summary.tsv"
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

# Task 4: Merge per-sample report
task MergeSampleReport {
  input {
    File class_counts
    File annotation_summary
    String sample
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
    echo "### SAMPLE: ~{sample}" > ~{sample}.report.tsv
    echo "" >> ~{sample}.report.tsv
    cat ~{class_counts} >> ~{sample}.report.tsv
    echo "" >> ~{sample}.report.tsv
    cat ~{annotation_summary} >> ~{sample}.report.tsv
  >>>

  output {
    File report = "~{sample}.report.tsv"
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


task CountVariantsPerSample {
  input {
    File vcf_gz
    File vcf_tbi
    String sample
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
    sample = "~{sample}"

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

    with open("~{sample}.variant_summary.tsv", "w") as out:
        out.write("SAMPLE\tCATEGORY\tVARIANT_COUNT\tUNIQUE_POS_COUNT\n")
        for c in categories:
            out.write(
                f"{sample}\t{c}\t{counts[c]}\t{len(positions[c])}\n"
            )
    EOF
  >>>

  output {
    File count_table = "~{sample}.variant_summary.tsv"
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



