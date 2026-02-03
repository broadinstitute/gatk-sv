version 1.0

import "Structs.wdl"

workflow SampleVariantReport {
  input {
    File vcf_gz
    File vcf_idx
    File sample_list
    String docker_image

    RuntimeAttr? runtime_attr_split
    RuntimeAttr? runtime_attr_count
    RuntimeAttr? runtime_attr_annot
    RuntimeAttr? runtime_attr_merge
  }

  scatter (sample in read_lines(sample_list)) {

    call SplitSampleNonRef {
      input:
        vcf_gz = vcf_gz,
        vcf_idx = vcf_idx,
        sample = sample,
        docker_image = docker_image,
        runtime_attr_override = runtime_attr_split
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

# =====================================================
# Task 1: Split per-sample non-ref VCF
# =====================================================
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

# =====================================================
# Task 2: Count variant size classes + PREDICTED_LOF
# =====================================================
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

  command <<<
    set -euo pipefail

    bcftools view ~{sample_vcf} -Ov | \
    awk '
    BEGIN {
      FS=OFS="\t"
    }
    !/^#/ {
      ref=$4
      split($5, alts, ",")
      info=$8

      lof=""
      if (info ~ /PREDICTED_LOF=/) {
        match(info,/PREDICTED_LOF=([^;]+)/,m)
        lof=m[1]
      }

      for (i in alts) {
        r=length(ref); a=length(alts[i])

        if (r==1 && a==1) cat="SNV"
        else if (r>a && r-a<=49) cat="DEL_1_49"
        else if (a>r && a-r<=49) cat="INS_1_49"
        else if (r>a) cat="DEL_GT49"
        else if (a>r) cat="INS_GT49"
        else continue

        total[cat]++

        if (lof!="") {
          lof_count[cat]++
          n=split(lof, g, ",")
          for (j=1;j<=n;j++) lof_gene[cat,g[j]]=1
        }
      }
    }
    END {
      print "CATEGORY\tTOTAL\tPREDICTED_LOF\tLOF_GENES"
      for (c in total) {
        genes=""
        for (k in lof_gene) {
          split(k,a,SUBSEP)
          if (a[1]==c) genes=(genes==""?a[2]:genes","a[2])
        }
        print c,total[c]+0,lof_count[c]+0,genes
      }
    }' > ~{sample}.class_counts.tsv
  >>>

  output {
    File class_table = "~{sample}.class_counts.tsv"
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

# =====================================================
# Task 3: VEP consequence + coding-disruptive genes
# =====================================================
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

  command <<<
    set -euo pipefail

    bcftools view ~{sample_vcf} -Ov | \
    awk '
    BEGIN { FS=OFS="\t" }
    !/^#/ {
      if ($8 ~ /vep=/) {
        match($8,/vep=([^;]+)/,m)
        n=split(m[1],csqs,",")
        for (i=1;i<=n;i++) {
          split(csqs[i],f,"|")
          cons=f[2]; gene=f[4]
          vep[cons]++
          if (cons ~ /frameshift|stop_gained|splice_(acceptor|donor)|start_lost/)
            gene_set[cons,gene]=1
        }
      }
    }
    END {
      print "CONSEQUENCE\tCOUNT\tGENES"
      for (k in vep) {
        genes=""
        for (g in gene_set) {
          split(g,a,SUBSEP)
          if (a[1]==k) genes=(genes==""?a[2]:genes","a[2])
        }
        print k,vep[k],genes
      }
    }' > ~{sample}.annotation.tsv
  >>>

  output {
    File annotation_table = "~{sample}.annotation.tsv"
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

# =====================================================
# Task 4: Merge per-sample report
# =====================================================
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