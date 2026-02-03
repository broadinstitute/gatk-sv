version 1.0

############################
# Workflow
############################

workflow SplitAndAnnotateVcf {
  input {
    File multi_sample_vcf
    File sample_list
    File annotation1
    File annotation2
    String docker_image

    RuntimeAttr? runtime_attr_override
  }

  call ReadSampleList {
    input:
      sample_list = sample_list,
      docker_image = docker_image,
      runtime_attr_override = runtime_attr_override
  }

  scatter (sample in ReadSampleList.samples) {

    call SplitNonRefSampleVcf {
      input:
        input_vcf = multi_sample_vcf,
        sample = sample,
        docker_image = docker_image,
        runtime_attr_override = runtime_attr_override
    }

    call AnnotateWithAnno1 {
      input:
        input_vcf = SplitNonRefSampleVcf.out_vcf,
        annotation1 = annotation1,
        docker_image = docker_image,
        runtime_attr_override = runtime_attr_override
    }

    call AnnotateWithAnno2 {
      input:
        input_vcf = AnnotateWithAnno1.out_vcf,
        annotation2 = annotation2,
        docker_image = docker_image,
        runtime_attr_override = runtime_attr_override
    }
  }

  output {
    Array[File] annotated_vcfs = AnnotateWithAnno2.out_vcf
  }
}

############################
# Struct
############################

struct RuntimeAttr {
  Int cpu_cores
  Int mem_gb
  Int disk_gb
  Int boot_disk_gb
  Int preemptible_tries
  Int max_retries
}

############################
# Tasks
############################

task ReadSampleList {
  input {
    File sample_list
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  command <<<
    set -euo pipefail
    cat ~{sample_list}
  >>>

  output {
    Array[String] samples = read_lines(stdout())
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 4,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

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

---

task SplitNonRefSampleVcf {
  input {
    File input_vcf
    String sample
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  command <<<
    set -euo pipefail

    bcftools view \
      -s ~{sample} \
      -c 1 \
      -Oz \
      -o ~{sample}.nonref.vcf.gz \
      ~{input_vcf}

    bcftools index ~{sample}.nonref.vcf.gz
  >>>

  output {
    File out_vcf = "~{sample}.nonref.vcf.gz"
    File out_vcf_index = "~{sample}.nonref.vcf.gz.csi"
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 4,
    disk_gb: ceil(10 + size(input_vcf, "GB") * 2),
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

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

---

task AnnotateWithAnno1 {
  input {
    File input_vcf
    File annotation1
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  command <<<
    set -euo pipefail

    python3 << 'PYTHON'
import gzip

def open_text(f):
    return gzip.open(f, "rt") if f.endswith(".gz") else open(f)

ann = {}
with open_text("~{annotation1}") as f:
    for line in f:
        fields = line.rstrip().split("\t")
        if len(fields) >= 6:
            ann[fields[4]] = fields[5]

with gzip.open("anno1.vcf.gz", "wt") as out:
    with open_text("~{input_vcf}") as fin:
        for line in fin:
            if line.startswith("#"):
                out.write(line)
                continue
            fields = line.rstrip().split("\t")
            vid = fields[2]
            if vid in ann:
                fields[7] += f";VEP={ann[vid]}"
            out.write("\t".join(fields) + "\n")
PYTHON

    bcftools index anno1.vcf.gz
  >>>

  output {
    File out_vcf = "anno1.vcf.gz"
    File out_vcf_index = "anno1.vcf.gz.csi"
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 4,
    disk_gb: ceil(10 + size(input_vcf, "GB") * 2),
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

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

---

task AnnotateWithAnno2 {
  input {
    File input_vcf
    File annotation2
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  command <<<
    set -euo pipefail

    python3 << 'PYTHON'
import gzip

def open_text(f):
    return gzip.open(f, "rt") if f.endswith(".gz") else open(f)

with open_text("~{annotation2}") as f:
    header = f.readline().rstrip().split("\t")
    keys = header[5:]
    ann = {}
    for line in f:
        fields = line.rstrip().split("\t")
        ann[fields[4]] = dict(zip(keys, fields[5:]))

with gzip.open("final.vcf.gz", "wt") as out:
    with open_text("~{input_vcf}") as fin:
        for line in fin:
            if line.startswith("#"):
                out.write(line)
                continue
            fields = line.rstrip().split("\t")
            vid = fields[2]
            if vid in ann:
                for k, v in ann[vid].items():
                    if v not in (".", ""):
                        fields[7] += f";{k}={v}"
            out.write("\t".join(fields) + "\n")
PYTHON

    bcftools index final.vcf.gz
  >>>

  output {
    File out_vcf = "final.vcf.gz"
    File out_vcf_index = "final.vcf.gz.csi"
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 4,
    disk_gb: ceil(10 + size(input_vcf, "GB") * 2),
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

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