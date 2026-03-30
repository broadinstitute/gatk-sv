version 1.0

import "Structs.wdl"
import "TasksMakeCohortVcf.wdl" as MiniTasks

workflow AnnotateBiallelicCNCounts {

  input {
    Array[File] vcfs
    String prefix

    File par_bed
    File ploidy_table

    String sv_pipeline_docker
    String sv_base_mini_docker

    RuntimeAttr? runtime_attr_annotate
    RuntimeAttr? runtime_attr_concat
  }

  scatter (vcf in vcfs) {
    call AnnotateCNCounts {
      input:
        vcf = vcf,
        par_bed = par_bed,
        ploidy_table = ploidy_table,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_annotate
    }
  }

  call MiniTasks.ConcatVcfs {
    input:
      vcfs = AnnotateCNCounts.annotated_vcf,
      naive = true,
      sites_only = true,
      outfile_prefix = "~{prefix}.cn_counts",
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_concat
  }

  output {
    File annotated_vcf = ConcatVcfs.concat_vcf
    File annotated_vcf_idx = ConcatVcfs.concat_vcf_idx
  }
}

task AnnotateCNCounts {
  input {
    File vcf
    File par_bed
    File ploidy_table
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  String output_file = basename(vcf, ".vcf.gz") + ".cn_counts.vcf.gz"

  Float input_size = size(vcf, "GB")
  RuntimeAttr runtime_default = object {
    mem_gb: 3.75,
    disk_gb: ceil(10.0 + input_size * 2.0),
    cpu_cores: 1,
    preemptible_tries: 3,
    max_retries: 1,
    boot_disk_gb: 10
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, runtime_default])

  runtime {
    memory: select_first([runtime_attr.mem_gb, runtime_default.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, runtime_default.disk_gb]) + " HDD"
    cpu: select_first([runtime_attr.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_attr.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euo pipefail

    python3 <<CODE
import pysam
from collections import defaultdict


PLOIDY = dict()
PAR = defaultdict(list)


with open("~{par_bed}", 'r') as inp:
    for line in inp:
        chrom, start, end, label = line.strip("\n").split("\t")
        PAR[chrom].append((int(start), int(end)))  # PAR[chrom] = [(start1, end1), (start2, end2)]


with open("~{ploidy_table}", 'r') as inp:
    first = True
    header = []
    for line in inp:
        fields = line.strip("\n").split("\t")
        if first:
            header = fields
            first = False
            continue
        PLOIDY[fields[0]] = dict()  # first column = sample
        for i in range(1, len(fields)):
            PLOIDY[fields[0]][header[i]] = int(fields[i])  # PLOIDY[sample][contig] = ploidy


with pysam.VariantFile("~{vcf}") as vcf_in:
    out_header = vcf_in.header
    out_header.add_line('##INFO=<ID=RD_CN_ESTIMATED_AF,Number=1,Type=Float,Description="Estimated AF from RD_CN for CNVs">')
    with pysam.VariantFile("~{output_file}", 'w', header=out_header) as vcf_out:
        for record in vcf_in:
            svtype = record.info.get("SVTYPE", None)
            if svtype in ("DEL", "DUP", "CNV"):
                ac = 0
                an = 0
                par = False
                if record.chrom == "chrX" or record.chrom == "chrY":
                    for start, end in PAR[record.chrom]:
                        if float(min(record.stop, end) - max(record.pos, start))/(record.stop - record.pos) > 0.5:
                            par = True  # if variant is >50% in PAR
                for sample in record.samples:
                    rd_cn = record.samples[sample].get("RD_CN", None)
                    ecn = PLOIDY[sample][record.chrom]
                    if par and ecn == 1:
                        ecn = 2  # males expected CN in PAR on chrX = 2; no variants contained within par on chrY in gnomAD
                    if rd_cn is None:
                        continue
                    else:
                        an += ecn
                    if rd_cn == ecn - 1 or rd_cn == ecn + 1:
                        ac += 1
                    elif rd_cn < ecn - 1 or rd_cn > ecn + 1:
                        ac += 2

                if an > 0:
                    record.info["RD_CN_ESTIMATED_AF"] = float(ac)/float(an)
            vcf_out.write(record)


CODE

    tabix -p vcf "~{output_file}"
  >>>

  output {
    File annotated_vcf = output_file
    File annotated_vcf_idx = output_file + ".tbi"
  }
}
