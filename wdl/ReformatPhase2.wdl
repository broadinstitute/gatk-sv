version 1.0

import "Structs.wdl"

workflow ReformatPhase2 {
  input {
    File pesr_vcf
    File depth_vcf
    String prefix
    String sv_pipeline_docker
  }

  call ReformatVcf as ReformatPesr {
    input:
      vcf=pesr_vcf,
      vcf_idx="~{pesr_vcf}.tbi",
      prefix="~{prefix}.pesr",
      sv_pipeline_docker=sv_pipeline_docker
  }

  call ReformatVcf as ReformatDepth {
    input:
      vcf=depth_vcf,
      vcf_idx="~{depth_vcf}.tbi",
      prefix="~{prefix}.depth",
      sv_pipeline_docker=sv_pipeline_docker
  }

  output {
    File reformatted_pesr_vcf = ReformatPesr.reformatted_vcf
    File reformatted_pesr_vcf_idx = ReformatPesr.reformatted_vcf_idx

    File reformatted_depth_vcf = ReformatDepth.reformatted_vcf
    File reformatted_depth_vcf_idx = ReformatDepth.reformatted_vcf_idx

    File original_pesr_vcf = pesr_vcf
    File original_depth_vcf = depth_vcf
  }
}

task ReformatVcf {
  input {
    File vcf
    File vcf_idx
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 3.75,
    disk_gb: ceil(10 + 3 * size(vcf, "GB")),
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File reformatted_vcf = "${prefix}.reformatted.vcf.gz"
    File reformatted_vcf_idx = "~{prefix}.reformatted.vcf.gz.tbi"
  }
  command <<<

    set -euo pipefail

    python3 <<CODE
import gzip
import pysam


bnd_end_dict = dict()
with gzip.open("~{vcf}", 'rt') as f, open("intermediate1.vcf", 'w') as out:
    for line in f:
        if line.startswith('#'):
            out.write(line)
            continue
        columns = line.split('\t', 8)
        alt = columns[4]
        # replace breakend alt alleles with symbolic <BND> allele
        if "]" in alt or "[" in alt:
            columns[4] = "<BND>"
        vid = columns[2]
        info = columns[7]
        # save END for BND/CTX in dictionary
        if 'SVTYPE=BND' in info or 'SVTYPE=CTX' in info:
            info_tokens = info.split(';')
            end_field_list = [x for x in info_tokens if x.startswith("END=")]
            if len(end_field_list) > 0:
                end = int(end_field_list[0].replace("END=", ""))
            else:
                # Special case where END and POS happen to be equal
                end = int(columns[1])
            bnd_end_dict[vid] = end
        out.write("\t".join(columns))

# set BND/CTX END2 to ends from dictionary if not present
with pysam.VariantFile("intermediate1.vcf", 'r') as vcf:
    header = vcf.header
    header.add_line('##INFO=<ID=END2,Number=1,Type=Integer,Description="End position of the structural variant on CHR2">')
    with pysam.VariantFile("intermediate2.vcf.gz", 'w', header=header) as out:
        for record in vcf:
            svtype = record.info.get('SVTYPE', None)
            if (svtype == 'BND' or svtype == 'CTX') and 'END2' not in record.info:
                record.info['END2'] = bnd_end_dict[record.id] if bnd_end_dict is not None \
                    else record.info.get('END2', record.stop)
            if svtype == 'BND' and 'CHR2' not in record.info:
                record.info['CHR2'] = record.chrom
            out.write(record)

CODE

    bcftools annotate -x INFO/MEMBERS intermediate2.vcf.gz -O z -o ~{prefix}.reformatted.vcf.gz
    tabix ~{prefix}.reformatted.vcf.gz

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
