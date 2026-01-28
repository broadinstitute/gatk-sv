version 1.0

workflow AnnotateAfForLrVcf {
  input {
    File vcf_gz
    String prefix
    Array[String] contig_list
    String docker_image
    String sv_pipeline_base_docker
    RuntimeAttr? runtime_attr_detect_contigs
    RuntimeAttr? runtime_attr_process_contig
    RuntimeAttr? runtime_attr_merge_vcfs
  }

  scatter (c in contig_list) {
    call ProcessContig {
      input:
        vcf_gz = vcf_gz,
        contig = c,
        docker_image = sv_pipeline_base_docker,
        runtime_attr_override = runtime_attr_process_contig
    }
  }

  call MergeVcfs {
    input:
      vcfs = ProcessContig.out_vcf,
      prefix = prefix,
      docker_image = sv_pipeline_base_docker,
      runtime_attr_override = runtime_attr_merge_vcfs
  }

  output {
    File final_vcf = MergeVcfs.merged_vcf
    File final_vcf_index = MergeVcfs.merged_vcf_index
  }
}

task ProcessContig {
  input {
    File vcf_gz
    String contig
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  # Define runtime attributes
  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 15,
    disk_gb: ceil(10 + size(vcf_gz, "GB") * 2),
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail

    # Extract contig
    zcat ~{vcf_gz} | \
    awk -v c="~{contig}" 'BEGIN{OFS="\t"} /^#/ {print; next} $1==c {print}' | \
    gzip > ~{contig}.vcf.gz

    # Run AC/AN/AF Python script
    python3 <<'PYTHON' ~{contig}.vcf.gz ~{contig}.out.vcf.gz
import sys
import gzip

vcf_in = sys.argv[1]
vcf_out = sys.argv[2]

def parse_gt(gt):
    if gt in ('.', './.', '.|.'):
        return None
    sep = '/' if '/' in gt else '|'
    alleles = gt.split(sep)
    return [int(a) for a in alleles if a!="." ]

with gzip.open(vcf_in, 'rt') as fin, gzip.open(vcf_out, 'wt') as fout:
    for line in fin:
        if line.startswith('##'):
            fout.write(line)
            continue
        if line.startswith('#CHROM'):
            fout.write('##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count">\n')
            fout.write('##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles">\n')
            fout.write('##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">\n')
            fout.write(line)
            continue
        fields = line.rstrip().split('\t')
        ref, alts = fields[3], fields[4]
        info = fields[7]
        format_fields = fields[8].split(':')
        samples = fields[9:]
        alt_alleles = alts.split(',')
        n_alt = len(alt_alleles)
        AC = [0]*n_alt
        AN = 0
        gt_index = format_fields.index('GT')
        for sample in samples:
            sample_fields = sample.split(':')
            gt = sample_fields[gt_index]
            alleles = parse_gt(gt)
            if alleles is None: continue
            for a in alleles:
                if a==0:
                    AN += 1
                elif 1 <= a <= n_alt:
                    AN += 1
                    AC[a-1] += 1
        AF = [ac/AN if AN>0 else 0 for ac in AC]
        info_items = [x for x in info.split(';') if not x.startswith(('AC=','AN=','AF='))]
        info_items.append(f"AC={','.join(map(str,AC))}")
        info_items.append(f"AN={AN}")
        info_items.append(f"AF={','.join(f'{x:.6g}' for x in AF)}")
        fields[7]=';'.join(info_items)
        fout.write('\t'.join(fields)+'\n')
PYTHON
  >>>

  output {
    File out_vcf = "~{contig}.out.vcf.gz"
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


task MergeVcfs {
  input {
    Array[File] vcfs
    String prefix
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  # Define runtime attributes
  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 15,
    disk_gb: ceil(10 + size(vcfs, "GB") * 2),
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail
    bcftools concat -a -Oz -o ~{prefix}.vcf.gz ~{sep=' ' vcfs}
    bcftools index ~{prefix}.vcf.gz
  >>>

  output {
    File merged_vcf = "~{prefix}.vcf.gz"
    File merged_vcf_index = "~{prefix}.vcf.gz.csi"
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
