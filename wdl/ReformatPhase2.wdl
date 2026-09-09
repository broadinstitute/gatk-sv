version 1.0

import "Structs.wdl"

workflow ReformatPhase2 {
  input {
    File clustered_manta_vcf
    File clustered_manta_vcf_index

    File clustered_depth_vcf
    File clustered_depth_vcf_index

    File clustered_scramble_vcf
    File clustered_scramble_vcf_index

    File clustered_wham_vcf
    File clustered_wham_vcf_index

    String prefix
    String sv_pipeline_docker
  }

  Array[File] vcfs = [clustered_manta_vcf, clustered_depth_vcf, clustered_scramble_vcf, clustered_wham_vcf]
  Array[File] vcf_indexes = [clustered_manta_vcf_index, clustered_depth_vcf_index, clustered_scramble_vcf_index, clustered_wham_vcf_index]
  Array[String] labels = ['manta', 'depth', 'scramble', 'wham']

  scatter (i in range(length(vcfs))) {
    call ReformatVcf {
      input:
        vcf=vcfs[i],
        vcf_idx=vcf_indexes[i],
        prefix="~{prefix}.~{labels[i]}",
        sv_pipeline_docker=sv_pipeline_docker
    }
  }

  output {
    File reformatted_manta_vcf = ReformatVcf.reformatted_vcf[0]
    File reformatted_manta_vcf_idx = ReformatVcf.reformatted_vcf_idx[0]

    File reformatted_depth_vcf = ReformatVcf.reformatted_vcf[1]
    File reformatted_depth_vcf_idx = ReformatVcf.reformatted_vcf_idx[1]

    File reformatted_scramble_vcf = ReformatVcf.reformatted_vcf[2]
    File reformatted_scramble_vcf_idx = ReformatVcf.reformatted_vcf_idx[2]

    File reformatted_wham_vcf = ReformatVcf.reformatted_vcf[3]
    File reformatted_wham_vcf_idx = ReformatVcf.reformatted_vcf_idx[3]

    File original_manta_vcf = clustered_manta_vcf
    File original_manta_vcf_index = clustered_manta_vcf_index
    File original_depth_vcf = clustered_depth_vcf
    File original_depth_vcf_index = clustered_depth_vcf_index
    File original_scramble_vcf = clustered_scramble_vcf
    File original_scramble_vcf_index = clustered_scramble_vcf_index
    File original_wham_vcf = clustered_wham_vcf
    File original_wham_vcf_index = clustered_wham_vcf_index
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
    with pysam.VariantFile("~{prefix}.reformatted.vcf.gz", 'w', header=header) as out:
        for record in vcf:
            svtype = record.info.get('SVTYPE', None)
            if (svtype == 'BND' or svtype == 'CTX'):
                if 'END2' not in record.info:
                    record.info['END2'] = bnd_end_dict[record.id] if bnd_end_dict is not None \
                        else record.info.get('END2', record.stop)
                record.stop = record.pos
            if svtype == 'BND' and 'CHR2' not in record.info:
                record.info['CHR2'] = record.chrom
            out.write(record)

CODE

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
