version 1.0

import "Structs.wdl"
import "ShardedAnnotateVcf.wdl" as sharded_annotate_vcf
import "TasksMakeCohortVcf.wdl" as MiniTasks
import "FilterGenotypes.wdl" as filter


workflow PreprocessForFederation {
  input {
    Array[File] vcfs
    File contig_list
    String prefix

    File sample_keep_list
    File? sample_id_rename_map

    File par_bed
    File sample_pop_assignments
    File ped_file
    File ploidy_table
    Int sv_per_shard

    Boolean is_gnomad
    Boolean is_aou

    String header_drop_fields

    String sv_pipeline_docker
    String sv_base_mini_docker
  }

  Array[String] contigs = read_lines(contig_list)

  scatter (i in range(length(contigs))) {
    call sharded_annotate_vcf.ShardedAnnotateVcf {
      input:
        vcf = vcfs[i],
        vcf_idx = vcfs[i] + ".tbi",
        prefix = prefix,
        contig = contigs[i],
        par_bed = par_bed,
        ped_file = ped_file,
        ploidy_table = ploidy_table,
        sv_per_shard = sv_per_shard,
        sample_keep_list = sample_keep_list,
        sample_id_rename_map = sample_id_rename_map,
        sample_pop_assignments=sample_pop_assignments,
        sv_pipeline_docker=sv_pipeline_docker,
        sv_base_mini_docker=sv_base_mini_docker
      }
  }

  # ShardedAnnotateVcf.sharded_annotated_vcf is is an Array[Array[File]] with one inner Array[File] of shards per contig
  Array[File] vcfs_for_concatenation = flatten(ShardedAnnotateVcf.sharded_annotated_vcf)
  Array[File] vcf_idxs_for_concatenation = flatten(ShardedAnnotateVcf.sharded_annotated_vcf_idx)

  call MiniTasks.ConcatVcfs {
    input:
      vcfs = vcfs_for_concatenation,
      vcfs_idx=vcf_idxs_for_concatenation,
      naive = true,
      sites_only = true,
      outfile_prefix = "~{prefix}.annotated",
      sv_base_mini_docker = sv_base_mini_docker
  }

  if (is_gnomad) {
    call MoveInsToIntervals {
      input:
        vcf = ConcatVcfs.concat_vcf,
        vcf_idx = ConcatVcfs.concat_vcf_idx,
        sv_pipeline_docker = sv_pipeline_docker
    }

    call FixEnds {
      input:
        vcf = MoveInsToIntervals.ins_intervals_vcf,
        vcf_idx = MoveInsToIntervals.ins_intervals_vcf_idx,
        sv_pipeline_docker = sv_pipeline_docker
    }
  }

  if (is_aou) {
    call MoveMultiallelic {
      input:
        vcf = ConcatVcfs.concat_vcf,
        vcf_idx = ConcatVcfs.concat_vcf_idx,
        sv_pipeline_docker = sv_pipeline_docker
    }
  }

  call filter.SanitizeHeader {
    input:
      vcf = select_first([MoveMultiallelic.multiallelic_info_vcf, FixEnds.fix_ends_vcf, ConcatVcfs.concat_vcf]),
      vcf_index = select_first([MoveMultiallelic.multiallelic_info_vcf_idx, FixEnds.fix_ends_vcf_idx, ConcatVcfs.concat_vcf_idx]),
      drop_fields = header_drop_fields,
      prefix = "~{prefix}.sanitized",
      sv_pipeline_docker = sv_pipeline_docker
  }


  output {
    File preprocessed_vcf = SanitizeHeader.out
    File preprocessed_vcf_idx = SanitizeHeader.out_index
  }
}


task MoveMultiallelic {
  input {
    File vcf
    File vcf_idx
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  String output_file = basename(vcf, ".vcf.gz") + ".multiallelic_info.vcf.gz"

  RuntimeAttr runtime_default = object {
    mem_gb: 3.75,
    disk_gb: ceil(10.0 + size(vcf, "GB") * 2.0),
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

with pysam.VariantFile("~{vcf}", 'r') as vcf:
  header = vcf.header
  header.add_line("##INFO=<ID=MULTIALLELIC,Number=0,Type=Flag,Description=\"Multiallelic site\">")
  with pysam.VariantFile("~{output_file}", 'w', header=header) as out:
    for record in vcf:
      if 'MULTIALLELIC' in record.filter:
        filters = [x for x in record.filter if x != "MULTIALLELIC"]
        record.filter.clear()
        for filt in filters:
          record.filter.add(filt)
        record.info['MULTIALLELIC'] = True
      out.write(record)
CODE

    tabix -p vcf "~{output_file}"
  >>>

  output {
    File multiallelic_info_vcf = output_file
    File multiallelic_info_vcf_idx = output_file + ".tbi"
  }
}





task FixEnds {
  input {
    File vcf
    File vcf_idx
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  String output_file = basename(vcf, ".vcf.gz") + ".ins_to_intervals.vcf.gz"

  RuntimeAttr runtime_default = object {
    mem_gb: 3.75,
    disk_gb: ceil(10.0 + size(vcf, "GB") * 2.0),
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
import gzip


def _parse_bnd_ends(vcf_path):
    """
    Since pysam automatically changes invalid END fields (i.e. when less than the start position), they must
    be parsed manually.

    Parameters
    ----------
    vcf_path: Text
        input vcf path

    Returns
    -------
    header: Dict[Text, int]
        map from variant ID to END position
    """
    bnd_end_dict = dict()
    with gzip.open(vcf_path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            columns = line.split('\t', 8)
            vid = columns[2]
            info = columns[7]
            if 'SVTYPE=BND' not in info and 'SVTYPE=CTX' not in info and 'SVTYPE=CPX' not in info:
                continue
            info_tokens = info.split(';')
            end_field_list = [x for x in info_tokens if x.startswith("END=")]
            if len(end_field_list) > 0:
                end = int(end_field_list[0].replace("END=", ""))
            else:
                # Special case where END and POS happen to be equal
                end = int(columns[1])
            bnd_end_dict[vid] = end
    return bnd_end_dict


bnd_end_dict = _parse_bnd_ends("~{vcf}")
with pysam.VariantFile("~{vcf}", 'r') as inp:
  with pysam.VariantFile("~{output_file}", 'w', header=inp.header) as out:
    for record in inp:

      if record.info['SVTYPE'] in ['BND', 'CTX'] and 'END2' not in record.info:
        if record.chrom != record.info['CHR2']:
          print(f"Missing END2 but CHR2 != CHROM: {record.id}")  # this would be problematic so print
        record.info['END2'] = bnd_end_dict[record.id]  # only update if missing - END2 already set for most BNDs here. use non-pysam END to be safe
        if bnd_end_dict[record.id] < record.start:
          record.stop = record.start + 1  # just in case
      elif record.info['SVTYPE'] == 'CPX':
        if bnd_end_dict[record.id] < record.start:
          print(f"END<POS but has END2: {record.id}, {record.info['CPX_TYPE']}")  # if just insertion-type CPX then let it get updated by pysam
        if 'CHR2' in record.info:
          record.info.pop('CHR2')
        if 'END2' in record.info:
          record.info.pop('END2')
      out.write(record)

CODE

    tabix -p vcf "~{output_file}"
  >>>

  output {
    File fix_ends_vcf = output_file
    File fix_ends_vcf_idx = output_file + ".tbi"
  }
}


task MoveInsToIntervals {
  input {
    File vcf
    File vcf_idx
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  String output_file = basename(vcf, ".vcf.gz") + ".ins_to_intervals.vcf.gz"

  RuntimeAttr runtime_default = object {
    mem_gb: 3.75,
    disk_gb: ceil(10.0 + size(vcf, "GB") * 2.0),
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
with pysam.VariantFile("~{vcf}", 'r') as vcf:
  with pysam.VariantFile("~{output_file}", 'w', header=vcf.header) as out:
    for record in vcf:
      if 'CPX_TYPE' in record.info and record.info['CPX_TYPE'] == 'INS_iDEL':
        lst_intervals = list(record.info['CPX_INTERVALS'])
        lst_intervals.append(record.info['SOURCE'])
        record.info.pop('SOURCE')
        record.info['CPX_INTERVALS'] = tuple(lst_intervals)
      out.write(record)
CODE

    tabix -p vcf "~{output_file}"
  >>>

  output {
    File ins_intervals_vcf = output_file
    File ins_intervals_vcf_idx = output_file + ".tbi"
  }
}


