version 1.0

import "PESRClustering.wdl" as pesr
import "TasksClusterBatch.wdl" as tasks

workflow ClusterTloc {
  input {
    String prefix

    Array[File] manta_tloc_vcfs  # >1 batches' worth of manta_tloc VCFs

    Float max_af

    File ped_file

    # Reference
    File contig_list
    File reference_fasta
    File reference_fasta_fai
    File reference_dict
    String? chr_x
    String? chr_y
    File cytobands

    # PESR-based variant clustering
    Int? pesr_min_size
    File pesr_exclude_intervals
    Float pesr_interval_overlap
    Int pesr_breakend_window
    String? pesr_clustering_algorithm

    String gatk_docker
    String sv_base_mini_docker
    String sv_pipeline_docker
    String linux_docker

    Float? java_mem_fraction

    RuntimeAttr? runtime_attr_tar_files
    RuntimeAttr? runtime_attr_create_ploidy
    RuntimeAttr? runtime_attr_prepare_pesr_vcfs
    RuntimeAttr? runtime_attr_svcluster_manta_tloc
    RuntimeAttr? runtime_override_concat_vcfs_pesr
    RuntimeAttr? runtime_attr_gatk_to_svtk_vcf_pesr
    RuntimeAttr? runtime_attr_exclude_intervals_pesr
    RuntimeAttr? runtime_attr_select_rare_label_arms
  }

  call PreprocessAndTarTlocVcfs {
    input:
      files=manta_tloc_vcfs,
      prefix="~{prefix}.clustered.manta_tloc",
      linux_docker=linux_docker,
      runtime_attr_override=runtime_attr_tar_files
  }

  # TODO : properly set allosome ploidy, which creates problems in RDTest for allosomes at the moment
  call tasks.CreatePloidyTableFromPed {
    input:
      ped_file=ped_file,
      contig_list=contig_list,
      retain_female_chr_y=true,
      chr_x=chr_x,
      chr_y=chr_y,
      output_prefix="~{prefix}.ploidy",
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_attr_create_ploidy
  }

  call pesr.ClusterPESR {
    input:
      vcf_tar=PreprocessAndTarTlocVcfs.out,
      ploidy_table=CreatePloidyTableFromPed.out,
      batch=prefix,
      caller="manta_tloc",
      min_size=select_first([pesr_min_size, 50]),
      exclude_intervals=pesr_exclude_intervals,
      contig_list=contig_list,
      pesr_interval_overlap=pesr_interval_overlap,
      pesr_breakend_window=pesr_breakend_window,
      clustering_algorithm=pesr_clustering_algorithm,
      reference_fasta=reference_fasta,
      reference_fasta_fai=reference_fasta_fai,
      reference_dict=reference_dict,
      java_mem_fraction=java_mem_fraction,
      gatk_docker=gatk_docker,
      sv_base_mini_docker=sv_base_mini_docker,
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_prepare_pesr_vcfs=runtime_attr_prepare_pesr_vcfs,
      runtime_attr_svcluster=runtime_attr_svcluster_manta_tloc,
      runtime_override_concat_vcfs_pesr=runtime_override_concat_vcfs_pesr,
      runtime_attr_gatk_to_svtk_vcf=runtime_attr_gatk_to_svtk_vcf_pesr,
      runtime_attr_exclude_intervals_pesr=runtime_attr_exclude_intervals_pesr
  }

  call SelectRareAndLabelArms {
    input:
      prefix=prefix,
      clustered_manta_tloc_vcf=ClusterPESR.clustered_vcf,
      clustered_manta_tloc_vcf_index=ClusterPESR.clustered_vcf_index,
      cytobands=cytobands,
      max_af=max_af,
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_attr_select_rare_label_arms

  }

  output {
    File? tloc_vcf = SelectRareAndLabelArms.rare_witharms_vcf
    File? tloc_vcf_index = SelectRareAndLabelArms.rare_witharms_vcf_index
  }
}


task PreprocessAndTarTlocVcfs {
  input {
    String prefix
    Array[File] files
    String linux_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 1.0,
                               disk_gb: ceil(10 + 2 * size(files, "GB")),
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail
    mkdir to_tar
    while read file; do
      # extract CTX events from manta_tloc VCFs before clustering
      bcftools view -i "INFO/SVTYPE=='CTX" $file -O z -o to_tar/"$(basename $file)"
    done < ~{write_lines(files)}
    # tar individual sample vcfs for clustering
    tar czf ~{prefix}.tar.gz -C to_tar/ .
  >>>

  output {
    File out = "~{prefix}.tar.gz"
  }

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: linux_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}


task SelectRareAndLabelArms {
  input {
    String prefix
    File clustered_manta_tloc_vcf
    File clustered_manta_tloc_vcf_index
    Float max_af
    File cytobands
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 3.75,
                               disk_gb: ceil(10 + 2 * size(clustered_manta_tloc_vcf, "GB")),
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail
    # Select rare CTX only: apply max AF filter
    bcftools +fill-tags ~{clustered_manta_tloc_vcf} -- -t AC,AN,AF \
      | bcftools view --no-update -i 'AF<~{max_af}' -Oz -o ~{prefix}.rare.clustered_manta_tloc.vcf.gz

    python <<CODE
import pysam
import gzip

def _parse_bnd_ends(vcf_path):
  # read in ENDs without pysam to avoid overwriting immediately
  bnd_end_dict = dict()
  with gzip.open(vcf_path, 'rt') as f:
    for line in f:
      if line.startswith('#'):
        continue
      columns = line.split('\t', 8)
      vid = columns[2]
      info = columns[7]
      if 'SVTYPE=BND' not in info and 'SVTYPE=CTX' not in info:
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

def get_arms(record, cytobands):
  regionA = '{0}:{1}-{1}'.format(record.chrom, record.pos)
  regionB = '{0}:{1}-{1}'.format(record.info['CHR2'], record.info['END2'])

  def _get_arm(region):
    print(region)
    arm = next(cytobands.fetch(region))
    return arm.split()[3][0]

  return _get_arm(regionA), _get_arm(regionB)

bnd_end_dict = _parse_bnd_ends("~{prefix}.rare.clustered_manta_tloc.vcf.gz")
cytobands = pysam.TabixFile("~{cytobands}")
vcf = pysam.VariantFile("~{prefix}.rare.clustered_manta_tloc.vcf.gz")
vcf.header.add_line('##INFO=<ID=END2,Number=1,Type=Integer,Description="Second position">')
outvcf = pysam.VariantFile("~{prefix}.witharms.rare.clustered_manta_tloc.vcf.gz", 'w', header=vcf.header)
for record in vcf:
  # Set END2
  record.info['END2'] = bnd_end_dict[record.id]
  armA, armB = get_arms(record, cytobands)
  # get CTX arm information
  if armA == armB:
    record.info['CPX_TYPE'] = "CTX_PP/QQ"
  else:
    record.info['CPX_TYPE'] = "CTX_PQ/QP"
  outvcf.write(record)

CODE

    tabix ~{prefix}.witharms.rare.clustered_manta_tloc.vcf.gz
  >>>

  output {
    File rare_witharms_vcf = "~{prefix}.witharms.rare.clustered_manta_tloc.vcf.gz"
    File rare_witharms_vcf_index = "~{prefix}.witharms.rare.clustered_manta_tloc.vcf.gz.tbi"
  }

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
