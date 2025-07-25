version 1.0

workflow BedToVcfPavPerGenome {
  input {
    String sample_id

    File ref_fasta

    File snv_bed_hap1
    File ins_bed_hap1
    File del_bed_hap1

    File snv_bed_hap2
    File ins_bed_hap2
    File del_bed_hap2

    File bed2vcf_snv
    File bed2vcf_ins
    File bed2vcf_del
    File merge_script
  }

  call BedToVCF as snv_hap1{
    input:
      bed_file=snv_bed_hap1,
      script=bed2vcf_snv,
      ref_fasta=ref_fasta
  }

  call BedToVCF as ins_hap1 {
    input:
      bed_file=ins_bed_hap1,
      script=bed2vcf_ins,
      ref_fasta=ref_fasta
  }

  call BedToVCF as del_hap1 {
    input:
      bed_file=del_bed_hap1,
      script=bed2vcf_del,
      ref_fasta=ref_fasta
  }

  call BedToVCF as snv_hap2 {
    input:
      bed_file=snv_bed_hap2,
      script=bed2vcf_snv,
      ref_fasta=ref_fasta
  }

  call BedToVCF as ins_hap2 {
    input:
      bed_file=ins_bed_hap2,
      script=bed2vcf_ins,
      ref_fasta=ref_fasta
  }

  call BedToVCF as del_hap2 {
    input:
      bed_file=del_bed_hap2,
      script=bed2vcf_del,
      ref_fasta=ref_fasta
  }

  call ConcatVCF as concat_hap1{
    input:
      vcf_1=snv_hap1.output_vcf,
      vcf_2=ins_hap1.output_vcf,
      vcf_3=del_hap1.output_vcf,

      vcf_idx_1=snv_hap1.output_vcf_idx,
      vcf_idx_2=ins_hap1.output_vcf_idx,
      vcf_idx_3=del_hap1.output_vcf_idx,

      output_name="~{sample_id}_hap1.vcf.gz"
  }

  call ConcatVCF as concat_hap2 {
    input:
      vcf_1=snv_hap2.output_vcf,
      vcf_2=ins_hap2.output_vcf,
      vcf_3=del_hap2.output_vcf,

      vcf_idx_1=snv_hap2.output_vcf_idx,
      vcf_idx_2=ins_hap2.output_vcf_idx,
      vcf_idx_3=del_hap2.output_vcf_idx,

      output_name="~{sample_id}_hap2.vcf.gz"
  }

  call MergeHaps {
    input:
      hap1_vcf=concat_hap1.output_vcf,
      hap1_vcf_idx=concat_hap1.output_vcf_idx,

      hap2_vcf=concat_hap2.output_vcf,
      hap2_vcf_idx=concat_hap2.output_vcf_idx,

      merge_script=merge_script,
      output_name="~{sample_id}.vcf.gz"
  }


  output{
        File output_vcf = MergeHaps.merged_vcf
        File output_vcf_idx = MergeHaps.merged_vcf_idx
  }
}

task BedToVCF {
  input {
    File bed_file
    File script
    File ref_fasta
  }

  String prefix = basename(bed_file, ".bed.gz")
  command {
    python ~{script} ~{bed_file} ~{ref_fasta}
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 10 + ceil(size(bed_file, "GiB")*5),
    disk_gb: 15 + ceil(size(bed_file, "GiB")*5),
    boot_disk_gb: 10,
    preemptible_tries: 1,
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

  output {
    File output_vcf = "~{prefix}.vcf.gz"
    File output_vcf_idx =  "~{prefix}.vcf.gz.tbi"
  }
}

task ConcatVCF {
  input {
    File vcf_1
    File vcf_2
    File vcf_3
    File vcf_idx_1
    File vcf_idx_2
    File vcf_idx_3
    String output_name
  }

  command {
    bcftools concat -a -Oz -o ~{output_name} ~{vcf1} ~{vcf2} ~{vcf3}
    tabix -p vcf ~{output_name}
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 10 + ceil(size(vcf1, "GiB")*5),
    disk_gb: 15 + ceil(size(vcf1, "GiB")*5),
    boot_disk_gb: 10,
    preemptible_tries: 1,
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

  output {
    File output_vcf = output_name
    File output_vcf_idx = output_name + ".tbi"
  }
}

task MergeHaps {
  input {
    File hap1_vcf
    File hap2_vcf
    File merge_script
    String output_name
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  String prefix = basename(output_name, ".vcf.gz")
  command {
    python ~{merge_script} ~{hap1_vcf} ~{hap2_vcf} ~{output_name}
    bcftools sort -Oz -o "~{prefix}.sorted.vcf.gz" ~{output_name}
    tabix -p vcf "~{prefix}.sorted.vcf.gz"

  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 10 + ceil(size(hap1_vcf, "GiB")*5),
    disk_gb: 15 + ceil(size(hap1_vcf, "GiB")*5),
    boot_disk_gb: 10,
    preemptible_tries: 1,
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

  output {
    File merged_vcf = "~{prefix}.sorted.vcf.gz"
    File merged_vcf_idx = "~{prefix}.sorted.vcf.gz.tbi"
  }

}