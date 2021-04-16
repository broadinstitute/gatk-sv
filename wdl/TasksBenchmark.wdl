version 1.0

import "Structs.wdl"

# Merge shards after VaPoR
task ConcatVaPoR {
  input {
    Array[File] shard_bed_files
    Array[File] shard_plots
    String prefix
    Boolean? index_output
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  Boolean call_tabix = select_first([index_output, true])
  String output_file="~{prefix}.bed.gz"

  # when filtering/sorting/etc, memory usage will likely go up (much of the data will have to
  # be held in memory or disk while working, potentially in a form that takes up more space)
  Float input_size = size(shard_bed_files, "GB")
  Float compression_factor = 5.0
  Float base_disk_gb = 5.0
  Float base_mem_gb = 2.0
  RuntimeAttr runtime_default = object {
                                  mem_gb: base_mem_gb + compression_factor * input_size,
                                  disk_gb: ceil(base_disk_gb + input_size * (2.0 + compression_factor)),
                                  cpu_cores: 1,
                                  preemptible_tries: 3,
                                  max_retries: 1,
                                  boot_disk_gb: 10
                                }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_base_mini_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -eu

    # note head -n1 stops reading early and sends SIGPIPE to zcat,
    # so setting pipefail here would result in early termination

    # no more early stopping
    set -o pipefail

    while read SPLIT; do
    zcat $SPLIT | tail -n+2
      done < ~{write_lines(shard_bed_files)} \
      | sort -Vk1,1 -k2,2n -k3,3n \
      | bgzip -c \
      > ~{output_file}

    if ~{call_tabix}; then
      tabix -f -p bed ~{output_file}
    else
      touch ~{output_file}.tbi
    fi

    mkdir ~{prefix}.plots
    while read SPLIT; do
      tar zxvf $SPLIT -C ~{prefix}.plots/
    done < ~{write_lines(shard_plots)}

    tar -czf ~{prefix}.plots.tar.gz ~{prefix}.plots/
  >>>

  output {
    File merged_bed_file = output_file
    File merged_bed_plot = "~{prefix}.plots.tar.gz"
  }
}

#localize a specific contig of a bam/cram file
task LocalizeCram {
  input{
    String contig
    File ref_fasta
    File ref_fai
    File ref_dict
    String bam_or_cram_file
    String bam_or_cram_index
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 15, 
    disk_gb: 40,
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  Float mem_gb = select_first([runtime_attr.mem_gb, default_attr.mem_gb])
  Int java_mem_mb = ceil(mem_gb * 1000 * 0.8)

  output{
    File local_bam = "~{contig}.bam"
    File local_bai = "~{contig}.bam.bai"
  }

  command <<<
    set -Eeuo pipefail
    
    java -Xmx~{java_mem_mb}M -jar ${GATK_JAR}  PrintReads \
      -I ~{bam_or_cram_file} \
      -L ~{contig} \
      -O ~{contig}.bam \
      -R ~{ref_fasta}
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

task LocalizeCramRequestPay{
  input{
    String contig
    File ref_fasta
    File ref_fai
    File ref_dict
    String project_id
    String bam_or_cram_file
    String bam_or_cram_index
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 3.75, 
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  Float mem_gb = select_first([runtime_attr.mem_gb, default_attr.mem_gb])
  Int java_mem_mb = ceil(mem_gb * 1000 * 0.8)

  output{
    File local_bam = "~{contig}.bam"
    File local_bai = "~{contig}.bam.bai"
  }

  command <<<
    set -Eeuo pipefail
    
    java -Xmx~{java_mem_mb}M -jar ${GATK_JAR}  PrintReads \
      -I ~{bam_or_cram_file} \
      -L ~{contig} \
      -O ~{contig}.bam \
      -R ~{ref_fasta} \
      --gcs-project-for-requester-pays ~{project_id}
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

#extract specific contig from vcf
task SplitBed{
  input{
    String contig
    File? bed_file
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 3.75, 
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output{
    File contig_bed = "~{contig}.bed"
  }

  command <<<
    if [[ ~{bed_file} == *.gz ]] ;  then
      zcat ~{bed_file} | awk '{if ($1=="~{contig}") print}'  > ~{contig}.bed
    else
      awk '{if ($1=="~{contig}") print}' ~{bed_file} > ~{contig}.bed
    fi
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
task SplitVcf{
  input{
    String contig
    File? vcf_file
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 3.75, 
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output{
    File contig_vcf = "~{contig}.vcf.gz"
    File contig_vcf_index = "~{contig}.vcf.gz.tbi"
  }

  command <<<
    if [[ ~{vcf_file} == *.gz ]] ;  then
      tabix -f -p vcf ~{vcf_file}
      tabix -h ~{vcf_file} ~{contig} | bgzip > ~{contig}.vcf.gz
      tabix -p vcf ~{contig}.vcf.gz
    else
      bgzip ~{vcf_file}
      tabix -f -p vcf ~{vcf_file}.gz
      tabix -h ~{vcf_file}.gz ~{contig} | bgzip > ~{contig}.vcf.gz
      tabix -p vcf ~{contig}.vcf.gz
    fi
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

task vcf2bed{
  input{
    File vcf
    File? vcf_index
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 10, 
    disk_gb: 100,
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  String filename = basename(vcf, ".vcf.gz")

  output {
    File bed = "${filename}.bed"
  }

  command <<<

    set -Eeuo pipefail
    
    svtk vcf2bed -i SVTYPE -i SVLEN ~{vcf} tmp1.bed
    
    cat \
      <(awk '{if ($5=="DEL") print}' tmp1.bed | cut -f1-5)  \
      <(awk '{if ($5=="DUP") print}' tmp1.bed | cut -f1-5) \
      <(awk '{if ($5=="INV") print}' tmp1.bed | cut -f1-5) \
      > ~{filename}.bed

    paste -d '_' \
    <(awk '{if ($5=="INS") print}' tmp1.bed | cut -f1-5) \
    <(awk '{if ($5=="INS") print}' tmp1.bed | cut -f8) \
    >> ~{filename}.bed

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








