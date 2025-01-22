version 1.0

import "Structs.wdl"

workflow Vapor {
  input {
    String prefix
    File bam_or_cram_file
    File bam_or_cram_index
    File bed_file
    String sample_id

    Boolean save_plots  # Control whether plots are final output

    File ref_fasta
    File ref_fai
    File ref_dict
    File contigs

    String vapor_docker
    String sv_base_mini_docker
    String sv_pipeline_docker

    RuntimeAttr? runtime_attr_vapor
    RuntimeAttr? runtime_attr_bcf2vcf
    RuntimeAttr? runtime_attr_vcf2bed
    RuntimeAttr? runtime_attr_split_vcf
    RuntimeAttr? runtime_attr_concat_beds
    RuntimeAttr? runtime_attr_LocalizeCram

    File? NONE_FILE_ # Create a null file - do not use this input
  }

  scatter (contig in read_lines(contigs)) {

    call PreprocessBedForVapor {
      input:
        prefix = "~{prefix}.~{contig}.preprocess",
        contig = contig,
        sample_to_extract = sample_id,
        bed_file = bed_file,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_split_vcf
    }

    call RunVaporWithCram {
      input:
        prefix = "~{prefix}.~{contig}",
        contig = contig,
        bam_or_cram_file = bam_or_cram_file,
        bam_or_cram_index = bam_or_cram_index,
        bed = PreprocessBedForVapor.contig_bed,
        ref_fasta = ref_fasta,
        ref_fai = ref_fai,
        ref_dict = ref_dict,
        vapor_docker = vapor_docker,
        runtime_attr_override = runtime_attr_vapor
    }
  }

  call ConcatVapor {
    input:
      shard_bed_files = RunVaporWithCram.vapor,
      shard_plots = RunVaporWithCram.vapor_plot,
      prefix = prefix,
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_concat_beds
  }

  output {
    File vapor_bed = ConcatVapor.merged_bed_file
    File? vapor_plots = if save_plots then ConcatVapor.merged_bed_plot else NONE_FILE_
  }
}

# extract specific contig from BED, and sites for sample if provided, and add SVLEN to INS if header contains SVLEN column
task PreprocessBedForVapor {
  input {
    String prefix
    String contig
    String? sample_to_extract
    File bed_file  # first 5 columns must be chrom, start, end, name, svtype (or Vapor description). if >5 columns, use header or assume samples is 6th. Need header & SVLEN column unless already appended to INS descriptions
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 3.75,
                               disk_gb: 10,
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File contig_bed = "~{prefix}.bed"
  }

  command <<<
    set -euo pipefail
    python /opt/sv-pipeline/scripts/preprocess_bed_for_vapor.py \
      --contig ~{contig} \
      --bed-in ~{bed_file} \
      --bed-out ~{prefix}.bed \
      ~{"-s " + sample_to_extract}
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

task RunVaporWithCram {
  input {
    String prefix
    String contig
    String bam_or_cram_file
    String bam_or_cram_index
    File bed
    File ref_fasta
    File ref_fai
    File ref_dict
    String vapor_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 15,
    disk_gb: 30,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File vapor = "~{prefix}.~{contig}.vapor.gz"
    File vapor_plot = "~{prefix}.~{contig}.tar.gz"
  }

  command <<<

    set -Eeuo pipefail

    # localize cram files
    export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
    samtools view -h -T ~{ref_fasta} -o ~{contig}.bam ~{bam_or_cram_file} ~{contig}
    samtools index ~{contig}.bam

    # run vapor
    mkdir ~{prefix}.~{contig}

    vapor bed \
      --sv-input ~{bed} \
      --output-path ~{prefix}.~{contig} \
      --output-file ~{prefix}.~{contig}.vapor \
      --reference ~{ref_fasta} \
      --PB-supp 0 \
      --pacbio-input ~{contig}.bam

    tar -czf ~{prefix}.~{contig}.tar.gz ~{prefix}.~{contig}
    bgzip ~{prefix}.~{contig}.vapor
  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: vapor_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

# Merge shards after Vapor
task ConcatVapor {
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
  RuntimeAttr runtime_default = object {
                                  mem_gb: 2.0 + compression_factor * input_size,
                                  disk_gb: ceil(10.0 + input_size * (2.0 + compression_factor)),
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

    zcat ~{shard_bed_files[0]} | head -n1 > header.txt
    # note head -n1 stops reading early and sends SIGPIPE to zcat,
    # so setting pipefail here would result in early termination

    # no more early stopping
    set -o pipefail

    while read SPLIT; do
      zcat $SPLIT | tail -n+2
    done < ~{write_lines(shard_bed_files)} \
      | sort -Vk1,1 -k2,2n -k3,3n \
      | cat header.txt - \
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