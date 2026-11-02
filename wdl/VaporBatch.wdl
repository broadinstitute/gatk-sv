version 1.0

import "Structs.wdl"
import "Vapor.wdl" as vapor_bed

workflow VaporBatch {
  input {
    Array[String] samples
    Array[File] bam_or_cram_files
    Array[File] bam_or_cram_indexes
    File bed_file  # Multi-sample bed file, generated with svtk vcf2bed (e.g. as in MainVcfQc)

    File ref_fasta
    File ref_fai
    File ref_dict
    File contig_list

    String vapor_docker
    String sv_base_mini_docker
    String sv_pipeline_docker

    Boolean save_plots = false

    RuntimeAttr? runtime_attr_split_bed_by_contig
    RuntimeAttr? runtime_attr_subset_sample
    RuntimeAttr? runtime_attr_vcf_to_bed
    RuntimeAttr? runtime_attr_vapor
    RuntimeAttr? runtime_attr_split_vcf
    RuntimeAttr? runtime_attr_concat_beds
  }

  # Split the multi-sample bed by contig once, so each per-sample Vapor shard
  # only localizes the records for its own contig.
  call SplitBedByContig {
    input:
      bed_file = bed_file,
      contig_list = contig_list,
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_split_bed_by_contig
  }

  scatter (i in range(length(bam_or_cram_files))) {
    call vapor_bed.Vapor {
      input:
        bam_or_cram_file = bam_or_cram_files[i],
        bam_or_cram_index = bam_or_cram_indexes[i],
        bed_files = SplitBedByContig.contig_beds,
        sample_id = samples[i],
        ref_fasta = ref_fasta,
        ref_fai = ref_fai,
        ref_dict = ref_dict,
        contig_list = contig_list,
        save_plots = save_plots,
        vapor_docker = vapor_docker,
        sv_base_mini_docker = sv_base_mini_docker,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_vapor = runtime_attr_vapor,
        runtime_attr_subset_sample = runtime_attr_subset_sample,
        runtime_attr_vcf_to_bed = runtime_attr_vcf_to_bed,
        runtime_attr_split_vcf = runtime_attr_split_vcf,
        runtime_attr_concat_beds = runtime_attr_concat_beds
    }
  }
  output {
    Array[File] vapor_batch_beds = Vapor.vapor_bed
    Array[File?] vapor_batch_plots = Vapor.vapor_plots
  }
}

# Split a (optionally gzipped) bed into one gzipped bed per contig, in the order of contig_list.
# Header lines (starting with '#') are copied to every output.
task SplitBedByContig {
  input {
    File bed_file
    File contig_list
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 3.75,
                               disk_gb: 10 + ceil(size(bed_file, "GiB") * 12),
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    Array[File] contig_beds = read_lines("contig_beds.list")
  }

  command <<<
    set -euo pipefail

    mkdir -p split
    if [[ ~{bed_file} == *.gz ]]; then
      zcat ~{bed_file} > in.bed
    else
      cp ~{bed_file} in.bed
    fi

    grep '^#' in.bed > header.txt || true

    i=0
    while read -r contig; do
      [[ -z "${contig}" ]] && continue
      out="split/${i}.${contig}.bed"
      cp header.txt "${out}"
      i=$((i+1))
    done < ~{contig_list}

    # One pass over the bed, appending each record to its contig's file
    awk -F'\t' -v OFS='\t' '
      FNR==NR { idx[$1] = FNR-1; next }
      /^#/ { next }
      ($1 in idx) { print >> ("split/" idx[$1] "." $1 ".bed") }
    ' <(grep -v '^[[:space:]]*$' ~{contig_list}) in.bed

    i=0
    : > contig_beds.list
    while read -r contig; do
      [[ -z "${contig}" ]] && continue
      gzip "split/${i}.${contig}.bed"
      echo "split/${i}.${contig}.bed.gz" >> contig_beds.list
      i=$((i+1))
    done < ~{contig_list}
  >>>

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_base_mini_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}
