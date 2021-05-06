## This WDL pipeline implements SV calling with DELLY2 by Tobias Rausch (https://github.com/dellytools/delly)
##
## Requirements/expectations :
## - Human whole-genome pair-end sequencing data in mapped BAM format

version 1.0

import "Structs.wdl"

workflow Delly {
  # Run Delly SV detection algorithm on whole genomes in array of
  #   bam or cram files.
  input {
    File bam_or_cram_file
    File? bam_or_cram_index
    String sample_id
    File reference_fasta
    File? reference_index
    File exclude_intervals_file
    Array[String]? sv_types
    Int? read_pairs
    String sv_base_mini_docker
    String delly_docker
    RuntimeAttr? runtime_attr_delly
    RuntimeAttr? runtime_attr_gather
  }
    
  parameter_meta {
    bam_or_cram_file: ".bam or .cram file to search for SVs. crams are preferable because they localise faster and use less disk."
    bam_or_cram_index: "[optional] Index for bam_or_cram_file. If not specified, index is assumed to be at bam_file_path + '.bai' or cram_file_path + '.crai'"
    sample_id: "sample name. Outputs will be sample_name + 'manta.vcf.gz' and sample_name + 'manta.vcf.gz.tbi'"
    reference_fasta: ".fasta file with reference used to align bam or cram file"
    reference_index: "[optional] If omitted, the WDL will look for an index by appending .fai to the .fasta file"
    exclude_intervals_file: "text file with lines specifying genomic intervals where SVs should not be called.Each line in (tab-separated) format: 'contig   begin   end interval-type' or 'contig'"
    sv_types: "[optional] array of event types for Delly to search for. Defaults to ['DEL', 'DUP', 'INV']."
    read_pairs: "[optional] Value of READ_PAIRS obtained from CollectInsertSizeMetrics, used for optimal estimate of VM memory needs."
  }
  
  meta {
      author: "Ted Brookings"
      email: "tbrookin@broadinstitute.org"
  }

  Array[String] event_types = select_first([sv_types, ["DEL", "DUP", "INV"]])

  scatter (event_type in event_types) {
    call RunDelly {
      input:
        bam_or_cram_file = bam_or_cram_file,
        bam_or_cram_index = bam_or_cram_index,
        sample_id = sample_id,
        reference_fasta = reference_fasta,
        reference_index = reference_index,
        exclude_intervals_file = exclude_intervals_file,
        event_type = event_type,
        read_pairs = read_pairs,
        delly_docker = delly_docker,
        runtime_attr_override = runtime_attr_delly
    }
  }

  call GatherBCFs {
    input:
      bcfs = RunDelly.bcf,
      sample_id = sample_id,
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_gather
  }

  output {
    File vcf = GatherBCFs.vcf
    File index = GatherBCFs.index
  }
}

task RunDelly {
  input {
    File bam_or_cram_file
    File? bam_or_cram_index
    String sample_id
    File reference_fasta
    File? reference_index
    File exclude_intervals_file
    String event_type
    Int? read_pairs
    String delly_docker
    RuntimeAttr? runtime_attr_override
  }

  Boolean is_bam = basename(bam_or_cram_file, ".bam") + ".bam" == basename(bam_or_cram_file)
  String index_path = if is_bam then bam_or_cram_file + ".bai" else bam_or_cram_file + ".crai"
  File bam_or_cram_index_file = select_first([bam_or_cram_index, index_path])
  File ref_index = select_first([reference_index, reference_fasta + ".fai"])

  # ensure there's sufficient disk space
  Float disk_overhead = 10.0
  Float bam_or_cram_size = size(bam_or_cram_file, "GiB")
  Float bam_or_cram_index_size = size(bam_or_cram_index_file, "GiB")
  Float ref_size = size(reference_fasta, "GiB")
  Float ref_index_size = size(ref_index, "GiB")
  Float exclude_list_size = size(exclude_intervals_file, "GiB")
  Int vm_disk_size = ceil(bam_or_cram_size + bam_or_cram_index_size + ref_size + ref_index_size + exclude_list_size + disk_overhead)
  
  # ensure there's sufficient memory. These coefficients are obtained
  # via linear regression to test data, assuming uncertainty in memory
  # usage proportional to bam or cram size, and adding 3 standard
  # deviations to best fit estimate of memory usage.
  Float mem_per_bam_size = 0.03937
  Float mem_bam_offset = 4.4239
  Float mem_per_cram_size = 0.08579
  Float mem_cram_offset = 4.8633
  Float mem_per_read_pairs = "9.745e-9"
  Float mem_read_pairs_offset = 1.947
  Float mem_size_gb =
      if defined(read_pairs) then
          mem_read_pairs_offset + mem_per_read_pairs * select_first([read_pairs])
      else if is_bam then
          mem_per_bam_size * bam_or_cram_size + mem_bam_offset
      else
          mem_per_cram_size * bam_or_cram_size + mem_cram_offset

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: mem_size_gb, 
    disk_gb: vm_disk_size,
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File bcf = "${sample_id}.delly.${event_type}.bcf"
  }
  command <<<

    set -Eeuo pipefail
    
    BCF="~{sample_id}.delly.~{event_type}.bcf"
    echo "Running delly on event_type=~{event_type}, output=$BCF"
    delly call \
      -t ~{event_type} \
      -g "~{reference_fasta}" \
      -x "~{exclude_intervals_file}" \
      -o "$BCF" \
      -n \
      "~{bam_or_cram_file}"
    
  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: delly_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

}

task GatherBCFs {
  input {
    Array[File] bcfs
    String sample_id
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  File first_bcf = bcfs[0]
  Int num_bcfs = length(bcfs)

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 1.7, 
    disk_gb: 20,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File vcf = "${sample_id}.delly.vcf.gz"
    File index = "${sample_id}.delly.vcf.gz.tbi"
  }
  command <<<

    set -Eeuo pipefail

    BCFS="~{sep='\n' bcfs}"
    echo "Extracting vcf header"
    bcftools view ~{first_bcf} | grep "^#" > header.vcf
    VCFS="header.vcf"
    for ((BCF_IND=1; BCF_IND <= ~{num_bcfs}; BCF_IND++)); do
      echo "Processing BCF $BCF_IND"
      BCF=$(echo "$BCFS" | sed "$BCF_IND""q;d")
      VCF=$(echo "$BCF" | sed -e 's/.bcf$/.vcf/')
      echo "Converting $BCF to $VCF, skipping header"
      bcftools view  $BCF | grep -v "^#" > $VCF
      VCFS="$VCFS $VCF"
    done

    VCF_OUT="~{sample_id}.delly.vcf.gz"
    echo "Concatenating vcfs into $VCF_OUT"
    cat $VCFS \
      | vcf-sort -c \
      | bcftools reheader -s <(echo "~{sample_id}") \
      | bgzip -c \
      > $VCF_OUT
    echo "Indexing $VCF_OUT"
    tabix "$VCF_OUT"
    
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

