version 1.0

import "Structs.wdl"

# Run Whamg SV detection algorithm on whole genome in bam or cram file
# If include_list is provided, run whamg on explicitly included subset of genome,
#   concatenate, then fix output vcf headers to contain good contig length data.

workflow Whamg {
  input {
    File bam_or_cram_file
    File? bam_or_cram_index
    String sample_id
    File reference_fasta
    File? reference_index
    File include_bed_file
    File primary_contigs_list
    Int? pf_reads_improper_pairs
    Float? pct_exc_total
    String wham_docker
    RuntimeAttr? runtime_attr_wham
  }

  parameter_meta {
    bam_or_cram_file: ".bam or .cram file to search for SVs. bams are preferable, crams will be converted to bams."
    bam_or_cram_index: "[optional] associated index file. If omitted, the WDL will look for an index file by appending .bai/.crai to the .bam/.cram file."
    sample_id: "sample name. Outputs will be sample_name + '.wham.vcf.gz' and sample_name + '.wham.vcf.gz.tbi'."
    reference_fasta: ".fasta file with reference used to align bam or cram file."
    reference_index: "[optional] reference index file. If omitted, the WDL will look for an index by appending .fai to the .fasta file."
    include_bed_file: "bed file with regions that wham will evaluate."
    primary_contigs_list: "list of canonical contigs in the reference, for estimating read stats."
    pf_reads_improper_pairs: "[optional] value of PF_READS_IMPROPER_PAIRS obtained from CollectAlignmentSummaryMetrics, used for optimal estimate of VM memory needs."
    pct_exc_total: "[optional] value of PCT_EXC_TOTAL obtained from CollectWgsMetrics, used for optimal estimate of VM memory needs."
  }

  Boolean is_bam = basename(bam_or_cram_file, ".bam") + ".bam" == basename(bam_or_cram_file)
  File reference_index_ = select_first([reference_index, reference_fasta + ".fai"])

  if (is_bam) {
    call RunWhamgOnBam {
      input:
        bam_file = bam_or_cram_file,
        bam_index = select_first([bam_or_cram_index, bam_or_cram_file + ".bai"]),
        reference_fasta = reference_fasta,
        reference_index = reference_index_,
        sample_id = sample_id,
        include_bed_file = include_bed_file,
        primary_contigs_list = primary_contigs_list,
        pf_reads_improper_pairs = pf_reads_improper_pairs,
        pct_exc_total = pct_exc_total,
        wham_docker = wham_docker,
        runtime_attr_override = runtime_attr_wham
    }
  }
  if (!is_bam) {
    call RunWhamgOnCram {
      input:
        cram_file = bam_or_cram_file,
        cram_index = select_first([bam_or_cram_index, bam_or_cram_file + ".crai"]),
        reference_fasta = reference_fasta,
        reference_index = reference_index_,
        sample_id = sample_id,
        include_bed_file = include_bed_file,
        primary_contigs_list = primary_contigs_list,
        pf_reads_improper_pairs = pf_reads_improper_pairs,
        pct_exc_total = pct_exc_total,
        wham_docker = wham_docker,
        runtime_attr_override = runtime_attr_wham
    }
  }

  output {
    File index = select_first([RunWhamgOnBam.index, RunWhamgOnCram.index])
    File vcf = select_first([RunWhamgOnBam.vcf, RunWhamgOnCram.vcf])
  }
}

task RunWhamgOnBam {
  input {
    File bam_file
    File bam_index
    File reference_fasta
    File reference_index
    String sample_id
    File include_bed_file
    File primary_contigs_list
    Int? pf_reads_improper_pairs
    Float? pct_exc_total
    String wham_docker
    RuntimeAttr? runtime_attr_override
  }

  Array[String] chr_list = read_lines(primary_contigs_list)

  # Calculate default disk size
  Float bam_size = size(bam_file, "GiB")
  Int vm_disk_size = ceil(bam_size + size(reference_fasta, "GiB") + 5.0)

  # Ensure there's sufficient memory. Use picard metrics if
  # available, otherwise estimate based on bam size
  Float mem_per_bam_size = 0.0
  Float mem_bam_offset = 4.9
  Float mem_per_improper_pairs = "7.872e-8"
  Float mem_per_exc = 1.511
  Float mem_metrics_offset = 2.639
  Float mem_size_gb =
    if defined(pf_reads_improper_pairs) && defined(pct_exc_total) then
      mem_metrics_offset + mem_per_improper_pairs * select_first([pf_reads_improper_pairs]) + mem_per_exc * select_first([pct_exc_total])
    else
      mem_per_bam_size * bam_size + mem_bam_offset

  RuntimeAttr default_attr = object {
    mem_gb: mem_size_gb,
    cpu_cores: 4,
    disk_gb: vm_disk_size,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  Int cpu_cores = select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])

  output {
    File vcf = "${sample_id}.wham.vcf.gz"
    File index = "${sample_id}.wham.vcf.gz.tbi"
  }

  command <<<
    set -euo pipefail

    # print some info that may be useful for debugging
    df -h
    echo "whamg $(./whamg 2>&1 | grep Version)"

    # ensure that index files are present in appropriate locations
    if [ ! -e "~{bam_file}.bai" ]; then
      ln -s "~{bam_index}" "~{bam_file}.bai"
    fi
    if [ ! -e '~{basename(bam_file,".bam")}.bai' ]; then
      ln -s "~{bam_index}" '~{basename(bam_file,".bam")}.bai'
    fi
    if [ ! -e "~{reference_fasta}.fai" ]; then
      ln -s "~{reference_index}" "~{reference_fasta}.fai"
    fi

    # run whamg on all specified intervals
    mkdir tmpVcfs
    cd tmpVcfs
    awk 'BEGIN{FS=OFS="\t"}{printf("%07d\t%s\n",NR,$1":"$2"-"$3)}' ~{include_bed_file} |\
      while read -r line interval; do
        vcfFile="$line.wham.vcf.gz"
        whamg \
            -c "~{sep="," chr_list}" \
            -x ~{cpu_cores} \
            -a "~{reference_fasta}" \
            -f "~{bam_file}" \
            -r $interval \
          | bgzip -c > $vcfFile
        bcftools index -t $vcfFile
      done

    # We need to update both the VCF sample ID and the TAGS INFO field in the WHAM output VCFs.
    # WHAM uses both to store the sample identifier, and by default uses the SM identifier from the BAM file.
    # We need to update both to reflect the potentially-renamed sample identifier used by the pipeline (sample_id) --
    # svtk standardize_vcf uses the TAGS field to identify the sample for WHAM VCFs.
    ls -1 *.wham.vcf.gz > vcf.list
    bcftools concat -a -Ov -f vcf.list | \
      sed -e 's/^#CHROM\t.*/#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t~{sample_id}/' -e 's/;TAGS=[^;]*;/;TAGS=~{sample_id};/' | \
      bgzip -c > "../~{sample_id}.wham.vcf.gz"
    cd ..
    bcftools index -t "~{sample_id}.wham.vcf.gz"

    df -h
    ls -l
  >>>
  runtime {
    cpu: cpu_cores
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " SSD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: wham_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

task RunWhamgOnCram {
  input {
    File cram_file
    File cram_index
    File reference_fasta
    File reference_index
    String sample_id
    File include_bed_file
    File primary_contigs_list
    Int? pf_reads_improper_pairs
    Float? pct_exc_total
    String wham_docker
    RuntimeAttr? runtime_attr_override
  }

  parameter_meta {
    cram_file: {
      localization_optional: true
    }
    cram_index: {
      localization_optional: true
    }
  }

  Array[String] chr_list = read_lines(primary_contigs_list)

  # Calculate default disk size
  Float cram_size = size(cram_file, "GiB")
  Float bam_size = 3.5 * cram_size
  Int vm_disk_size = ceil(cram_size + bam_size + size(reference_fasta, "GiB") + 5.0)

  # Ensure there's sufficient memory. Use picard metrics if
  # available, otherwise estimate based on bam size
  Float mem_per_bam_size = 0.0
  Float mem_bam_offset = 5.4
  Float mem_per_improper_pairs = "7.872e-8"
  Float mem_per_exc = 1.511
  Float mem_metrics_offset = 2.639
  Float mem_size_gb =
    if defined(pf_reads_improper_pairs) && defined(pct_exc_total) then
      mem_metrics_offset + mem_per_improper_pairs * select_first([pf_reads_improper_pairs]) + mem_per_exc * select_first([pct_exc_total])
    else
      mem_per_bam_size * bam_size + mem_bam_offset

  RuntimeAttr default_attr = object {
    mem_gb: mem_size_gb,
    cpu_cores: 4,
    disk_gb: vm_disk_size,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  Int cpu_cores = select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])

  output {
    File vcf = "${sample_id}.wham.vcf.gz"
    File index = "${sample_id}.wham.vcf.gz.tbi"
  }
  command <<<
    set -euo pipefail

    # print some info that may be useful for debugging
    df -h
    echo "whamg $(whamg 2>&1 | grep Version)"

    # necessary for getting permission to read from google bucket directly
    export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`

    # covert cram to bam
    samtools view -b1@ ~{cpu_cores} -T "~{reference_fasta}" "~{cram_file}" > sample.bam

    # index bam file
    samtools index -@ ~{cpu_cores} sample.bam

    # ensure that index files are present in appropriate locations
    ln -s sample.bam.bai sample.bai
    if [ ! -e "~{reference_fasta}.fai" ]; then
      ln -s "~{reference_index}" "~{reference_fasta}.fai"
    fi

    # run whamg on all specified intervals
    mkdir tmpVcfs
    cd tmpVcfs
    awk 'BEGIN{FS=OFS="\t"}{printf("%07d\t%s\n",NR,$1":"$2"-"$3)}' ~{include_bed_file} |\
      while read -r line interval; do
        vcfFile="$line.wham.vcf.gz"
        whamg \
            -c "~{sep="," chr_list}" \
            -x ~{cpu_cores} \
            -a "~{reference_fasta}" \
            -f ../sample.bam \
            -r $interval \
          | bgzip -c > $vcfFile
        bcftools index -t $vcfFile
      done

    # We need to update both the VCF sample ID and the TAGS INFO field in the WHAM output VCFs.
    # WHAM uses both to store the sample identifier, and by default uses the SM identifier from the BAM file.
    # We need to update both to reflect the potentially-renamed sample identifier used by the pipeline (sample_id) --
    # svtk standardize_vcf uses the TAGS field to identify the sample for WHAM VCFs.
    ls -1 *.wham.vcf.gz > vcf.list
    bcftools concat -a -Ov -f vcf.list | \
      sed -e 's/^#CHROM\t.*/#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t~{sample_id}/' -e 's/;TAGS=[^;]*;/;TAGS=~{sample_id};/' | \
      bgzip -c > "../~{sample_id}.wham.vcf.gz"
    cd ..
    bcftools index -t "~{sample_id}.wham.vcf.gz"

    df -h
    ls -l
  >>>
  runtime {
    cpu: cpu_cores
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " SSD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: wham_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

