version 1.0

import "Structs.wdl"
import "CramToBam.wdl" as ctb

# Run Whamg SV detection algorithm on whole genome in bam or cram
#   file, or if include_list is provided, run whamg on explicitly included
#   subset of genome and concatenate.
# Then fix output vcf headers to contain good contig length data.

workflow Whamg {
  input {
    File bam_or_cram_file
    File? bam_or_cram_index
    String sample_id
    File reference_fasta
    File? reference_index
    File? include_bed_file
    File chr_file
    Int? pf_reads_improper_pairs
    Float? pct_exc_total
    String samtools_cloud_docker
    String wham_docker
    RuntimeAttr? runtime_attr_includelist
    RuntimeAttr? runtime_attr_cram_to_bam
    RuntimeAttr? runtime_attr_wham
  }
    
  parameter_meta {
    bam_or_cram_file: ".bam or .cram file to search for SVs. bams are preferable, crams will be converted to bams."
    bam_or_cram_index: "[optional] associated index file. If omitted, the WDL will look for an index file by appending .bai/.crai to the .bam/.cram file"
    sample_id: "sample name. Outputs will be sample_name + '.vcf.gz' and sample_name + '.vcf.gz.tbi'"
    reference_fasta: ".fasta file with reference used to align bam or cram file"
    reference_index: "[optional] reference index file. If omitted, the WDL will look for an index by appending .fai to the .fasta file"
    chr_file: "text file with newline-separated list of contigs that whamg will use to estimate template size. Typically you will want only primary contigs."
    include_bed_file: "[optional] bed file with intervals where whamg should make calls. If omitted, whamg will run on whole genome."
    pf_reads_improper_pairs: "[optional] Value of PF_READS_IMPROPER_PAIRS obtained from CollectAlignmentSummaryMetrics, used for optimal estimate of VM memory needs."
    pct_exc_total: "[optional] Value of PCT_EXC_TOTAL obtained from CollectWgsMetrics, used for optimal estimate of VM memory needs."
  }
  meta {
      author: "Ted Brookings"
      email: "tbrookin@broadinstitute.org"
  }

  if (basename(bam_or_cram_file, ".bam") + ".bam" != basename(bam_or_cram_file)) {
    call ctb.RunCramToBam as RunCramToBam {
      input:
        cram_file = bam_or_cram_file,
        reference_fasta = reference_fasta,
        reference_index = reference_index,
        samtools_cloud_docker = samtools_cloud_docker,
        runtime_attr_override = runtime_attr_cram_to_bam
    }
  }

  File bam_file = select_first([RunCramToBam.bam_file, bam_or_cram_file])
  File bam_index = select_first([RunCramToBam.bam_index, bam_or_cram_index, bam_or_cram_file + ".bai"])
  
  # decide whether to use includelist version or baseline
  Boolean use_include_list = defined(include_bed_file)

  if (use_include_list) {
    call RunWhamgIncludelist {
      input:
        bam_file = bam_file,
        bam_index = bam_index,
        reference_fasta = reference_fasta,
        reference_index = reference_index,
        sample_id = sample_id,
        chr_file = chr_file,
        include_bed_file = select_first([include_bed_file]),
        pf_reads_improper_pairs = pf_reads_improper_pairs,
        pct_exc_total = pct_exc_total,
        wham_docker = wham_docker,
        runtime_attr_override = runtime_attr_includelist
    }
  } # else
  if (!use_include_list) {
    call RunWhamg {
      input:
        bam_file = bam_file,
        bam_index = bam_index,
        reference_fasta = reference_fasta,
        reference_index = reference_index,
        sample_id = sample_id,
        chr_file = chr_file,
        pf_reads_improper_pairs = pf_reads_improper_pairs,
        pct_exc_total = pct_exc_total,
        wham_docker = wham_docker,
        runtime_attr_override = runtime_attr_wham
    }
  }

  output {
    File index = select_first([RunWhamg.index, RunWhamgIncludelist.index])
    File vcf = select_first([RunWhamg.vcf, RunWhamgIncludelist.vcf])
  }
}

task RunWhamg {
  input {
    File bam_file
    File? bam_index
    File reference_fasta
    File? reference_index
    String sample_id
    File chr_file
    Int? pf_reads_improper_pairs
    Float? pct_exc_total
    String wham_docker
    RuntimeAttr? runtime_attr_override
  }

  File bam_index_file = select_first([bam_index, bam_file + ".bai"])
  File reference_index_file = select_first([reference_index, reference_fasta + ".fai"])
  
  Int num_cpu = if defined(runtime_attr_override) then select_first([select_first([runtime_attr_override]).cpu_cores, 4]) else 4
  
  # Ensure there's sufficient disk space.
  Float disk_overhead = 5.0
  Float bam_size = size(bam_file, "GiB")
  Float bam_index_size = size(bam_index_file, "GiB")
  Float ref_size = size(reference_fasta, "GiB")
  Float ref_index_size = size(reference_index_file, "GiB")
  Float chr_size = size(chr_file, "GiB")
  Int vm_disk_size = ceil(bam_size + bam_index_size + ref_size + ref_index_size + chr_size + disk_overhead)

  # Ensure there's sufficient memory. Use picard metrics if
  # available, otherwise estimate based on bam size
  Float mem_per_bam_size = 0.1479
  Float mem_bam_offset = 18.96
  Float mem_per_improper_pairs = "4.401e-7"
  Float mem_per_exc = 20.51
  Float mem_metrics_offset = 4.570
  Float mem_size_gb =
    if defined(pf_reads_improper_pairs) && defined(pct_exc_total) then
      mem_metrics_offset + mem_per_improper_pairs * select_first([pf_reads_improper_pairs]) + mem_per_exc * select_first([pct_exc_total])
    else
      mem_per_bam_size * bam_size + mem_bam_offset

  Array[String] chr_list = read_lines(chr_file)

  RuntimeAttr default_attr = object {
    mem_gb: mem_size_gb, 
    disk_gb: vm_disk_size,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File vcf = "${sample_id}.wham.vcf.gz"
    File index = "${sample_id}.wham.vcf.gz.tbi"
  }
  command <<<

    set -euo pipefail
    # print some info that may be useful for debugging
    df -h
    cat /sys/fs/cgroup/memory/memory.stat
    echo "whamg $(./whamg 2>&1 >/dev/null | grep Version)"
    
    # ensure that index files are present in appropriate locations
    if [ ! -e "~{bam_file}.bai" ]; then
        ln -s "~{bam_index_file}" "~{bam_file}.bai"
    fi
    REPLACED_VERSION=$(echo "~{bam_file}" | sed 's/.bam$/.bai/g')
    if [ ! -e "$REPLACED_VERSION" ]; then
        ln -s "~{bam_index_file}" "$REPLACED_VERSION"
    fi
    if [ ! -e "~{reference_fasta}.fai" ]; then
        ln -s "~{reference_index_file}" "~{reference_fasta}.fai"
    fi
        
    # We need to update both the VCF sample ID and the TAGS INFO field in the WHAM output VCFs.
    # WHAM uses both to store the sample identifier, and by default uses the SM identifier from the BAM file.
    # We need to update both to reflect the potentially-renamed sample identifier used by the pipeline (sample_id) --
    # svtk standardize_vcf uses the TAGS field to identify the sample for WHAM VCFs.

    VCF_FILE_FIXED_HEADER="~{sample_id}.wham.vcf.gz"
    VCF_FILE_BAD_TAGS="~{sample_id}.wham_bad_header_bad_tags.vcf.gz"
    VCF_FILE_BAD_HEADER="~{sample_id}.wham_bad_header.vcf.gz"

    echo "Invoking whamg"

    whamg -c "~{sep=","  chr_list}" -x ~{num_cpu} -a ~{reference_fasta} -f ~{bam_file} \
      | bgzip -c > "$VCF_FILE_BAD_TAGS"
    tabix -p vcf "$VCF_FILE_BAD_TAGS"

    # write out a an annotation table with the new sample ID for each variant record.
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t~{sample_id}\n' $VCF_FILE_BAD_TAGS | bgzip -c > tags_annotation_file.tsv.gz
    tabix -f -s1 -b2 -e2 tags_annotation_file.tsv.gz
    # Use bcftools annotate to re-write the TAGS field in each variant based on the annotation table.
    bcftools annotate -a tags_annotation_file.tsv.gz -c CHROM,POS,REF,ALT,INFO/TAGS  $VCF_FILE_BAD_TAGS | bgzip -c > $VCF_FILE_BAD_HEADER
    tabix -p vcf "$VCF_FILE_BAD_HEADER"

    echo "Adding contigs with length to vcf header"
    # get the existing header
    OLD_HEADER=$(bcftools view -h "$VCF_FILE_BAD_HEADER" | grep -v '##contig=')
    # create new header lines with the contig lengths
    CONTIGS_PATTERN="^$(cat ~{chr_file} | paste -sd "," - | sed "s/,/\\\t|^/g")\t"
    CONTIGS_HEADER=$(grep $CONTIGS_PATTERN -P ~{reference_index_file} | awk '{print "##contig=<ID=" $1 ",length=" $2 ">"}')
    # Create a new header with
    #   -all but last line of old header, followed by newline
    #   -contig header lines, followed by newline
    #   -last line of old header
    # Replace old header with new header
    echo "Replacing header"
    bcftools reheader \
        -h <(echo "$OLD_HEADER" | sed \$d ; echo "$CONTIGS_HEADER" ; echo "$OLD_HEADER" | tail -n 1) \
        -s <(echo "~{sample_id}") \
        "$VCF_FILE_BAD_HEADER" > $VCF_FILE_FIXED_HEADER
    
    echo "Indexing vcf"
    tabix "$VCF_FILE_FIXED_HEADER"
    
    echo "finished RunWhamg"
    
  >>>
  runtime {
    cpu: num_cpu
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: wham_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

}

task RunWhamgIncludelist {
  input {
    File bam_file
    File? bam_index
    File reference_fasta
    File? reference_index
    String sample_id
    File chr_file
    File include_bed_file
    Int? pf_reads_improper_pairs
    Float? pct_exc_total
    String wham_docker
    RuntimeAttr? runtime_attr_override
  }

  File bam_index_file = select_first([bam_index, bam_file + ".bai"])
  File reference_index_file = select_first([reference_index, reference_fasta + ".fai"])
  
  Int num_cpu = if defined(runtime_attr_override) then select_first([select_first([runtime_attr_override]).cpu_cores, 4]) else 4
  
  # Ensure there's sufficient disk space.
  Float disk_overhead = 5.0
  Float bam_size = size(bam_file, "GiB")
  Float bam_index_size = size(bam_index_file, "GiB")
  Float ref_size = size(reference_fasta, "GiB")
  Float ref_index_size = size(reference_index_file, "GiB")
  Float chr_size = size(chr_file, "GiB")
  Float include_list_size = size(include_bed_file, "GiB")
  Int vm_disk_size = ceil(bam_size + bam_index_size + ref_size + ref_index_size + chr_size + include_list_size + disk_overhead)

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
  

  Array[String] chr_list = read_lines(chr_file)
  Array[String] good_intervals = read_lines(include_bed_file)

  RuntimeAttr default_attr = object {
    mem_gb: mem_size_gb, 
    disk_gb: vm_disk_size,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File vcf = "${sample_id}.wham.vcf.gz"
    File index = "${sample_id}.wham.vcf.gz.tbi"
  }
  command <<<

    set -euo pipefail
    # print some info that may be useful for debugging
    df -h
    cat /sys/fs/cgroup/memory/memory.stat
    echo "whamg $(./whamg 2>&1 >/dev/null | grep Version)"
    
    # ensure that index files are present in appropriate locations
    if [ ! -e "~{bam_file}.bai" ]; then
      ln -s "~{bam_index_file}" "~{bam_file}.bai"
    fi
    REPLACED_VERSION=$(echo "~{bam_file}" | sed 's/.bam$/.bai/g')
    if [ ! -e "$REPLACED_VERSION" ]; then
      ln -s "~{bam_index_file}" "$REPLACED_VERSION"
    fi
    if [ ! -e "~{reference_fasta}.fai" ]; then
      ln -s "~{reference_index_file}" "~{reference_fasta}.fai"
    fi
    
    # run whamg on all good intervals
    GOOD_INTERVALS=$(echo "~{sep='\n' good_intervals}" \
      | awk '{print $1 ":" $2 "-" $3}')
    VCFS=""
    for GOOD_INTERVAL in $GOOD_INTERVALS; do
      echo "Analyzing $GOOD_INTERVAL"
      NEW_VCF="$GOOD_INTERVAL.wham.vcf.gz"
      # run wham, gzip results
      whamg \
        -c "~{sep=","  chr_list}" \
        -x ~{num_cpu} \
        -a "~{reference_fasta}"\
        -f "~{bam_file}" \
        -r "$GOOD_INTERVAL" \
        | bgzip -c > "$NEW_VCF"
      # index results vcf
      tabix "$NEW_VCF"
      # add vcf file name to list of vcfs, separated by spaces
      if [ -z "$VCFS" ]; then
        VCFS="$NEW_VCF"
      else
        VCFS="$VCFS $NEW_VCF"
      fi
    done
    
    # We need to update both the VCF sample ID and the TAGS INFO field in the WHAM output VCFs.
    # WHAM uses both to store the sample identifier, and by default uses the SM identifier from the BAM file.
    # We need to update both to reflect the potentially-renamed sample identifier used by the pipeline (sample_id) --
    # svtk standardize_vcf uses the TAGS field to identify the sample for WHAM VCFs.

    VCF_FILE_FIXED_HEADER="~{sample_id}.wham.vcf.gz"
    VCF_FILE_BAD_TAGS="~{sample_id}.wham_bad_header_bad_tags.vcf.gz"
    VCF_FILE_BAD_HEADER="~{sample_id}.wham_bad_header.vcf.gz"

    # concatenate resulting vcfs into final vcf
    echo "Concatenating results"
    bcftools concat -a -O z -o "$VCF_FILE_BAD_TAGS" $VCFS
    tabix -p vcf "$VCF_FILE_BAD_TAGS"

    # write out a an annotation table with the new sample ID for each variant record.
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t~{sample_id}\n' $VCF_FILE_BAD_TAGS | bgzip -c > tags_annotation_file.tsv.gz
    tabix -f -s1 -b2 -e2 tags_annotation_file.tsv.gz
    # Use bcftools annotate to re-write the TAGS field in each variant based on the annotation table.
    bcftools annotate -a tags_annotation_file.tsv.gz -c CHROM,POS,REF,ALT,INFO/TAGS  $VCF_FILE_BAD_TAGS | bgzip -c > $VCF_FILE_BAD_HEADER
    tabix -p vcf "$VCF_FILE_BAD_HEADER"

    echo "Getting existing header"
    # get the existing header
    OLD_HEADER=$(bcftools view -h "$VCF_FILE_BAD_HEADER" | grep -v '##contig=')
    # create new header lines with the contig lengths
    echo "Making grep pattern to extract contigs from reference index"
    CONTIGS_PATTERN="^$(cat ~{chr_file} | paste -sd "," - | sed "s/,/\\\t|^/g")\t"
    echo "Adding contigs with length to vcf header"        
    CONTIGS_HEADER=$(grep "$CONTIGS_PATTERN" -P "~{reference_index_file}" | awk '{print "##contig=<ID=" $1 ",length=" $2 ">"}')
    # Create a new header with
    #   -all but last line of old header, followed by newline
    #   -contig header lines, followed by newline
    #   -last line of old header
    # Replace old header with new header
    echo "Replacing header"
    echo "$OLD_HEADER" | grep -v "^#CHROM" > new_header.txt
    echo "$CONTIGS_HEADER" >> new_header.txt
    echo "$OLD_HEADER" | grep "^#CHROM" >> new_header.txt
    bcftools reheader \
      -h new_header.txt \
      -s <(echo "~{sample_id}") \
      "$VCF_FILE_BAD_HEADER" \
      > "$VCF_FILE_FIXED_HEADER"
    
    echo "Indexing vcf"
    tabix "$VCF_FILE_FIXED_HEADER"
    
    echo "finished RunWhamg"
    
  >>>
  runtime {
    cpu: num_cpu
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: wham_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

}

