## This WDL pipeline implements SV calling of mobile elements with Scramble by Rebecca I. Torene.
## https://github.com/GeneDx/scramble
##
## Requirements/expectations :
## - Human whole-genome pair-end sequencing data in mapped BAM format

version 1.0

import "Structs.wdl"

workflow Scramble {
  input {
    File bam_or_cram_file
    File bam_or_cram_index
    String sample_name
    File reference_fasta
    File reference_index
    File regions_list
    String scramble_docker
    Int? part2_threads
    RuntimeAttr? runtime_attr_scramble_part1
    RuntimeAttr? runtime_attr_scramble_part2
  }
    
  parameter_meta {
    bam_or_cram_file: "A .bam or .cram file to search for SVs. crams are preferable because they localise faster and use less disk."
    bam_or_cram_index: "Index for bam_or_cram_file."
    sample_name: "A sample name. Outputs will be sample_name+'.scramble.insertions.vcf.gz' and sample_name+'.scramble.deletions.vcf.gz'."
    reference_fasta: "A FASTA file with the reference used to align bam or cram file."
    detect_deletions: "Run deletion detection as well as mobile element insertion."
  }
  
  call ScramblePart1 {
    input:
      bam_or_cram_file = bam_or_cram_file,
      bam_or_cram_index = bam_or_cram_index,
      sample_name = sample_name,
      regions_list = regions_list,
      reference_fasta = reference_fasta,
      reference_index = reference_index,
      scramble_docker = scramble_docker,
      runtime_attr_override = runtime_attr_scramble_part1
  }

  call ScramblePart2 {
    input:
      clusters_file = ScramblePart1.clusters_file,
      sample_name = sample_name,
      reference_fasta = reference_fasta,
      threads = part2_threads,
      scramble_docker = scramble_docker,
      runtime_attr_override = runtime_attr_scramble_part2
  }

  output {
    File vcf = ScramblePart2.vcf
    File index = ScramblePart2.index
  }
}

task ScramblePart1 {
  input {
    File bam_or_cram_file
    File bam_or_cram_index
    String sample_name
    File regions_list
    File reference_fasta
    File reference_index
    String scramble_docker
    RuntimeAttr? runtime_attr_override
  }

  Int disk_size_gb = ceil(size(bam_or_cram_file,"GiB") + size(reference_fasta,"GiB")*1.5 + 50)

  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 2.0,
                               disk_gb: disk_size_gb,
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File clusters_file = "~{sample_name}.scramble_clusters.tsv.gz"
  }
  command <<<
    set -euo pipefail

    # Identify clusters of split reads
    while read region; do
      time /app/scramble-gatk-sv/cluster_identifier/src/build/cluster_identifier -l -r "${region}" -t ~{reference_fasta} ~{bam_or_cram_file} \
        | gzip >> ~{sample_name}.scramble_clusters.tsv.gz
    done < ~{regions_list}
  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: scramble_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

task ScramblePart2 {
  input {
    File clusters_file
    String sample_name
    File reference_fasta
    String scramble_docker
    Int percent_align_cutoff = 70
    Int threads = 7  # Number of threads
    RuntimeAttr? runtime_attr_override
  }

  Int disk_size_gb = ceil(10*size(clusters_file,"GiB") + size(reference_fasta,"GiB") + 10)

  RuntimeAttr default_attr = object {
    cpu_cores: 8,
    mem_gb: 12.0,
    disk_gb: disk_size_gb,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File vcf = "~{sample_name}.scramble.vcf.gz"
    File index = "~{sample_name}.scramble.vcf.gz.tbi"
  }
  command <<<
    set -euo pipefail

    xDir=$PWD
    clusterFile=$xDir/clusters
    scrambleDir="/app/scramble-gatk-sv"
    meiRef=$scrambleDir/cluster_analysis/resources/MEI_consensus_seqs.fa

    # create a blast db from the reference
    cat ~{reference_fasta} | makeblastdb -in - -parse_seqids -title ref -dbtype nucl -out ref

    gunzip -c ~{clusters_file} > $clusterFile

    # Produce ${clusterFile}_MEIs.txt
    Rscript --vanilla $scrambleDir/cluster_analysis/bin/SCRAMble.R --out-name $clusterFile \
            --cluster-file $clusterFile --install-dir $scrambleDir/cluster_analysis/bin \
            --mei-refs $meiRef --ref $xDir/ref --no-vcf --eval-meis --cores ~{threads} \
            --pct-align ~{percent_align_cutoff}

    # create a header for the output vcf
    echo \
    '##fileformat=VCFv4.3
    ##reference=~{reference_fasta}
    ##source=scramble' > tmp.vcf

    grep '^>' $meiRef | awk \
    '{mei=toupper(substr($0,2)); if (mei=="L1") mei="LINE1"
      print "##ALT=<ID=INS:ME:" mei ",Description=\"" mei " element insertion\">"}' >> tmp.vcf

    echo \
    '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
    ##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
    ##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">
    ##INFO=<ID=ALGORITHMS,Number=.,Type=String,Description="Source algorithms">
    ##INFO=<ID=STRANDS,Number=1,Type=String,Description="Breakpoint strandedness [++,+-,-+,--]">
    ##INFO=<ID=CHR2,Number=1,Type=String,Description="Chromosome for END coordinate">
    ##INFO=<ID=MEI_START,Number=1,Type=Integer,Description="Start of alignment to canonical MEI sequence">
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' >> tmp.vcf

    blastdbcmd -db ref -entry all -outfmt '##contig=<ID=%a,length=%l>' >> tmp.vcf
    echo "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	~{sample_name}" >> tmp.vcf

    # use awk to write the first part of an awk script that initializes an awk associative array
    # mapping MEI name onto its consensus sequence length
    awk \
    'BEGIN { FS=OFS="\t"; print "BEGIN{" }
     /^>/  {if ( seq != "" ) print "seqLen[\"" seq "\"]=" len; seq = substr($0,2); len = 0}
     !/^>/ {len += length($0)}
     END   {if ( seq != "" ) print "seqLen[\"" seq "\"]=" len; print "}"}' $meiRef > awkScript.awk

    # write the rest of the awk script that transforms the contents of the *_MEIs.txt files into a VCF
    echo \
    'BEGIN{ FS=OFS="\t" }
    { if(FNR<2)next
      split($1,loc,":")
      start=loc[2]+1
      end=start+1
      len=seqLen[$2]-$10
      mei=toupper($2); if (mei=="L1") mei="LINE1"
      print loc[1],start,".","N","<INS:ME:" mei ">",int($6),"PASS",\
            "END=" end ";SVTYPE=INS;SVLEN=" len ";MEI_START=" $10 ";STRANDS=+-;CHR2=" loc[1] ";ALGORITHMS=scramble",\
            "GT","0/1" }' >> awkScript.awk

    # transform the MEI descriptions into VCF lines
    awk -f awkScript.awk ${clusterFile}_MEIs.txt >> tmp.vcf

    # sort and index the output VCF
    bcftools sort -Oz <tmp.vcf >"~{sample_name}.scramble.vcf.gz"
    bcftools index -ft "~{sample_name}.scramble.vcf.gz"
  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: scramble_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}
