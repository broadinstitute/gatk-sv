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
    String sv_pipeline_docker
    String scramble_docker
    Int? part2_threads
    RuntimeAttr? runtime_attr_scramble_part1
    RuntimeAttr? runtime_attr_scramble_part2
    RuntimeAttr? runtime_attr_scramble_make_vcf
  }

  parameter_meta {
    bam_or_cram_file: "A .bam or .cram file to search for SVs. crams are preferable because they localise faster and use less disk."
    bam_or_cram_index: "Index for bam_or_cram_file."
    sample_name: "A sample name. Outputs will be sample_name+'.scramble.insertions.vcf.gz' and sample_name+'.scramble.deletions.vcf.gz'."
    reference_fasta: "A FASTA file with the reference used to align bam or cram file."
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

  call MakeScrambleVcf {
    input:
      scramble_table = ScramblePart2.table,
      sample_name = sample_name,
      reference_fasta = reference_fasta,
      reference_index = reference_index,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_scramble_make_vcf
  }

  output {
    File vcf = MakeScrambleVcf.vcf
    File index = MakeScrambleVcf.vcf_index
    File clusters = ScramblePart1.clusters_file
    File table = ScramblePart2.table
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
    File table = "~{sample_name}.scramble.tsv.gz"
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

    # Save raw outputs
    mv ${clusterFile}_MEIs.txt ~{sample_name}.scramble.tsv
    gzip ~{sample_name}.scramble.tsv
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

task MakeScrambleVcf {
  input {
    File scramble_table
    File reference_fasta
    File reference_index
    String sample_name
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr runtime_default = object {
                                  mem_gb: 0.9,
                                  disk_gb: 10,
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
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euo pipefail
    # Sort by position

    zcat ~{scramble_table} \
      | awk -F '\t' -v OFS='\t' '{print "chrom","pos","end",$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16; exit 0;}' \
      > header.txt

    zcat ~{scramble_table} \
      | sed 1d \
      | tr ':' '\t' \
      | awk -F'\t' -v OFS='\t' '{print $1,$2,$2+1,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17}' \
      | sort -k1,1V -k2,2n \
      | cat header.txt - \
      | bgzip \
      > scramble.bed.gz

    python /opt/sv-pipeline/scripts/make_scramble_vcf.py \
      --table ~{scramble_table} \
      --sample ~{sample_name} \
      --reference ~{reference_fasta} \
      --out unsorted.vcf.gz
    bcftools sort unsorted.vcf.gz -Oz -o ~{sample_name}.scramble.vcf.gz
    tabix ~{sample_name}.scramble.vcf.gz
  >>>
  output {
    File vcf = "~{sample_name}.scramble.vcf.gz"
    File vcf_index = "~{sample_name}.scramble.vcf.gz.tbi"
  }
}