## This WDL pipeline implements SV calling of mobile elements with Scramble by Rebecca I. Torene.
## https://github.com/GeneDx/scramble
##
## Requirements/expectations :
## - Human whole-genome pair-end sequencing data in mapped BAM or CRAM format

version 1.0

import "Structs.wdl"

workflow Scramble {
  input {
    File bam_or_cram_file
    File bam_or_cram_index
    File original_bam_or_cram_file
    File original_bam_or_cram_index
    File counts_file
    File input_vcf
    String sample_name
    File reference_fasta
    File reference_index
    File regions_list

    # Critical parameter for sensitivity/specificity
    # Recommended values for aligners:
    #   BWA-MEM: 90
    #   DRAGEN-3.7.8: 60
    Int alignment_score_cutoff

    Float min_clipped_reads_fraction = 0.22
    Int percent_align_cutoff = 70
    File mei_bed
    File? scramble_vcf_script
    String? make_scramble_vcf_args
    String sv_pipeline_docker
    String scramble_docker
    Int? part2_threads
    RuntimeAttr? runtime_attr_scramble_part1
    RuntimeAttr? runtime_attr_scramble_part2
    RuntimeAttr? runtime_attr_scramble_make_vcf
  }

  call ScramblePart1 {
    input:
      bam_or_cram_file = bam_or_cram_file,
      bam_or_cram_index = bam_or_cram_index,
      counts_file = counts_file,
      sample_name = sample_name,
      regions_list = regions_list,
      reference_fasta = reference_fasta,
      reference_index = reference_index,
      min_clipped_reads_fraction = min_clipped_reads_fraction,
      scramble_docker = scramble_docker,
      runtime_attr_override = runtime_attr_scramble_part1
  }

  call ScramblePart2 {
    input:
      clusters_file = ScramblePart1.clusters_file,
      sample_name = sample_name,
      reference_fasta = reference_fasta,
      percent_align_cutoff = percent_align_cutoff,
      alignment_score_cutoff = alignment_score_cutoff,
      min_clipped_reads = ScramblePart1.min_clipped_reads,
      threads = part2_threads,
      scramble_docker = scramble_docker,
      runtime_attr_override = runtime_attr_scramble_part2
  }

  call MakeScrambleVcf {
    input:
      scramble_table = ScramblePart2.table,
      original_bam_or_cram_file = original_bam_or_cram_file,
      original_bam_or_cram_index = original_bam_or_cram_index,
      input_vcf = input_vcf,
      sample_name = sample_name,
      reference_fasta = reference_fasta,
      reference_index = reference_index,
      mei_bed = mei_bed,
      scramble_vcf_script = scramble_vcf_script,
      script_args = make_scramble_vcf_args,
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
    File counts_file
    String sample_name
    File regions_list
    File reference_fasta
    File reference_index
    Float min_clipped_reads_fraction
    String scramble_docker
    RuntimeAttr? runtime_attr_override
  }

  Int disk_size_gb = ceil(size(bam_or_cram_file,"GiB") + size(reference_fasta,"GiB")*1.5 + 50)

  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 3.0,
                               disk_gb: disk_size_gb,
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File clusters_file = "~{sample_name}.scramble_clusters.tsv.gz"
    Int min_clipped_reads = read_int("cutoff.txt")
  }
  command <<<
    set -euo pipefail

    # Calibrate clipped reads cutoff based on median coverage
    zcat ~{counts_file} \
      | awk '$0!~"@"' \
      | sed 1d \
      | awk 'NR % 100 == 0' \
      | cut -f4 \
      | Rscript -e "cat(round(~{min_clipped_reads_fraction}*median(data.matrix(read.csv(file(\"stdin\"))))))" \
      > cutoff.txt
    MIN_CLIPPED_READS=$(cat cutoff.txt)
    echo "MIN_CLIPPED_READS: ${MIN_CLIPPED_READS}"

    # Identify clusters of split reads
    while read region; do
      time /app/scramble-gatk-sv/cluster_identifier/src/build/cluster_identifier -l -s ${MIN_CLIPPED_READS} -r "${region}" -t ~{reference_fasta} ~{bam_or_cram_file} \
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
    Int min_clipped_reads
    String scramble_docker
    Int percent_align_cutoff
    Int alignment_score_cutoff
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
            --pct-align ~{percent_align_cutoff} -n ~{min_clipped_reads} --mei-score ~{alignment_score_cutoff}

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
    File original_bam_or_cram_file
    File original_bam_or_cram_index
    File input_vcf
    File reference_fasta
    File reference_index
    File mei_bed
    String sample_name
    File? scramble_vcf_script
    String? script_args
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Int disk_size_gb = ceil(size(original_bam_or_cram_file, "GiB") + size(reference_fasta, "GiB") + 10)

  RuntimeAttr runtime_default = object {
                                  mem_gb: 4,
                                  disk_gb: disk_size_gb,
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
    set -euxo pipefail
    python ~{default="/opt/sv-pipeline/scripts/make_scramble_vcf.py" scramble_vcf_script} \
      --table ~{scramble_table} \
      --input-vcf ~{input_vcf} \
      --alignments-file ~{original_bam_or_cram_file} \
      --sample ~{sample_name} \
      --reference ~{reference_fasta} \
      --mei-bed ~{mei_bed} \
      --out unsorted.vcf.gz \
      ~{script_args}
    bcftools sort unsorted.vcf.gz -Oz -o ~{sample_name}.scramble.vcf.gz
    tabix ~{sample_name}.scramble.vcf.gz
  >>>
  output {
    File vcf = "~{sample_name}.scramble.vcf.gz"
    File vcf_index = "~{sample_name}.scramble.vcf.gz.tbi"
  }
}