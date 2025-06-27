version 1.0

import "Structs.wdl"
import "CollectCoverage.wdl" as cov
import "CollectSVEvidence.wdl" as coev
import "Manta.wdl" as manta
import "MELT.wdl" as melt
import "Scramble.wdl" as scramble
import "Whamg.wdl" as wham
import "GatherSampleEvidenceMetrics.wdl" as metrics

# Runs selected tools on BAM/CRAM files

workflow GatherSampleEvidence {
  input {
    File bam_or_cram_file
    File? bam_or_cram_index

    # Note: raw and "safe" CRAM/BAM IDs can be generated with GetSampleID
    String sample_id

    Boolean? is_dragen_3_7_8
    # Optionally pass in DRAGEN-SV VCF if available
    File? dragen_vcf

    # Evidence collection flags
    Boolean collect_coverage = true
    Boolean collect_pesr = true

    # Google Cloud Platform (GCP) and Azure users can safely enable and get a slight cost improvement.
    # Users running on shared file systems should NOT enable this feature.
    # Enabling this option when running the workflow on GCP or Azure will result in moving the **localized**
    # files from one path to another on the VM, without impacting the files in their source persistent location.
    # It will lead to using a slightly smaller disk size and running faster, thereby providing a slight cost improvement.
    # However, when run on shared file systems (e.g., HPC), it will, by default, create a copy of the
    # input files, and all subsequent operations will run on the deep copy of the input file.
    Boolean move_bam_or_cram_files = false

    # Localize reads parameters
    # set to true on default, skips localize_reads if set to false
    Boolean run_localize_reads = true

    # Common parameters
    File primary_contigs_list
    File reference_fasta
    File reference_index    # Index (.fai), must be in same dir as fasta
    File reference_dict     # Dictionary (.dict), must be in same dir as fasta
    String? reference_version   # Either "38" or "19"

    # Reference bwa index files, only required for alignments with Dragen 3.7.8
    File? reference_bwa_alt
    File? reference_bwa_amb
    File? reference_bwa_ann
    File? reference_bwa_bwt
    File? reference_bwa_pac
    File? reference_bwa_sa

    # Coverage collection inputs
    File preprocessed_intervals
    Float? mem_gb_for_collect_counts
    Int? disk_space_gb_for_collect_counts

    # Manta inputs
    File manta_region_bed
    File manta_region_bed_index
    Float? manta_jobs_per_cpu
    Int? manta_mem_gb_per_job

    # PESR inputs
    File sd_locs_vcf

    # Melt inputs
    File? melt_standard_vcf_header # required if run_melt True
    File? melt_metrics_intervals
    Float? insert_size
    Int? read_length
    Float? coverage
    File? metrics_intervals
    Float? pct_chimeras
    Float? total_reads
    Int? pf_reads_improper_pairs

    # Scramble inputs
    File mei_bed
    Int? scramble_alignment_score_cutoff
    Int? scramble_percent_align_cutoff
    Float? scramble_min_clipped_reads_fraction
    Int? scramble_part2_threads
    File? scramble_vcf_script

    # Required if running Scramble but not running Manta
    File? manta_vcf_input
    File? manta_vcf_index_input

    # Wham inputs
    File wham_include_list_bed_file

    # Module metrics parameters
    # Run module metrics workflow at the end - on by default
    Boolean run_module_metrics = true
    File? primary_contigs_fai # required if run_module_metrics = true
    File? baseline_manta_vcf # baseline files are optional for metrics workflow
    File? baseline_wham_vcf
    File? baseline_melt_vcf
    File? baseline_scramble_vcf

    # Docker
    String sv_pipeline_docker
    String sv_base_mini_docker
    String samtools_cloud_docker
    String? manta_docker
    String? melt_docker
    String? scramble_docker
    String? wham_docker
    String gatk_docker
    String? gatk_docker_pesr_override
    String genomes_in_the_cloud_docker
    String cloud_sdk_docker

    # Runtime configuration overrides
    RuntimeAttr? runtime_attr_localize_reads
    RuntimeAttr? runtime_attr_split_cram
    RuntimeAttr? runtime_attr_concat_bam
    RuntimeAttr? runtime_attr_manta
    RuntimeAttr? runtime_attr_melt_coverage
    RuntimeAttr? runtime_attr_melt_metrics
    RuntimeAttr? runtime_attr_melt
    RuntimeAttr? runtime_attr_scramble_part1
    RuntimeAttr? runtime_attr_scramble_part2
    RuntimeAttr? runtime_attr_scramble_make_vcf
    RuntimeAttr? runtime_attr_realign_soft_clips
    RuntimeAttr? runtime_attr_scramble_part1_realigned
    RuntimeAttr? runtime_attr_scramble_part2_realigned
    RuntimeAttr? runtime_attr_scramble_make_vcf_realigned
    RuntimeAttr? runtime_attr_pesr
    RuntimeAttr? runtime_attr_wham

    # Never assign these values! (workaround until None type is implemented)
    Float? NONE_FLOAT_
    Int? NONE_INT_
    File? NONE_FILE_
  }

  Boolean run_manta = defined(manta_docker)
  Boolean run_melt = defined(melt_docker)
  Boolean run_scramble = defined(scramble_docker)
  Boolean run_wham = defined(wham_docker)

  Boolean is_bam_ = basename(bam_or_cram_file, ".bam") + ".bam" == basename(bam_or_cram_file)
  String index_ext_ = if is_bam_ then ".bai" else ".crai"
  File bam_or_cram_index_ = select_first([bam_or_cram_index, bam_or_cram_file + index_ext_])

  # move the reads nearby -- handles requester_pays and makes cross-region transfers just once
  if (run_localize_reads) {
    call LocalizeReads {
    input:
      reads_path = bam_or_cram_file,
      reads_index = bam_or_cram_index_,
      move_files = move_bam_or_cram_files,
      runtime_attr_override = runtime_attr_localize_reads
    }
  }

  File reads_file_ = select_first([LocalizeReads.output_file, bam_or_cram_file])
  File reads_index_ = select_first([LocalizeReads.output_index, bam_or_cram_index_])

  if (collect_coverage || run_melt || run_scramble) {
    call cov.CollectCounts {
      input:
        intervals = preprocessed_intervals,
        cram_or_bam = reads_file_,
        cram_or_bam_idx = reads_index_,
        sample_id = sample_id,
        ref_fasta = reference_fasta,
        ref_fasta_fai = reference_index,
        ref_fasta_dict = reference_dict,
        gatk_docker = gatk_docker,
        mem_gb = mem_gb_for_collect_counts,
        disk_space_gb = disk_space_gb_for_collect_counts,
        disabled_read_filters = ["MappingQualityReadFilter"]
    }
  }

  if (run_manta) {
    call manta.Manta {
      input:
        bam_or_cram_file = reads_file_,
        bam_or_cram_index = reads_index_,
        sample_id = sample_id,
        reference_fasta = reference_fasta,
        reference_index = reference_index,
        region_bed = manta_region_bed,
        region_bed_index = manta_region_bed_index,
        jobs_per_cpu = manta_jobs_per_cpu,
        mem_gb_per_job = manta_mem_gb_per_job,
        manta_docker = select_first([manta_docker]),
        runtime_attr_override = runtime_attr_manta
    }
  }

  if (collect_pesr) {
    call coev.CollectSVEvidence {
      input:
        bam_or_cram_file = reads_file_,
        bam_or_cram_index = reads_index_,
        sample_id = sample_id,
        reference_fasta = reference_fasta,
        reference_index = reference_index,
        reference_dict = reference_dict,
        sd_locs_vcf = sd_locs_vcf,
        primary_contigs_list = primary_contigs_list,
        preprocessed_intervals = preprocessed_intervals,
        gatk_docker = select_first([gatk_docker_pesr_override, gatk_docker]),
        runtime_attr_override = runtime_attr_pesr
    }
  }

  if (run_melt) {
    call melt.MELT {
      input:
        bam_or_cram_file = reads_file_,
        bam_or_cram_index = reads_index_,
        counts_file = select_first([CollectCounts.counts]),
        sample_id = sample_id,
        reference_fasta = reference_fasta,
        reference_index = reference_index,
        reference_dict = reference_dict,
        reference_version = reference_version,
        melt_standard_vcf_header = select_first([melt_standard_vcf_header]),
        insert_size = insert_size,
        read_length = read_length,
        coverage = coverage,
        wgs_metrics_intervals = melt_metrics_intervals,
        pct_chimeras = pct_chimeras,
        total_reads = total_reads,
        pf_reads_improper_pairs = pf_reads_improper_pairs,
        runtime_attr_coverage = runtime_attr_melt_coverage,
        runtime_attr_metrics = runtime_attr_melt_metrics,
        samtools_cloud_docker = samtools_cloud_docker,
        gatk_docker = gatk_docker,
        genomes_in_the_cloud_docker = genomes_in_the_cloud_docker,
        melt_docker = select_first([melt_docker]),
        runtime_attr_melt = runtime_attr_melt
    }
  }

  if (run_scramble) {
    if (!defined(is_dragen_3_7_8)) {
      # check if the reads were aligned with dragen 3.7.8
      call CheckAligner {
        input:
          reads_path = reads_file_,
          reads_index = reads_index_,
          reference_fasta = if is_bam_ then NONE_FILE_ else reference_fasta,
          reference_index = if is_bam_ then NONE_FILE_ else reference_index,
          reference_dict = if is_bam_ then NONE_FILE_ else reference_dict,
          sample_id = sample_id,
          gatk_docker = gatk_docker,
          runtime_attr_override = runtime_attr_localize_reads
      }
    }
    Boolean is_dragen_3_7_8_ = (defined(CheckAligner.is_dragen_3_7_8) && select_first([CheckAligner.is_dragen_3_7_8]) > 0)
      || (!defined(CheckAligner.is_dragen_3_7_8) && select_first([is_dragen_3_7_8]))

    # Significant modulator of precision/sensitivity
    Int auto_alignment_score_cutoff = if is_dragen_3_7_8_ then 60 else 90
    Int scramble_alignment_score_cutoff_ = select_first([scramble_alignment_score_cutoff, auto_alignment_score_cutoff])

    call scramble.Scramble {
      input:
        bam_or_cram_file = reads_file_,
        bam_or_cram_index = reads_index_,
        original_bam_or_cram_file = reads_file_,
        original_bam_or_cram_index = reads_index_,
        counts_file = select_first([CollectCounts.counts]),
        input_vcf = select_first([dragen_vcf, Manta.vcf, manta_vcf_input]),
        sample_name = sample_id,
        reference_fasta = reference_fasta,
        reference_index = reference_index,
        mei_bed = mei_bed,
        regions_list = primary_contigs_list,
        alignment_score_cutoff = scramble_alignment_score_cutoff_,
        percent_align_cutoff = scramble_percent_align_cutoff,
        min_clipped_reads_fraction = scramble_min_clipped_reads_fraction,
        part2_threads = scramble_part2_threads,
        scramble_vcf_script = scramble_vcf_script,
        scramble_docker = select_first([scramble_docker]),
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_scramble_part1 = runtime_attr_scramble_part1,
        runtime_attr_scramble_part2 = runtime_attr_scramble_part2,
        runtime_attr_scramble_make_vcf = runtime_attr_scramble_make_vcf
    }
    if (is_dragen_3_7_8_) {
      # addresses bug in dragmap where some reads are incorrectly soft-clipped
      call RealignSoftClippedReads {
        input:
          reads_path = reads_file_,
          reads_index = reads_index_,
          scramble_table = Scramble.table,
          is_bam = is_bam_,
          sample_id = sample_id,
          reference_fasta = reference_fasta,
          reference_index = reference_index,
          reference_bwa_alt = select_first([reference_bwa_alt]),
          reference_bwa_amb = select_first([reference_bwa_amb]),
          reference_bwa_ann = select_first([reference_bwa_ann]),
          reference_bwa_bwt = select_first([reference_bwa_bwt]),
          reference_bwa_pac = select_first([reference_bwa_pac]),
          reference_bwa_sa = select_first([reference_bwa_sa]),
          sv_base_mini_docker = sv_base_mini_docker,
          runtime_attr_override = runtime_attr_realign_soft_clips
      }
      call scramble.Scramble as ScrambleRealigned {
        input:
          bam_or_cram_file = RealignSoftClippedReads.out,
          bam_or_cram_index = RealignSoftClippedReads.out_index,
          original_bam_or_cram_file = reads_file_,
          original_bam_or_cram_index = reads_index_,
          counts_file = select_first([CollectCounts.counts]),
          input_vcf = select_first([dragen_vcf, Manta.vcf, manta_vcf_input]),
          sample_name = sample_id,
          reference_fasta = reference_fasta,
          reference_index = reference_index,
          mei_bed = mei_bed,
          regions_list = primary_contigs_list,
          alignment_score_cutoff = scramble_alignment_score_cutoff_,
          percent_align_cutoff = scramble_percent_align_cutoff,
          min_clipped_reads_fraction = scramble_min_clipped_reads_fraction,
          part2_threads = scramble_part2_threads,
          scramble_vcf_script = scramble_vcf_script,
          scramble_docker = select_first([scramble_docker]),
          sv_pipeline_docker = sv_pipeline_docker,
          runtime_attr_scramble_part1 = runtime_attr_scramble_part1_realigned,
          runtime_attr_scramble_part2 = runtime_attr_scramble_part2_realigned,
          runtime_attr_scramble_make_vcf = runtime_attr_scramble_make_vcf_realigned
      }
    }
  }

  if (run_wham) {
    call wham.Whamg {
      input:
        bam_or_cram_file = reads_file_,
        bam_or_cram_index = reads_index_,
        sample_id = sample_id,
        reference_fasta = reference_fasta,
        reference_index = reference_index,
        include_bed_file = wham_include_list_bed_file,
        primary_contigs_list = primary_contigs_list,
        wham_docker = select_first([wham_docker]),
        runtime_attr_wham = runtime_attr_wham
    }
  }

  if (run_module_metrics) {
    call metrics.GatherSampleEvidenceMetrics {
      input:
        sample = sample_id,
        coverage_counts = CollectCounts.counts,
        pesr_disc = CollectSVEvidence.disc_out,
        pesr_split = CollectSVEvidence.split_out,
        manta_vcf = Manta.vcf,
        melt_vcf = MELT.vcf,
        scramble_vcf = Scramble.vcf,
        wham_vcf = Whamg.vcf,
        baseline_manta_vcf = baseline_manta_vcf,
        baseline_melt_vcf = baseline_melt_vcf,
        baseline_scramble_vcf = baseline_scramble_vcf,
        baseline_wham_vcf = baseline_wham_vcf,
        contig_list = primary_contigs_list,
        contig_index = select_first([primary_contigs_fai]),
        sv_pipeline_docker = sv_pipeline_docker
    }
  }

  output {
    File? coverage_counts = CollectCounts.counts

    File? manta_vcf = if run_manta then select_first([Manta.vcf]) else manta_vcf_input
    File? manta_index = if run_manta then select_first([Manta.index]) else manta_vcf_index_input

    File? melt_vcf = MELT.vcf
    File? melt_index = MELT.index
    Float? melt_coverage = MELT.coverage_out
    Int? melt_read_length = MELT.read_length_out
    Float? melt_insert_size = MELT.insert_size_out

    File? scramble_vcf = if run_scramble then select_first([ScrambleRealigned.vcf, Scramble.vcf]) else NONE_FILE_
    File? scramble_index = if run_scramble then select_first([ScrambleRealigned.index, Scramble.index]) else NONE_FILE_
    File? scramble_clusters = if run_scramble then select_first([ScrambleRealigned.clusters, Scramble.clusters]) else NONE_FILE_
    File? scramble_table = if run_scramble then select_first([ScrambleRealigned.table, Scramble.table]) else NONE_FILE_

    File? pesr_disc = CollectSVEvidence.disc_out
    File? pesr_disc_index = CollectSVEvidence.disc_out_index
    File? pesr_split = CollectSVEvidence.split_out
    File? pesr_split_index = CollectSVEvidence.split_out_index
    File? pesr_sd = CollectSVEvidence.sd_out
    File? pesr_sd_index = CollectSVEvidence.sd_out_index

    File? wham_vcf = Whamg.vcf
    File? wham_index = Whamg.index

    Array[File]? sample_metrics_files = GatherSampleEvidenceMetrics.sample_metrics_files
  }
}


task LocalizeReads {
  input {
    File reads_path
    File reads_index
    Boolean move_files = false
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = if move_files then size(reads_path, "GB") * 2 else size(reads_path, "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 3.75,
                                  disk_gb: ceil(50.0 + input_size),
                                  cpu_cores: 2,
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
    docker: "ubuntu:18.04"
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  Int disk_size = ceil(50 + size(reads_path, "GB"))

  command {
    set -exuo pipefail

    # When this pipeline is run on an HPC, moving files could lead to
    # moving the files from their original source, compared to moving
    # them from one directory of the VM to another when run on Cloud.
    # Therefore, to avoid moving files unexpectedly, we provide both
    # options for moving and copying, and set the copy as default.
    # Note that, when copying the files, the task can be slower depending
    # on the file size and IO performance and will need additional disk
    # space, hence it will be more expensive to run.

    if ~{move_files}; then
      mv ~{reads_path} $(basename ~{reads_path})
      mv ~{reads_index} $(basename ~{reads_index})
    else
      cp ~{reads_path} $(basename ~{reads_path})
      cp ~{reads_index} $(basename ~{reads_index})
    fi
  }
  output {
    File output_file = basename(reads_path)
    File output_index = basename(reads_index)
  }
}


task CheckAligner {
  input {
    File reads_path
    File reads_index
    File? reference_fasta
    File? reference_index
    File? reference_dict
    String sample_id

    String gatk_docker
    RuntimeAttr? runtime_attr_override
  }

  parameter_meta {
    reads_path: {
                 localization_optional: true
               }
  }

  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 1.0,
                               disk_gb: 10,
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File header = "~{sample_id}.header.sam"
    Int is_dragen_3_7_8 = read_int("is_dragen_3_7_8.txt")
  }
  command <<<
    set -euo pipefail

    gatk PrintReadsHeader \
      -I ~{reads_path} \
      --read-index ~{reads_index} \
      -O ~{sample_id}.header.sam \
      ~{"-R " + reference_fasta}

    awk '$0~"@PG" && $0~"ID: DRAGEN SW build" && $0~"VN: 05.021.604.3.7.8"' ~{sample_id}.header.sam \
      | wc -l \
      > is_dragen_3_7_8.txt
  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: gatk_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

task RealignSoftClippedReads {
  input {
    File reads_path
    File reads_index

    File scramble_table

    Boolean is_bam
    String sample_id
    File reference_fasta
    File reference_index

    File reference_bwa_alt
    File reference_bwa_amb
    File reference_bwa_ann
    File reference_bwa_bwt
    File reference_bwa_pac
    File reference_bwa_sa

    String sv_base_mini_docker
    Boolean use_ssd = true
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr runtime_default = object {
                               cpu_cores: 4,
                               mem_gb: 20,
                               disk_gb: ceil(10 + size(reads_path, "GB") * 2 + size([reference_bwa_bwt, reference_bwa_pac, reference_bwa_sa], "GB")),
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  String disk_type = if use_ssd then "SSD" else "HDD"

  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GiB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} ~{disk_type}"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_base_mini_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  output {
    File out = "~{sample_id}.realign_soft_clipped_reads.bam"
    File out_index = "~{sample_id}.realign_soft_clipped_reads.bam.bai"
  }
  command <<<
    set -exuo pipefail
    # Get insertion intervals
    zcat ~{scramble_table} \
      | sed 1d \
      | cut -f1 \
      | tr ':' '\t' \
      | awk -F'\t' -v OFS='\t' '{print $1,$2,$2+1}' \
      | sort -k1,1V -k2,2n \
      | bedtools slop -i - -g ~{reference_index} -b 150 \
      | bedtools merge \
      > intervals.bed
    mkdir tmp/
    samtools view --header-only ~{reads_path} > header.sam
    N_CORES=$(nproc)
    time samtools view --no-header \
      -T ~{reference_fasta} \
      -ML intervals.bed \
      ~{reads_path} \
      | awk -F'\t' -v OFS='\t' '$6~"S"' \
      | sort -u \
      | cat header.sam - \
      | samtools fastq \
      > reads.fastq
    bwa mem -H header.sam -K 100000000 -v 3 -t ${N_CORES} -Y ~{reference_fasta} reads.fastq \
      | samtools sort -T tmp \
      | samtools view -1 -h -O BAM -o ~{sample_id}.realign_soft_clipped_reads.bam
    samtools index -@${N_CORES} ~{sample_id}.realign_soft_clipped_reads.bam
  >>>
}
