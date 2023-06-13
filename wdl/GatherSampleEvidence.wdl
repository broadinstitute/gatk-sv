version 1.0

import "Structs.wdl"
import "CollectCoverage.wdl" as cov
import "CollectSVEvidence.wdl" as coev
import "CramToBam.ReviseBase.wdl" as rb
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

    # Evidence collection flags
    Boolean collect_coverage = true
    Boolean collect_pesr = true

    # Convert ambiguous bases (e.g. K, S, Y, etc.) to N
    # Only use if encountering errors (expensive!)
    Boolean revise_base = false

    # Common parameters
    File primary_contigs_list
    File reference_fasta
    File reference_index    # Index (.fai), must be in same dir as fasta
    File reference_dict     # Dictionary (.dict), must be in same dir as fasta
    String? reference_version   # Either "38" or "19"

    # Coverage collection inputs
    File preprocessed_intervals
    Float? mem_gb_for_collect_counts
    Int? disk_space_gb_for_collect_counts

    # Manta inputs
    File manta_region_bed
    File? manta_region_bed_index
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

    # Wham inputs
    File wham_include_list_bed_file

    # Module metrics parameters
    # Run module metrics workflow at the end - on by default
    Boolean run_module_metrics = true
    File? primary_contigs_fai # required if run_module_metrics = true
    String? sv_pipeline_base_docker  # required if run_module_metrics = true
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
    RuntimeAttr? runtime_attr_revise_base
    RuntimeAttr? runtime_attr_concat_bam
    RuntimeAttr? runtime_attr_manta
    RuntimeAttr? runtime_attr_melt_coverage
    RuntimeAttr? runtime_attr_melt_metrics
    RuntimeAttr? runtime_attr_melt
    RuntimeAttr? runtime_attr_scramble
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
  call LocalizeReads {
    input:
      reads_path = bam_or_cram_file,
      reads_index = bam_or_cram_index_,
      runtime_attr_override = runtime_attr_localize_reads
  }

  if (revise_base) {
    call rb.CramToBamReviseBase {
      input:
        cram_file = LocalizeReads.output_file,
        cram_index = LocalizeReads.output_index,
        reference_fasta = reference_fasta,
        reference_index = reference_index,
        contiglist = select_first([primary_contigs_fai]),
        samtools_cloud_docker = samtools_cloud_docker,
        runtime_attr_split_cram = runtime_attr_split_cram,
        runtime_attr_revise_base = runtime_attr_revise_base,
        runtime_attr_concat_bam = runtime_attr_concat_bam
    }
  }

  File reads_file_ = select_first([CramToBamReviseBase.bam_file, LocalizeReads.output_file])
  File reads_index_ = select_first([CramToBamReviseBase.bam_index, LocalizeReads.output_index])

  if (collect_coverage || run_melt || run_module_metrics) {
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

  if (collect_pesr || run_module_metrics) {
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
    call scramble.Scramble {
      input:
        bam_or_cram_file = reads_file_,
        bam_or_cram_index = reads_index_,
        sample_name = sample_id,
        reference_fasta = reference_fasta,
        detect_deletions = false,
        scramble_docker = select_first([scramble_docker]),
        runtime_attr_scramble = runtime_attr_scramble
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
        sv_pipeline_base_docker = select_first([sv_pipeline_base_docker])
    }
  }

  output {
    File? coverage_counts = CollectCounts.counts

    File? manta_vcf = Manta.vcf
    File? manta_index = Manta.index

    File? melt_vcf = MELT.vcf
    File? melt_index = MELT.index
    Float? melt_coverage = MELT.coverage_out
    Int? melt_read_length = MELT.read_length_out
    Float? melt_insert_size = MELT.insert_size_out

    File? scramble_vcf = Scramble.vcf
    File? scramble_index = Scramble.index

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
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size(reads_path, "GB")
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
    ln -s ~{reads_path}
    ln -s ~{reads_index}
  }
  output {
    File output_file = basename(reads_path)
    File output_index = basename(reads_index)
  }
}
