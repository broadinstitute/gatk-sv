version 1.0

import "Module00a.wdl" as m00a

workflow Module00aBatch {
  input {
    Array[File] bam_or_cram_files
    Array[File]? bam_or_cram_indexes
    Array[String] sample_ids

    # Use only for crams in requester pays buckets
    Boolean requester_pays_crams = false

    # Caller flags
    Boolean collect_coverage = true
    Boolean collect_pesr = true

    # If true, any intermediate BAM files will be deleted after the algorithms have completed.
    # NOTE: If the workflow (ie any algorithm) fails, the bam will NOT be deleted.
    Boolean delete_intermediate_bam = false

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

    # Delly inputs
    File? delly_blacklist_intervals_file  # Required if run_delly True
    Array[String]? delly_sv_types

    # Manta inputs
    File manta_region_bed
    File? manta_region_bed_index
    Float? manta_jobs_per_cpu
    Int? manta_mem_gb_per_job

    # Melt inputs
    File? melt_standard_vcf_header # required if use_melt True
    File? melt_metrics_intervals
    Array[Float]? insert_size
    Array[Int]? read_length
    Array[Float]? coverage
    File? metrics_intervals
    Array[Float]? pct_chimeras
    Array[Float]? total_reads
    Array[Int]? pf_reads_improper_pairs

    # Wham inputs
    File wham_whitelist_bed_file

    # Docker
    String sv_pipeline_docker
    String sv_base_mini_docker
    String samtools_cloud_docker
    String? delly_docker
    String? manta_docker
    String? melt_docker
    String? wham_docker
    String gatk_docker
    String? gatk_docker_pesr_override
    String genomes_in_the_cloud_docker

    # Runtime configuration overrides
    RuntimeAttr? runtime_attr_merge_vcfs
    RuntimeAttr? runtime_attr_baf_sample
    RuntimeAttr? runtime_attr_cram_to_bam
    RuntimeAttr? runtime_attr_delly
    RuntimeAttr? runtime_attr_delly_gather
    RuntimeAttr? runtime_attr_manta
    RuntimeAttr? runtime_attr_melt_coverage
    RuntimeAttr? runtime_attr_melt_metrics
    RuntimeAttr? runtime_attr_melt
    RuntimeAttr? runtime_attr_pesr
    RuntimeAttr? runtime_attr_wham
    RuntimeAttr? runtime_attr_wham_whitelist

    File? NONE_FILE_
    String? NONE_STRING_
  }

  scatter (i in range(length(bam_or_cram_files))) {
    call m00a.Module00a {
      input:
        bam_or_cram_file = bam_or_cram_files[i],
        bam_or_cram_index = if defined(bam_or_cram_indexes) then select_first([bam_or_cram_indexes])[i] else NONE_FILE_,
        sample_id = sample_ids[i],
        requester_pays_crams = requester_pays_crams,
        collect_coverage = collect_coverage,
        collect_pesr = collect_pesr,
        delete_intermediate_bam = delete_intermediate_bam,
        primary_contigs_list = primary_contigs_list,
        reference_fasta = reference_fasta,
        reference_index = reference_index,
        reference_dict = reference_dict,
        reference_version = reference_version,
        preprocessed_intervals = preprocessed_intervals,
        mem_gb_for_collect_counts = mem_gb_for_collect_counts,
        disk_space_gb_for_collect_counts = disk_space_gb_for_collect_counts,
        delly_blacklist_intervals_file = delly_blacklist_intervals_file,
        delly_sv_types = delly_sv_types,
        manta_region_bed = manta_region_bed,
        manta_region_bed_index = manta_region_bed_index,
        manta_jobs_per_cpu = manta_jobs_per_cpu,
        manta_mem_gb_per_job = manta_mem_gb_per_job,
        melt_standard_vcf_header = melt_standard_vcf_header,
        melt_metrics_intervals = melt_metrics_intervals,
        insert_size = insert_size,
        read_length = read_length,
        coverage = coverage,
        metrics_intervals = metrics_intervals,
        pct_chimeras = pct_chimeras,
        total_reads = total_reads,
        pf_reads_improper_pairs = pf_reads_improper_pairs,
        wham_whitelist_bed_file = wham_whitelist_bed_file,
        sv_pipeline_docker = sv_pipeline_docker,
        sv_base_mini_docker = sv_base_mini_docker,
        samtools_cloud_docker = samtools_cloud_docker,
        delly_docker = delly_docker,
        manta_docker = manta_docker,
        melt_docker = melt_docker,
        wham_docker = wham_docker,
        gatk_docker = gatk_docker,
        gatk_docker_pesr_override = gatk_docker_pesr_override,
        genomes_in_the_cloud_docker = genomes_in_the_cloud_docker,
        runtime_attr_merge_vcfs = runtime_attr_merge_vcfs,
        runtime_attr_cram_to_bam = runtime_attr_cram_to_bam,
        runtime_attr_delly = runtime_attr_delly,
        runtime_attr_delly_gather = runtime_attr_delly_gather,
        runtime_attr_manta = runtime_attr_manta,
        runtime_attr_melt_coverage = runtime_attr_melt_coverage,
        runtime_attr_melt_metrics = runtime_attr_melt_metrics,
        runtime_attr_melt = runtime_attr_melt,
        runtime_attr_pesr = runtime_attr_pesr,
        runtime_attr_wham = runtime_attr_wham,
        runtime_attr_wham_whitelist = runtime_attr_wham_whitelist
    }
  }

  output {
    Array[File?] coverage_counts = Module00a.coverage_counts

    Array[File?] delly_vcf = Module00a.delly_vcf
    Array[File?] delly_index = Module00a.delly_index

    Array[File?] manta_vcf = Module00a.manta_vcf
    Array[File?] manta_index = Module00a.manta_index

    Array[File?] melt_vcf = Module00a.melt_vcf
    Array[File?] melt_index = Module00a.melt_index
    Array[Float?] melt_coverage = Module00a.melt_coverage
    Array[Int?] melt_read_length = Module00a.melt_read_length
    Array[Float?] melt_insert_size = Module00a.melt_insert_size

    Array[File?] pesr_disc = Module00a.pesr_disc
    Array[File?] pesr_disc_index = Module00a.pesr_disc_index
    Array[File?] pesr_split = Module00a.pesr_split
    Array[File?] pesr_split_index = Module00a.pesr_split_index

    Array[File?] wham_vcf = Module00a.wham_vcf
    Array[File?] wham_index = Module00a.wham_index
  }
}