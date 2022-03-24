version 1.0

import "GatherSampleEvidence.wdl" as sampleevidence
import "TestUtils.wdl" as tu

workflow GatherSampleEvidenceBatch {
  input {
    Array[File] bam_or_cram_files
    Array[File]? bam_or_cram_indexes
    Array[String] sample_ids

    # Use only for crams in requester pays buckets
    Boolean requester_pays_crams = false

    # Use to revise Y, R, W, S, K, M, D, H, V, B, X bases in BAM to N. Use only if providing a CRAM file as input 
    # May be more expensive - use only if necessary
    Boolean revise_base_cram_to_bam = false
    File? primary_contigs_fai # required if using revise_base_cram_to_bam (or if run_module_metrics = true)

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
    File? delly_exclude_intervals_file  # Required if run_delly True
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
    File wham_include_list_bed_file

    # Module metrics parameters
    # Run module metrics workflow at the end - on by default
    Boolean? run_module_metrics
    String? batch  # required if run_module_metrics = true
    String? sv_pipeline_base_docker  # required if run_module_metrics = true
    String? linux_docker  # required if run_module_metrics = true
    File? baseline_delly_vcf  # baseline files are optional for metrics workflow
    File? baseline_manta_vcf
    File? baseline_wham_vcf
    File? baseline_melt_vcf
    File? baseline_scramble_vcf

    # Docker
    String sv_pipeline_docker
    String sv_base_mini_docker
    String samtools_cloud_docker
    String? delly_docker
    String? manta_docker
    String? melt_docker
    String? scramble_docker
    String? wham_docker
    String gatk_docker
    String? gatk_docker_pesr_override
    String genomes_in_the_cloud_docker
    String cloud_sdk_docker

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
    RuntimeAttr? runtime_attr_scramble
    RuntimeAttr? runtime_attr_pesr
    RuntimeAttr? runtime_attr_wham
    RuntimeAttr? runtime_attr_wham_include_list
    RuntimeAttr? runtime_attr_ReviseBaseInBam
    RuntimeAttr? runtime_attr_ConcatBam
    RuntimeAttr? runtime_attr_cat_metrics

    File? NONE_FILE_
    String? NONE_STRING_
    Float? NONE_FLOAT_
    Int? NONE_INT_
  }

  Boolean run_module_metrics_ = if defined(run_module_metrics) then select_first([run_module_metrics]) else true

  scatter (i in range(length(bam_or_cram_files))) {
    call sampleevidence.GatherSampleEvidence {
      input:
        bam_or_cram_file = bam_or_cram_files[i],
        bam_or_cram_index = if defined(bam_or_cram_indexes) then select_first([bam_or_cram_indexes])[i] else NONE_FILE_,
        sample_id = sample_ids[i],
        requester_pays_crams = requester_pays_crams,
        revise_base_cram_to_bam = revise_base_cram_to_bam,
        collect_coverage = collect_coverage,
        collect_pesr = collect_pesr,
        delete_intermediate_bam = delete_intermediate_bam,
        primary_contigs_list = primary_contigs_list,
        primary_contigs_fai = primary_contigs_fai,
        reference_fasta = reference_fasta,
        reference_index = reference_index,
        reference_dict = reference_dict,
        reference_version = reference_version,
        preprocessed_intervals = preprocessed_intervals,
        mem_gb_for_collect_counts = mem_gb_for_collect_counts,
        disk_space_gb_for_collect_counts = disk_space_gb_for_collect_counts,
        delly_exclude_intervals_file = delly_exclude_intervals_file,
        delly_sv_types = delly_sv_types,
        manta_region_bed = manta_region_bed,
        manta_region_bed_index = manta_region_bed_index,
        manta_jobs_per_cpu = manta_jobs_per_cpu,
        manta_mem_gb_per_job = manta_mem_gb_per_job,
        melt_standard_vcf_header = melt_standard_vcf_header,
        melt_metrics_intervals = melt_metrics_intervals,
        insert_size = if defined(insert_size) then select_first([insert_size])[i] else NONE_FLOAT_,
        read_length = if defined(read_length) then select_first([read_length])[i] else NONE_INT_,
        coverage = if defined(coverage) then select_first([coverage])[i] else NONE_FLOAT_,
        metrics_intervals = metrics_intervals,
        pct_chimeras = if defined(pct_chimeras) then select_first([pct_chimeras])[i] else NONE_FLOAT_,
        total_reads = if defined(total_reads) then select_first([total_reads])[i] else NONE_FLOAT_,
        pf_reads_improper_pairs = if defined(pf_reads_improper_pairs) then select_first([pf_reads_improper_pairs])[i] else NONE_INT_,
        wham_include_list_bed_file = wham_include_list_bed_file,
        run_module_metrics = run_module_metrics_,
        sv_pipeline_base_docker = sv_pipeline_base_docker,
        baseline_delly_vcf = baseline_delly_vcf,
        baseline_manta_vcf = baseline_manta_vcf,
        baseline_melt_vcf = baseline_melt_vcf,
        baseline_scramble_vcf = baseline_scramble_vcf,
        baseline_wham_vcf = baseline_wham_vcf,
        sv_pipeline_docker = sv_pipeline_docker,
        sv_base_mini_docker = sv_base_mini_docker,
        samtools_cloud_docker = samtools_cloud_docker,
        delly_docker = delly_docker,
        manta_docker = manta_docker,
        melt_docker = melt_docker,
        scramble_docker = scramble_docker,
        wham_docker = wham_docker,
        gatk_docker = gatk_docker,
        gatk_docker_pesr_override = gatk_docker_pesr_override,
        genomes_in_the_cloud_docker = genomes_in_the_cloud_docker,
        cloud_sdk_docker = cloud_sdk_docker,
        runtime_attr_merge_vcfs = runtime_attr_merge_vcfs,
        runtime_attr_cram_to_bam = runtime_attr_cram_to_bam,
        runtime_attr_delly = runtime_attr_delly,
        runtime_attr_delly_gather = runtime_attr_delly_gather,
        runtime_attr_manta = runtime_attr_manta,
        runtime_attr_melt_coverage = runtime_attr_melt_coverage,
        runtime_attr_melt_metrics = runtime_attr_melt_metrics,
        runtime_attr_melt = runtime_attr_melt,
        runtime_attr_scramble = runtime_attr_scramble,
        runtime_attr_pesr = runtime_attr_pesr,
        runtime_attr_wham = runtime_attr_wham,
        runtime_attr_wham_include_list = runtime_attr_wham_include_list,
        runtime_attr_ReviseBaseInBam = runtime_attr_ReviseBaseInBam,
        runtime_attr_ConcatBam = runtime_attr_ConcatBam
    }
  }

  if (run_module_metrics_) {
    Array[Array[File]] sample_metrics_files_ = select_all(GatherSampleEvidence.sample_metrics_files)
    call tu.CatMetrics {
      input:
        prefix = "GatherSampleEvidence." + select_first([batch]),
        metric_files = flatten(sample_metrics_files_),
        linux_docker = select_first([linux_docker]),
        runtime_attr_override = runtime_attr_cat_metrics
    }
  }

  output {
    Array[File?] coverage_counts = GatherSampleEvidence.coverage_counts

    Array[File?] delly_vcf = GatherSampleEvidence.delly_vcf
    Array[File?] delly_index = GatherSampleEvidence.delly_index

    Array[File?] manta_vcf = GatherSampleEvidence.manta_vcf
    Array[File?] manta_index = GatherSampleEvidence.manta_index

    Array[File?] melt_vcf = GatherSampleEvidence.melt_vcf
    Array[File?] melt_index = GatherSampleEvidence.melt_index
    Array[Float?] melt_coverage = GatherSampleEvidence.melt_coverage
    Array[Int?] melt_read_length = GatherSampleEvidence.melt_read_length
    Array[Float?] melt_insert_size = GatherSampleEvidence.melt_insert_size

    Array[File?] scramble_vcf = GatherSampleEvidence.scramble_vcf
    Array[File?] scramble_index = GatherSampleEvidence.scramble_index

    Array[File?] pesr_disc = GatherSampleEvidence.pesr_disc
    Array[File?] pesr_disc_index = GatherSampleEvidence.pesr_disc_index
    Array[File?] pesr_split = GatherSampleEvidence.pesr_split
    Array[File?] pesr_split_index = GatherSampleEvidence.pesr_split_index

    Array[File?] wham_vcf = GatherSampleEvidence.wham_vcf
    Array[File?] wham_index = GatherSampleEvidence.wham_index

    File? metrics_file_sampleevidence = CatMetrics.out
  }
}
