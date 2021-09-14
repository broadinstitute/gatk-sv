version 1.0

import "Structs.wdl"
import "CollectCoverage.wdl" as cov
import "CramToBam.wdl" as ctb
import "CramToBam.ReviseBase.wdl" as ctb_revise
import "Delly.wdl" as delly
import "Manta.wdl" as manta
import "MELT.wdl" as melt
import "GatherSampleEvidenceMetrics.wdl" as metrics
import "PESRCollection.wdl" as pesr
import "Whamg.wdl" as wham

# Runs selected tools on BAM/CRAM files

workflow GatherSampleEvidence {
  input {
    File bam_or_cram_file
    File? bam_or_cram_index

    # Use only for crams in requester pays buckets
    Boolean requester_pays_crams = false

    # Use to revise Y, R, W, S, K, M, D, H, V, B, X bases in BAM to N. Use only if providing a CRAM file as input 
    # May be more expensive - use only if necessary
    Boolean revise_base_cram_to_bam = false
    File? primary_contigs_fai # required if using revise_base_cram_to_bam (or if run_module_metrics = true)

    # Note: raw and "safe" CRAM/BAM IDs can be generated with GetSampleID
    String sample_id

    # Evidence collection flags
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
    Boolean? run_module_metrics
    String? sv_pipeline_base_docker  # required if run_module_metrics = true
    File? baseline_delly_vcf  # baseline files are optional for metrics workflow
    File? baseline_manta_vcf
    File? baseline_wham_vcf
    File? baseline_melt_vcf

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
    RuntimeAttr? runtime_attr_pesr
    RuntimeAttr? runtime_attr_wham
    RuntimeAttr? runtime_attr_wham_include_list
    RuntimeAttr? runtime_attr_ReviseBaseInBam
    RuntimeAttr? runtime_attr_ConcatBam

    # Never assign these values! (workaround until None type is implemented)
    Float? NONE_FLOAT_
    Int? NONE_INT_
    File? NONE_FILE_
  }

  Boolean run_delly = defined(delly_docker)
  Boolean run_manta = defined(manta_docker)
  Boolean run_melt = defined(melt_docker)
  Boolean run_wham = defined(wham_docker)

  Boolean is_bam_ = basename(bam_or_cram_file, ".bam") + ".bam" == basename(bam_or_cram_file)
  String index_ext_ = if is_bam_ then ".bai" else ".crai"
  File bam_or_cram_index_ = if defined(bam_or_cram_index) then select_first([bam_or_cram_index]) else bam_or_cram_file + index_ext_

  # Convert to BAM if we have a CRAM
  if (!is_bam_) {
    if (!revise_base_cram_to_bam) {
      call ctb.CramToBam {
        input:
          cram_file = bam_or_cram_file,
          reference_fasta = reference_fasta,
          reference_index = reference_index,
          requester_pays = requester_pays_crams,
          samtools_cloud_docker = samtools_cloud_docker,
          runtime_attr_override = runtime_attr_cram_to_bam
      }
    }
    if (revise_base_cram_to_bam) {
      call ctb_revise.CramToBamReviseBase {
        input:
          cram_file = bam_or_cram_file,
          cram_index = bam_or_cram_index_,
          reference_fasta = reference_fasta,
          reference_index = reference_index,
          requester_pays = requester_pays_crams,
          contiglist = select_first([primary_contigs_fai]),
          samtools_cloud_docker = samtools_cloud_docker,
          runtime_attr_override = runtime_attr_cram_to_bam,
          runtime_attr_ReviseBaseInBam = runtime_attr_ReviseBaseInBam,
          runtime_attr_ConcatBam = runtime_attr_ConcatBam
      }
    }
  }

  File bam_file_ = select_first([CramToBam.bam_file, CramToBamReviseBase.bam_file, bam_or_cram_file])
  File bam_index_ = select_first([CramToBam.bam_index, CramToBamReviseBase.bam_index, bam_or_cram_index_])

  if (collect_coverage) {
    call cov.CollectCounts {
      input:
        intervals = preprocessed_intervals,
        bam = bam_file_,
        bam_idx = bam_index_,
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

  if (run_delly) {
    call delly.Delly {
      input:
        bam_or_cram_file = bam_file_,
        bam_or_cram_index = bam_index_,
        sample_id = sample_id,
        reference_fasta = reference_fasta,
        reference_index = reference_index,
        exclude_intervals_file = select_first([delly_exclude_intervals_file]),
        sv_types = delly_sv_types,
        sv_base_mini_docker = sv_base_mini_docker,
        delly_docker = select_first([delly_docker]),
        runtime_attr_delly = runtime_attr_delly,
        runtime_attr_gather = runtime_attr_delly_gather
    }
  }

  if (run_manta) {
    call manta.Manta {
      input:
        bam_or_cram_file = bam_file_,
        bam_or_cram_index = bam_index_,
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
    call pesr.PESRCollection {
      input:
        cram = bam_file_,
        cram_index = bam_index_,
        sample_id = sample_id,
        reference_fasta = reference_fasta,
        reference_index = reference_index,
        reference_dict = reference_dict,
        gatk_docker = select_first([gatk_docker_pesr_override, gatk_docker]),
        runtime_attr_override = runtime_attr_pesr
    }
  }

  if (run_melt) {
    call melt.MELT {
      input:
        bam_or_cram_file = bam_file_,
        bam_or_cram_index = bam_index_,
        counts_file = select_first([CollectCounts.counts]),
        sample_id = sample_id,
        reference_fasta = reference_fasta,
        reference_index = reference_index,
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

  if (run_wham) {
    call wham.Whamg {
      input:
        bam_or_cram_file = bam_file_,
        bam_or_cram_index = bam_index_,
        sample_id = sample_id,
        reference_fasta = reference_fasta,
        reference_index = reference_index,
        include_bed_file = wham_include_list_bed_file,
        chr_file = primary_contigs_list,
        samtools_cloud_docker = samtools_cloud_docker,
        wham_docker = select_first([wham_docker]),
        runtime_attr_includelist = runtime_attr_wham_include_list,
        runtime_attr_wham = runtime_attr_wham
    }
  }

  # Avoid storage costs
  if (!is_bam_) {
    if (delete_intermediate_bam) {
      Array[File] ctb_dummy = select_all([CollectCounts.counts, Delly.vcf, Manta.vcf, PESRCollection.disc_out, PESRCollection.split_out, MELT.vcf, Whamg.vcf])
      call DeleteIntermediateFiles {
        input:
          intermediates = select_all([CramToBam.bam_file, MELT.filtered_bam]),
          dummy = ctb_dummy,
          cloud_sdk_docker = cloud_sdk_docker
      }
    }
  }

  Boolean run_module_metrics_ = if defined(run_module_metrics) then select_first([run_module_metrics]) else true
  if (run_module_metrics_) {
    call metrics.GatherSampleEvidenceMetrics {
      input:
        sample = sample_id,
        coverage_counts = CollectCounts.counts,
        pesr_disc = PESRCollection.disc_out,
        pesr_split = PESRCollection.split_out,
        delly_vcf = Delly.vcf,
        manta_vcf = Manta.vcf,
        melt_vcf = MELT.vcf,
        wham_vcf = Whamg.vcf,
        baseline_delly_vcf = baseline_delly_vcf,
        baseline_manta_vcf = baseline_manta_vcf,
        baseline_melt_vcf = baseline_melt_vcf,
        baseline_wham_vcf = baseline_wham_vcf,
        contig_list = primary_contigs_list,
        contig_index = select_first([primary_contigs_fai]),
        sv_pipeline_base_docker = select_first([sv_pipeline_base_docker])
    }
  }

  output {
    File? coverage_counts = CollectCounts.counts

    File? delly_vcf = Delly.vcf
    File? delly_index = Delly.index

    File? manta_vcf = Manta.vcf
    File? manta_index = Manta.index

    File? melt_vcf = MELT.vcf
    File? melt_index = MELT.index
    Float? melt_coverage = MELT.coverage_out
    Int? melt_read_length = MELT.read_length_out
    Float? melt_insert_size = MELT.insert_size_out

    File? pesr_disc = PESRCollection.disc_out
    File? pesr_disc_index = PESRCollection.disc_out_index
    File? pesr_split = PESRCollection.split_out
    File? pesr_split_index = PESRCollection.split_out_index

    File? wham_vcf = Whamg.vcf
    File? wham_index = Whamg.index

    Array[File]? sample_metrics_files = GatherSampleEvidenceMetrics.sample_metrics_files
  }
}

task GetBamID {
  input {
    File bam_file
    String samtools_cloud_docker
  }

  parameter_meta {
    bam_file: {
      localization_optional: true
    }
  }

  output {
    String out = read_lines("sample.txt")[0]
  }
  command <<<
    set -euo pipefail
    export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
    samtools view -H ~{bam_file} \
      | grep "^@RG" \
      | awk -F "\t" '{for(i=2;i<=NF;i++){if($i~/^SM:/){a=$i}} print substr(a,4)}' \
      | sort \
      | uniq \
      > sample.txt
    NUM_IDS=$(wc -l < sample.txt)
    if [[ ${NUM_IDS} -eq 0 ]]; then
      echo "No sample IDs were found in the BAM header"
      exit 1
    fi
    if [[ ${NUM_IDS} -gt 1 ]]; then
      echo "Multiple sample IDs were found in the BAM header"
      exit 1
    fi
  >>>
  runtime {
    docker: samtools_cloud_docker
    memory: "1 GB"
    cpu: "1"
    disks: "local-disk 10 HDD"
    preemptible: "3"
    maxRetries: "1"
  }
}

task InternalSampleID {
  input {
    String external_sample_id
    Int hash_length
    String sv_pipeline_docker
  }

  output {
    String out = read_lines("external_id.txt")[0]
  }
  command <<<
    set -euo pipefail
    HASH=$(echo -n "~{external_sample_id}" | openssl sha1 | awk '{print substr($2,0,~{hash_length})}')
    SAFE_ID=$(echo -n "~{external_sample_id}" | sed 's/[^a-zA-Z0-9]/_/g')
    echo "__${SAFE_ID}__${HASH}" > external_id.txt
  >>>
  runtime {
    docker: sv_pipeline_docker
    memory: "1 GB"
    cpu: "1"
    disks: "local-disk 10 HDD"
    preemptible: "3"
    maxRetries: "1"
  }
}

task DeleteIntermediateFiles {
  input {
    Array[File] intermediates
    Array[File]? dummy # Pass in outputs that must be complete before cleanup
    String cloud_sdk_docker # Cloud provider SDK docker image. For GCP use "google/cloud-sdk" and for AWS use "amazon/aws-cli"
  }
  parameter_meta {
    intermediates: {
      localization_optional: true
    }
    dummy: {
      localization_optional: true
    }
  }

  command <<<
    {
      gsutil rm -I < ~{write_lines(intermediates)}
    } || {
      while read line; do
        aws s3 rm $line
      done < ${write_lines(intermediates)}
    }
  >>>
  runtime {
    docker: cloud_sdk_docker
    memory: "1 GB"
    cpu: "1"
    disks: "local-disk 10 HDD"
    preemptible: "3"
    maxRetries: "1"
  }
}
