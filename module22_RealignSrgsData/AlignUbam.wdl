version 1.0

import "Structs.wdl"

struct ReferenceFasta {
  File ref_fasta
  File ref_fasta_index
  File ref_dict
  File ref_alt
  File ref_amb
  File ref_ann
  File ref_bwt
  File ref_pac
  File ref_sa
}

workflow AlignUbam {
  input {
    # Parallel arrays with matched indices
    Array[File] input_bams
    Array[File?] fastq1s
    Array[File?] fastq2s

    String sample_name
    ReferenceFasta reference_fasta

    # BQSR Known Sites Inputs
    File dbSNP_vcf
    File dbSNP_vcf_index
    Array[File] known_indels_sites_VCFs
    Array[File] known_indels_sites_indices

    # Contamination Site Inputs
    File contamination_sites_ud
    File contamination_sites_bed
    File contamination_sites_mu

    String bwa_commandline = "bwa mem -K 100000000 -v 3 -t 16 -Y $bash_ref_fasta"
    Int compression_level = 2
    Boolean hard_clip_reads = false
    Boolean unmap_contaminant_reads = true
    Boolean allow_empty_ref_alt = false
    Int preemptible_tries = 3

    RuntimeAttr? runtime_attr_sam_to_fastq_override
    RuntimeAttr? runtime_attr_bwa_mem_override
    RuntimeAttr? runtime_attr_sort_sam_override
    RuntimeAttr? runtime_attr_merge_alignment_override
    RuntimeAttr? runtime_attr_mark_duplicates_override
  }

  # Process each readgroup in parallel
  scatter (i in range(length(input_bams))) {
    File input_bam = input_bams[i]
    File? fastq1   = fastq1s[i]
    File? fastq2   = fastq2s[i]

    Boolean has_paired_fastqs = defined(fastq1) && defined(fastq2)

    if (!has_paired_fastqs) {
      call SamToFastq {
        input:
          input_bam             = input_bam,
          preemptible_tries     = preemptible_tries,
          runtime_attr_override = runtime_attr_sam_to_fastq_override
      }
    }

    Array[File] bwa_input_fastqs = if has_paired_fastqs then select_all([fastq1, fastq2]) else select_all([SamToFastq.fastq])

    call BwaMem {
      input:
        input_bam             = input_bam,
        input_fastqs          = bwa_input_fastqs,
        bwa_commandline       = if has_paired_fastqs then bwa_commandline else bwa_commandline + " -p",
        reference_fasta       = reference_fasta,
        preemptible_tries     = preemptible_tries,
        allow_empty_ref_alt   = allow_empty_ref_alt,
        runtime_attr_override = runtime_attr_bwa_mem_override
    }

    call SortSamByQueryName as SortAlignedSam {
      input:
        input_sam_or_bam      = BwaMem.aligned_sam,
        preemptible_tries     = preemptible_tries,
        runtime_attr_override = runtime_attr_sort_sam_override
    }

    call MergeAlignment {
      input:
        input_bam               = input_bam,
        aligned_sam             = SortAlignedSam.sorted_bam,
        bwa_stderr_log          = BwaMem.bwa_stderr_log,
        bwa_version             = BwaMem.bwa_version,
        bwa_commandline         = bwa_commandline,
        reference_fasta         = reference_fasta,
        compression_level       = compression_level,
        hard_clip_reads         = hard_clip_reads,
        unmap_contaminant_reads = unmap_contaminant_reads,
        preemptible_tries       = preemptible_tries,
        runtime_attr_override   = runtime_attr_merge_alignment_override
    }

    call CollectUnsortedReadgroupBamQualityMetrics {
      input:
        input_bam         = MergeAlignment.output_bam,
        reference_fasta   = reference_fasta,
        preemptible_tries = preemptible_tries
    }
  }

  # Aggregate merged readgroup BAMs and Mark Duplicates
  call MarkDuplicates {
    input:
      input_bams            = MergeAlignment.output_bam,
      output_bam_basename  = sample_name,
      metrics_filename     = sample_name + ".duplicate_metrics",
      compression_level    = compression_level,
      preemptible_tries    = preemptible_tries,
      runtime_attr_override = runtime_attr_mark_duplicates_override
  }

  # Coordinate Sort of the sample-level duplicate-marked BAM
  call SortSampleBam {
    input:
      input_bam         = MarkDuplicates.output_bam,
      output_basename   = sample_name + ".aligned.duplicates_marked.sorted",
      preemptible_tries = preemptible_tries,
      runtime_attr_override = runtime_attr_sort_sam_override
  }

  # Create sequence groupings for scattered BaseRecalibrator
  call CreateSequenceGroupingTSV {
    input:
      ref_dict          = reference_fasta.ref_dict,
      preemptible_tries = preemptible_tries
  }

  # Estimate level of cross-sample contamination
  call CheckContamination {
    input:
      input_bam               = SortSampleBam.output_bam,
      input_bam_index         = SortSampleBam.output_bai,
      contamination_sites_ud  = contamination_sites_ud,
      contamination_sites_bed = contamination_sites_bed,
      contamination_sites_mu  = contamination_sites_mu,
      reference_fasta         = reference_fasta,
      output_prefix           = sample_name + ".preBqsr",
      preemptible_tries       = preemptible_tries
  }

  # BQSR Scattered over sequence groupings
  scatter (subgroup in CreateSequenceGroupingTSV.sequence_grouping) {
    call BaseRecalibrator {
      input:
        input_bam                   = SortSampleBam.output_bam,
        input_bam_index             = SortSampleBam.output_bai,
        recalibration_report_filename = sample_name + ".recal_data.csv",
        sequence_group_interval     = subgroup,
        dbSNP_vcf                   = dbSNP_vcf,
        dbSNP_vcf_index             = dbSNP_vcf_index,
        known_indels_sites_VCFs     = known_indels_sites_VCFs,
        known_indels_sites_indices  = known_indels_sites_indices,
        reference_fasta             = reference_fasta,
        preemptible_tries           = preemptible_tries
    }
  }

  # Gather scattered BQSR recalibration reports
  call GatherBqsrReports {
    input:
      input_bqsr_reports = BaseRecalibrator.recalibration_report,
      output_filename    = sample_name + ".recal_data.csv",
      preemptible_tries  = preemptible_tries
  }

  # Apply BQSR recalibration scattered over sequence groupings
  scatter (subgroup in CreateSequenceGroupingTSV.sequence_grouping) {
    call ApplyBQSR {
      input:
        input_bam               = SortSampleBam.output_bam,
        input_bam_index         = SortSampleBam.output_bai,
        output_bam_basename     = sample_name + ".recalibrated",
        recalibration_report    = GatherBqsrReports.output_bqsr_report,
        sequence_group_interval = subgroup,
        reference_fasta         = reference_fasta,
        compression_level       = compression_level,
        preemptible_tries       = preemptible_tries
    }
  }

  # Gather recalibrated BAM files into a single final BAM
  call GatherBamFiles {
    input:
      input_bams          = ApplyBQSR.recalibrated_bam,
      output_bam_basename = sample_name + ".aligned.duplicates_marked.recalibrated",
      compression_level   = compression_level,
      preemptible_tries   = preemptible_tries
  }

  output {
    File final_bam                          = GatherBamFiles.output_bam
    File final_bai                          = GatherBamFiles.output_bai
    File duplicate_metrics                  = MarkDuplicates.duplicate_metrics
    File contamination_metrics              = CheckContamination.contamination_metrics
    File bqsr_report                        = GatherBqsrReports.output_bqsr_report
    Array[File] readgroup_aligned_bams      = MergeAlignment.output_bam
    Array[File] readgroup_quality_metrics   = CollectUnsortedReadgroupBamQualityMetrics.quality_metrics
  }
}

# --- Tasks ---

task SamToFastq {
  input {
    File input_bam
    Int preemptible_tries
    RuntimeAttr? runtime_attr_override
  }

  Float unmapped_bam_size = size(input_bam, "GiB")
  Int disk_gb_default = ceil(unmapped_bam_size + (4.0 * unmapped_bam_size) + 20)
  Float mem_gb_default = 4.0

  RuntimeAttr runtime_attr_str_to_fastq_default = object {
    cpu_cores:          2,
    mem_gb:             mem_gb_default,
    disk_gb:            disk_gb_default,
    boot_disk_gb:       15,
    preemptible_tries:  preemptible_tries,
    max_retries:        1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, runtime_attr_str_to_fastq_default])

  String output_bam_basename = output_basename(input_bam, ".bam")

  command <<<
    set -o pipefail
    set -e

    java -Xms~{ceil(select_first([runtime_attr.mem_gb, mem_gb_default]) * 800)}m \
      -Xmx~{ceil(select_first([runtime_attr.mem_gb, mem_gb_default]) * 800)}m \
      -jar /usr/gitc/picard.jar \
      SamToFastq \
      INPUT=~{input_bam} \
      FASTQ=~{output_bam_basename}.interleaved.fastq \
      INTERLEAVE=true \
      NON_PF=true
  >>>

  runtime {
    docker:             "us.gcr.io/broad-gotc-prod/samtools-picard-bwa:1.0.2-0.7.15-2.26.10-1643840748"
    cpu:                select_first([runtime_attr.cpu_cores, runtime_attr_str_to_fastq_default.cpu_cores])
    memory:             select_first([runtime_attr.mem_gb, runtime_attr_str_to_fastq_default.mem_gb]) + " GiB"
    disks:              "local-disk " + select_first([runtime_attr.disk_gb, runtime_attr_str_to_fastq_default.disk_gb]) + " HDD"
    bootDiskSizeGb:     select_first([runtime_attr.boot_disk_gb, runtime_attr_str_to_fastq_default.boot_disk_gb])
    preemptible:        select_first([runtime_attr.preemptible_tries, runtime_attr_str_to_fastq_default.preemptible_tries])
    maxRetries:         select_first([runtime_attr.max_retries, runtime_attr_str_to_fastq_default.max_retries])
  }

  output {
    File fastq = "~{output_bam_basename}.interleaved.fastq"
  }
}

task BwaMem {
  input {
    File input_bam
    Array[File] input_fastqs
    String bwa_commandline
    String output_bam_basename
    ReferenceFasta reference_fasta
    Int preemptible_tries
    Boolean allow_empty_ref_alt = false
    RuntimeAttr? runtime_attr_override
  }

  Float fastq_size = size(input_fastqs, "GiB")
  Float ref_size = size(reference_fasta.ref_fasta, "GiB") + size(reference_fasta.ref_fasta_index, "GiB") + size(reference_fasta.ref_dict, "GiB")
  Float bwa_ref_size = ref_size + size(reference_fasta.ref_alt, "GiB") + size(reference_fasta.ref_amb, "GiB") + size(reference_fasta.ref_ann, "GiB") + size(reference_fasta.ref_bwt, "GiB") + size(reference_fasta.ref_pac, "GiB") + size(reference_fasta.ref_sa, "GiB")
  Int disk_gb_default = ceil(fastq_size + bwa_ref_size + (3.0 * fastq_size) + 20)
  Float mem_gb_scaled = (bwa_ref_size * 2.0) + 6.0
  Float mem_gb_default = if mem_gb_scaled > 14.0 then mem_gb_scaled else 14.0

  RuntimeAttr runtime_attr_bwa_mem_default = object {
    cpu_cores:          16,
    mem_gb:             mem_gb_default,
    disk_gb:            disk_gb_default,
    boot_disk_gb:       15,
    preemptible_tries:  preemptible_tries,
    max_retries:        1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, runtime_attr_bwa_mem_default])

  String output_bam_basename = output_basename(input_bam, ".bam")

  command <<<
    BWA_VERSION=$(/usr/gitc/bwa 2>&1 | grep -e '^Version' | sed 's/Version: //')
    set -o pipefail
    set -e

    if [ -z "${BWA_VERSION}" ]; then exit 1; fi
    echo "${BWA_VERSION}" > bwa_version.txt

    bash_ref_fasta=~{reference_fasta.ref_fasta}

    if [ -s ~{reference_fasta.ref_alt} ] || ~{allow_empty_ref_alt}; then
      /usr/gitc/~{bwa_commandline} \
        ~{sep=' ' input_fastqs} \
        > ~{output_bam_basename}.aligned.unsorted.bwa.sam \
        2> >(tee ~{output_bam_basename}.bwa.stderr.log >&2)
    else
      echo "ref_alt input is empty or not provided." >&2
      exit 1
    fi
  >>>

  runtime {
    docker:             "us.gcr.io/broad-gotc-prod/samtools-picard-bwa:1.0.2-0.7.15-2.26.10-1643840748"
    cpu:                select_first([runtime_attr.cpu_cores, runtime_attr_bwa_mem_default.cpu_cores])
    memory:             select_first([runtime_attr.mem_gb, runtime_attr_bwa_mem_default.mem_gb]) + " GiB"
    disks:              "local-disk " + select_first([runtime_attr.disk_gb, runtime_attr_bwa_mem_default.disk_gb]) + " HDD"
    bootDiskSizeGb:     select_first([runtime_attr.boot_disk_gb, runtime_attr_bwa_mem_default.boot_disk_gb])
    preemptible:        select_first([runtime_attr.preemptible_tries, runtime_attr_bwa_mem_default.preemptible_tries])
    maxRetries:         select_first([runtime_attr.max_retries, runtime_attr_bwa_mem_default.max_retries])
  }

  output {
    File aligned_sam    = "~{output_bam_basename}.aligned.unsorted.bwa.sam"
    File bwa_stderr_log = "~{output_bam_basename}.bwa.stderr.log"
    String bwa_version  = read_string("bwa_version.txt")
  }
}

task SortSamByQueryName {
  input {
    File input_sam_or_bam
    Int preemptible_tries
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size(input_sam_or_bam, "GiB")
  Int disk_gb_default = ceil(input_size + (2.5 * input_size) + 20)
  Float mem_gb_default = 8.0

  RuntimeAttr runtime_attr_sort_sam_default = object {
    cpu_cores:          2,
    mem_gb:             mem_gb_default,
    disk_gb:            disk_gb_default,
    boot_disk_gb:       15,
    preemptible_tries:  preemptible_tries,
    max_retries:        1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, runtime_attr_sort_sam_default])

  String output_basename = output_basename(input_sam_or_bam)
  command <<<
    set -o pipefail
    set -e

    java -Xms~{ceil(select_first([runtime_attr.mem_gb, mem_gb_default]) * 800)}m \
      -Xmx~{ceil(select_first([runtime_attr.mem_gb, mem_gb_default]) * 800)}m \
      -jar /usr/gitc/picard.jar \
      SortSam \
      INPUT=~{input_sam_or_bam} \
      OUTPUT=~{output_basename}.query_sorted.bam \
      SORT_ORDER=queryname
  >>>

  runtime {
    docker:             "us.gcr.io/broad-gotc-prod/samtools-picard-bwa:1.0.2-0.7.15-2.26.10-1643840748"
    cpu:                select_first([runtime_attr.cpu_cores, runtime_attr_sort_sam_default.cpu_cores])
    memory:             select_first([runtime_attr.mem_gb, runtime_attr_sort_sam_default.mem_gb]) + " GiB"
    disks:              "local-disk " + select_first([runtime_attr.disk_gb, runtime_attr_sort_sam_default.disk_gb]) + " HDD"
    bootDiskSizeGb:     select_first([runtime_attr.boot_disk_gb, runtime_attr_sort_sam_default.boot_disk_gb])
    preemptible:        select_first([runtime_attr.preemptible_tries, runtime_attr_sort_sam_default.preemptible_tries])
    maxRetries:         select_first([runtime_attr.max_retries, runtime_attr_sort_sam_default.max_retries])
  }

  output {
    File sorted_bam = "~{output_basename}.query_sorted.bam"
  }
}

task MergeAlignment {
  input {
    File input_bam
    File aligned_sam
    File bwa_stderr_log
    String bwa_version
    String bwa_commandline
    ReferenceFasta reference_fasta
    Int compression_level
    Int preemptible_tries
    Boolean hard_clip_reads = false
    Boolean unmap_contaminant_reads = true
    RuntimeAttr? runtime_attr_override
  }

  Float unmapped_bam_size = size(input_bam, "GiB")
  Float aligned_sam_size = size(aligned_sam, "GiB")
  Float ref_size = size(reference_fasta.ref_fasta, "GiB") + size(reference_fasta.ref_fasta_index, "GiB") + size(reference_fasta.ref_dict, "GiB")
  Int disk_gb_default = ceil(unmapped_bam_size + aligned_sam_size + ref_size + (2.5 * unmapped_bam_size) + 20)
  Float mem_gb_scaled = 4.0 + (aligned_sam_size / 50.0)
  Float mem_gb_default = if mem_gb_scaled > 5.0 then mem_gb_scaled else 5.0

  RuntimeAttr runtime_attr_merge_alignment_default = object {
    cpu_cores:          2,
    mem_gb:             mem_gb_default,
    disk_gb:            disk_gb_default,
    boot_disk_gb:       15,
    preemptible_tries:  preemptible_tries,
    max_retries:        1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, runtime_attr_merge_alignment_default])

  String output_bam_basename = basename(input_bam, ".bam")
  command <<<
    set -o pipefail
    set -e

    java -Dsamjdk.compression_level=~{compression_level} \
      -Xms~{ceil(select_first([runtime_attr.mem_gb, mem_gb_default]) * 800)}m \
      -Xmx~{ceil(select_first([runtime_attr.mem_gb, mem_gb_default]) * 800)}m \
      -jar /usr/gitc/picard.jar \
      MergeBamAlignment \
      VALIDATION_STRINGENCY=SILENT \
      EXPECTED_ORIENTATIONS=FR \
      ATTRIBUTES_TO_RETAIN=X0 \
      ATTRIBUTES_TO_REMOVE=NM \
      ATTRIBUTES_TO_REMOVE=MD \
      ALIGNED_BAM=~{aligned_sam} \
      UNMAPPED_BAM=~{input_bam} \
      OUTPUT=~{output_bam_basename}.bam \
      REFERENCE_SEQUENCE=~{reference_fasta.ref_fasta} \
      SORT_ORDER="unsorted" \
      IS_BISULFITE_SEQUENCE=false \
      ALIGNED_READS_ONLY=false \
      CLIP_ADAPTERS=false \
      ~{true='CLIP_OVERLAPPING_READS=true' false='' hard_clip_reads} \
      ~{true='CLIP_OVERLAPPING_READS_OPERATOR=H' false='' hard_clip_reads} \
      MAX_RECORDS_IN_RAM=2000000 \
      ADD_MATE_CIGAR=true \
      MAX_INSERTIONS_OR_DELETIONS=-1 \
      PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
      PROGRAM_RECORD_ID="bwamem" \
      PROGRAM_GROUP_VERSION="~{bwa_version}" \
      PROGRAM_GROUP_COMMAND_LINE="~{bwa_commandline}" \
      PROGRAM_GROUP_NAME="bwamem" \
      UNMAPPED_READ_STRATEGY=COPY_TO_TAG \
      ALIGNER_PROPER_PAIR_FLAGS=true \
      UNMAP_CONTAMINANT_READS=~{unmap_contaminant_reads} \
      ADD_PG_TAG_TO_READS=false
  >>>

  runtime {
    docker:             "us.gcr.io/broad-gotc-prod/samtools-picard-bwa:1.0.2-0.7.15-2.26.10-1643840748"
    cpu:                select_first([runtime_attr.cpu_cores, runtime_attr_merge_alignment_default.cpu_cores])
    memory:             select_first([runtime_attr.mem_gb, runtime_attr_merge_alignment_default.mem_gb]) + " GiB"
    disks:              "local-disk " + select_first([runtime_attr.disk_gb, runtime_attr_merge_alignment_default.disk_gb]) + " HDD"
    bootDiskSizeGb:     select_first([runtime_attr.boot_disk_gb, runtime_attr_merge_alignment_default.boot_disk_gb])
    preemptible:        select_first([runtime_attr.preemptible_tries, runtime_attr_merge_alignment_default.preemptible_tries])
    maxRetries:         select_first([runtime_attr.max_retries, runtime_attr_merge_alignment_default.max_retries])
  }

  output {
    File output_bam = "~{output_bam_basename}.bam"
  }
}

task CollectUnsortedReadgroupBamQualityMetrics {
  input {
    File input_bam
    ReferenceFasta reference_fasta
    Int preemptible_tries
  }

  String output_basename = basename(input_bam, ".bam")
  command <<<
    set -e
    java -Xms2000m -jar /usr/gitc/picard.jar \
      CollectMultipleMetrics \
      INPUT=~{input_bam} \
      OUTPUT=~{output_basename} \
      REFERENCE_SEQUENCE=~{reference_fasta.ref_fasta} \
      ASSUME_SORTED=true \
      PROGRAM=null \
      PROGRAM=QualityScoreDistribution \
      PROGRAM=MeanQualityByCycle \
      METRIC_ACCUMULATION_LEVEL=null \
      METRIC_ACCUMULATION_LEVEL=READ_GROUP
  >>>

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/samtools-picard-bwa:1.0.2-0.7.15-2.26.10-1643840748"
    memory: "3 GiB"
    cpu: 1
    disks: "local-disk 50 HDD"
    preemptible: preemptible_tries
  }

  output {
    File quality_metrics = "~{output_basename}.quality_distribution_metrics"
  }
}

task MarkDuplicates {
  input {
    Array[File] input_bams
    String output_bam_basename
    String metrics_filename
    Int compression_level
    Int preemptible_tries
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size(input_bams, "GiB")
  Int disk_gb_default = ceil(input_size * 3.0 + 30)
  Int default_mem_gb = ceil(input_size * 2.0 + 20)

  RuntimeAttr runtime_attr_md_default = object {
    cpu_cores:          2,
    mem_gb:             14.0,
    disk_gb:            disk_gb_default,
    boot_disk_gb:       15,
    preemptible_tries:  preemptible_tries,
    max_retries:        1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, runtime_attr_md_default])

  command <<<
    set -e
    java -Dsamjdk.compression_level=~{compression_level} \
      -Xms~{ceil(select_first([runtime_attr.mem_gb, default_mem_gb]) * 800)}m \
      -jar /usr/gitc/picard.jar \
      MarkDuplicates \
      INPUT=~{sep=' INPUT=' input_bams} \
      OUTPUT=~{output_bam_basename}.bam \
      METRICS_FILE=~{metrics_filename} \
      VALIDATION_STRINGENCY=SILENT \
      OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
      ASSUME_SORT_ORDER=queryname

  >>>

  runtime {
    docker:             "us.gcr.io/broad-gotc-prod/samtools-picard-bwa:1.0.2-0.7.15-2.26.10-1643840748"
    cpu:                select_first([runtime_attr.cpu_cores, runtime_attr_md_default.cpu_cores])
    memory:             select_first([runtime_attr.mem_gb, runtime_attr_md_default.mem_gb]) + " GiB"
    disks:              "local-disk " + select_first([runtime_attr.disk_gb, runtime_attr_md_default.disk_gb]) + " HDD"
    bootDiskSizeGb:     select_first([runtime_attr.boot_disk_gb, runtime_attr_md_default.boot_disk_gb])
    preemptible:        select_first([runtime_attr.preemptible_tries, runtime_attr_md_default.preemptible_tries])
    maxRetries:         select_first([runtime_attr.max_retries, runtime_attr_md_default.max_retries])
  }

  output {
    File output_bam        = "~{output_bam_basename}.bam"
    File duplicate_metrics = "~{metrics_filename}"
  }
}

task SortSampleBam {
  input {
    File input_bam
    String output_basename
    Int preemptible_tries
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size(input_bam, "GiB")
  Int disk_gb_default = ceil(input_size * 3.0 + 30)
  Int default_mem_gb = ceil(input_size * 2.0 + 20)

  RuntimeAttr runtime_attr_sort_default = object {
    cpu_cores:          2,
    mem_gb:             8.0,
    disk_gb:            disk_gb_default,
    boot_disk_gb:       15,
    preemptible_tries:  preemptible_tries,
    max_retries:        1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, runtime_attr_sort_default])

  command <<<
    set -e
    java -Xms~{ceil(select_first([runtime_attr.mem_gb, default_mem_gb]) * 800)}m \
      -jar /usr/gitc/picard.jar \
      SortSam \
      INPUT=~{input_bam} \
      OUTPUT=~{output_basename}.bam \
      SORT_ORDER=coordinate \
      CREATE_INDEX=true
  >>>

  runtime {
    docker:             "us.gcr.io/broad-gotc-prod/samtools-picard-bwa:1.0.2-0.7.15-2.26.10-1643840748"
    cpu:                select_first([runtime_attr.cpu_cores, runtime_attr_sort_default.cpu_cores])
    memory:             select_first([runtime_attr.mem_gb, runtime_attr_sort_default.mem_gb]) + " GiB"
    disks:              "local-disk " + select_first([runtime_attr.disk_gb, runtime_attr_sort_default.disk_gb]) + " HDD"
    bootDiskSizeGb:     select_first([runtime_attr.boot_disk_gb, runtime_attr_sort_default.boot_disk_gb])
    preemptible:        select_first([runtime_attr.preemptible_tries, runtime_attr_sort_default.preemptible_tries])
    maxRetries:         select_first([runtime_attr.max_retries, runtime_attr_sort_default.max_retries])
  }

  output {
    File output_bam = "~{output_basename}.bam"
    File output_bai = "~{output_basename}.bai"
  }
}

task CreateSequenceGroupingTSV {
  input {
    File ref_dict
    Int preemptible_tries
  }

  command <<<
    python3 <<SCRIPT
    with open("~{ref_dict}", "r") as ref_dict_file:
      sequence_tuple_list = []
      longest_sequence = 0
      for line in ref_dict_file:
        if line.startswith("@SQ"):
          line_split = line.split("\t")
          for token in line_split:
            if token.startswith("SN"):
              sequence_name = token.split(":")[1]
            if token.startswith("LN"):
              sequence_length = int(token.split(":")[1])
          sequence_tuple_list.append((sequence_name, sequence_length))
          if sequence_length > longest_sequence:
            longest_sequence = sequence_length

      sequence_tuple_list.sort(key=lambda x: x[1], reverse=True)

      tsv_group_list = []
      tsv_sub_group = []
      tsv_sub_group_old = []
      current_group_size = 0

      for sequence_tuple in sequence_tuple_list:
        if current_group_size + sequence_tuple[1] <= longest_sequence:
          tsv_sub_group.append(sequence_tuple[0])
          current_group_size += sequence_tuple[1]
        else:
          tsv_group_list.append(tsv_sub_group)
          tsv_sub_group = [sequence_tuple[0]]
          current_group_size = sequence_tuple[1]

      if tsv_sub_group:
        tsv_group_list.append(tsv_sub_group)

      with open("sequence_grouping.txt", "w") as sequence_grouping_file:
        for group in tsv_group_list:
          sequence_grouping_file.write("\t".join(group) + "\n")

      with open("sequence_grouping_with_unmapped.txt", "w") as sequence_grouping_with_unmapped_file:
        for group in tsv_group_list:
          sequence_grouping_with_unmapped_file.write("\t".join(group) + "\n")
        sequence_grouping_with_unmapped_file.write("unmapped\n")
    SCRIPT
  >>>

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/samtools-picard-bwa:1.0.2-0.7.15-2.26.10-1643840748"
    memory: "2 GiB"
    cpu: 1
    preemptible: preemptible_tries
  }

  output {
    Array[Array[String]] sequence_grouping               = read_tsv("sequence_grouping.txt")
    Array[Array[String]] sequence_grouping_with_unmapped = read_tsv("sequence_grouping_with_unmapped.txt")
  }
}

task CheckContamination {
  input {
    File input_bam
    File input_bam_index
    File contamination_sites_ud
    File contamination_sites_bed
    File contamination_sites_mu
    ReferenceFasta reference_fasta
    String output_prefix
    Int preemptible_tries
  }

  command <<<
    set -e
    /usr/gitc/verifyBamID2 \
      --DisableSanityCheck \
      --UDPath ~{contamination_sites_ud} \
      --BedPath ~{contamination_sites_bed} \
      --MeanPath ~{contamination_sites_mu} \
      --Reference ~{reference_fasta.ref_fasta} \
      --BamFile ~{input_bam} \
      --Output ~{output_prefix}
  >>>

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/samtools-picard-bwa:1.0.2-0.7.15-2.26.10-1643840748"
    memory: "4 GiB"
    cpu: 1
    disks: "local-disk 50 HDD"
    preemptible: preemptible_tries
  }

  output {
    File contamination_metrics = "~{output_prefix}.selfSM"
  }
}

task BaseRecalibrator {
  input {
    File input_bam
    File input_bam_index
    String recalibration_report_filename
    Array[String] sequence_group_interval
    File dbSNP_vcf
    File dbSNP_vcf_index
    Array[File] known_indels_sites_VCFs
    Array[File] known_indels_sites_indices
    ReferenceFasta reference_fasta
    Int preemptible_tries
  }

  command <<<
    set -e
    gatk --java-options "-Xms3000m" \
      BaseRecalibrator \
      -R ~{reference_fasta.ref_fasta} \
      -I ~{input_bam} \
      --use-original-qualities \
      -O ~{recalibration_report_filename} \
      --known-sites ~{dbSNP_vcf} \
      ~{sep=' ' prefix('--known-sites ', known_indels_sites_VCFs)} \
      -L ~{sep=' -L ' sequence_group_interval}
  >>>

  runtime {
    docker: "us.gcr.io/broad-gatc-prod/genomes-in-the-cloud:2.4.3-1564013544"
    memory: "4 GiB"
    cpu: 1
    disks: "local-disk 100 HDD"
    preemptible: preemptible_tries
  }

  output {
    File recalibration_report = "~{recalibration_report_filename}"
  }
}

task GatherBqsrReports {
  input {
    Array[File] input_bqsr_reports
    String output_filename
    Int preemptible_tries
  }

  command <<<
    set -e
    gatk --java-options "-Xms3000m" \
      GatherBQSRReports \
      -I ~{sep=' -I ' input_bqsr_reports} \
      -O ~{output_filename}
  >>>

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.3-1564013544"
    memory: "3.5 GiB"
    cpu: 1
    disks: "local-disk 30 HDD"
    preemptible: preemptible_tries
  }

  output {
    File output_bqsr_report = "~{output_filename}"
  }
}

task ApplyBQSR {
  input {
    File input_bam
    File input_bam_index
    String output_bam_basename
    File recalibration_report
    Array[String] sequence_group_interval
    ReferenceFasta reference_fasta
    Int compression_level
    Int preemptible_tries
  }

  command <<<
    set -e
    gatk --java-options "-Xms3000m" \
      ApplyBQSR \
      -R ~{reference_fasta.ref_fasta} \
      -I ~{input_bam} \
      -O ~{output_bam_basename}.bam \
      -bqsr ~{recalibration_report} \
      --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 \
      --add-output-sam-program-record \
      --use-original-qualities \
      -L ~{sep=' -L ' sequence_group_interval}
  >>>

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.3-1564013544"
    memory: "4 GiB"
    cpu: 1
    disks: "local-disk 100 HDD"
    preemptible: preemptible_tries
  }

  output {
    File recalibrated_bam = "~{output_bam_basename}.bam"
  }
}

task GatherBamFiles {
  input {
    Array[File] input_bams
    String output_bam_basename
    Int compression_level
    Int preemptible_tries
  }

  command <<<
    set -e
    java -Dsamjdk.compression_level=~{compression_level} -Xms2000m \
      -jar /usr/gitc/picard.jar \
      GatherBamFiles \
      INPUT=~{sep=' INPUT=' input_bams} \
      OUTPUT=~{output_bam_basename}.bam \
      CREATE_INDEX=true
  >>>

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/samtools-picard-bwa:1.0.2-0.7.15-2.26.10-1643840748"
    memory: "3 GiB"
    cpu: 1
    disks: "local-disk 100 HDD"
    preemptible: preemptible_tries
  }

  output {
    File output_bam = "~{output_bam_basename}.bam"
    File output_bai = "~{output_bam_basename}.bai"
  }
}