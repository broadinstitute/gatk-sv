version 1.0

## Copyright (structure adapted from Broad Institute WARP, tasks/wdl/Alignment.wdl,
## tag WholeGenomeReprocessing_v3.3.7:
## https://github.com/broadinstitute/warp/blob/WholeGenomeReprocessing_v3.3.7/tasks/wdl/Alignment.wdl

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
    File input_bam
    String output_bam_basename
    String bwa_commandline
    ReferenceFasta reference_fasta
    Int compression_level = 2
    Boolean hard_clip_reads = false
    Boolean unmap_contaminant_reads = true
    Boolean allow_empty_ref_alt = false
    Int preemptible_tries = 3
  }

  call SamToFastq {
    input:
      input_bam           = input_bam,
      output_bam_basename = output_bam_basename,
      preemptible_tries   = preemptible_tries
  }

  call BwaMem {
    input:
      input_fastq          = SamToFastq.fastq,
      output_bam_basename  = output_bam_basename,
      bwa_commandline      = bwa_commandline,
      reference_fasta      = reference_fasta,
      preemptible_tries    = preemptible_tries,
      allow_empty_ref_alt  = allow_empty_ref_alt
  }

  call SortSamByQueryName {
    input:
      aligned_sam         = BwaMem.aligned_sam,
      output_bam_basename = output_bam_basename,
      preemptible_tries   = preemptible_tries
  }

  call MergeAlignment {
    input:
      input_bam               = input_bam,
      aligned_sam             = SortSamByQueryName.sorted_aligned_sam,
      bwa_stderr_log          = BwaMem.bwa_stderr_log,
      bwa_version             = BwaMem.bwa_version,
      bwa_commandline         = bwa_commandline,
      output_bam_basename     = output_bam_basename,
      reference_fasta         = reference_fasta,
      compression_level       = compression_level,
      hard_clip_reads         = hard_clip_reads,
      unmap_contaminant_reads = unmap_contaminant_reads,
      preemptible_tries       = preemptible_tries
  }

  output {
    File output_bam         = MergeAlignment.output_bam
    File fastq              = SamToFastq.fastq
    File aligned_sam        = BwaMem.aligned_sam
    File sorted_aligned_sam = SortSamByQueryName.sorted_aligned_sam
    File bwa_stderr_log     = BwaMem.bwa_stderr_log
  }
}

# Step 1: uBAM -> interleaved FASTQ.
task SamToFastq {
  input {
    File input_bam
    String output_bam_basename
    Int preemptible_tries
  }

  Float unmapped_bam_size = size(input_bam, "GiB")
  Float disk_multiplier = 4.0
  Int disk_size = ceil(unmapped_bam_size + (disk_multiplier * unmapped_bam_size) + 20)

  command <<<
    set -o pipefail
    set -e

    java -Xms1000m -Xmx1000m -jar /usr/gitc/picard.jar \
      SamToFastq \
      INPUT=~{input_bam} \
      FASTQ=~{output_bam_basename}.interleaved.fastq \
      INTERLEAVE=true \
      NON_PF=true
  >>>

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/samtools-picard-bwa:1.0.2-0.7.15-2.26.10-1643840748"
    preemptible: preemptible_tries
    memory: "3 GiB"
    cpu: "2"
    disks: "local-disk " + disk_size + " HDD"
  }

  output {
    File fastq = "~{output_bam_basename}.interleaved.fastq"
  }
}

# Step 2: interleaved FASTQ -> bwa mem
task BwaMem {
  input {
    File input_fastq
    String bwa_commandline
    String output_bam_basename
    ReferenceFasta reference_fasta
    Int preemptible_tries
    Boolean allow_empty_ref_alt = false
  }

  Float fastq_size = size(input_fastq, "GiB")
  Float ref_size = size(reference_fasta.ref_fasta, "GiB") + size(reference_fasta.ref_fasta_index, "GiB") + size(reference_fasta.ref_dict, "GiB")
  Float bwa_ref_size = ref_size + size(reference_fasta.ref_alt, "GiB") + size(reference_fasta.ref_amb, "GiB") + size(reference_fasta.ref_ann, "GiB") + size(reference_fasta.ref_bwt, "GiB") + size(reference_fasta.ref_pac, "GiB") + size(reference_fasta.ref_sa, "GiB")

  Float disk_multiplier = 3.0
  Int disk_size = ceil(fastq_size + bwa_ref_size + (disk_multiplier * fastq_size) + 20)

  command <<<
    set -o pipefail
    set -e

    BWA_VERSION=$(/usr/gitc/bwa 2>&1 | \
      grep -e '^Version' | \
      sed 's/Version: //')

    if [ -z "${BWA_VERSION}" ]; then
      exit 1
    fi

    echo "${BWA_VERSION}" > bwa_version.txt

    if [ -s ~{reference_fasta.ref_alt} ] || ~{allow_empty_ref_alt}; then
      /usr/gitc/~{bwa_commandline} ~{input_fastq} \
        > ~{output_bam_basename}.aligned.unsorted.bwa.sam \
        2> >(tee ~{output_bam_basename}.bwa.stderr.log >&2)

      if ~{!allow_empty_ref_alt}; then
        grep -m1 "read .* ALT contigs" ~{output_bam_basename}.bwa.stderr.log | \
          grep -v "read 0 ALT contigs"
      fi
    else
      echo "ref_alt input is empty or not provided." >&2
      exit 1
    fi
  >>>

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/samtools-picard-bwa:1.0.2-0.7.15-2.26.10-1643840748"
    preemptible: preemptible_tries
    memory: "14 GiB"
    cpu: "16"
    disks: "local-disk " + disk_size + " HDD"
  }

  output {
    File aligned_sam    = "~{output_bam_basename}.aligned.unsorted.bwa.sam"
    File bwa_stderr_log = "~{output_bam_basename}.bwa.stderr.log"
    String bwa_version  = read_string("bwa_version.txt")
  }
}

# Step 3: Explicit queryname sort on aligned SAM using samtools
task SortSamByQueryName {
  input {
    File aligned_sam
    String output_bam_basename
    Int preemptible_tries
  }

  Float aligned_sam_size = size(aligned_sam, "GiB")
  Float disk_multiplier = 2.5
  Int disk_size = ceil(aligned_sam_size + (disk_multiplier * aligned_sam_size) + 20)

  command <<<
    set -o pipefail
    set -e

    samtools sort -n -O SAM -@ 4 \
      -o ~{output_bam_basename}.aligned.query_sorted.sam \
      ~{aligned_sam}
  >>>

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/samtools-picard-bwa:1.0.2-0.7.15-2.26.10-1643840748"
    preemptible: preemptible_tries
    memory: "8 GiB"
    cpu: "4"
    disks: "local-disk " + disk_size + " HDD"
  }

  output {
    File sorted_aligned_sam = "~{output_bam_basename}.aligned.query_sorted.sam"
  }
}

# Step 4: Merge alignment back into uBAM metadata
task MergeAlignment {
  input {
    File input_bam
    File aligned_sam
    File bwa_stderr_log
    String bwa_version
    String bwa_commandline
    String output_bam_basename
    ReferenceFasta reference_fasta
    Int compression_level
    Int preemptible_tries
    Boolean hard_clip_reads = false
    Boolean unmap_contaminant_reads = true
  }

  Float unmapped_bam_size = size(input_bam, "GiB")
  Float aligned_sam_size = size(aligned_sam, "GiB")
  Float ref_size = size(reference_fasta.ref_fasta, "GiB") + size(reference_fasta.ref_fasta_index, "GiB") + size(reference_fasta.ref_dict, "GiB")
  Float disk_multiplier = 2.5
  Int disk_size = ceil(unmapped_bam_size + aligned_sam_size + ref_size + (disk_multiplier * unmapped_bam_size) + 20)

  command <<<
    set -o pipefail
    set -e

    java -Dsamjdk.compression_level=~{compression_level} -Xms1000m -Xmx1000m -jar /usr/gitc/picard.jar \
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
    docker: "us.gcr.io/broad-gotc-prod/samtools-picard-bwa:1.0.2-0.7.15-2.26.10-1643840748"
    preemptible: preemptible_tries
    memory: "5 GiB"
    cpu: "2"
    disks: "local-disk " + disk_size + " HDD"
  }

  output {
    File output_bam = "~{output_bam_basename}.bam"
  }
}
