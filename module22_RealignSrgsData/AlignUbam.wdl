version 1.0

## Copyright (structure adapted from Broad Institute WARP, tasks/wdl/Alignment.wdl,
## tag WholeGenomeReprocessing_v3.3.7:
## https://github.com/broadinstitute/warp/blob/WholeGenomeReprocessing_v3.3.7/tasks/wdl/Alignment.wdl

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
        File input_bam
        String output_bam_basename
        ReferenceFasta reference_fasta

        File? fastq1
        File? fastq2

        String bwa_commandline = if (defined(fastq1) && defined(fastq2))
            then "bwa mem -K 100000000 -v 3 -t 16 -Y $bash_ref_fasta"
            else "bwa mem -K 100000000 -p -v 3 -t 16 -Y $bash_ref_fasta"

        Int compression_level = 2
        Boolean hard_clip_reads = false
        Boolean unmap_contaminant_reads = true
        Boolean allow_empty_ref_alt = false
        Int preemptible_tries = 3

        RuntimeAttr? runtime_attr_sam_to_fastq_override
        RuntimeAttr? runtime_attr_bwa_mem_override
        RuntimeAttr? runtime_attr_sort_sam_override
        RuntimeAttr? runtime_attr_merge_alignment_override
    }

    Boolean has_paired_fastqs = defined(fastq1) && defined(fastq2)

    if (!has_paired_fastqs) {
        call SamToFastq {
            input:
                input_bam             = input_bam,
                output_bam_basename   = output_bam_basename,
                preemptible_tries     = preemptible_tries,
                runtime_attr_override = runtime_attr_sam_to_fastq_override
        }
    }

    Array[File] bwa_input_fastqs = if has_paired_fastqs
        then select_all([fastq1, fastq2])
        else select_all([SamToFastq.fastq])

    call BwaMem {
        input:
            input_fastqs          = bwa_input_fastqs,
            output_bam_basename   = output_bam_basename,
            bwa_commandline       = bwa_commandline,
            reference_fasta       = reference_fasta,
            preemptible_tries     = preemptible_tries,
            allow_empty_ref_alt   = allow_empty_ref_alt,
            runtime_attr_override = runtime_attr_bwa_mem_override
    }

    # Parallel Picard Sort Step 1: Queryname sort the aligned SAM from BwaMem
    call SortSamByQueryName as SortAlignedSam {
        input:
            input_sam_or_bam      = BwaMem.aligned_sam,
            output_basename       = output_bam_basename + ".aligned",
            preemptible_tries     = preemptible_tries,
            runtime_attr_override = runtime_attr_sort_sam_override
    }

    # Parallel Picard Sort Step 2: Queryname sort the input uBAM using the exact same Picard comparator
    call SortSamByQueryName as SortUnmappedBam {
        input:
            input_sam_or_bam      = input_bam,
            output_basename       = output_bam_basename + ".unmapped",
            preemptible_tries     = preemptible_tries,
            runtime_attr_override = runtime_attr_sort_sam_override
    }

    # Step 3: Merge the matched, Picard-queryname-sorted alignment and uBAM files
    call MergeAlignment {
        input:
            input_bam               = SortUnmappedBam.sorted_bam,
            aligned_sam             = SortAlignedSam.sorted_bam,
            bwa_stderr_log          = BwaMem.bwa_stderr_log,
            bwa_version             = BwaMem.bwa_version,
            bwa_commandline         = bwa_commandline,
            output_bam_basename     = output_bam_basename,
            reference_fasta         = reference_fasta,
            compression_level       = compression_level,
            hard_clip_reads         = hard_clip_reads,
            unmap_contaminant_reads = unmap_contaminant_reads,
            preemptible_tries       = preemptible_tries,
            runtime_attr_override   = runtime_attr_merge_alignment_override
    }

    output {
        File output_bam            = MergeAlignment.output_bam
        File? fastq                = SamToFastq.fastq
        File aligned_sam           = BwaMem.aligned_sam
        File sorted_aligned_bam    = SortAlignedSam.sorted_bam
        File sorted_unmapped_bam   = SortUnmappedBam.sorted_bam
        File bwa_stderr_log        = BwaMem.bwa_stderr_log
    }
}

# Task 1: uBAM -> interleaved FASTQ
task SamToFastq {
    input {
        File input_bam
        String output_bam_basename
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

# Task 2: interleaved FASTQ -> BWA-MEM alignment
task BwaMem {
    input {
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

    command <<<
        BWA_VERSION=$(/usr/gitc/bwa 2>&1 | \
            grep -e '^Version' | \
            sed 's/Version: //')

        set -o pipefail
        set -e

        if [ -z "${BWA_VERSION}" ]; then
            exit 1
        fi
        echo "${BWA_VERSION}" > bwa_version.txt

        bash_ref_fasta=~{reference_fasta.ref_fasta}

        if [ -s ~{reference_fasta.ref_alt} ] || ~{allow_empty_ref_alt}; then
            /usr/gitc/~{bwa_commandline} ~{sep=' ' input_fastqs} \
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

# Task 3: Picard SortSam task (called in parallel for both aligned SAM and unmapped BAM)
task SortSamByQueryName {
    input {
        File input_sam_or_bam
        String output_basename
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

# Task 4: Picard MergeBamAlignment on strictly matched queryname-sorted streams
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
