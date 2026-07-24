version 1.0

## Copyright (structure adapted from Broad Institute WARP, tasks/wdl/Alignment.wdl,
## tag WholeGenomeReprocessing_v3.3.7:
## https://github.com/broadinstitute/warp/blob/WholeGenomeReprocessing_v3.3.7/tasks/wdl/Alignment.wdl
##
## This WDL aligns a single unmapped BAM (uBAM) with bwa mem and merges the alignment
## back into the uBAM's metadata (read groups, tags, etc.) using Picard MergeBamAlignment,
## producing an aligned, unsorted BAM.
##
## DIFFERENCE FROM THE ORIGINAL WARP TASK:
## The original task pipes bwa mem's SAM output directly into MergeBamAlignment via
## ALIGNED_BAM=/dev/stdin in a single shell pipeline. That streaming pattern can trigger a
## known Picard/htsjdk bug in SamAlignmentMerger.getDictionaryForMergedBam(), which throws:
##   "Do not use this function to merge dictionaries with different sequences in them...
##    Found [] and [chr1, chr2, ...]"
## even after MergeBamAlignment has successfully read and merged every alignment record.
## This has been reported independently by multiple users running this exact WARP-style
## pipeline (Biostars, GATK support forum, Drop-seq GitHub issues), including on modern
## GATK4 releases, so it appears to be a persistent upstream limitation rather than
## something specific to any one dataset.
##
## To avoid it, this WDL writes bwa mem's output to an actual file first, then runs
## MergeBamAlignment as a second, separate command reading that completed file (this is
## the fix that resolved the issue for other users who hit the same bug). This costs a
## bit of extra local disk space and I/O for the intermediate aligned SAM/BAM, which is
## reflected in the disk_size calculation below.
##
## These steps are split into three separate tasks (SamToFastq, BwaMem, MergeAlignment)
## rather than one combined task. Each intermediate file - the interleaved FASTQ and the
## aligned SAM - is a first-class WDL output you can pull up on its own: check read counts
## in the FASTQ, run samtools quickcheck/flagstat/view -H on the aligned SAM, etc. - before
## the next step ever touches it. It also means Cromwell retries and caches each step
## independently, so a failure downstream doesn't force earlier steps to rerun.
##
## LICENSING:
## Released under the WDL source code license (BSD-3), matching the license terms of the
## original WARP task this is adapted from. The programs it calls (bwa, Picard, samtools)
## may be subject to different licenses - check the docker image for details.

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

    call MergeAlignment {
        input:
            input_bam               = input_bam,
            aligned_sam             = BwaMem.aligned_sam,
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
        File output_bam     = MergeAlignment.output_bam
        File fastq          = SamToFastq.fastq
        File aligned_sam    = BwaMem.aligned_sam
        File bwa_stderr_log = BwaMem.bwa_stderr_log
    }
}

# Step 1: uBAM -> interleaved FASTQ. Its output (fastq) is a normal task output, so you
# can QC it directly - read counts, spot-check read names/pairing - before bwa mem
# ever touches it.
task SamToFastq {
    input {
        File input_bam
        String output_bam_basename
        Int preemptible_tries
    }

    Float unmapped_bam_size = size(input_bam, "GiB")

    # Interleaved FASTQ (names + sequence + quals, uncompressed) runs larger than the
    # source BAM; leave generous headroom plus spillage room.
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

# Step 2: interleaved FASTQ -> bwa mem, writing the aligned output to a real file rather
# than piping it into MergeBamAlignment. Its output (aligned_sam) is a normal task
# output, so you can QC it directly - samtools quickcheck/flagstat, check @SQ lines are
# complete, compare read counts against the FASTQ - before MergeAlignment ever touches it.
task BwaMem {
    input {
        File input_fastq
        String bwa_commandline
        String output_bam_basename

        # reference_fasta.ref_alt is the .alt file from bwa-kit
        # (https://github.com/lh3/bwa/tree/master/bwakit),
        # listing the reference contigs that are "alternative".
        ReferenceFasta reference_fasta

        Int preemptible_tries
        Boolean allow_empty_ref_alt = false
    }

    Float fastq_size = size(input_fastq, "GiB")
    Float ref_size = size(reference_fasta.ref_fasta, "GiB") + size(reference_fasta.ref_fasta_index, "GiB") + size(reference_fasta.ref_dict, "GiB")
    Float bwa_ref_size = ref_size + size(reference_fasta.ref_alt, "GiB") + size(reference_fasta.ref_amb, "GiB") + size(reference_fasta.ref_ann, "GiB") + size(reference_fasta.ref_bwt, "GiB") + size(reference_fasta.ref_pac, "GiB") + size(reference_fasta.ref_sa, "GiB")

    # The aligned SAM output is uncompressed text, roughly on par with the FASTQ input
    # plus alignment tags, so give it generous headroom.
    Float disk_multiplier = 3.0
    Int disk_size = ceil(fastq_size + bwa_ref_size + (disk_multiplier * fastq_size) + 20)

    command <<<
        set -o pipefail
        set -e

        # This is done before "set -o pipefail" applies to later commands because bwa
        # prints its version banner with rc=1, and we don't want that rc to trip pipefail.
        BWA_VERSION=$(/usr/gitc/bwa 2>&1 | \
            grep -e '^Version' | \
            sed 's/Version: //')

        if [ -z "${BWA_VERSION}" ]; then
            exit 1
        fi
        echo "${BWA_VERSION}" > bwa_version.txt

        # if reference_fasta.ref_alt has data in it or allow_empty_ref_alt is set
        if [ -s ~{reference_fasta.ref_alt} ] || ~{allow_empty_ref_alt}; then

            /usr/gitc/~{bwa_commandline} ~{input_fastq} \
                > ~{output_bam_basename}.aligned.unsorted.bwa.sam \
                2> >(tee ~{output_bam_basename}.bwa.stderr.log >&2)

            if ~{!allow_empty_ref_alt}; then
                grep -m1 "read .* ALT contigs" ~{output_bam_basename}.bwa.stderr.log | \
                    grep -v "read 0 ALT contigs"
            fi

        # else reference_fasta.ref_alt is empty or could not be found
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

# Step 3: merge the completed alignment file back into the unmapped BAM's metadata
# (read groups, tags, etc.) via Picard MergeBamAlignment. Reading ALIGNED_BAM from a
# finished file here - rather than /dev/stdin in the same pipeline as bwa mem - avoids a
# known Picard/htsjdk bug in SamAlignmentMerger.getDictionaryForMergedBam() that throws
# "Do not use this function to merge dictionaries..." even after successfully reading
# every alignment record. That bug has been reported independently by multiple users
# running this same style of pipeline (Biostars, GATK support forum, Drop-seq GitHub
# issues), including on modern GATK4 releases.
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
