version 1.0

workflow ConvertPairedFastQsToUnmappedBamWf {

  input {
    String sample_name

    # Provide EITHER fastq_1+fastq_2 OR ubam
    File? fastq_1
    File? fastq_2
    File? ubam

    # Required only when FASTQs are provided
    Array[String]? readgroup_names
    Array[String]? library_names
    Array[String]? platform_units
    Array[String]? run_dates
    Array[String]? platform_names
    Array[String]? sequencing_centers
    File? ref_dict

    String samtools_docker = "staphb/samtools:1.17"
    String gatk_docker     = "broadinstitute/gatk:4.5.0.0"

    Int additional_disk_space_gb = 50
    Int machine_mem_gb           = 7
    Int preemptible_attempts     = 1
  }

  # FASTQ path: convert paired FASTQs to uBAM
  if (defined(fastq_1) && defined(fastq_2)) {
    call PairedFastQsToUnmappedBAM {
      input:
        sample_name              = sample_name,
        fastq_1                  = select_first([fastq_1]),
        fastq_2                  = select_first([fastq_2]),
        readgroup_names          = select_first([readgroup_names]),
        library_names            = select_first([library_names]),
        platform_units           = select_first([platform_units]),
        run_dates                = select_first([run_dates]),
        platform_names           = select_first([platform_names]),
        sequencing_centers       = select_first([sequencing_centers]),
        ref_dict                 = select_first([ref_dict]),
        docker                   = samtools_docker,
        additional_disk_space_gb = additional_disk_space_gb,
        machine_mem_gb           = machine_mem_gb,
        preemptible_attempts     = preemptible_attempts
    }
  }

  # Always sort — either the uBAM from the FASTQ path or the provided uBAM
  File bam_to_sort = select_first([PairedFastQsToUnmappedBAM.output_unmapped_bam, ubam])

  call SortSam {
    input:
      input_bam                = bam_to_sort,
      sample_name              = sample_name,
      docker                   = gatk_docker,
      machine_mem_gb           = machine_mem_gb,
      additional_disk_space_gb = additional_disk_space_gb,
      preemptible_attempts     = preemptible_attempts
  }

  output {
    File output_unmapped_bam = SortSam.sorted_bam
  }
}

task PairedFastQsToUnmappedBAM {
  input {
    String        sample_name
    File          fastq_1
    File          fastq_2
    Array[String] readgroup_names
    Array[String] library_names
    Array[String] platform_units
    Array[String] run_dates
    Array[String] platform_names
    Array[String] sequencing_centers

    File   ref_dict

    Int    additional_disk_space_gb = 50
    Int    machine_mem_gb           = 7
    Int    preemptible_attempts     = 1
    String docker
  }

  Int disk_space_gb = ceil((size(fastq_1, "GB") + size(fastq_2, "GB")) * 4) + additional_disk_space_gb

  command <<<
    set -euo pipefail

    # Read arrays from files to avoid shell quoting issues
    mapfile -t RG_NAMES    < ~{write_lines(readgroup_names)}
    mapfile -t LIB_NAMES   < ~{write_lines(library_names)}
    mapfile -t PLAT_UNITS  < ~{write_lines(platform_units)}
    mapfile -t RUN_DATES   < ~{write_lines(run_dates)}
    mapfile -t PLAT_NAMES  < ~{write_lines(platform_names)}
    mapfile -t SEQ_CENTERS < ~{write_lines(sequencing_centers)}

    n_rg=${#RG_NAMES[@]}

    # Validate all arrays are the same length
    for arr_name in LIB_NAMES PLAT_UNITS RUN_DATES PLAT_NAMES SEQ_CENTERS; do
      declare -n arr="$arr_name"
      n=${#arr[@]}
      if [ "$n" -ne "$n_rg" ]; then
        echo "ERROR: $arr_name length ($n) does not match readgroup_names length ($n_rg)" >&2
        exit 1
      fi
    done

    # Build header with @SQ lines from ref_dict.
    # Picard 2.26.x MergeBamAlignment requires the uBAM sequence dictionary
    # to match the reference — an empty dictionary causes a merge failure.
    printf "@HD\tVN:1.6\tSO:queryname\n" > new_header.sam
    grep "^@SQ" ~{ref_dict} >> new_header.sam

    for (( i=0; i<n_rg; i++ )); do
      printf "@RG\tID:%s\tSM:%s\tLB:%s\tPU:%s\tDT:%s\tPL:%s\tCN:%s\n" \
        "${RG_NAMES[$i]}" \
        "~{sample_name}" \
        "${LIB_NAMES[$i]}" \
        "${PLAT_UNITS[$i]}" \
        "${RUN_DATES[$i]}" \
        "${PLAT_NAMES[$i]}" \
        "${SEQ_CENTERS[$i]}" \
        >> new_header.sam
    done

    echo "Header to be applied:"
    cat new_header.sam
    echo "RG lines: $(grep -c "^@RG" new_header.sam)"
    echo "SQ lines: $(grep -c "^@SQ" new_header.sam)"

    FIRST_RG=$(printf "ID:%s\tSM:%s\tLB:%s\tPU:%s\tDT:%s\tPL:%s\tCN:%s" \
      "${RG_NAMES[0]}" \
      "~{sample_name}" \
      "${LIB_NAMES[0]}" \
      "${PLAT_UNITS[0]}" \
      "${RUN_DATES[0]}" \
      "${PLAT_NAMES[0]}" \
      "${SEQ_CENTERS[0]}")

    samtools import \
      -1 ~{fastq_1} \
      -2 ~{fastq_2} \
      -r "$FIRST_RG" \
      -O BAM \
      -o raw.unmapped.bam

    samtools reheader new_header.sam raw.unmapped.bam > ~{sample_name}.unmapped.bam

    echo "Done. Verifying output BAM header:"
    samtools view -H ~{sample_name}.unmapped.bam | grep "^@RG"
    echo "SQ lines: $(samtools view -H ~{sample_name}.unmapped.bam | grep -c "^@SQ")"
    echo "Flagstat:"
    samtools flagstat ~{sample_name}.unmapped.bam
  >>>

  output {
    File output_unmapped_bam = "~{sample_name}.unmapped.bam"
  }

  runtime {
    docker:      docker
    memory:      machine_mem_gb + " GB"
    disks:       "local-disk " + disk_space_gb + " HDD"
    preemptible: preemptible_attempts
    maxRetries:  preemptible_attempts
  }
}

task SortSam {
  input {
    File   input_bam
    String sample_name
    String docker
    Int    machine_mem_gb           = 7
    Int    additional_disk_space_gb = 50
    Int    preemptible_attempts     = 1
  }

  Int command_mem_gb = machine_mem_gb - 1
  Int disk_space_gb  = ceil(size(input_bam, "GB") * 4) + additional_disk_space_gb

  command <<<
    set -euo pipefail

    mkdir -p tmp

    gatk --java-options "-Xmx~{command_mem_gb}g -Djava.io.tmpdir=./tmp" \
      SortSam \
      -I ~{input_bam} \
      -O ~{sample_name}.query_sorted.unmapped.bam \
      --SORT_ORDER queryname \
      --TMP_DIR ./tmp

    echo "Sort complete. Flagstat:"
    samtools flagstat ~{sample_name}.query_sorted.unmapped.bam
  >>>

  output {
    File sorted_bam = "~{sample_name}.query_sorted.unmapped.bam"
  }

  runtime {
    docker:      docker
    memory:      machine_mem_gb + " GB"
    disks:       "local-disk " + disk_space_gb + " HDD"
    preemptible: preemptible_attempts
    maxRetries:  preemptible_attempts
  }
}
