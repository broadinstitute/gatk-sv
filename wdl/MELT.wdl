## This WDL pipeline runs MELT for mobile element insertion/deletion discovery
## Note: runs patched v2.0.5, obtained directly from the MELT authors
##
## Requirements/expectations :
## - Human whole-genome sequencing data in mapped BAM format
## - BAM index
## - Reference assembly FASTA
## - Approximate haploid nucleotide coverage of library (optional)
## - Library read length and average insert size (optional)

version 1.0

import "Structs.wdl"
import "CollectCoverage.wdl" as cov

workflow MELT {
  input {
    File bam_or_cram_file
    File? bam_or_cram_index
    File counts_file
    String sample_id
    File reference_fasta
    String? reference_index
    File reference_dict
    String? reference_version
    File melt_standard_vcf_header

    File? wgs_metrics_intervals
    File? multiple_metrics_intervals

    Float? insert_size
    Int? read_length
    Float? coverage

    Float? pct_chimeras
    Float? total_reads
    Int? pf_reads_improper_pairs

    String samtools_cloud_docker
    String gatk_docker
    String genomes_in_the_cloud_docker
    String melt_docker
    RuntimeAttr? runtime_attr_coverage
    RuntimeAttr? runtime_attr_metrics
    RuntimeAttr? runtime_attr_cram_to_bam
    RuntimeAttr? runtime_attr_melt
  }

  parameter_meta {
      bam_or_cram_file: ".bam or .cram file to search for SVs. bams are preferable, crams will be converted to bams."
      bam_or_cram_index: "[optional] associated index file. If omitted, the WDL will look for an index file by appending .bai/.crai to the .bam/.cram file"
      sample_id: "sample name. Outputs will be sample_name + '.vcf.gz' and sample_name + '.vcf.gz.tbi'"
      reference_fasta: ".fasta file with reference used to align bam or cram file"
      reference_index: "[optional] reference index file. If omitted, the WDL will look for an index by appending .fai to the .fasta file"
      reference_version: "[optional] Number of genome reference used (19 or 38). If omitted, defaults to 38."
      melt_standard_vcf_header: "vcf header used to correct MELT output and ensure proper vcf header format."
      coverage: "[optional] Value of MEAN_COVERAGE obtained from CollectWgsMetrics, needed to run MELT. If omitted, a task will be created to calculate it."
      read_length: "[optional] Value of MEAN_READ_LENGTH obtained from CollectAlignmentSummaryMetrics, needed to run MELT. If omitted, a task will be created to calculate it."
      insert_size: "[optional] Value of MEAN_INSERT_SIZE obtained from CollectInsertSizeMetrics, needed to run MELT. If omitted, a task will be created to calculate it."

      pf_reads_improper_pairs: "[optional] Value of PF_READS_IMPROPER_PAIRS obtained from CollectAlignmentSummaryMetrics, used for optimal estimate of VM memory and disk needs. If omitted, estimates based on coverage and insert size are used."
      pct_exc_total: "[optional] Value of PCT_EXC_TOTAL obtained from CollectWgsMetrics, used for optimal estimate of VM memory and disk needs. If omitted, estimates based on coverage and insert size are used."
    }

  Boolean is_bam = basename(bam_or_cram_file, ".bam") + ".bam" == basename(bam_or_cram_file)
  File bam_or_cram_index_file = select_first([bam_or_cram_index, if is_bam then bam_or_cram_file + ".bai" else bam_or_cram_file + ".crai"])
  File reference_index_file = select_first([reference_index, reference_fasta + ".fai"])

  Boolean have_multiple_metrics = defined(insert_size) && defined(read_length)
     && defined(pct_chimeras) && defined(total_reads) && defined(pf_reads_improper_pairs)
  if (!have_multiple_metrics) {
    call GetMultipleMetrics {
      input:
        bam_or_cram_file = bam_or_cram_file,
        bam_or_cram_index = bam_or_cram_index_file,
        reference_fasta = reference_fasta,
        reference_index = reference_index_file,
        intervals = multiple_metrics_intervals,
        gatk_docker = gatk_docker,
        runtime_attr_override = runtime_attr_metrics
    }

    # form map with properties:
      # TOTAL_READS, PF_READS, PCT_PF_READS, PF_NOISE_READS, PF_READS_ALIGNED, PCT_PF_READS_ALIGNED,
        # PF_ALIGNED_BASES, PF_HQ_ALIGNED_READS, PF_HQ_ALIGNED_BASES, PF_HQ_ALIGNED_Q20_BASES,
        # PF_HQ_MEDIAN_MISMATCHES, PF_MISMATCH_RATE, PF_HQ_ERROR_RATE, PF_INDEL_RATE, MEAN_READ_LENGTH,
        # READS_ALIGNED_IN_PAIRS, PCT_READS_ALIGNED_IN_PAIRS, PF_READS_IMPROPER_PAIRS, PCT_PF_READS_IMPROPER_PAIRS,
        # BAD_CYCLES, STRAND_BALANCE, PCT_CHIMERAS, PCT_ADAPTER
      # MEDIAN_INSERT_SIZE, MEDIAN_ABSOLUTE_DEVIATION, MIN_INSERT_SIZE, MAX_INSERT_SIZE, MEAN_INSERT_SIZE,
        # STANDARD_DEVIATION, READ_PAIRS, PAIR_ORIENTATION, WIDTH_OF_10_PERCENT, WIDTH_OF_20_PERCENT,
        # WIDTH_OF_30_PERCENT, WIDTH_OF_40_PERCENT, WIDTH_OF_50_PERCENT, WIDTH_OF_60_PERCENT, WIDTH_OF_70_PERCENT,
        # WIDTH_OF_80_PERCENT, WIDTH_OF_90_PERCENT, WIDTH_OF_99_PERCENT
      # WINDOW_SIZE, TOTAL_CLUSTERS, ALIGNED_READS, AT_DROPOUT, GC_DROPOUT, GC_NC_0_19, GC_NC_20_39,
          # GC_NC_40_59, GC_NC_60_79, GC_NC_80_100
    Map[String, String] multiple_metrics_map = read_map(GetMultipleMetrics.metrics_table_file)

    Float calculated_insert_size = multiple_metrics_map["MEAN_INSERT_SIZE"]
    Int calculated_read_length = round(multiple_metrics_map["MEAN_READ_LENGTH"])
    Float calculated_pct_chimeras = multiple_metrics_map["PCT_CHIMERAS"]
    Float calculated_total_reads = multiple_metrics_map["TOTAL_READS"]
    Int calculated_pf_reads_improper_pairs = multiple_metrics_map["PF_READS_IMPROPER_PAIRS"]
  }

  Boolean have_wgs_metrics = defined(coverage)
  if (!have_wgs_metrics) {
    call GetWgsMetrics {
      input:
        bam_or_cram_file = bam_or_cram_file,
        bam_or_cram_index = bam_or_cram_index_file,
        read_length = select_first([read_length, calculated_read_length]),
        reference_fasta = reference_fasta,
        reference_index = reference_index_file,
        intervals = wgs_metrics_intervals,
        genomes_in_the_cloud_docker = genomes_in_the_cloud_docker,
        runtime_attr_override = runtime_attr_coverage
    }
    # form map with available properties:
    # GENOME_TERRITORY, MEAN_COVERAGE, SD_COVERAGE, MEDIAN_COVERAGE, MAD_COVERAGE, PCT_EXC_MAPQ, PCT_EXC_DUPE,
    # PCT_EXC_UNPAIRED, PCT_EXC_BASEQ, PCT_EXC_OVERLAP, PCT_EXC_CAPPED, PCT_EXC_TOTAL, PCT_1X, PCT_5X, PCT_10X,
    # PCT_15X, PCT_20X, PCT_25X, PCT_30X, PCT_40X, PCT_50X, PCT_60X, PCT_70X, PCT_80X, PCT_90X, PCT_100X,
    # HET_SNP_SENSITIVITY, HET_SNP_Q
    Map[String, String] wgs_metrics_map = read_map(GetWgsMetrics.metrics_file)

    Float calculated_coverage = wgs_metrics_map["MEAN_COVERAGE"]
  }

  Float melt_coverage = select_first([coverage, calculated_coverage])
  Int melt_read_length = select_first([read_length, calculated_read_length])
  Float melt_insert_size = select_first([insert_size, calculated_insert_size])
  Float melt_pct_chimeras = select_first([pct_chimeras, calculated_pct_chimeras])
  Float melt_total_reads = select_first([total_reads, calculated_total_reads])
  Int melt_pf_reads_improper_pairs = select_first([pf_reads_improper_pairs, calculated_pf_reads_improper_pairs])

  call RunMELT {
    input:
      bam_or_cram_file = bam_or_cram_file,
      bam_or_cram_index = bam_or_cram_index_file,
      sample_id = sample_id,
      reference_fasta = reference_fasta,
      reference_index = reference_index_file,
      reference_dict = reference_dict,
      reference_version = reference_version,
      counts_file = counts_file,
      melt_standard_vcf_header = melt_standard_vcf_header,
      coverage = melt_coverage,
      read_length = melt_read_length,
      insert_size = melt_insert_size,
      pct_chimeras = melt_pct_chimeras,
      total_reads = melt_total_reads,
      pf_reads_improper_pairs = melt_pf_reads_improper_pairs,
      melt_docker = melt_docker,
      runtime_attr_override = runtime_attr_melt
  }

  output {
    File vcf = RunMELT.vcf
    File index = RunMELT.index
    Float coverage_out = melt_coverage
    Int read_length_out = melt_read_length
    Float insert_size_out = melt_insert_size
  }
}

task GetMultipleMetrics {
  input {
    File bam_or_cram_file
    File bam_or_cram_index
    File reference_fasta
    File reference_index
    File? intervals
    String gatk_docker
    RuntimeAttr? runtime_attr_override
  }

  String metrics_base = "multiple_metrics"
  String alignment_summary_filename = metrics_base + ".alignment_summary_metrics"
  String insert_size_filename = metrics_base + ".insert_size_metrics"
  String gc_bias_filename = metrics_base + ".gc_bias.summary_metrics"
  String metrics_table_filename=metrics_base + "_table.tsv"

  Float bam_or_cram_size = size(bam_or_cram_file, "GiB")
  Float disk_overhead = 20.0
  Float ref_size = size(reference_fasta, "GiB")
  Int vm_disk_size = ceil(bam_or_cram_size + ref_size + disk_overhead)

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 7.0,
    disk_gb: vm_disk_size,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  Int cpu_cores = select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
  Float mem_use_gb = select_first([runtime_attr.mem_gb, cpu_cores * 6.0])
  Float java_mem_pad_gb = cpu_cores * 1.0
  Float mem_size_gb = mem_use_gb + java_mem_pad_gb

  runtime {
    cpu: cpu_cores
    memory: mem_use_gb + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: gatk_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
  Int java_mem_mb = round(1024 * (select_first([runtime_attr.mem_gb, default_attr.mem_gb]) - java_mem_pad_gb))

  command <<<

    set -Eeuo pipefail

    gatk --java-options -Xmx~{java_mem_mb}m CollectMultipleMetrics \
      -I "~{bam_or_cram_file}" \
      -O "~{metrics_base}" \
      -R "~{reference_fasta}" \
      ~{"--INTERVALS " + intervals} \
      --ASSUME_SORTED true \
      --PROGRAM null \
      --PROGRAM CollectAlignmentSummaryMetrics \
      --PROGRAM CollectInsertSizeMetrics \
      --PROGRAM CollectSequencingArtifactMetrics \
      --PROGRAM CollectGcBiasMetrics \
      --PROGRAM QualityScoreDistribution \
      --METRIC_ACCUMULATION_LEVEL null \
      --METRIC_ACCUMULATION_LEVEL SAMPLE

    function transpose_table() {
          cat \
          | awk ' {
              for (col = 1; col <= NF; ++col) {
                table[NR, col] = $col
              }
              if(NF > num_cols) {
                num_cols = NF
              }
            } END {
              for (row = 1; row <= num_cols; ++row) {
                printf "%s", table[1, row]
                for (col = 2; col <= NR; ++col) {
                  printf "\t%s", table[col, row]
                }
                printf "\n"
              }
            }'
        }

    # get alignment summary metrics for all PAIR from sample,
    # and transpose to get one
    #    property_name -tab- value
    # per line. This is readable by cromwell read_map()
    # Then remove the CATEGORY, LIBRARY, SAMPLE, and READ_GROUP properties
    # to get table with properties:
    #   TOTAL_READS, PF_READS, PCT_PF_READS, PF_NOISE_READS, PF_READS_ALIGNED, PCT_PF_READS_ALIGNED,
    #   PF_ALIGNED_BASES, PF_HQ_ALIGNED_READS, PF_HQ_ALIGNED_BASES, PF_HQ_ALIGNED_Q20_BASES,
    #   PF_HQ_MEDIAN_MISMATCHES, PF_MISMATCH_RATE, PF_HQ_ERROR_RATE, PF_INDEL_RATE, MEAN_READ_LENGTH,
    #   READS_ALIGNED_IN_PAIRS, PCT_READS_ALIGNED_IN_PAIRS, PF_READS_IMPROPER_PAIRS, PCT_PF_READS_IMPROPER_PAIRS,
    #   BAD_CYCLES, STRAND_BALANCE, PCT_CHIMERAS, PCT_ADAPTER
    grep -A4 "^## METRICS CLASS" "~{alignment_summary_filename}" \
      | grep -E "^CATEGORY\s|^PAIR\s" \
      | transpose_table \
      | grep -Ev "^CATEGORY\s|^LIBRARY\s|^SAMPLE\s|^READ_GROUP\s" \
      > "~{metrics_table_filename}"

    # get insert size metrics from sample,
    # and transpose to format compatible with cromwell read_map
    # to get table with properties:
    #   MEDIAN_INSERT_SIZE, MEDIAN_ABSOLUTE_DEVIATION, MIN_INSERT_SIZE, MAX_INSERT_SIZE, MEAN_INSERT_SIZE,
    #   STANDARD_DEVIATION, READ_PAIRS, PAIR_ORIENTATION, WIDTH_OF_10_PERCENT, WIDTH_OF_20_PERCENT,
    #   WIDTH_OF_30_PERCENT, WIDTH_OF_40_PERCENT, WIDTH_OF_50_PERCENT, WIDTH_OF_60_PERCENT, WIDTH_OF_70_PERCENT,
    #   WIDTH_OF_80_PERCENT, WIDTH_OF_90_PERCENT, WIDTH_OF_99_PERCENT
    grep -A1 MEDIAN_INSERT_SIZE "~{insert_size_filename}"  \
      | transpose_table \
      | grep -Ev "^LIBRARY\s|^SAMPLE\s|^READ_GROUP\s" \
      >> "~{metrics_table_filename}"

    # get gc bias summary metrics from sample,
    # and transpose to format compatible with cromwell read_map
    # to get table with properties:
    #   WINDOW_SIZE, TOTAL_CLUSTERS, ALIGNED_READS, AT_DROPOUT, GC_DROPOUT, GC_NC_0_19, GC_NC_20_39,
    #   GC_NC_40_59, GC_NC_60_79, GC_NC_80_100
    grep -A1 "^ACCUMULATION_LEVEL" "~{gc_bias_filename}" \
      | transpose_table \
      | grep -Ev "^ACCUMULATION_LEVEL\s|^READS_USED\s|^LIBRARY\s|^SAMPLE\s|^READ_GROUP\s" \
      >> "~{metrics_table_filename}"
  >>>

  output {
    Array[File] metrics_files = glob("~{metrics_base}.*")
    File metrics_table_file = metrics_table_filename
  }
}

task GetWgsMetrics {
  input {
    File bam_or_cram_file
    File bam_or_cram_index
    Int read_length
    File reference_fasta
    File reference_index
    File? intervals
    Boolean use_fast_algorithm = true
    String genomes_in_the_cloud_docker
    RuntimeAttr? runtime_attr_override
  }

  String metrics_file_name = "wgs_metrics.txt"

  Int num_cpu = if defined(runtime_attr_override) then select_first([select_first([runtime_attr_override]).cpu_cores, 1]) else 1
  Float mem_size_gb = num_cpu * 8.0
  Int java_mem_mb = round(mem_size_gb * 0.8 * 1000)
  Float bam_or_cram_size = size(bam_or_cram_file, "GiB")
  Float disk_overhead = 20.0
  Float ref_size = size(reference_fasta, "GiB")
  Int vm_disk_size = ceil(bam_or_cram_size + ref_size + disk_overhead)

  RuntimeAttr default_attr = object {
    cpu_cores: num_cpu,
    mem_gb: mem_size_gb,
    disk_gb: vm_disk_size,
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: genomes_in_the_cloud_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

  command <<<

    set -Eeuo pipefail

    # calculate coverage
    java -Xms~{java_mem_mb}m -jar /usr/gitc/picard.jar \
      CollectWgsMetrics \
      INPUT=~{bam_or_cram_file} \
      VALIDATION_STRINGENCY=SILENT \
      REFERENCE_SEQUENCE=~{reference_fasta} \
      READ_LENGTH=~{read_length} \
      INCLUDE_BQ_HISTOGRAM=true \
      ~{if defined(intervals) then "INTERVALS=~{intervals}" else ""} \
      OUTPUT="raw_~{metrics_file_name}" \
      USE_FAST_ALGORITHM=~{use_fast_algorithm}

    function transpose_table() {
      cat \
      | awk ' {
          for (col = 1; col <= NF; ++col) {
            table[NR, col] = $col
          }
          if(NF > num_cols) {
            num_cols = NF
          }
        } END {
          for (row = 1; row <= num_cols; ++row) {
            printf "%s", table[1, row]
            for (col = 2; col <= NR; ++col) {
              printf "\t%s", table[col, row]
            }
            printf "\n"
          }
        }'
    }

    # get the two output lines corresponding to the calculated properties header and values,
    # and transpose to get one
    #    property_name -tab- value
    # per line. This is readable by cromwell read_map()
    grep -A1 MEAN_COVERAGE "raw_~{metrics_file_name}" \
      | transpose_table \
      > "~{metrics_file_name}"
  >>>

  output {
    File metrics_file = metrics_file_name
  }
}

task RunMELT {
  input {
    File bam_or_cram_file
    File bam_or_cram_index
    String sample_id
    File reference_fasta
    File reference_index
    File reference_dict
    String reference_version = "38"
    File counts_file
    File melt_standard_vcf_header

    # picard metrics needed for MELT operation
    Float coverage
    Int read_length
    Float insert_size

    # needed for VM sizing
    Float pct_chimeras
    Float total_reads
    Int pf_reads_improper_pairs

    # needed for filtering out high coverage intervals
    Int interval_padding = 100
    Int threshold = 1000

    String melt_docker
    RuntimeAttr? runtime_attr_override
  }

  parameter_meta {
    bam_or_cram_file: {
      localization_optional: true
    }
    bam_or_cram_index: {
      localization_optional: true
    }
  }


  # Estimate the chances of successfully running preemptively. Data collected by spec-ops suggests a failure rate of
  # roughly 0.1% per hour, with some wiggles and a rapid increase at around 500 minutes.
  # Running with preemptible = 1 is better than non-preemptively when the failure rate + the preemptible cost fraction
  # is less than 100%, (the cost fraction is 20%, so this means failure rate < 80%). If preemption events were
  # independent, optimal strategy would be to use preemptible = inf below this cutoff, however we suspect they are not
  # independent, so we use preemptible = 1 immediately above the cutoff, and set it to 3 well below the cutoff.
  # The estimated runtime is a binlinear regression based on the sv benchmark set + prenatal samples.
  Float breakeven_runtime_hours = 8.0
  Float risk_multiple_tries_hours = 6.0
  Float runtime_hours_per_pct_chimeras = 617.0
  Float runtime_hours_per_total_reads = "7.178e-9"
  Float runtime_hours_offset = -5.685
  Float runtime_hours_estimated =
    runtime_hours_offset
    + runtime_hours_per_pct_chimeras * pct_chimeras
    + runtime_hours_per_total_reads * total_reads
  Int preemptible_tries =
    if runtime_hours_estimated > breakeven_runtime_hours then 0
    else if runtime_hours_estimated > risk_multiple_tries_hours then 1
    else 3


  # Ensure there's sufficient memory. Estimate using extra metrics
  Float mem_per_pct_chimeras = 50.69
  Float mem_per_improper_pairs = "1.451e-8"
  Float mem_offset = 6.833
  Float mem_size_gb =
    mem_offset + mem_per_pct_chimeras * pct_chimeras + mem_per_improper_pairs * pf_reads_improper_pairs
  Float java_mem_fraction = 0.85
  Int java_mem_mb = round(mem_size_gb * java_mem_fraction * 1000)

  # Ensure there's sufficient disk space. Estimate using extra metrics
  Float disk_per_improper_pairs = "3.901e-7"
  Float disk_per_bam_size = 0.3088
  Float disk_per_cram_size = 3.5 * disk_per_bam_size
  Boolean is_bam = basename(bam_or_cram_file, ".bam") + ".bam" == basename(bam_or_cram_file)
  Float input_size = size(bam_or_cram_file, "GiB")
  Float disk_per_input_size = if is_bam then disk_per_bam_size else disk_per_cram_size
  Float disk__offset = 32.57
  Float disk_overhead =
    disk__offset + disk_per_improper_pairs * pf_reads_improper_pairs + disk_per_input_size * input_size

  Float filtered_bam_size = input_size * (if is_bam then 1.0 else 3.5)
  Float ref_size = size(reference_fasta, "GiB")
  Int vm_disk_size = ceil(filtered_bam_size + ref_size + disk_overhead)


  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: mem_size_gb,
    disk_gb: vm_disk_size,
    boot_disk_gb: 10,
    preemptible_tries: preemptible_tries,
    max_retries: 1,
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File vcf = "${sample_id}.melt.vcf.gz"
    File index = "${sample_id}.melt.vcf.gz.tbi"
  }
  command <<<

    set -Eeuo pipefail

    zgrep -v -E "^@|^CONTIG|^#" ~{counts_file} | \
      awk 'BEGIN{FS=OFS="\t"}{if ($4 > ~{threshold}) {print $1,$2,$3}}' > highCountIntervals.bed

    # create a bam file from the input bam or cram, filtering out high coverage intervals
    # and ambiguous bases
    java -Xmx~{java_mem_mb}m -jar /opt/gatk.jar PrintReads \
          -XL highCountIntervals.bed \
          --interval-exclusion-padding ~{interval_padding} \
          -I ~{bam_or_cram_file} \
          -R ~{reference_fasta} \
          -O /dev/stdout | \
      samtools view -h - | \
      awk 'BEGIN{FS=OFS="\t"}{gsub(/[BDHVRYKMSW]/,"N",$10);print}' | \
      samtools view -b1 - > ~{sample_id}.bam

    samtools index ~{sample_id}.bam

    # these locations should be stable
    MELT_DIR="/MELT"

    # `cromwell_root` and `cromwell-executions` are the **default**
    # root directory for cromwell deployments on GCP and Azure respectively.
    # The last option, i.e., $PWD, is the fall back option if the
    # default root directory of the cromwell instance was changed,
    # however, it is not a reliable option as it fails on GCP.
    # It does not seem Cromwell sets a runtime environment variable
    # exposing the configured value of the root directory,
    # which could have provided a portable solution for this.
    # The following solution works with the Cromwell deployments
    # we are currently using on GCP and Azure.
    if [ -d "/cromwell_root" ]; then
      CROMWELL_ROOT="/cromwell_root"
    elif [ -d "/cromwell-executions" ]; then
      CROMWELL_ROOT="/cromwell-executions"
    else
      CROMWELL_ROOT="$PWD"
    fi

    # these locations may vary based on MELT version number, so find them:
    MELT_ROOT=$(find "$MELT_DIR" -name "MELT.jar" | xargs -n1 dirname)
    MELT_SCRIPT=$(ls "$MELT_DIR/run_MELT"*.sh)

    # call MELT
    "$MELT_SCRIPT" \
      ~{sample_id}.bam \
      "~{reference_fasta}" \
      ~{coverage} \
      ~{read_length} \
      ~{insert_size} \
      "$MELT_ROOT" \
      "$CROMWELL_ROOT" \
      ~{reference_version}

    cat "~{melt_standard_vcf_header}" \
        <(echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t~{sample_id}") \
        <(grep -v "^#" SVA.final_comp.vcf) \
        <(grep -v "^#" LINE1.final_comp.vcf) \
        <(grep -v "^#" ALU.final_comp.vcf) \
      | sed -e "2i##fileDate=$(date +'%Y%m%d')" -e "s/No Difference/No_Difference/" \
      | bcftools sort -Oz - > "~{sample_id}.melt.vcf.gz"
    bcftools index -t "~{sample_id}.melt.vcf.gz"

    df -h
    ls -l
  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: melt_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}
