version 1.0

# Standalone monolithic 00a, which may be faster on HPC/local backends
# Written only for BAM inputs to generate counts, PE/SR, manta, melt, and wham outputs

workflow GatherSampleEvidenceHPC {
  input {
    Array[File] bam_files
    Array[File] bam_indexes
    Array[String] sample_ids
  }

  scatter(i in range(length(bam_files))) {
    call RunCallers {
      input:
        bam_file = bam_files[i],
        bam_index = bam_indexes[i],
        sample_id = sample_ids[i]
    }
  }

  output {
    Array[File] coverage_counts = RunCallers.coverage_counts

    Array[File] manta_vcf = RunCallers.manta_vcf
    Array[File] manta_index = RunCallers.manta_index

    Array[File] melt_vcf = RunCallers.melt_vcf
    Array[File] melt_index = RunCallers.melt_index
    Array[File] wgs_metrics_file = RunCallers.wgs_metrics_file

    File split_out = RunCallers.split_out
    File split_out_index = RunCallers.split_out_index
    File disc_out = RunCallers.disc_out
    File disc_out_index = RunCallers.disc_out_index
    File sd_out = RunCallers.sd_out
    File sd_out_index = RunCallers.sd_out_index

    Array[File] wham_vcf = RunCallers.wham_vcf
    Array[File] wham_index = RunCallers.wham_index
  }
}

task RunCallers {
  input {

    File bam_file
    File bam_index
    String sample_id

    # Resources - do not exceed requested VM size!
    Int num_cpu
    Float mem_gb

    # Common parameters
    File primary_contigs_list
    File reference_fasta
    File reference_index    # Index (.fai), must be in same dir as fasta
    File reference_dict     # Dictionary (.dict), must be in same dir as fasta
    String reference_version   # Either "38" or "19"

    # SD options
    File sd_locs_vcf
    Int site_depth_min_mapq = 6
    Int site_depth_min_baseq = 10

    # Coverage collection inputs
    File preprocessed_intervals

    # Manta inputs
    File manta_region_bed
    File manta_region_bed_index

    # Melt inputs
    File melt_standard_vcf_header

    # Wham inputs
    File wham_include_list_bed_file

    # Docker
    String sv_pipeline_callers_docker

  }

  Int command_mem_mb = ceil(0.8 * mem_gb * 1000)

  String metrics_file_name = "wgs_metrics.txt"
  Boolean fast_algorithm = true

  String metrics_base = "multiple_metrics"
  String alignment_summary_filename = metrics_base + ".alignment_summary_metrics"
  String insert_size_filename = metrics_base + ".insert_size_metrics"
  String gc_bias_filename = metrics_base + ".gc_bias.summary_metrics"
  String metrics_table_filename=metrics_base + "_table.tsv"

  Int interval_padding = 100
  Int threshold = 1000

  # select number of jobs (threads) to run
  Float jobs_per_cpu_use = 1.3
  Int manta_num_jobs = round(num_cpu * jobs_per_cpu_use)

  Array[String] good_intervals = read_lines(wham_include_list_bed_file)

  Array[String] chr_list = read_lines(primary_contigs_list)

  command <<<

    set -Eeuxo pipefail

    ###################################################
    ###### CollectCounts
    ###################################################

    java -Xmx~{command_mem_mb}m -jar /opt/gatk.jar CollectReadCounts \
      -L ~{preprocessed_intervals} \
      --input ~{bam_file} \
      --read-index ~{bam_index} \
      --reference ~{reference_fasta} \
      --format TSV \
      --interval-merging-rule OVERLAPPING_ONLY \
      --output ~{sample_id}.counts.tsv \
      --disable-read-filter MappingQualityReadFilter

    sed -ri "s/@RG\tID:GATKCopyNumber\tSM:.+/@RG\tID:GATKCopyNumber\tSM:~{sample_id}/g" ~{sample_id}.counts.tsv
    bgzip ~{sample_id}.counts.tsv

    ###################################################
    ###### CollectSVEvidence
    ###################################################

    java -Xmx~{command_mem_mb}m -jar /opt/gatk.jar CollectSVEvidence \
      -I "~{bam_file}" \
      --sample-name "~{sample_id}" \
      -F "~{sd_locs_vcf}" \
      -SR "~{sample_id}.sr.txt.gz" \
      -PE "~{sample_id}.pe.txt.gz" \
      -SD "~{sample_id}.sd.txt.gz" \
      --site-depth-min-mapq "~{site_depth_min_mapq}" \
      --site-depth-min-baseq "~{site_depth_min_baseq}" \
      ~{"-R " + reference_fasta} \
      ~{"-L " + primary_contigs_list}

    ###################################################
    ###### MELT
    ###################################################

    ###### GetWgsMetrics

    # calculate coverage
    java -Xms~{command_mem_mb}m -jar /opt/picard.jar \
      CollectWgsMetrics \
      INPUT=~{bam_file} \
      VALIDATION_STRINGENCY=SILENT \
      REFERENCE_SEQUENCE=~{reference_fasta} \
      INCLUDE_BQ_HISTOGRAM=true \
      OUTPUT="raw_~{metrics_file_name}" \
      USE_FAST_ALGORITHM=~{fast_algorithm}

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

    COVERAGE=$(fgrep -w "MEAN_COVERAGE" "~{metrics_file_name}" | cut -f2)
    echo "COVERAGE=${COVERAGE}"

    ###### FilterHighCoverageIntervals

    zgrep -v -E "^@|^CONTIG" ~{sample_id}.counts.tsv.gz | awk -F "\t" -v OFS="\t" '{if ($4 > ~{threshold}) {print $1,$2,$3}}' > highCountIntervals.bed

    java -Xmx~{command_mem_mb}m -jar /opt/gatk.jar PrintReads \
      -XL highCountIntervals.bed \
      --interval-exclusion-padding ~{interval_padding} \
      -I ~{bam_file} \
      -O "filtered.bam"


    ###### GetMultipleMetrics

    java -Xmx~{command_mem_mb}m -jar /opt/gatk.jar CollectMultipleMetrics \
      -I "~{bam_file}" \
      -O "~{metrics_base}" \
      -R "~{reference_fasta}" \
      --ASSUME_SORTED true \
      --PROGRAM null \
      --PROGRAM CollectAlignmentSummaryMetrics \
      --PROGRAM CollectInsertSizeMetrics \
      --PROGRAM CollectSequencingArtifactMetrics \
      --PROGRAM CollectGcBiasMetrics \
      --PROGRAM QualityScoreDistribution \
      --METRIC_ACCUMULATION_LEVEL null \
      --METRIC_ACCUMULATION_LEVEL SAMPLE

    function transpose_table2() {
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
      | transpose_table2 \
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
      | transpose_table2 \
      | grep -Ev "^LIBRARY\s|^SAMPLE\s|^READ_GROUP\s" \
      >> "~{metrics_table_filename}"

    # get gc bias summary metrics from sample,
    # and transpose to format compatible with cromwell read_map
    # to get table with properties:
    #   WINDOW_SIZE, TOTAL_CLUSTERS, ALIGNED_READS, AT_DROPOUT, GC_DROPOUT, GC_NC_0_19, GC_NC_20_39,
    #   GC_NC_40_59, GC_NC_60_79, GC_NC_80_100
    grep -A1 "^ACCUMULATION_LEVEL" "~{gc_bias_filename}" \
      | transpose_table2 \
      | grep -Ev "^ACCUMULATION_LEVEL\s|^READS_USED\s|^LIBRARY\s|^SAMPLE\s|^READ_GROUP\s" \
      >> "~{metrics_table_filename}"


    READ_LENGTH=$(fgrep -w "MEAN_READ_LENGTH" "~{metrics_table_filename}" | cut -f2)
    echo "READ_LENGTH=${READ_LENGTH}"
    INSERT_SIZE=$(fgrep -w "MEAN_INSERT_SIZE" "~{metrics_table_filename}" | cut -f2)
    echo "INSERT_SIZE=${INSERT_SIZE}"

    ###### RunMELT

    # MELT expects the BAM index to have extension ".bam.bai"
    mv filtered.bam ~{sample_id}.bam
    mv filtered.bai ~{sample_id}.bam.bai

    # these locations should be stable
    MELT_DIR="/MELT"
    CROMWELL_ROOT=$(pwd)

    # these locations may vary based on MELT version number, so find them:
    MELT_ROOT=$(find "$MELT_DIR" -name "MELT.jar" | xargs -n1 dirname)
    MELT_SCRIPT=$(ls "$MELT_DIR/run_MELT"*.sh)

    # call MELT
    "$MELT_SCRIPT" \
      "~{sample_id}.bam" \
      "~{reference_fasta}" \
      $COVERAGE \
      $READ_LENGTH \
      $INSERT_SIZE \
      "$MELT_ROOT" \
      "$CROMWELL_ROOT" \
      ~{reference_version}

    # combine different mobile element VCFs into a single sample VCF
    # then sort into position order and compress with bgzip
    cat SVA.final_comp.vcf | grep "^#" > "~{sample_id}.header.txt"
    cat SVA.final_comp.vcf | grep -v "^#" > "~{sample_id}.sva.vcf"
    cat LINE1.final_comp.vcf | grep -v "^#"> "~{sample_id}.line1.vcf"
    cat ALU.final_comp.vcf | grep -v "^#"> "~{sample_id}.alu.vcf"
    cat ~{sample_id}.header.txt ~{sample_id}.sva.vcf \
      ~{sample_id}.line1.vcf ~{sample_id}.alu.vcf \
      | vcf-sort -c | bgzip -c > ~{sample_id}.melt.vcf.gz

    # fix three known problems with MELT output vcf:
    # 1) header sometimes doesn't have all chromosomes, and misses
    #    some INF annotations;
    # 2) sample name in the last column of header is problematic;
    # 3) some INFO values contain space " "
    fix_melt() {
      vcf_file=$1
      vcf_text=$(bgzip -cd $vcf_file)
      # replace basename of bam file in header with sample id
      grep '^#CHR' <<<"$vcf_text" | sed 's|~{basename(bam_file, ".bam")}|~{sample_id}|g' > temp_fix_header.txt
      # remove space from INFO values
      grep -v '^#' <<<"$vcf_text" | sed 's/No Difference/No_Difference/g' >> temp_fix_header.txt
      FILEDATE=$(grep -F 'fileDate=' <<<"$vcf_text")
      cat "~{melt_standard_vcf_header}" temp_fix_header.txt \
      | sed "2i$FILEDATE" | bgzip -c > "$vcf_file"
    }
    # apply the fix
    fix_melt ~{sample_id}.melt.vcf.gz

    # rename sample id
    mv ~{sample_id}.melt.vcf.gz temp.vcf.gz
    bcftools reheader -s <(echo "~{sample_id}") temp.vcf.gz > ~{sample_id}.melt.vcf.gz
    rm temp.vcf.gz

    # index vcf
    tabix -p vcf "~{sample_id}.melt.vcf.gz"

    rm ~{sample_id}.bam

    ###################################################
    ###### MANTA
    ###################################################

    # if a preemptible instance restarts and runWorkflow.py already
    # exists, manta will throw an error
    if [ -f ./runWorkflow.py ]; then
    rm ./runWorkflow.py
    fi

    ln -s ~{bam_file} sample.bam
    ln -s ~{bam_index} sample.bam.bai

    # prepare the analysis job
    /usr/local/bin/manta/bin/configManta.py \
    --bam sample.bam \
    --referenceFasta ~{reference_fasta} \
    --runDir . \
    --callRegions ~{manta_region_bed}

    # always tell manta there are 2 GiB per job, otherwise it will
    # scale back the requested number of jobs, even if they won't
    # need that much memory
    ./runWorkflow.py \
    --mode local \
    --jobs ~{manta_num_jobs} \
    --memGb $((~{manta_num_jobs} * 2))

    # inversion conversion, then compression and index
    python2 /usr/local/bin/manta/libexec/convertInversion.py \
    /usr/local/bin/samtools \
    ~{reference_fasta} \
    results/variants/diploidSV.vcf.gz \
    | bcftools reheader -s <(echo "~{sample_id}") \
    > diploidSV.vcf

    bgzip -c diploidSV.vcf > ~{sample_id}.manta.vcf.gz
    tabix -p vcf ~{sample_id}.manta.vcf.gz

    ###################################################
    ###### WHAM
    ###################################################

    # print some info that may be useful for debugging
    df -h
    echo "whamg $(./whamg 2>&1 >/dev/null | grep Version)"

    # ensure that index files are present in appropriate locations
    if [ ! -e "~{bam_file}.bai" ]; then
      ln -s "~{bam_index}" "~{bam_file}.bai"
    fi
    REPLACED_VERSION=$(echo "~{bam_file}" | sed 's/.bam$/.bai/g')
    if [ ! -e "$REPLACED_VERSION" ]; then
      ln -s "~{bam_index}" "$REPLACED_VERSION"
    fi
    if [ ! -e "~{reference_fasta}.fai" ]; then
      ln -s "~{reference_index}" "~{reference_fasta}.fai"
    fi

    # run whamg on all good intervals
    GOOD_INTERVALS=$(echo "~{sep='\n' good_intervals}" \
      | awk '{print $1 ":" $2 "-" $3}')
    VCFS=""
    for GOOD_INTERVAL in $GOOD_INTERVALS; do
    echo "Analyzing $GOOD_INTERVAL"
    NEW_VCF="$GOOD_INTERVAL.wham.vcf.gz"
    # run wham, gzip results
    whamg \
      -c "~{sep=","  chr_list}" \
      -x ~{num_cpu} \
      -a "~{reference_fasta}"\
      -f "~{bam_file}" \
      -r "$GOOD_INTERVAL" \
      | bgzip -c > "$NEW_VCF"
    # index results vcf
    tabix "$NEW_VCF"
    # add vcf file name to list of vcfs, separated by spaces
    if [ -z "$VCFS" ]; then
      VCFS="$NEW_VCF"
    else
      VCFS="$VCFS $NEW_VCF"
    fi
    done

    # We need to update both the VCF sample ID and the TAGS INFO field in the WHAM output VCFs.
    # WHAM uses both to store the sample identifier, and by default uses the SM identifier from the BAM file.
    # We need to update both to reflect the potentially-renamed sample identifier used by the pipeline (sample_id) --
    # svtk standardize_vcf uses the TAGS field to identify the sample for WHAM VCFs.

    VCF_FILE_FIXED_HEADER="~{sample_id}.wham.vcf.gz"
    VCF_FILE_BAD_TAGS="~{sample_id}.wham_bad_header_bad_tags.vcf.gz"
    VCF_FILE_BAD_HEADER="~{sample_id}.wham_bad_header.vcf.gz"

    # concatenate resulting vcfs into final vcf
    echo "Concatenating results"
    bcftools concat -a -O z -o "$VCF_FILE_BAD_TAGS" $VCFS
    tabix -p vcf "$VCF_FILE_BAD_TAGS"

    # write out a an annotation table with the new sample ID for each variant record.
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t~{sample_id}\n' $VCF_FILE_BAD_TAGS | bgzip -c > tags_annotation_file.tsv.gz
    tabix -f -s1 -b2 -e2 tags_annotation_file.tsv.gz
    # Use bcftools annotate to re-write the TAGS field in each variant based on the annotation table.
    bcftools annotate -a tags_annotation_file.tsv.gz -c CHROM,POS,REF,ALT,INFO/TAGS  $VCF_FILE_BAD_TAGS | bgzip -c > $VCF_FILE_BAD_HEADER
    tabix -p vcf "$VCF_FILE_BAD_HEADER"

    echo "Getting existing header"
    # get the existing header
    OLD_HEADER=$(bcftools view -h "$VCF_FILE_BAD_HEADER" | grep -v '##contig=')
    # create new header lines with the contig lengths
    echo "Making grep pattern to extract contigs from reference index"
    CONTIGS_PATTERN="^$(cat ~{primary_contigs_list} | paste -sd "," - | sed "s/,/\\\t|^/g")\t"
    echo "Adding contigs with length to vcf header"
    CONTIGS_HEADER=$(grep "$CONTIGS_PATTERN" -P "~{reference_index}" | awk '{print "##contig=<ID=" $1 ",length=" $2 ">"}')
    # Create a new header with
    #   -all but last line of old header, followed by newline
    #   -contig header lines, followed by newline
    #   -last line of old header
    # Replace old header with new header
    echo "Replacing header"
    echo "$OLD_HEADER" | grep -v "^#CHROM" > new_header.txt
    echo "$CONTIGS_HEADER" >> new_header.txt
    echo "$OLD_HEADER" | grep "^#CHROM" >> new_header.txt
    bcftools reheader \
      -h new_header.txt \
      -s <(echo "~{sample_id}") \
      "$VCF_FILE_BAD_HEADER" \
      > "$VCF_FILE_FIXED_HEADER"

    echo "Indexing vcf"
    tabix "$VCF_FILE_FIXED_HEADER"

    echo "finished RunWhamg"

  >>>
  output {
    File coverage_counts = "~{sample_id}.counts.tsv.gz"

    File manta_vcf = "~{sample_id}.manta.vcf.gz"
    File manta_index = "~{sample_id}.manta.vcf.gz.tbi"

    File melt_vcf = "~{sample_id}.melt.vcf.gz"
    File melt_index = "~{sample_id}.melt.vcf.gz.tbi"
    File wgs_metrics_file = "~{metrics_file_name}"
    File multiple_metrics_file = "~{metrics_table_filename}"

    File pesr_disc = "~{sample_id}.disc.txt.gz"
    File pesr_disc_index = "~{sample_id}.disc.txt.gz.tbi"
    File pesr_split = "~{sample_id}.split.txt.gz"
    File pesr_split_index = "~{sample_id}.split.txt.gz.tbi"

    File wham_vcf = "~{sample_id}.wham.vcf.gz"
    File wham_index = "~{sample_id}.wham.vcf.gz.tbi"

    File split_out = "~{sample_id}.sr.txt.gz"
    File split_out_index = "~{sample_id}.sr.txt.gz.tbi"
    File disc_out = "~{sample_id}.pe.txt.gz"
    File disc_out_index = "~{sample_id}.pe.txt.gz.tbi"
    File sd_out = "~{sample_id}.sd.txt.gz"
    File sd_out_index = "~{sample_id}.sd.txt.gz.tbi"
  }
  runtime {
    docker: sv_pipeline_callers_docker
  }
}