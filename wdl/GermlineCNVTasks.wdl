version 1.0

import "Structs.wdl"

task AnnotateIntervals {
    input {
      File intervals
      File ref_fasta
      File ref_fasta_fai
      File ref_fasta_dict
      File? mappability_track_bed
      File? mappability_track_bed_idx
      File? segmental_duplication_track_bed
      File? segmental_duplication_track_bed_idx
      Int? feature_query_lookahead
      File? gatk4_jar_override
      String gatk_docker
      RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(select_all([intervals, ref_fasta, mappability_track_bed, segmental_duplication_track_bed]), "GB")

    RuntimeAttr default_attr = object {
      cpu_cores: 1,
      mem_gb: 1.7,
      disk_gb: ceil(10.0 + input_size),
      boot_disk_gb: 10,
      preemptible_tries: 3,
      max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    Float mem_gb = select_first([runtime_attr.mem_gb, default_attr.mem_gb])
    Int command_mem_mb = ceil(mem_gb * 1000 - 500)

    # Determine output filename
    String base_filename = basename(intervals, ".interval_list")

    command <<<
        set -euo pipefail
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk4_jar_override}

        gatk --java-options "-Xmx~{command_mem_mb}m" AnnotateIntervals \
            -L ~{intervals} \
            --reference ~{ref_fasta} \
            ~{"--mappability-track " + mappability_track_bed} \
            ~{"--segmental-duplication-track " + segmental_duplication_track_bed} \
            --feature-query-lookahead ~{default=1000000 feature_query_lookahead} \
            --interval-merging-rule OVERLAPPING_ONLY \
            --output ~{base_filename}.annotated.tsv
    >>>
    runtime {
      cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
      memory: mem_gb + " GiB"
      disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
      bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
      docker: gatk_docker
      preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
      maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }

    output {
        File annotated_intervals = "~{base_filename}.annotated.tsv"
    }
}

task FilterIntervals {
    input {
      File intervals
      File? exclude_intervals
      File? annotated_intervals
      Array[File] read_count_files
      Float? minimum_gc_content
      Float? maximum_gc_content
      Float? minimum_mappability
      Float? maximum_mappability
      Float? minimum_segmental_duplication_content
      Float? maximum_segmental_duplication_content
      Int? low_count_filter_count_threshold
      Float? low_count_filter_percentage_of_samples
      Float? extreme_count_filter_minimum_percentile
      Float? extreme_count_filter_maximum_percentile
      Float? extreme_count_filter_percentage_of_samples
      File? gatk4_jar_override
      String gatk_docker
      RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
      cpu_cores: 1,
      mem_gb: 6.0,
      disk_gb: 10,
      boot_disk_gb: 10,
      preemptible_tries: 3,
      max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    Float mem_gb = select_first([runtime_attr.mem_gb, default_attr.mem_gb])
    Int command_mem_mb = ceil(mem_gb * 1000 - 500)

    # Determine output filename
    String base_filename = basename(intervals, ".interval_list")

    command <<<
        set -euo pipefail
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk4_jar_override}

        read_count_files_list=~{write_lines(read_count_files)}
        grep gz$ $read_count_files_list | xargs -l1 -P0 gunzip
        sed 's/\.gz$//' $read_count_files_list | \
            awk '{print "--input "$0}' > read_count_files.args

        gatk --java-options "-Xmx~{command_mem_mb}m" FilterIntervals \
            -L ~{intervals} \
            ~{"-XL " + exclude_intervals} \
            ~{"--annotated-intervals " + annotated_intervals} \
            --arguments_file read_count_files.args \
            --minimum-gc-content ~{default="0.1" minimum_gc_content} \
            --maximum-gc-content ~{default="0.9" maximum_gc_content} \
            --minimum-mappability ~{default="0.9" minimum_mappability} \
            --maximum-mappability ~{default="1.0" maximum_mappability} \
            --minimum-segmental-duplication-content ~{default="0.0" minimum_segmental_duplication_content} \
            --maximum-segmental-duplication-content ~{default="0.5" maximum_segmental_duplication_content} \
            --low-count-filter-count-threshold ~{default="5" low_count_filter_count_threshold} \
            --low-count-filter-percentage-of-samples ~{default="90.0" low_count_filter_percentage_of_samples} \
            --extreme-count-filter-minimum-percentile ~{default="1.0" extreme_count_filter_minimum_percentile} \
            --extreme-count-filter-maximum-percentile ~{default="99.0" extreme_count_filter_maximum_percentile} \
            --extreme-count-filter-percentage-of-samples ~{default="90.0" extreme_count_filter_percentage_of_samples} \
            --interval-merging-rule OVERLAPPING_ONLY \
            --output ~{base_filename}.filtered.interval_list
    >>>
    runtime {
      cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
      memory: mem_gb + " GiB"
      disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
      bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
      docker: gatk_docker
      preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
      maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }

    output {
        File filtered_intervals = "~{base_filename}.filtered.interval_list"
    }
}

task ScatterIntervals {
    input{
      File interval_list
      Int num_intervals_per_scatter
      String? output_dir
      File? gatk4_jar_override
      String gatk_docker
      RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
      cpu_cores: 1,
      mem_gb: 2.0,
      disk_gb: 10,
      boot_disk_gb: 10,
      preemptible_tries: 3,
      max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    Float mem_gb = select_first([runtime_attr.mem_gb, default_attr.mem_gb])
    Int command_mem_mb = ceil(mem_gb * 1000 - 500)

    # If optional output_dir not specified, use "out";
    String output_dir_ = select_first([output_dir, "out"])

    String base_filename = basename(interval_list, ".interval_list")

    command <<<
        set -euo pipefail
        mkdir ~{output_dir_}
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk4_jar_override}

        {
            >&2 echo "Attempting to run IntervalListTools..."
            gatk --java-options "-Xmx~{command_mem_mb}m" IntervalListTools \
                --INPUT ~{interval_list} \
                --SUBDIVISION_MODE INTERVAL_COUNT \
                --SCATTER_CONTENT ~{num_intervals_per_scatter} \
                --OUTPUT ~{output_dir_} &&
            # output files are named output_dir_/temp_0001_of_N/scattered.interval_list, etc. (N = num_intervals_per_scatter);
            # we rename them as output_dir_/base_filename.scattered.0000.interval_list, etc.
            ls ~{output_dir_}/*/scattered.interval_list | \
                cat -n | \
                while read n filename; do mv $filename ~{output_dir_}/~{base_filename}.scattered.$(printf "%04d" $n).interval_list; done
        } || {
            # if only a single shard is required, then we can just rename the original interval list
            >&2 echo "IntervalListTools failed because only a single shard is required. Copying original interval list..."
            cp ~{interval_list} ~{output_dir_}/~{base_filename}.scattered.1.interval_list
        }
    >>>
    runtime {
      cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
      memory: mem_gb + " GiB"
      disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
      bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
      docker: gatk_docker
      preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
      maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }

    output {
        Array[File] scattered_interval_lists = glob("~{output_dir_}/~{base_filename}.scattered.*.interval_list")
    }
}

task ExplodePloidyCalls {
    input {
      File contig_ploidy_calls_tar
      Array[String] samples
      String linux_docker
      RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
      cpu_cores: 1,
      mem_gb: 1.5,
      disk_gb: 10,
      boot_disk_gb: 10,
      preemptible_tries: 3,
      max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    Int num_samples = length(samples)
    String out_dir = "calls_renamed"

    command <<<
      set -euo pipefail

      # Extract ploidy calls
      mkdir calls
      tar xzf ~{contig_ploidy_calls_tar} -C calls/

      # Archive call files by sample, renaming so they will be glob'd in order
      sample_ids=(~{sep=" " samples})
      for (( i=0; i<~{num_samples}; i++ ))
      do
        sample_id=${sample_ids[$i]}
        sample_no=`printf %06d $i`
        tar -czf sample_${sample_no}.${sample_id}.contig_ploidy_calls.tar.gz -C calls/SAMPLE_${i} .
      done
    >>>
    runtime {
      cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
      memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
      disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
      bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
      docker: linux_docker
      preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
      maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }

    output {
        Array[File] sample_contig_ploidy_calls_tar = glob("sample_*.contig_ploidy_calls.tar.gz")
    }
}

task BundlePostprocessingInvariants {
    input {
        Array[File] calls_tars
        Array[File] model_tars
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_override
    }

    Int disk_gb_base = 10
    Float files_size_gb = size(calls_tars, "GiB") + size(model_tars, "GiB")
    # Usage factor takes into account that files sometimes sit on disk in both compressed and decompressed form
    Float usage_factor = 10.0
    Int disk_gb = disk_gb_base + ceil(usage_factor * files_size_gb)

    RuntimeAttr default_attr = object {
      cpu_cores: 1,
      mem_gb: 2.0,
      disk_gb: disk_gb,
      boot_disk_gb: 10,
      preemptible_tries: 3,
      max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    command <<<
        set -euo pipefail
        mkdir -p out

        calls_files_tar_list=~{write_lines(calls_tars)}
        model_files_tar_list=~{write_lines(model_tars)}

        cat $calls_files_tar_list | sort -V > calls_files_tar_list.sorted
        cat $model_files_tar_list | sort -V > model_files_tar_list.sorted

        paste calls_files_tar_list.sorted model_files_tar_list.sorted |\
              awk '{print (NR-1)"\t"$0}' > file_pairs.sorted
        OIFS=$IFS
        IFS=$'\t'
        NUM_PAIRS=$(($(wc -l < file_pairs.sorted)))
        NUM_COMPLETED=0
        while read index calls_tar model_tar; do
            printf "%s\t %d/%d pairs uncompressed\n" "$(date)" $NUM_COMPLETED $NUM_PAIRS
            mkdir -p out/CALLS_$index
            mkdir -p out/MODEL_$index
            tar xzf $calls_tar -C out/CALLS_$index
            tar xzf $model_tar -C out/MODEL_$index
            rm $calls_tar $model_tar
            ((++NUM_COMPLETED))
        done < file_pairs.sorted
        printf "%s\t %d/%d pairs uncompressed\n" "$(date)" $NUM_COMPLETED $NUM_PAIRS
        IFS=$OIFS

        printf "Recompressing all pairs\n"
        tar c -C out . | gzip -1 > case-gcnv-postprocessing-invariants.tar.gz
        rm -Rf out
    >>>
    runtime {
      cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
      memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
      disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
      bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
      docker: sv_base_mini_docker
      preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
      maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }

    output {
        File bundle_tar = "case-gcnv-postprocessing-invariants.tar.gz"
    }
}

task BundledPostprocessGermlineCNVCalls {
    input {
        File invariants_tar
        String entity_id
        File contig_ploidy_calls_tar
        Array[String]? allosomal_contigs
        Int ref_copy_number_autosomal_contigs
        Int sample_index
        File? gatk4_jar_override
        String gatk_docker
        RuntimeAttr? runtime_attr_override
    }

    Int disk_gb_baseline = 10
    Float disk_usage_factor = 10.0
    Int disk_gb = disk_gb_baseline + ceil(disk_usage_factor * size([invariants_tar, contig_ploidy_calls_tar, gatk4_jar_override], "GiB"))

    RuntimeAttr default_attr = object {
      cpu_cores: 1,
      mem_gb: 8.5,
      disk_gb: disk_gb,
      boot_disk_gb: 10,
      preemptible_tries: 3,
      max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    Float mem_gb = select_first([runtime_attr.mem_gb, default_attr.mem_gb])
    Int command_mem_mb = ceil(mem_gb * 1000 - 500)

    String genotyped_intervals_vcf_filename = "genotyped-intervals-~{entity_id}.vcf.gz"
    String genotyped_segments_vcf_filename = "genotyped-segments-~{entity_id}.vcf.gz"
    Boolean allosomal_contigs_specified = defined(allosomal_contigs)

    command <<<
        set -euo pipefail
        
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk4_jar_override}

        # untar calls to CALLS_0, CALLS_1, etc directories and build the command line
        # also copy over shard config and interval files
        time tar xzf ~{invariants_tar}
        rm ~{invariants_tar}
        number_of_shards=`find . -name 'CALLS_*' | wc -l` 

        touch calls_and_model_args.txt
        for i in $(seq 0 `expr $number_of_shards - 1`); do 
            echo "--calls-shard-path CALLS_$i" >> calls_and_model_args.txt
            echo "--model-shard-path MODEL_$i" >> calls_and_model_args.txt
        done

        mkdir -p extracted-contig-ploidy-calls
        tar xzf ~{contig_ploidy_calls_tar} -C extracted-contig-ploidy-calls
        rm ~{contig_ploidy_calls_tar}

        allosomal_contigs_args="--allosomal-contig ~{sep=" --allosomal-contig " allosomal_contigs}"

        time gatk --java-options "-Xmx~{command_mem_mb}m" PostprocessGermlineCNVCalls \
             --arguments_file calls_and_model_args.txt \
            ~{true="$allosomal_contigs_args" false="" allosomal_contigs_specified} \
            --autosomal-ref-copy-number ~{ref_copy_number_autosomal_contigs} \
            --contig-ploidy-calls extracted-contig-ploidy-calls \
            --sample-index ~{sample_index} \
            --output-genotyped-intervals ~{genotyped_intervals_vcf_filename} \
            --output-genotyped-segments ~{genotyped_segments_vcf_filename}

        rm -Rf extracted-contig-ploidy-calls
    >>>
    runtime {
      cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
      memory: mem_gb + " GiB"
      disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
      bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
      docker: gatk_docker
      preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
      maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }

    output {
        File genotyped_intervals_vcf = genotyped_intervals_vcf_filename
        File genotyped_segments_vcf = genotyped_segments_vcf_filename
    }
}


task PostprocessGermlineCNVCalls {
    input {
        String entity_id
        Array[File] gcnv_calls_tars
        Array[File] gcnv_model_tars
        Array[File] calling_configs
        Array[File] denoising_configs
        Array[File] gcnvkernel_version
        Array[File] sharded_interval_lists
        File contig_ploidy_calls_tar
        Array[String]? allosomal_contigs
        Int ref_copy_number_autosomal_contigs
        Int sample_index
        File? gatk4_jar_override
        String gatk_docker
        RuntimeAttr? runtime_attr_override
    }

    Int disk_gb_baseline = 10
    Float disk_usage_factor = 10.0
    Int disk_gb = disk_gb_baseline + ceil(disk_usage_factor * (
        size(gcnv_calls_tars, "GiB") + size(gcnv_model_tars, "GiB") + size(calling_configs, "GiB") +
        size(denoising_configs, "GiB") + size(gcnvkernel_version, "GiB") +
        size(sharded_interval_lists, "GiB") + size(contig_ploidy_calls_tar, "GiB")
    ))

    RuntimeAttr default_attr = object {
      cpu_cores: 1,
      mem_gb: 12,
      disk_gb: disk_gb,
      boot_disk_gb: 10,
      preemptible_tries: 3,
      max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    Float mem_gb = select_first([runtime_attr.mem_gb, default_attr.mem_gb])
    Int command_mem_mb = ceil(mem_gb * 1000 * 0.6)

    String genotyped_intervals_vcf_filename = "genotyped-intervals-~{entity_id}.vcf.gz"
    String genotyped_segments_vcf_filename = "genotyped-segments-~{entity_id}.vcf.gz"
    String denoised_copy_ratios_filename ="denoised_copy_ratios-~{entity_id}.tsv"
    Boolean allosomal_contigs_specified = defined(allosomal_contigs)

    command <<<
        set -euo pipefail

        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk4_jar_override}

        # untar calls to CALLS_0, CALLS_1, etc directories and build the command line
        # also copy over shard config and interval files
        touch calls_and_model_args.txt
        CALL_INDEX=0
        paste ~{write_lines(gcnv_calls_tars)} ~{write_lines(calling_configs)} ~{write_lines(denoising_configs)} \
              ~{write_lines(gcnvkernel_version)} ~{write_lines(sharded_interval_lists)} ~{write_lines(gcnv_model_tars)} \
            | while read GCNV_CALLS_TAR CALLING_CONFIG DENOISING_CONFIG GCNVKERNEL_VERSION SHARDED_INTERVAL_LIST GCNV_MODEL_TAR; do
                CALL_DIR="CALLS_$CALL_INDEX"
                mkdir -p $CALL_DIR/SAMPLE_~{sample_index}
                tar xzf $GCNV_CALLS_TAR -C $CALL_DIR/SAMPLE_~{sample_index}
                ln -rs $CALLING_CONFIG $CALL_DIR
                ln -rs $DENOISING_CONFIG $CALL_DIR
                ln -rs $GCNVKERNEL_VERSION $CALL_DIR
                ln -rs $SHARDED_INTERVAL_LIST $CALL_DIR
                echo "--calls-shard-path $CALL_DIR" >> calls_and_model_args.txt

                MODEL_DIR="MODEL_$CALL_INDEX"
                mkdir $MODEL_DIR
                tar xzf $GCNV_MODEL_TAR -C "$MODEL_DIR"
                echo "--model-shard-path $MODEL_DIR" >> calls_and_model_args.txt
                ((++CALL_INDEX))
              done


        mkdir -p contig-ploidy-calls
        tar xzf ~{contig_ploidy_calls_tar} -C contig-ploidy-calls
        rm ~{contig_ploidy_calls_tar}

        allosomal_contigs_args="--allosomal-contig ~{sep=" --allosomal-contig " allosomal_contigs}"

        time gatk --java-options "-Xmx~{command_mem_mb}m" PostprocessGermlineCNVCalls \
            --arguments_file calls_and_model_args.txt \
            ~{true="$allosomal_contigs_args" false="" allosomal_contigs_specified} \
            --autosomal-ref-copy-number ~{ref_copy_number_autosomal_contigs} \
            --contig-ploidy-calls contig-ploidy-calls \
            --sample-index ~{sample_index} \
            --output-genotyped-intervals ~{genotyped_intervals_vcf_filename} \
            --output-genotyped-segments ~{genotyped_segments_vcf_filename} \
            --output-denoised-copy-ratios ~{denoised_copy_ratios_filename}
    >>>
    runtime {
      cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
      memory: mem_gb + " GiB"
      disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
      bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
      docker: gatk_docker
      preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
      maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }

    output {
        File genotyped_intervals_vcf = genotyped_intervals_vcf_filename
        File genotyped_segments_vcf = genotyped_segments_vcf_filename
        File denoised_copy_ratios = denoised_copy_ratios_filename
    }
}
