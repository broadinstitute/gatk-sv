version 1.0

import "Structs.wdl"

task SVCluster {
    input {
        # Either vcfs of vcfs_tar should be provided
        Array[File] vcfs = []  # Can't use optional because of write_lines() call
        File? vcfs_tar

        File ploidy_table
        String output_prefix

        String? contig

        Boolean? fast_mode
        Boolean? omit_members
        Boolean? enable_cnv
        Boolean? default_no_call

        String? algorithm
        String? insertion_length_summary_strategy
        String? breakpoint_summary_strategy
        String? alt_allele_summary_strategy

        Float? defrag_padding_fraction
        Float? defrag_sample_overlap

        Float? depth_sample_overlap
        Float? depth_interval_overlap
        Float? depth_size_similarity
        Int? depth_breakend_window
        Float? mixed_sample_overlap
        Float? mixed_interval_overlap
        Float? mixed_size_similarity
        Int? mixed_breakend_window
        Float? pesr_sample_overlap
        Float? pesr_interval_overlap
        Float? pesr_size_similarity
        Int? pesr_breakend_window

        File reference_fasta
        File reference_fasta_fai
        File reference_dict

        Float? java_mem_fraction
        String? additional_args
        String? variant_prefix

        String gatk_docker
        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        vcfs: {
                  localization_optional: true
              }
    }

    RuntimeAttr default_attr = object {
                                   cpu_cores: 1,
                                   mem_gb: 3.75,
                                   disk_gb: ceil(10 + size(vcfs, "GB") * 3 + size(reference_fasta, "GB")),
                                   boot_disk_gb: 10,
                                   preemptible_tries: 3,
                                   max_retries: 1
                               }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output {
        File out = "~{output_prefix}.vcf.gz"
        File out_index = "~{output_prefix}.vcf.gz.tbi"
    }
    command <<<
        set -euxo pipefail

        function getJavaMem() {
            # get JVM memory in MiB by getting total memory from /proc/meminfo
            # and multiplying by java_mem_fraction
            cat /proc/meminfo \
                | awk -v MEM_FIELD="$1" '{
                    f[substr($1, 1, length($1)-1)] = $2
                } END {
                    printf "%dM", f[MEM_FIELD] * ~{default="0.85" java_mem_fraction} / 1024
                }'
        }
        JVM_MAX_MEM=$(getJavaMem MemTotal)
        echo "JVM memory: $JVM_MAX_MEM"

        if ~{length(vcfs) > 0}; then
            awk '{print "-V "$0}' ~{write_lines(vcfs)} > arguments.txt
        elif ~{defined(vcfs_tar)}; then
            mkdir vcfs
            tar xzf ~{vcfs_tar} -C vcfs/
            ls vcfs/*.vcf.gz | awk '{print "-V "$0}' > arguments.txt
        else
            echo "ERROR: neither vcfs nor vcfs_tar was provided"
            exit 1
        fi

        gatk --java-options "-Xmx${JVM_MAX_MEM}" SVCluster \
            --arguments_file arguments.txt \
            --output ~{output_prefix}.vcf.gz \
            --ploidy-table ~{ploidy_table} \
            --reference ~{reference_fasta} \
            ~{"-L " + contig} \
            ~{true="--fast-mode" false="" fast_mode} \
            ~{true="--enable-cnv" false="" enable_cnv} \
            ~{true="--omit-members" false="" omit_members} \
            ~{true="--default-no-call" false="" default_no_call} \
            ~{"--variant-prefix " + variant_prefix} \
            ~{"--algorithm " + algorithm} \
            ~{"--defrag-padding-fraction " + defrag_padding_fraction} \
            ~{"--defrag-sample-overlap " + defrag_sample_overlap} \
            ~{"--depth-sample-overlap " + depth_sample_overlap} \
            ~{"--depth-interval-overlap " + depth_interval_overlap} \
            ~{"--depth-size-similarity " + depth_size_similarity} \
            ~{"--depth-breakend-window " + depth_breakend_window} \
            ~{"--mixed-sample-overlap " + mixed_sample_overlap} \
            ~{"--mixed-interval-overlap " + mixed_interval_overlap} \
            ~{"--mixed-size-similarity " + mixed_size_similarity} \
            ~{"--mixed-breakend-window " + mixed_breakend_window} \
            ~{"--pesr-sample-overlap " + pesr_sample_overlap} \
            ~{"--pesr-interval-overlap " + pesr_interval_overlap} \
            ~{"--pesr-size-similarity " + pesr_size_similarity} \
            ~{"--pesr-breakend-window " + pesr_breakend_window} \
            ~{"--insertion-length-summary-strategy " + insertion_length_summary_strategy} \
            ~{"--breakpoint-summary-strategy " + breakpoint_summary_strategy} \
            ~{"--alt-allele-summary-strategy " + alt_allele_summary_strategy} \
            ~{additional_args}
    >>>
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: gatk_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task ExcludeIntervalsByEndpoints {
    input {
        File vcf
        File reference_fasta_fai
        File intervals
        File intervals_index
        String output_prefix
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
                                   cpu_cores: 1,
                                   mem_gb: 3.75,
                                   disk_gb: ceil(10 + size(vcf, "GB") * 2),
                                   boot_disk_gb: 10,
                                   preemptible_tries: 3,
                                   max_retries: 1
                               }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output {
        File out = "~{output_prefix}.vcf.gz"
        File out_index = "~{output_prefix}.vcf.gz.tbi"
    }
    command <<<
        set -euo pipefail
        cut -f1,2 ~{reference_fasta_fai} > genome.file
        bcftools query -f '%CHROM\t%POS\t%POS\t%ID\t%SVTYPE\n%CHROM\t%END\t%END\t%ID\t%SVTYPE\n%CHR2\t%END2\t%END2\t%ID\t%SVTYPE\n' ~{vcf} \
            | awk '$1!="."' \
            | sort -k1,1V -k2,2n -k3,3n \
            > ends.bed
        bedtools intersect -sorted -u -wa -g genome.file -wa -a ends.bed -b ~{intervals} | cut -f4 | sort | uniq \
            > excluded_vids.list
        bcftools view -i '%ID!=@excluded_vids.list' ~{vcf} -Oz -o ~{output_prefix}.vcf.gz
        tabix ~{output_prefix}.vcf.gz
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
}

task ExcludeIntervalsByIntervalOverlap {
    input {
        File vcf
        Float overlap_fraction
        File intervals
        File intervals_index
        File reference_fasta_fai
        String output_prefix
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
                                   cpu_cores: 1,
                                   mem_gb: 3.75,
                                   disk_gb: ceil(10 + size(vcf, "GB") * 2),
                                   boot_disk_gb: 10,
                                   preemptible_tries: 3,
                                   max_retries: 1
                               }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output {
        File out = "~{output_prefix}.vcf.gz"
        File out_index = "~{output_prefix}.vcf.gz.tbi"
    }
    command <<<
        set -euo pipefail
        cut -f1,2 ~{reference_fasta_fai} > genome.file
        bcftools query -f '%CHROM\t%POS\t%END\t%ID\t%SVTYPE\n' ~{vcf} > variants.bed
        bedtools coverage -sorted -g genome.file -f ~{overlap_fraction} -a variants.bed -b ~{intervals} \
            | awk -F"\t" '$6>0' \
            | cut -f4 \
            > excluded_vids.list
        bcftools view -i '%ID!=@excluded_vids.list' ~{vcf} -Oz -o ~{output_prefix}.vcf.gz
        tabix ~{output_prefix}.vcf.gz
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
}

task GatkToSvtkVcf {
    input {
        File vcf
        File? script
        String source
        File contig_list
        String? remove_infos
        String? remove_formats
        Boolean set_pass = false
        String output_prefix
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
                                   cpu_cores: 1,
                                   mem_gb: 3.75,
                                   disk_gb: ceil(10 + size(vcf, "GB") * 2),
                                   boot_disk_gb: 10,
                                   preemptible_tries: 3,
                                   max_retries: 1
                               }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output {
        File out = "~{output_prefix}.vcf.gz"
        File out_index = "~{output_prefix}.vcf.gz.tbi"
    }
    command <<<
        set -euo pipefail
        python ~{default="/opt/sv-pipeline/scripts/format_gatk_vcf_for_svtk.py" script} \
            --vcf ~{vcf} \
            --out ~{output_prefix}.vcf.gz \
            --source ~{source} \
            --contigs ~{contig_list} \
            ~{"--remove-infos " + remove_infos} \
            ~{"--remove-formats " + remove_formats} \
            ~{if set_pass then "--set-pass" else ""}
        tabix ~{output_prefix}.vcf.gz
    >>>
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: sv_pipeline_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task CNVBedToGatkVcf {
    input {
        File bed
        File? script
        File sample_list
        File contig_list
        File ploidy_table
        File? reference_fasta_fai
        String vid_prefix
        String output_prefix
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
                                   cpu_cores: 1,
                                   mem_gb: 3.75,
                                   disk_gb: ceil(10 + size(bed, "GB") * 2),
                                   boot_disk_gb: 10,
                                   preemptible_tries: 3,
                                   max_retries: 1
                               }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output {
        File out = "~{output_prefix}.vcf.gz"
        File out_index = "~{output_prefix}.vcf.gz.tbi"
    }
    command <<<
        set -euo pipefail
        python ~{default="/opt/sv-pipeline/scripts/convert_bed_to_gatk_vcf.py" script} \
            --bed ~{bed} \
            --out ~{output_prefix}.vcf.gz \
            --samples ~{sample_list} \
            --contigs ~{contig_list} \
            --vid-prefix ~{vid_prefix} \
            --ploidy-table ~{ploidy_table} \
            ~{"--fai " + reference_fasta_fai}
        tabix ~{output_prefix}.vcf.gz
    >>>
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: sv_pipeline_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task CreatePloidyTableFromPed {
    input {
        File ped_file
        File? script
        File contig_list
        Boolean retain_female_chr_y = false
        String? chr_x
        String? chr_y
        String output_prefix
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
                                   cpu_cores: 1,
                                   mem_gb: 3.75,
                                   disk_gb: 10,
                                   boot_disk_gb: 10,
                                   preemptible_tries: 3,
                                   max_retries: 1
                               }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    String output_file = if retain_female_chr_y then "~{output_prefix}.FEMALE_chrY_1.tsv" else "~{output_prefix}.tsv"

    output {
        File out = "~{output_file}"
    }
    command <<<
        set -euo pipefail
        python ~{default="/opt/sv-pipeline/scripts/ploidy_table_from_ped.py" script} \
            --ped ~{ped_file} \
            --out tmp.tsv \
            --contigs ~{contig_list} \
            ~{"--chr-x " + chr_x} \
            ~{"--chr-y " + chr_y}

        # TODO : For now we retain female Y genotypes for clustering
        if ~{retain_female_chr_y}; then
            sed -e 's/\t0/\t1/g' tmp.tsv > ~{output_file}
        else
            mv tmp.tsv ~{output_file}
        fi
    >>>
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: sv_pipeline_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task ScatterCompressedBedOmitHeaders {
    input {
        File bed
        String prefix
        Int records_per_shard
        Int n_digits = 6
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(bed, "GB")
    Float base_disk_gb = 10.0

    RuntimeAttr runtime_default = object {
                                      mem_gb: 3.75,
                                      disk_gb: ceil(base_disk_gb + input_size * 10.0),
                                      cpu_cores: 1,
                                      preemptible_tries: 3,
                                      max_retries: 1,
                                      boot_disk_gb: 10
                                  }
    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: sv_base_mini_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command <<<
        set -euo pipefail
        mkdir out
        gunzip -c ~{bed} \
            | awk '$0!~"#"' \
            | split -d -a ~{n_digits} -l ~{records_per_shard} - out/~{prefix}
        for file in out/~{prefix}*; do
            mv $file $file.bed
            bgzip $file.bed
        done
    >>>
    output {
        Array[File] out = glob("out/~{prefix}*.bed.gz")
    }
}

task TarFiles {
    input {
        Array[File] files
        String prefix
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(files, "GB")
    Float base_disk_gb = 10.0

    RuntimeAttr runtime_default = object {
                                      mem_gb: 3.75,
                                      disk_gb: ceil(base_disk_gb + input_size * 3.0),
                                      cpu_cores: 1,
                                      preemptible_tries: 3,
                                      max_retries: 1,
                                      boot_disk_gb: 10
                                  }
    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: sv_base_mini_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command <<<
        set -euo pipefail
        FILES="~{write_lines(files)}"
        # Note that all filenames must be unique
        awk 'a[$0]++{print "Duplicate file found:",$0; exit(1)}' $FILES
        # Avoids argument limits encountered with ls; uses hard links avoid issues with local file systems
        mkdir files
        cat $FILES | xargs -I{} -n1 ln {} files/
        tar czf ~{prefix}.tar.gz -h -C files/ .
    >>>
    output {
        File out = "~{prefix}.tar.gz"
    }
}

task TarFilesFromList {
    input {
        File file_list
        String prefix
        String sv_base_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr runtime_default = object {
                                      mem_gb: 3.75,
                                      disk_gb: 100,
                                      cpu_cores: 1,
                                      preemptible_tries: 3,
                                      max_retries: 1,
                                      boot_disk_gb: 10
                                  }
    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: sv_base_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command <<<
        set -euo pipefail
        mkdir files
        # Note that all filenames must be unique
        cat ~{file_list} | gsutil -m cp -I files/
        tar czf ~{prefix}.tar.gz -C files/ .
    >>>
    output {
        File out = "~{prefix}.tar.gz"
    }
}