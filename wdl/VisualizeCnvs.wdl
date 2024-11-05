version 1.0

import "Structs.wdl"

# Plots CNV depth profiles across batches

workflow VisualizeCnvs {
  input{
    # Note vcf will be faster
    File vcf_or_bed  # bed columns: chrom,start,end,name,svtype,samples
    String plot_prefix
    File ped_file
    Int min_size
    Int? variants_per_shard
    String rdtest_flags

    File sample_table # TSV with sample_id, sample_set_id
    Array[File] median_files
    Array[File] rd_files
    Array[File] rd_file_indicies

    String linux_docker
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_format_vcf_or_bed
    RuntimeAttr? runtime_attr_shard_variants
    RuntimeAttr? runtime_attr_subset_rd_matrices
    RuntimeAttr? runtime_attr_rdtest
    RuntimeAttr? runtime_attr_merge_plot_tars
  }

  call FormatVCFOrBED {
    input:
      vcf_or_bed = vcf_or_bed,
      min_size = min_size,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_format_vcf_or_bed
  }

  call ShardVariants {
    input:
      variants = FormatVCFOrBED.variants,
      variants_per_shard = variants_per_shard,
      linux_docker = linux_docker,
      runtime_attr_override = runtime_attr_shard_variants
  }

  scatter (i in range(length(rd_files))) {
    call SubsetRDMatrices {
      input:
        rd_file = rd_files[i],
        rd_file_index = rd_file_indicies[i],
        intervals = FormatVCFOrBED.intervals,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_subset_rd_matrices
    }
  }

  scatter (shard in ShardVariants.shards) {
    call RdTestPlot {
      input:
        variants = shard,
        rd_files = SubsetRDMatrices.rd_subset,
        rd_file_indicies = SubsetRDMatrices.rd_subset_index,
        median_files = median_files,
        sample_table = sample_table,
        ped_file = ped_file,
        plot_prefix = plot_prefix,
        sv_pipeline_docker = sv_pipeline_docker,
        flags = rdtest_flags,
        runtime_attr_override = runtime_attr_rdtest
    }
  }

  call MergePlotTars {
    input:
      plot_tars = RdTestPlot.plots,
      plot_prefix = plot_prefix,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_merge_plot_tars
  }

  output {
    File rdtest_plots = MergePlotTars.plots
  }
}

task FormatVCFOrBED {
  input {
    File vcf_or_bed
    Int min_size
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 4,
    disk_gb: ceil(size(vcf_or_bed, "GB") * 2) + 16,
    boot_disk_gb: 16,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    cat2() {
      if [[ "$1" == *.gz ]]; then
        zcat "$1"
      else
        cat "$1"
      fi
    }
    if [[ '~{vcf_or_bed}' == *.vcf.gz ]]; then
      bcftools view --include '(SVTYPE == "DEL" || SVTYPE == "DUP") && SVLEN >= ~{min_size}' '~{vcf_or_bed}' \
        | svtk vcf2bed stdin raw.bed
      awk -F'\t' -v OFS='\t' '!/^#/{print $1,$2,$3,$4,$6,$5}' raw.bed \
        | LC_ALL=C sort -k1,1 -k2,2n > cnvs.bed
      rm raw.bed
    elif [[ '~{vcf_or_bed}' == *.bed || '~{vcf_or_bed}' == *.bed.gz ]]; then
      cat2 '~{vcf_or_bed}' \
        | awk -F'\t' -v OFS='\t' '
          !/^#/ && $3 - $2 >= ~{min_size} && ($5 == "DEL" || $5 == "DUP") {print $1,$2,$3,$4,$6,$5}' \
        | LC_ALL=C sort -k1,1 -k2,2n > cnvs.bed
    else
      echo "Invalid extension for input calls. Must be .vcf.gz, .bed.gz, or .bed"
      exit 1
    fi

    bedtools merge -i cnvs.bed | cut -f1-3 > merged_intervals.bed
    gzip cnvs.bed
  >>>

  output {
    File variants = 'cnvs.bed.gz'
    File intervals = 'merged_intervals.bed'
  }
}

task ShardVariants {
  input {
    File variants
    Int variants_per_shard = 40
    String linux_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 1,
    disk_gb: ceil(size(variants, "GB") * 5) + 16,
    boot_disk_gb: 16,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: linux_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    variants='~{variants}'
    gzip -d "${variants}"
    split -l ~{variants_per_shard} "${variants%.gz}" cnvs_
  >>>

  output {
    Array[File] shards = glob("cnvs_*")
  }
}

task SubsetRDMatrices {
  input {
    File rd_file
    File rd_file_index
    File intervals
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 1,
    disk_gb: ceil(size(rd_file, "GB") * 2 + size(intervals, "GB")) + 16,
    boot_disk_gb: 16,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

  String rd_file_bn = basename(rd_file)

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    tabix --print-header --regions '~{intervals}' '~{rd_file}' \
      | bgzip > '~{rd_file_bn}'
    tabix --preset bed '~{rd_file_bn}'
  >>>

  output {
    File rd_subset = rd_file_bn
    File rd_subset_index = rd_file_bn + ".tbi"
  }
}

task RdTestPlot {
  input {
    File variants
    Array[File] rd_files
    Array[File] rd_file_indicies
    Array[File] median_files
    File sample_table
    File ped_file
    String plot_prefix
    String sv_pipeline_docker
    String flags
    RuntimeAttr? runtime_attr_override
  }

  Int variant_count = length(read_lines(variants))
  Float input_size = size(rd_files, "GB")
    + size(rd_file_indicies, "GB")
    + size(median_files, "GB")
    + size(sample_table, "GB")
    + size(ped_file, "GB")
    + size(variants, "GB")

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 8,
    disk_gb: ceil(input_size + variant_count * 0.01) + 16,
    boot_disk_gb: 16,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    mkdir median_files
    xargs -I '{}' mv '{}' median_files < '~{write_lines(median_files)}'

    mkdir rd_files
    xargs -I '{}' mv '{}' rd_files < '~{write_lines(rd_files)}'
    xargs -I '{}' mv '{}' rd_files < '~{write_lines(rd_file_indicies)}'

    mkdir rd_links
    while IFS=$'\t' read -r chr start end vid samples svtype; do
      : > variant.bed
      : > batch_ids.list
      : > merged_medians.bed

      printf '%s\t%s\t%s\t%s\t%s\t%s\n' "${chr}" "${start}" "${end}" "${vid}" "${samples}" "${svtype}" > variant.bed
      tr ',' '\n' <<< "${samples}" \
        | awk -F'\t' 'NR==FNR{a[$1]} NR>FNR && ($1 in a){print $2}' - '~{sample_table}' \
        | sort -u > batch_ids.list
      awk '{
        print ENVIRON["PWD"] "/rd_files/" $0 ".RD.txt.gz"
        print "rd_links/" $0 ".RD.txt.gz"
        print ENVIRON["PWD"] "/rd_files/" $0 ".RD.txt.gz.tbi"
        print "rd_links/" $0 ".RD.txt.gz.tbi"
      }' batch_ids.list \
        | xargs -L 2 ln -s
      awk '{print "median_files/" $0 "_medianCov.transposed.bed"}' batch_ids.list \
        | xargs paste > merged_medians.bed

      Rscript /opt/RdTest/RdTest.R \
        -b variant.bed \
        -n '~{plot_prefix}' \
        -c rd_links \
        -m merged_medians.bed \
        -p TRUE \
        ~{flags}

      find rd_links -type l -delete
    done < '~{variants}'

    mkdir ~{plot_prefix}_rd_plots
    mv *.jpg ~{plot_prefix}_rd_plots
    tar -czvf ~{plot_prefix}_rd_plots.tar.gz ~{plot_prefix}_rd_plots/
  >>>

  output {
    File plots = "${plot_prefix}_rd_plots.tar.gz"
  }
}

task MergePlotTars {
  input {
    Array[File] plot_tars
    String plot_prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 1,
    disk_gb: ceil(size(plot_tars, "GB") * 3) + 16,
    boot_disk_gb: 16,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    while read -r pt; do
      tar -xzf "${pt}"
    done < '~{write_lines(plot_tars)}'

    tar -czf ~{plot_prefix}_rd_plots.tar.gz ~{plot_prefix}_rd_plots
  >>>

  output {
    File plots = "${plot_prefix}_rd_plots.tar.gz"
  }
}
