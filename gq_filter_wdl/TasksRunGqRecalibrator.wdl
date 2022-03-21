version development

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

# Split a vcf into even chunks
#   -if contig is specified, only take records from that contig
#   -if n_shards is specified, use it as a maximum value for number of shards,
#    but use lower value if necessary to ensure that min_vars_per_shard is obeyed
#   -if vcf_idx is present, use it if it will speed up processing, otherwise ignore
task SplitVcf {
  input {
    File vcf
    File? vcf_idx
    String prefix
    Int? n_shards
    Int min_vars_per_shard
    String? contig
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  # Worst-case scenario for disk size: simultaneously storing uncompressed vcf and records vcf
  # Note that the header is multiplied by the number of shards, but it's so small as to be a rounding error unless
  # thousands of shards are produced (and even then, small compared to base_disk_gb)
  Float input_size = size(vcf, "GB")
  RuntimeAttr runtime_default = object {
    mem_gb: 2.0,
    disk_gb: ceil(10 + input_size * 30),
    cpu_cores: 1,
    preemptible_tries: 3,
    max_retries: 1,
    boot_disk_gb: 10
  }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: select_first([runtime_override.mem_gb, runtime_default.mem_gb]) + " GB"
    disks: "local-disk " + select_first([runtime_override.disk_gb, runtime_default.disk_gb]) + " HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_base_mini_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -eu -o pipefail

    # Uncompress vcf. If index is present and contig is specified, use tabix to extract desired contig
    ~{if defined(vcf_idx) && defined(contig) then "tabix -h ~{vcf} ~{contig}" else "zcat ~{vcf}"} \
      > uncompressed.vcf

    # Extract vcf header:
    #  search for first line not starting with '#', stop immediately,
    #  take everything up to that point, then remove last line.
    ONLY_HEADER=false
    grep -B9999999999 -m1 -Ev "^#" uncompressed.vcf | sed '$ d' > header.vcf \
      || ONLY_HEADER=true

    if $ONLY_HEADER; then
      # No records so nothing to split. Just move the original file to be named like a chunk
      bgzip -c uncompressed.vcf > "~{prefix}1.vcf.gz"
    else
      N_HEADER=$(wc -l < header.vcf)

      # Select the desired records.
      #  if contig is specified, then select only records from appropriate contig,
      #  otherwise take all lines after the end of the header
      tail -n+$((N_HEADER+1)) uncompressed.vcf~{if defined(contig) && !defined(vcf_idx) then
                                                ' | { grep -w "^' + contig + '" || true; }'
                                                else ''} \
        > records.vcf
      N_RECORDS=$(wc -l < records.vcf)
      rm uncompressed.vcf

      # specifying split -n instead of split -l produces more even splits
      N_CHUNKS=$((N_RECORDS / ~{min_vars_per_shard}))
      if [ $N_CHUNKS -gt 1 ]; then
        rm "~{vcf}"
        MAX_CHUNKS=~{if defined(n_shards) then n_shards else 0}
        if [ $MAX_CHUNKS -gt 0 ] && [ $MAX_CHUNKS -lt $N_CHUNKS ]; then
          N_CHUNKS=$MAX_CHUNKS
        fi

        # figure out how many digits we need in split suffixes
        N_DIGITS=${#N_CHUNKS}

        # produce headerless uncompressed files with VCF records, each
        # having an numeric suffix
        split -n l/$N_CHUNKS -a $N_DIGITS \
          --numeric-suffixes=$(printf "%0${N_DIGITS}d" 1) \
          records.vcf "~{prefix}"

        rm records.vcf

        # loop over split records
        find . -name "~{prefix}*" | while read VCF_RECORD; do
          # add header, compress, add .vcf.gz extension
          cat header.vcf $VCF_RECORD \
            | bgzip -c \
            > $VCF_RECORD.vcf.gz
          rm $VCF_RECORD
        done
      else
        # just one chunk, so use full records.vcf. add header, compress, and name like a chunk
        # use records.vcf in case vcf_idx not defined and uncompressed.vcf not result of tabix call
        cat header.vcf records.vcf | bgzip -c > "~{prefix}1.vcf.gz"
      fi
    fi
  >>>

  output {
    Array[File] vcf_shards = glob("~{prefix}*.vcf.gz")
  }
}

# Combine multiple sorted VCFs
task ConcatVcfs {
  input {
    Array[File] vcfs
    Array[File]? vcfs_idx
    Boolean allow_overlaps = false
    Boolean naive = false
    Boolean generate_index = true
    Boolean sites_only = false
    Boolean sort_vcf_list = false
    String? outfile_prefix
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  # when filtering/sorting/etc, memory usage will likely go up (much of the data will have to
  # be held in memory or disk while working, potentially in a form that takes up more space)
  RuntimeAttr runtime_default = object {
    mem_gb: 3.75,
    disk_gb: ceil(10 + size(vcfs, "GB") * 2),
    cpu_cores: 1,
    preemptible_tries: 3,
    max_retries: 1,
    boot_disk_gb: 10
  }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: select_first([runtime_override.mem_gb, runtime_default.mem_gb]) + " GB"
    disks: "local-disk " + select_first([runtime_override.disk_gb, runtime_default.disk_gb]) + " SSD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_base_mini_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  String outfile_name = outfile_prefix + ".vcf.gz"
  String allow_overlaps_flag = if allow_overlaps then "--allow-overlaps" else ""
  String naive_flag = if naive then "--naive" else ""
  String concat_output_type = if (sites_only) then "v" else "z"
  String sites_only_command = if (sites_only) then "| bcftools view --no-version -G -Oz" else ""
  String generate_index_command = if (generate_index) then "tabix ~{outfile_name}" else "touch ~{outfile_name}.tbi"

  command <<<
    set -euo pipefail
    VCFS="~{write_lines(vcfs)}"
    if ~{sort_vcf_list}; then
      cat $VCFS | awk -F '/' '{print $NF"\t"$0}' | sort -k1,1V | awk '{print $2}' > vcfs.list
    else
      cp $VCFS vcfs.list
    fi
    bcftools concat --no-version ~{allow_overlaps_flag} ~{naive_flag} -O~{concat_output_type} --file-list vcfs.list \
      ~{sites_only_command} \
      > ~{outfile_name}
    ~{generate_index_command}
  >>>

  output {
    File concat_vcf = outfile_name
    File concat_vcf_idx = outfile_name + ".tbi"
  }
}
