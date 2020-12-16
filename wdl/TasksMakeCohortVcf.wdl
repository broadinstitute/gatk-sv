version 1.0

import "Structs.wdl"

# use zcat to concatenate compressed files
# -replaces "combine" task in some workflows
# -if filter_command is omitted, input files will be concatenated as
#  usual
# -if filter_command is passed, it must be a valid bash command,
#  accepting the resulting file via pipe on stdin, and outputing the
#  desired file on stdout
task ZcatCompressedFiles {
  input {
    Array[File] shards
    String? outfile_name
    String? filter_command
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  String output_file_name = select_first([outfile_name, "output.txt.gz"])
  Boolean do_filter = defined(filter_command) && select_first([filter_command]) != ""

  # when filtering/sorting/etc, memory usage will likely go up (much of the data will have to
  # be held in memory or disk while working, potentially in a form that takes up more space)
  Float input_size = size(shards, "GB")
  Float compression_factor = 5.0
  Float base_disk_gb = 5.0
  RuntimeAttr runtime_default = object {
    mem_gb: 2.0,
    disk_gb: ceil(base_disk_gb + input_size * if do_filter then 2.0 + compression_factor else 2.0),
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

  command {
    set -eu -o pipefail

    while read SHARD; do
      if [ -n "$SHARD" ]; then
        zcat "$SHARD"
      fi
    done < ~{write_lines(shards)} \
      ~{if do_filter then "| " + select_first([filter_command]) else ""} \
      | bgzip -c \
      > "~{outfile_name}"
  }

  output {
    File outfile=output_file_name
  }
}

# concatenate uncompressed files
# -replaces "combine" task in some workflows
# -if filter_command is omitted, input files will be concatenated as
#  usual
# -if filter_command is passed, it must be a valid bash command,
#  accepting the resulting file via pipe on stdin, and outputing the
#  desired file on stdout
task CatUncompressedFiles {
  input {
    Array[File] shards
    String? outfile_name
    String? filter_command
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  String output_file_name = select_first([outfile_name, "output.txt"])
  Boolean do_filter = defined(filter_command) && select_first([filter_command]) != ""

  # when filtering/sorting/etc, memory usage will likely go up (much of the data will have to
  # be held in memory or disk while working, potentially in a form that takes up more space)
  Float input_size = size(shards, "GB")
  Float base_disk_gb = 5.0
  RuntimeAttr runtime_default = object {
    mem_gb: 2.0,
    disk_gb: ceil(base_disk_gb + input_size * (if do_filter then 3.0 else 2.0)),
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
    while read SHARD; do
      if [ -n "$SHARD" ]; then
        cat "$SHARD"
      fi
    done < ~{write_lines(shards)} \
      ~{if do_filter then "| " + select_first([filter_command]) else ""} \
      > ~{output_file_name}
  >>>
  
  output {
    File outfile=output_file_name
  }
}

task SortVcf {
  input {
    File vcf
    String? outfile_prefix
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  String outfile_name = outfile_prefix + ".vcf.gz"

  RuntimeAttr runtime_default = object {
                                  mem_gb: 3.75,
                                  disk_gb: ceil(10.0 +  size(vcf, "GB") * 20),
                                  cpu_cores: 1,
                                  preemptible_tries: 3,
                                  max_retries: 1,
                                  boot_disk_gb: 10
                                }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

  Float runtime_mem_gb = select_first([runtime_override.mem_gb, runtime_default.mem_gb])

  command <<<
    set -euo pipefail
    mkdir temp
    bcftools sort \
        --temp-dir temp \
        --output-type z \
        --output-file ~{outfile_name} \
        ~{vcf}
    tabix ~{outfile_name}
  >>>

  output {
    File out = outfile_name
    File out_index = outfile_name + ".tbi"
  }
  runtime {
    memory: runtime_mem_gb + " GiB"
    disks: "local-disk " + select_first([runtime_override.disk_gb, runtime_default.disk_gb]) + " HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_base_mini_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
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
    String? outfile_prefix
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  String outfile_name = outfile_prefix + ".vcf.gz"
  String allow_overlaps_flag = if allow_overlaps then "--allow-overlaps" else ""
  String naive_flag = if naive then "--naive" else ""

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

  command <<<
    set -euo pipefail
    VCFS="~{write_lines(vcfs)}"
    if ~{!defined(vcfs_idx)}; then
      cat ${VCFS} | xargs -n1 tabix
    fi
    bcftools concat --no-version ~{allow_overlaps_flag} ~{naive_flag} --output-type z --file-list ${VCFS} --output "~{outfile_name}"
    if ~{generate_index}; then
      tabix "~{outfile_name}"
    else
      touch ~{outfile_name}.tbi
    fi
  >>>

  output {
    File concat_vcf = outfile_name
    File concat_vcf_idx = outfile_name + ".tbi"
  }
}

# Merge shards after VCF stats collection
task ConcatBeds {
  input {
    Array[File] shard_bed_files
    String prefix
    Boolean? index_output
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  Boolean call_tabix = select_first([index_output, true])
  String output_file="~{prefix}.bed.gz"

  # when filtering/sorting/etc, memory usage will likely go up (much of the data will have to
  # be held in memory or disk while working, potentially in a form that takes up more space)
  Float input_size = size(shard_bed_files, "GB")
  Float compression_factor = 5.0
  Float base_disk_gb = 5.0
  Float base_mem_gb = 2.0
  RuntimeAttr runtime_default = object {
    mem_gb: 2.0,
    disk_gb: ceil(base_disk_gb + input_size * (2.0 + compression_factor)),
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
    set -eux

    # note head -n1 stops reading early and sends SIGPIPE to zcat,
    # so setting pipefail here would result in early termination
    zcat ~{shard_bed_files[0]} | head -n1 > header.txt

    # no more early stopping
    set -o pipefail

    while read SPLIT; do
      zcat $SPLIT
    done < ~{write_lines(shard_bed_files)} \
      | (grep -Ev "^#" || printf "") \
      | sort -Vk1,1 -k2,2n -k3,3n \
      | cat header.txt - \
      | bgzip -c \
      > ~{output_file}

    if ~{call_tabix}; then
      tabix -p bed ~{output_file}
    else
      touch ~{output_file}.tbi
    fi
  >>>

  output {
    File merged_bed_file = output_file
    File merged_bed_idx = output_file + ".tbi"
  }
}


# Task to merge VID lists across shards
task FilesToTarredFolder {
  input {
    Array[File] in_files
    String? folder_name
    String? tarball_prefix
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  String tar_folder_name = select_first([folder_name, "merged"])
  String outfile_name = select_first([tarball_prefix, tar_folder_name]) + ".tar.gz"

  # Since the input files are often/always compressed themselves, assume compression factor for tarring is 1.0
  Float input_size = size(in_files, "GB")
  Float base_disk_gb = 5.0
  RuntimeAttr runtime_default = object {
    mem_gb: 2.0,
    disk_gb: ceil(base_disk_gb + input_size * 2.0),
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
    # Create final output directory
    mkdir "~{tar_folder_name}"

    while read VID_LIST; do
      mv "$VID_LIST" "~{tar_folder_name}"
    done < ~{write_lines(in_files)}

    # Compress final output directory
    tar -czvf "~{outfile_name}" "~{tar_folder_name}"
  >>>

  output {
    File tarball = outfile_name
  }
}


#Create input file for per-batch genotyping of predicted CPX CNV intervals
task PasteFiles {
  input {
    Array[String] input_strings
    Array[File] input_files
    String outfile_name
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  # when filtering/sorting/etc, memory usage will likely go up (much of the data will have to
  # be held in memory or disk while working, potentially in a form that takes up more space)
  Float input_size = size(input_files, "GB")
  Float base_disk_gb = 5.0
  RuntimeAttr runtime_default = object {
    mem_gb: 2.0,
    disk_gb: ceil(base_disk_gb + input_size * 2.0),
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

    paste ~{sep=' ' input_files} \
      > ~{outfile_name}
  >>>

  output {
    File outfile = outfile_name
  }
}

# Select a subset of vcf records by passing a filter command
# records_filter must be a bcftools expression
task FilterVcf {
  input {
    File vcf
    File vcf_index
    String outfile_prefix
    String records_filter
    Boolean? use_ssd
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  String outfile_name = outfile_prefix + ".vcf.gz"
  String disk_type = if (defined(use_ssd) && select_first([use_ssd])) then "SSD" else "HDD"

  # when filtering/sorting/etc, memory usage will likely go up (much of the data will have to
  # be held in memory or disk while working, potentially in a form that takes up more space)
  RuntimeAttr runtime_default = object {
    mem_gb: 2.0,
    disk_gb: ceil(10.0 + size(vcf, "GB") * 2),
    cpu_cores: 1,
    preemptible_tries: 3,
    max_retries: 1,
    boot_disk_gb: 10
  }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: select_first([runtime_override.mem_gb, runtime_default.mem_gb]) + " GB"
    disks: "local-disk " + select_first([runtime_override.disk_gb, runtime_default.disk_gb]) + " " + disk_type
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_base_mini_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -eu
    bcftools view --no-version --no-update -i '~{records_filter}' -O z -o ~{outfile_name} ~{vcf}
    tabix ~{outfile_name}
  >>>

  output {
    File filtered_vcf = outfile_name
    File filtered_vcf_idx = outfile_name + ".tbi"
  }
}

# Find intersection of Variant IDs from vid_list with those present in vcf, return as filtered_vid_list
task SubsetVariantList {
  input {
    File vid_list
    File vcf
    String outfile_name
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  # when filtering/sorting/etc, memory usage will likely go up (much of the data will have to
  # be held in memory or disk while working, potentially in a form that takes up more space)
  Float vid_list_size = size(vid_list, "GB")
  Float vcf_size = size(vcf, "GB")
  RuntimeAttr runtime_default = object {
    mem_gb: 2.0,
    disk_gb: ceil(10.0 + vcf_size + vid_list_size * 2.0),
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
    #Get list of variant IDs present in VCF
    zcat ~{vcf} | (grep -vE "^#" || printf "") | cut -f3 > valid_vids.list
    #Restrict input variant ID list to valid VIDs
    (fgrep -wf valid_vids.list ~{vid_list} || printf "") > "~{outfile_name}"
  >>>

  output {
    File filtered_vid_list = outfile_name
  }
}


# evenly split text file into even chunks
#   if shuffle_file is set to true, shuffle the file before splitting (default = false)
task SplitUncompressed {
  input {
    File whole_file
    Int lines_per_shard
    String? shard_prefix
    Boolean? shuffle_file
    Int? random_seed
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  String split_prefix=select_first([shard_prefix, "shard_"])
  Boolean do_shuffle=select_first([shuffle_file, false])
  Int random_seed_ = if defined(random_seed) then select_first([random_seed]) else 0

  # when filtering/sorting/etc, memory usage will likely go up (much of the data will have to
  # be held in memory or disk while working, potentially in a form that takes up more space)
  Float input_size = size(whole_file, "GB")
  RuntimeAttr runtime_default = object {
    mem_gb: 2.0,
    disk_gb: ceil(10.0 + input_size * 2.0),
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
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -eu -o pipefail

    function get_seeded_random() {
      openssl enc -aes-256-ctr -pass pass:"$1" -nosalt </dev/zero 2>/dev/null
    }
    # note for if ~do_shuffle is true: shuf is faster than sort --random-sort, but
    # sort --random-sort will predictably fit in memory, making it a better choice for VMs
    ~{if do_shuffle then
        "sort --random-sort --random-source=<(get_seeded_random ~{random_seed_}) -o ~{whole_file} ~{whole_file}"
      else
        ""
      }

    N_LINES=$(wc -l < ~{whole_file})
    N_CHUNKS=$((N_LINES / ~{lines_per_shard}))
    if [ "$N_CHUNKS" -eq "0" ]; then N_CHUNKS=1; fi
    N_DIGITS=${#N_CHUNKS}

    split -d -a $N_DIGITS -n l/$N_CHUNKS \
      --numeric-suffixes=$(printf "%0${N_DIGITS}d" 1) \
      ~{whole_file} \
      ~{shard_prefix}

    # remove whole file if its name starts with split_prefix, to prevent including in glob
    if [[ "~{whole_file}" =~ ^"~{split_prefix}".* ]]; then
      rm -f "~{whole_file}"
    fi
  >>>

  output {
     Array[File] shards=glob("~{shard_prefix}*")
  }
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

#Update either SR bothside_pass or background_fail files
task UpdateSrList {
  input {
    File vcf
    File original_list
    String outfile
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  # when filtering/sorting/etc, memory usage will likely go up (much of the data will have to
  # be held in memory or disk while working, potentially in a form that takes up more space)
  Float input_size = size([vcf, original_list], "GiB")
  Float compression_factor = 5.0
  Float base_disk_gb = 10.0
  Float base_mem_gb = 2.0
  RuntimeAttr runtime_default = object {
    mem_gb: 3.75,
    disk_gb: ceil(base_disk_gb + input_size * (2.0 + 2.0 * compression_factor)),
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
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euxo pipefail

    ##append new ids to original list##
    svtk vcf2bed ~{vcf} int.bed -i MEMBERS

    ##remove header and match id one per line##
    awk '{if (NR>1) print $4 "\t" $NF}' int.bed \
      | awk -F'[,\t]' '{for(i=2; i<=NF; ++i) print $i "\t" $1 }' \
      | sort -k1,1\
      > newidlist.txt

    join -j 1 -t $'\t' <(awk '{print $NF "\t" $0}' ~{original_list} | sort -k1,1) newidlist.txt \
      | cut -f2-  \
      > ~{outfile}
  >>>

  output {
    File updated_list = outfile
  }
}

  
task ShardVids {
  input {
    File clustered_vcf
    String prefix
    Int records_per_shard
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size(clustered_vcf, "GiB")
  Float base_disk_gb = 10.0
  Float input_disk_scale = 1.0
  RuntimeAttr runtime_default = object {
                                  mem_gb: 2.0,
                                  disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
                                  cpu_cores: 1,
                                  preemptible_tries: 3,
                                  max_retries: 1,
                                  boot_disk_gb: 10
                                }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GiB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euo pipefail
    
    python3 <<CODE
    import sys
    import os
    import pysam

    vcf = pysam.VariantFile("~{clustered_vcf}")

    current_cluster = None
    current_cluster_vids = []
    current_shard = 0
    current_shard_size = 0
    shard_path_format = "~{prefix}.vids.shard_{}.list"
    shard_path = shard_path_format.format(current_shard)
    fout = open(shard_path, 'w')
    if fout is None:
      raise IOError("Could not open '{}'".format(shard_path))
      sys.exit(1)

    for record in vcf.fetch():
      cluster_id = record.info['CLUSTER']
      if cluster_id == current_cluster:
        current_cluster_vids.append(record.id)
      else:
        for vid in current_cluster_vids:
          fout.write(vid + '\n')
        current_shard_size += len(current_cluster_vids)
        if current_shard_size >= ~{records_per_shard}:
          current_shard += 1
          current_shard_size = 0
          fout.close()
          shard_path = shard_path_format.format(current_shard)
          fout = open(shard_path, 'w')
          if fout is None:
            raise IOError("Could not open '{}'".format(shard_path))
            sys.exit(1)
        current_cluster_vids = [record.id]
        current_cluster = cluster_id

    # Write last cluster
    for vid in current_cluster_vids:
      fout.write(vid + '\n')
    current_shard_size += len(current_cluster_vids)
    fout.close()

    # Delete trailing empty shard
    if current_shard > 0 and current_shard_size == 0:
      os.remove(shard_path)
    CODE
  >>>

  output {
    Array[File] out = glob("~{prefix}.vids.shard_*.list")
  }
}