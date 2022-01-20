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
  Boolean do_filter = defined(filter_command) && filter_command != ""

  # when filtering/sorting/etc, memory usage will likely go up (much of the data will have to
  # be held in memory or disk while working, potentially in a form that takes up more space)
  Float input_size = size(shards, "GB")
  Float compression_factor = 5.0
  Float base_disk_gb = 5.0
  Float base_mem_gb = 2.0
  RuntimeAttr runtime_default = object {
    mem_gb: base_mem_gb + if do_filter then compression_factor * input_size else 0.0,
    disk_gb: ceil(base_disk_gb + input_size * if do_filter then 2.0 + compression_factor else 2.0),
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
    File outfile=outfile_name
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
  Boolean do_filter = defined(filter_command) && filter_command != ""

  # when filtering/sorting/etc, memory usage will likely go up (much of the data will have to
  # be held in memory or disk while working, potentially in a form that takes up more space)
  Float input_size = size(shards, "GB")
  Float base_mem_gb = 2.0
  Float base_disk_gb = 5.0
  RuntimeAttr runtime_default = object {
    mem_gb: base_mem_gb + (if do_filter then input_size else 0.0),
    disk_gb: ceil(base_disk_gb + input_size * (if do_filter then 3.0 else 2.0)),
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

# Combine multiple sorted VCFs
task ConcatVcfs {
  input {
    Array[File] vcfs
    Array[File]? vcfs_idx
    Boolean merge_sort = false
    String? outfile_prefix
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  String outfile_name = outfile_prefix + ".vcf.gz"
  String merge_flag = if merge_sort then "--allow-overlaps" else ""

  # when filtering/sorting/etc, memory usage will likely go up (much of the data will have to
  # be held in memory or disk while working, potentially in a form that takes up more space)
  Float input_size = size(vcfs, "GB")
  Float compression_factor = 5.0
  Float base_disk_gb = 5.0
  Float base_mem_gb = 2.0
  RuntimeAttr runtime_default = object {
    mem_gb: base_mem_gb + compression_factor * input_size,
    disk_gb: ceil(base_disk_gb + input_size * (2.0 + compression_factor)),
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
    VCFS="~{write_lines(vcfs)}"
    if ~{!defined(vcfs_idx)}; then
      cat ${VCFS} | xargs -n1 tabix
    fi
    bcftools concat -a ~{merge_flag} --output-type z --file-list ${VCFS} --output ~{outfile_name}
    tabix -p vcf ~{outfile_name}
  >>>

  output {
    File concat_vcf = "~{outfile_name}"
    File concat_vcf_idx = "~{outfile_name}.tbi"
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

  Boolean call_tabix = select_first([index_output, false])
  String output_file="~{prefix}.bed.gz"

  # when filtering/sorting/etc, memory usage will likely go up (much of the data will have to
  # be held in memory or disk while working, potentially in a form that takes up more space)
  Float input_size = size(shard_bed_files, "GB")
  Float compression_factor = 5.0
  Float base_disk_gb = 5.0
  Float base_mem_gb = 2.0
  RuntimeAttr runtime_default = object {
    mem_gb: base_mem_gb + compression_factor * input_size,
    disk_gb: ceil(base_disk_gb + input_size * (2.0 + compression_factor)),
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

task ConcatGTGQs {
    input {
        Array[File] files
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_override
    }

    String sample = basename(files[0], ".gtgq.gz")

    # when filtering/sorting/etc, memory usage will likely go up (much of the data will have to
    # be held in memory or disk while working, potentially in a form that takes up more space)
    Float input_size = size(files, "GB")
    Float compression_factor = 5.0
    Float base_disk_gb = 5.0
    Float base_mem_gb = 2.0
    RuntimeAttr runtime_default = object {
      mem_gb: base_mem_gb + compression_factor * input_size,
      disk_gb: ceil(base_disk_gb + input_size * (2.0 + compression_factor)),
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
        set -eu
        zcat ~{files[0]} | head -1 > header.txt

        set -euo pipefail

        while read file; do
          zcat $file | tail -n+2
        done < ~{write_lines(files)} \
          | cat header.txt - \
          | bgzip \
        > ~{sample}.all_chr.gtgq.gz

    >>>
    output {
        File merged_gtgq_file = "~{sample}.all_chr.gtgq.gz"
    }
}

# Merge shards after VaPoR
task ConcatVaPoR {
  input {
    Array[File] shard_plots
    String prefix
    Boolean? index_output
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  Boolean call_tabix = select_first([index_output, true])
  String output_file="~{prefix}.bed.gz"

  # when filtering/sorting/etc, memory usage will likely go up (much of the data will have to
  # be held in memory or disk while working, potentially in a form that takes up more space)
  Float compression_factor = 5.0
  Float base_disk_gb = 5.0
  Float base_mem_gb = 2.0
  RuntimeAttr runtime_default = object {
    mem_gb: base_mem_gb,
    disk_gb: ceil(base_disk_gb  * (2.0 + compression_factor)),
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
    set -eu

    # note head -n1 stops reading early and sends SIGPIPE to zcat,
    # so setting pipefail here would result in early termination

    # no more early stopping
    set -o pipefail

    mkdir ~{prefix}.plots
    mkdir ~{prefix}.tmp_plots

    while read SPLIT; do
      tar zxvf $SPLIT -C ~{prefix}.tmp_plots/
      find -type f -name '~{prefix}.tmp_plots/*/*DEL.png' | xargs -n1 -I{} mv {} ~{prefix}.plots/
      find -type f -name '~{prefix}.tmp_plots/*/*INS.png' | xargs -n1 -I{} mv {} ~{prefix}.plots/
      find -type f -name '~{prefix}.tmp_plots/*/*.png' | xargs -n1 -I{} mv {} ~{prefix}.plots/

      #ls ~{prefix}.tmp_plots/*/*DEL.png | xargs -n1 -I{} mv {} ~{prefix}.plots/
      #ls ~{prefix}.tmp_plots/*/*INS.png | xargs -n1 -I{} mv {} ~{prefix}.plots/
      #ls ~{prefix}.tmp_plots/*/*.png | xargs -n1 -I{} mv {} ~{prefix}.plots/
      rm -r ~{prefix}.tmp_plots/*
    done < ~{write_lines(shard_plots)}

    tar -czf ~{prefix}.plots.tar.gz ~{prefix}.plots/
  >>>

  output {
    File merged_bed_plot = "~{prefix}.plots.tar.gz"
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
  Float base_mem_gb = 2.0
  RuntimeAttr runtime_default = object {
    mem_gb: base_mem_gb,
    disk_gb: ceil(base_disk_gb + input_size * 2.0),
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
  Float base_mem_gb = 2.0
  RuntimeAttr runtime_default = object {
    mem_gb: base_mem_gb,
    disk_gb: ceil(base_disk_gb + input_size * 2.0),
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
    set -eu -o pipefail

    paste ~{sep=' ' input_files} \
      > ~{outfile_name}
  >>>

  output {
    File outfile = outfile_name
  }
  }

# Select a subset of vcf records by passing a bash filter command
# records_filter must be a bash command accepting vcf records passed via
# pipe, and outputing the desired records to stdout
task FilterVcf {
  input {
    File vcf
    String outfile_prefix
    String records_filter
    Boolean? index_output
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  String outfile_name = outfile_prefix + ".vcf.gz"
  Boolean call_tabix = select_first([index_output, true])

  # when filtering/sorting/etc, memory usage will likely go up (much of the data will have to
  # be held in memory or disk while working, potentially in a form that takes up more space)
  Float input_size = size(vcf, "GB")
  Float compression_factor = 5.0
  Float base_disk_gb = 5.0
  Float base_mem_gb = 2.0
  RuntimeAttr runtime_default = object {
    mem_gb: base_mem_gb + compression_factor * input_size,
    disk_gb: ceil(base_disk_gb + input_size * 2.0 * (1 + compression_factor)),
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
    set -eu -o pipefail

    # uncompress vcf
    zcat ~{vcf} > uncompressed.vcf

    # Extract vcf header:
    #  search for first line not starting with '#', stop immediately,
    #  take everything up to that point, then remove last line.
    ONLY_HEADER=false
    grep -B9999999999 -m1 -Ev "^#" uncompressed.vcf | sed '$ d' > header.vcf \
      || ONLY_HEADER=true

    if $ONLY_HEADER; then
      # no records were found, so filter is trivial, just use original vcf
      mv ~{vcf} ~{outfile_name}
    else
      N_HEADER=$(wc -l < header.vcf)

      # Put filter inside subshell so that there is no pipefail if there are no matches
      # NOTE: this is dangerous in the event that the filter is buggy
      tail -n+$((N_HEADER+1)) uncompressed.vcf \
        | { ~{records_filter} || true; }\
        | cat header.vcf - \
        | vcf-sort \
        | bgzip -c \
        > "~{outfile_name}"
    fi

    if ~{call_tabix}; then
      tabix -p vcf -f "~{outfile_name}"
    else
      touch "~{outfile_name}.tbi"
    fi
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
  Float compression_factor = 5.0
  Float cut_factor = 2.0
  Float base_disk_gb = 5.0
  Float base_mem_gb = 2.0
  RuntimeAttr runtime_default = object {
    mem_gb: base_mem_gb + compression_factor / cut_factor * vcf_size,
    disk_gb: ceil(base_disk_gb + vid_list_size * 2.0 + vcf_size * (1 + compression_factor / cut_factor)),
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
  Float base_disk_gb = 5.0
  Float base_mem_gb = 2.0
  RuntimeAttr runtime_default = object {
    mem_gb: base_mem_gb,
    disk_gb: ceil(base_disk_gb + input_size * 2.0),
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


#localize a specific contig of a bam/cram file
task LocalizeCram{
  input{
    String contig
    File ref_fasta
    File ref_fai
    File ref_dict
    String? bam_or_cram_file
    String? bam_or_cram_index
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 15, 
    disk_gb: 40,
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  Float mem_gb = select_first([runtime_attr.mem_gb, default_attr.mem_gb])
  Int java_mem_mb = ceil(mem_gb * 1000 * 0.8)

  output{
    File local_bam = "~{contig}.bam"
    File local_bai = "~{contig}.bam.bai"
  }

  command <<<
    set -Eeuo pipefail
    
    java -Xmx~{java_mem_mb}M -jar ${GATK_JAR}  PrintReads \
    -I ~{bam_or_cram_file} \
    -L ~{contig} \
    -O ~{contig}.bam \
    -R ~{ref_fasta}

    samtools index ~{contig}.bam
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

task LocalizeCramRequestPay{
  input{
    String contig
    File ref_fasta
    File ref_fai
    File ref_dict
    String? project_id
    String? bam_or_cram_file
    String? bam_or_cram_index
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 3.75, 
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  Float mem_gb = select_first([runtime_attr.mem_gb, default_attr.mem_gb])
  Int java_mem_mb = ceil(mem_gb * 1000 * 0.8)

  output{
    File local_bam = "~{contig}.bam"
    File local_bai = "~{contig}.bam.bai"
  }

  command <<<
    set -Eeuo pipefail
    
    java -Xmx~{java_mem_mb}M -jar ${GATK_JAR}  PrintReads \
    -I ~{bam_or_cram_file} \
    -L ~{contig} \
    -O ~{contig}.bam \
    -R ~{ref_fasta} \
    --gcs-project-for-requester-pays ~{project_id}

    samtools index ~{contig}.bam
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

#extract specific contig from vcf
task SplitBed{
  input{
    String contig
    File? bed_file
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 3.75, 
    disk_gb: 15,
    boot_disk_gb: 15,
    preemptible_tries: 0,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output{
    File contig_bed = "~{contig}.bed.gz"
  }

  command <<<
    if [[ ~{bed_file} == *.gz ]] ;  then
      zcat ~{bed_file} | awk '{if ($1=="~{contig}") print}' | bgzip > ~{contig}.bed.gz
    else
      awk '{if ($1=="~{contig}") print}' ~{bed_file} | bgzip > ~{contig}.bed.gz
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

task SplitVcf{
  input{
    String contig
    File? vcf_file
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 3.75, 
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output{
    File contig_vcf = "~{contig}.vcf.gz"
    File contig_vcf_index = "~{contig}.vcf.gz.tbi"
  }

  command <<<
    if [[ ~{vcf_file} == *.gz ]] ;  then
      tabix -f -p vcf ~{vcf_file}
      tabix -h ~{vcf_file} ~{contig} | bgzip > ~{contig}.vcf.gz
      tabix -p vcf ~{contig}.vcf.gz
    else
      bgzip ~{vcf_file}
      tabix -f -p vcf ~{vcf_file}.gz
      tabix -h ~{vcf_file}.gz ~{contig} | bgzip > ~{contig}.vcf.gz
      tabix -p vcf ~{contig}.vcf.gz
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

task split_per_sample_vcf {
  input {
    File vcf
    String sample
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 7.5, 
    disk_gb: 20,
    boot_disk_gb: 20,
    preemptible_tries: 1,
    max_retries: 3
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail

    bcftools view -Oz -s ~{sample} ~{vcf} -o ~{sample}.vcf.gz

  >>>

  output {
    File vcf_file = "~{sample}.vcf.gz"
  }
  
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

task split_per_sample_gtgq {
  input {
    File clean_vcf
    File clean_vcf_idx
    Array[String] samples
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Int num_samples = length(samples)

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 7.5, 
    disk_gb: 20,
    boot_disk_gb: 20,
    preemptible_tries: 1,
    max_retries: 3
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail

    mkdir out
    #ulimit -S -n 2097152
    bcftools +split ~{clean_vcf} -Oz -o out/

    sample_ids=(~{sep=" " samples})
    for (( i=0; i<~{num_samples}; i++ ));
    do
      sample_id=${sample_ids[$i]}
      vcf=${sample_id}.vcf.gz
      sample_no=`printf %08d $i`
      echo -e "SVID\tGT\tCN\tCNQ\tEV\tGQ\tPE_GQ\tPE_GT\tRD_CN\tRD_GQ\tSR_GQ\tSR_GT" > out/gtgq_${sample_no}.${sample_id}.~{prefix}.gtgq
      zcat out/${vcf} | grep -v "#" | cut -f3,10 | sed -e "s/:/\t/g" | awk '{if ($2!="0/0") print}' >> out/gtgq_${sample_no}.${sample_id}.~{prefix}.gtgq
      bgzip out/gtgq_${sample_no}.${sample_id}.~{prefix}.gtgq
    done

  >>>

  output {
    Array[File] gtgq_file = glob("out/gtgq_*.gtgq.gz")
  }
  
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

task split_per_sample_gtgq_V2 {
  input {
    File clean_vcf
    File clean_vcf_idx
    Array[String] samples
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Int num_samples = length(samples)

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 7.5, 
    disk_gb: 20,
    boot_disk_gb: 20,
    preemptible_tries: 1,
    max_retries: 3
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail

    #mkdir out
    #ulimit -S -n 2097152
    #bcftools +split ~{clean_vcf} -Oz -o out/

    sample_ids=(~{sep=" " samples})
    for (( i=0; i<~{num_samples}; i++ ));
    do
      sample_id=${sample_ids[$i]}
      vcf=${sample_id}.vcf.gz
      sample_no=`printf %08d $i`
      echo "SVID\tGT\tCN\tCNQ\tEV\tGQ\tPE_GQ\tPE_GT\tRD_CN\tRD_GQ\tSR_GQ\tSR_GT" > gtgq_${sample_no}.${sample_id}.~{prefix}.gtgq
      bcftools view -s ${sample_id} | grep -v "#" | cut -f3,10 | sed -e "s/:/\t/g" >> gtgq_${sample_no}.${sample_id}.~{prefix}.gtgq
      bgzip gtgq_${sample_no}.${sample_id}.~{prefix}.gtgq
    done

  >>>

  output {
    Array[File] gtgq_file = glob("*.gtgq.gz")
  }
  
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

task split_per_sample_bed{
    input{
        File bed
        String sample
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
        }

    Float input_size = size(bed, 'GB')
    Float base_disk_gb = 5.0

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 1, 
        disk_gb: ceil(base_disk_gb + input_size*2),
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])


    output {
        File bed_file = "~{sample}.bed.gz"
    }
    command <<<

        set -Eeuo pipefail
        
        cat <(zcat ~{bed} | head -1) \
        <(zcat ~{bed} | grep ~{sample}) \
        | cut -f1-5,7- | bgzip > ~{sample}.bed.gz

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

task vcf2bed{
    input{
        String prefix
        String sample
        File vcf
        File? vcf_index
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 6,
        disk_gb: 25,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    String filename = basename(vcf, ".vcf.gz")

    output {
        File bed = "~{prefix}.bed.gz"
    }

    command <<<

        set -Eeuo pipefail
        
        gsutil cp ~{vcf} ./tmp.vcf.gz
        tabix -p vcf ./tmp.vcf.gz
        svtk vcf2bed -i SVTYPE -i SVLEN tmp.vcf.gz ~{prefix}.tmp

        head -1 ~{prefix}.tmp > ~{prefix}.bed
        awk '{if ($6=="~{sample}") print}' ~{prefix}.tmp >> ~{prefix}.bed
        bgzip ~{prefix}.bed
        
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

task vcf2bedSVSites{
    input{
        String prefix
        File vcf
        File? vcf_index
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 6,
        disk_gb: 25,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    String filename = basename(vcf, ".vcf.gz")

    output {
        File bed = "~{prefix}.bed"
    }

    command <<<

        set -Eeuo pipefail
        
        gsutil cp ~{vcf} ./tmp.vcf.gz
        tabix -p vcf ./tmp.vcf.gz
        svtk vcf2bed -i SVTYPE -i SVLEN -i EVIDENCE -i ALGORITHMS -i AF -i AC -i AN tmp.vcf.gz ~{prefix}.bed
        
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

task Bed2QueryAndRef{
    input{
        File bed
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 3, 
        disk_gb: 10,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output{
        File query = "${filebase}.query.gz"
        File ref = "${filebase}.ref.gz"
    }

    String filebase=basename(bed,".bed.gz")
    command <<<
        echo "#chroms tart end name SVTYPE SVLEN" | sed -e 's/ /\t/g' > ~{filebase}.query
        echo "#chrom start end VID svtype length AF samples" | sed -e 's/ /\t/g' > ~{filebase}.ref

        zcat ~{bed} | cut -f1-4,6,7 | grep -v "#" >> ~{filebase}.query
        zcat ~{bed} | cut -f1-4,6,7 | sed -e "s/$/\t0\t~{filebase}/" | grep -v "#" >> ~{filebase}.ref

        bgzip ~{filebase}.query
        bgzip ~{filebase}.ref

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

task BedComparison{
    input{
        File? query
        File? ref
        String prefix
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 3.75, 
        disk_gb: 10,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output{
        File comparison = "~{prefix}.bed"
    }

    command <<<
        bash /src/compare_callsets_V2.sh \
            -O ~{prefix}.bed -p ~{prefix} ~{query} ~{ref}
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

task InheritanceComparison{
    input{
        File? query_bed
        File? ref_bed
        String prefix
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 1, 
        disk_gb: 20,
        boot_disk_gb: 20,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output{
        File comparison = "~{prefix}.bed"
    }

    command <<<

        #for SVs other than DEL, DUP, CNV, require exact match between query and ref:
        zcat ~{ref_bed} | awk '{if ($5!="DEL" && $5!="DUP" && $5!="CNV" ) print}' > SVs_other_than_CNVs.bed
        grep -wf <(zcat ~{query_bed} | cut -f4 ) SVs_other_than_CNVs.bed > SVs_other_than_CNVs.match.bed

        #for CNVs under 5Kb, require exact match between query and ref:
        zcat ~{ref_bed} | awk '{if ($5=="DEL" || $5=="DUP" || $5=="CNV" ) print}' | awk '{if ($3-$2<5000) print}' > CNVs_under_5Kb.bed
        grep -wf <(zcat ~{query_bed} | cut -f4 ) CNVs_under_5Kb.bed > CNVs_under_5Kb.match.bed

        #for CNVs under 5Kb, require over50% reciprocal overlap with matched SV type
        zcat ~{ref_bed} | awk '{if ($5=="DEL" || $5=="DUP" || $5=="CNV" ) print}' | awk '{if (!$3-$2<5000) print}' > CNVs_over_5Kb.bed

        bedtools intersect -r -f .5 -wo -a CNVs_over_5Kb -b ~{query_bed} | awk '{if ($5==$11) print}' > bedtools_comparison_over5Kb.RO05.bed
        bedtools intersect -r -f .5 -wo -a CNVs_over_5Kb -b ~{query_bed} | awk '{if ($5=="CNV" && $11=="DEL") print}' >> bedtools_comparison_over5Kb.RO05.bed
        bedtools intersect -r -f .5 -wo -a CNVs_over_5Kb -b ~{query_bed} | awk '{if ($5=="CNV" && $11=="DUP") print}' >> bedtools_comparison_over5Kb.RO05.bed
        bedtools intersect -r -f .5 -wo -a CNVs_over_5Kb -b ~{query_bed} | awk '{if ($5=="DEL" && $11=="CNV") print}' >> bedtools_comparison_over5Kb.RO05.bed
        bedtools intersect -r -f .5 -wo -a CNVs_over_5Kb -b ~{query_bed} | awk '{if ($5=="DUP" && $11=="CNV") print}' >> bedtools_comparison_over5Kb.RO05.bed

        cat SVs_other_than_CNVs.match.bed CNVs_under_5Kb.match.bed bedtools_comparison_over5Kb.RO05.bed | cut -f1-5 > ~{prefix}.bed

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

task vcf2vapor{
    input{
        String prefix
        String sample
        File vcf
        File? vcf_index
        Int min_shard_size
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 6,
        disk_gb: 25,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    String filename = basename(vcf, ".vcf.gz")

    output {
        File bed = "~{prefix}.vapor.bed"
        Array[File] vapor_beds = glob("~{prefix}.vapor.split.*")
    }

    command <<<

        set -Eeuo pipefail
        
        gsutil cp ~{vcf} ./tmp.vcf.gz
        tabix -p vcf ./tmp.vcf.gz
        svtk vcf2bed -i SVTYPE -i SVLEN tmp.vcf.gz ~{prefix}.bed
        head -1 ~{prefix}.bed | cut -f1-4,7,8 > ~{prefix}.vapor.bed
        grep ~{sample} ~{prefix}.bed | grep -v "BND" | grep -v "CNV" | sed -e 's/CPX/INV/' | cut -f1-4,7,8 | sed -e 's/INS\t/INS_/' | cut -f1-5 >> ~{prefix}.vapor.bed
        split -l ~{min_shard_size} ~{prefix}.vapor.bed -a 6 ~{prefix}.vapor.split.
        
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

task bed2vapor{
    input{
        String prefix
        String sample
        File bed_file
        File? bed_index
        Int min_shard_size
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 3,
        disk_gb: 15,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output {
        File bed = "~{prefix}.vapor.bed"
        Array[File] vapor_beds = glob("~{prefix}.vapor.split.*")
    }

    command <<<

        set -Eeuo pipefail
         
        gsutil cp ~{bed_file} ./~{prefix}.bed.gz
        gunzip ~{prefix}.bed.gz
        grep -v "BND" ~{prefix}.bed | grep -v "CNV" | sed -e 's/CPX/INV/' | cut -f1-4,6,7 | sed -e 's/INS\t/INS_/' | cut -f1-5 >> ~{prefix}.vapor.bed
        split -l ~{min_shard_size} ~{prefix}.vapor.bed -a 6 ~{prefix}.vapor.split.

       
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

task Bcf2Vcf{
    input{
        String prefix
        String contig
        File bcf
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 3.75, 
        disk_gb: 10,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output{
        File vcf = "~{prefix}.~{contig}.duphold.vcf.gz"
    }

    command <<<
            set -Eeuo pipefail
            bcftools view ~{bcf} | bgzip > ~{prefix}.~{contig}.duphold.vcf.gz
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

task ExtracGTGQ{
    input{
        String prefix
        File vcf_file
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 7.5, 
        disk_gb: 10,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }
  
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output{
        File GQ_GT = "~{prefix}.SVID_gt.tsv"
    }

    command <<<
        zcat ~{vcf_file} | grep -v '#' > ~{prefix}.SVID_gt

        python <<CODE
        import os
        fin=open("~{prefix}.SVID_gt")
        svid_gt={}
        for line in fin:
          pin=line.strip().split()
          svid_gt[pin[2]]=[pin[9].split(':')[pin[8].split(':').index('GT')], 
                          pin[9].split(':')[pin[8].split(':').index('GQ')], 
                          pin[9].split(':')[pin[8].split(':').index('RD_CN')], 
                          pin[9].split(':')[pin[8].split(':').index('RD_GQ')], 
                          pin[9].split(':')[pin[8].split(':').index('PE_GT')], 
                          pin[9].split(':')[pin[8].split(':').index('PE_GQ')], 
                          pin[9].split(':')[pin[8].split(':').index('SR_GT')], 
                          pin[9].split(':')[pin[8].split(':').index('SR_GQ')]]
        fin.close()

        fo=open("~{prefix}.SVID_gt.tsv", 'w')
        print('\t'.join(['SVID','GT','GQ','RD_CN','RD_GQ','PE_GT','PE_GQ','SR_GT','SR_GQ']), file=fo)
        for i in svid_gt.keys():
          print('\t'.join([i]+svid_gt[i]), file=fo)
        fo.close()
        CODE

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

task RunRdPeSrAnnotation{
    input{
        String prefix
        File bed
        File? pe_matrix
        File? pe_index
        File? sr_matrix
        File? sr_index
        File? rd_matrix
        File? rd_index
        File ref_fasta
        File ref_fai
        File ref_dict
        String rdpesr_benchmark_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 15, 
        disk_gb: 20,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output {
        File cov = "~{filebase}.bed.Rd.gz"
        File cov_ri_flank = "~{filebase}.ri_flank.Rd.gz"
        File cov_le_flank = "~{filebase}.le_flank.Rd.gz"
        File pesr_anno = "~{filebase}.bed.PeSr.gz"
    }

    String filebase = basename(bed,".bed.gz")

    command <<<

        set -Eeuo pipefail

        Rscript /src/modify_bed_for_PE_SR_RD_labeling.R \
            -i ~{bed} \
            --le_bp ~{filebase}.le_bp \
            --ri_bp ~{filebase}.ri_bp \
            --le_flank ~{filebase}.le_flank \
            --ri_flank ~{filebase}.ri_flank

        zcat ~{rd_matrix} | grep -v '@' | grep -v CONTIG |bgzip >    bincov.tsv.gz
        Rscript /src/bincov_to_normCov.R -i bincov.tsv.gz
        bgzip normCov.tsv
        tabix -b 2 -e 2 normCov.tsv.gz

        zcat ~{bed} | cut -f1-4,7,8 > ~{filebase}.info
        python3 /src/add_RD_to_SVs.py ~{filebase}.info normCov.tsv.gz ~{filebase}.bed.Rd
        python3 /src/add_RD_to_SVs.py ~{filebase}.ri_flank normCov.tsv.gz ~{filebase}.ri_flank.Rd
        python3 /src/add_RD_to_SVs.py ~{filebase}.le_flank normCov.tsv.gz ~{filebase}.le_flank.Rd
        python3 /src/add_SR_PE_to_PB_INS.V2.py ~{filebase}.info ~{pe_matrix} ~{sr_matrix} ~{filebase}.bed.PeSr

        bgzip ~{filebase}.bed.Rd
        bgzip ~{filebase}.ri_flank.Rd
        bgzip ~{filebase}.le_flank.Rd
        bgzip ~{filebase}.bed.PeSr

    >>>
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: rdpesr_benchmark_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task RunGenomicContextAnnotation{
    input{
        File bed
        File? ref_SegDup
        File? ref_SimpRep
        File? ref_RepMask
        String prefix
        String rdpesr_benchmark_docker
        RuntimeAttr? runtime_attr_override
    }
    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 10, 
        disk_gb: 20,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    command <<<

        zcat ~{bed} | awk '{print $1,$2,$2,$4,$5}' | sed -e 's/ /\t/g' >    ~{prefix}.le_bp
        zcat ~{bed} | awk '{print $1,$3,$3,$4,$5}' | sed -e 's/ /\t/g' >    ~{prefix}.ri_bp
        bedtools coverage -a ~{prefix}.le_bp -b ~{ref_RepMask} | awk '{if ($9>0) print}'>    ~{prefix}.le_bp.vs.RM
        bedtools coverage -a ~{prefix}.le_bp -b ~{ref_SegDup}    | awk '{if ($9>0) print}'>    ~{prefix}.le_bp.vs.SD
        bedtools coverage -a ~{prefix}.le_bp -b ~{ref_SimpRep} | awk '{if ($9>0) print}'>    ~{prefix}.le_bp.vs.SR
        bedtools coverage -a ~{prefix}.ri_bp -b ~{ref_RepMask} | awk '{if ($9>0) print}'>    ~{prefix}.ri_bp.vs.RM
        bedtools coverage -a ~{prefix}.ri_bp -b ~{ref_SegDup}    | awk '{if ($9>0) print}'>    ~{prefix}.ri_bp.vs.SD
        bedtools coverage -a ~{prefix}.ri_bp -b ~{ref_SimpRep} | awk '{if ($9>0) print}'>    ~{prefix}.ri_bp.vs.SR


        Rscript /src/add_GC_anno_to_bed.R \
        -b ~{bed} \
        -o ~{prefix}.GC_anno.bed \
        --left_vs_SR    ~{prefix}.le_bp.vs.SR \
        --left_vs_SD    ~{prefix}.le_bp.vs.SD \
        --left_vs_RM    ~{prefix}.le_bp.vs.RM \
        --right_vs_SR ~{prefix}.ri_bp.vs.SR \
        --right_vs_SD ~{prefix}.ri_bp.vs.SD \
        --right_vs_RM ~{prefix}.ri_bp.vs.RM 
    >>>

    output{
        File anno_bed = "~{prefix}.GC_anno.bed"
    }

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: rdpesr_benchmark_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task ExtracAlgorithmEvidenceFilter{
  input{
    String prefix
    File vcf_file
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 3.75, 
        disk_gb: 10,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }
  
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output{
        File vcf_info = "~{prefix}.info"
    }

    command <<<
        zcat ~{vcf_file} | grep -v '##' | cut -f3,7  > ~{prefix}.SVID_filter
        svtk vcf2bed -i SVTYPE -i SVLEN -i ALGORITHMS -i EVIDENCE ~{vcf_file} ~{prefix}.bed
        paste  <(cut -f4,7-10 ~{prefix}.bed) \
               <(cut -f2 ~{prefix}.SVID_filter) \
               > ~{prefix}.info
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

