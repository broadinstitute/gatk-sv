version 1.0

import "Structs.wdl"

task GetSampleIdsFromVcf {
  input {
    File vcf
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  String sample_list = basename(vcf, ".vcf.gz") + ".samples.txt"

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 0.9,
    disk_gb: 2 + ceil(size(vcf, "GiB")),
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<

    set -eu
    bcftools query -l ~{vcf} > ~{sample_list}

  >>>

  output {
    File out_file = sample_list
    Array[String] out_array = read_lines(sample_list)
  }

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

task GetSampleIdsFromVcfArray {
  input {
    Array[File] vcfs
    String prefix
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 0.9,
                               disk_gb: 2 + ceil(size(vcfs, "GiB")),
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<

    set -eu
    touch ~{prefix}.txt
    while read VCF; do
      bcftools query -l $VCF >> ~{prefix}.txt
    done < ~{write_lines(vcfs)}

  >>>

  output {
    File out_file = "~{prefix}.txt"
    Array[String] out_array = read_lines("~{prefix}.txt")
  }

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

task GetSampleIdsFromVcfTar {
  input {
    File vcf_tar
    String prefix
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 0.9,
                               disk_gb: ceil(50 + 2 * size(vcf_tar, "GiB")),
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail
    mkdir vcfs
    tar xzf ~{vcf_tar} -C vcfs/
    ls vcfs/*.vcf.gz | xargs -n1 bcftools query -l | sort -u > ~{prefix}.txt
  >>>

  output {
    File out_file = "~{prefix}.txt"
    Array[String] out_array = read_lines("~{prefix}.txt")
  }

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

task CountSamples {
  input {
    File vcf
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 3.75,
                               disk_gb: 10 + ceil(size(vcf, "GiB")),
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -eu
    bcftools query -l ~{vcf} | wc -l > sample_count.txt
  >>>

  output {
    Int num_samples = read_int("sample_count.txt")
  }

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

task GetSampleIdsFromMedianCoverageFile {
  input {
    File median_file
    String name
    String linux_docker
    RuntimeAttr? runtime_attr_override
  }

  String sample_list = name + ".samples.txt"

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 0.9,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<

    set -euo pipefail
    head -1 ~{median_file} | sed -e 's/\t/\n/g' > ~{sample_list}

  >>>

  output {
    File out_file = sample_list
    Array[String] out_array = read_lines(sample_list)
  }

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: linux_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

task RunQC {
  input {
    String name
    File metrics
    File qc_definitions
    String sv_pipeline_docker
    Float mem_gib = 1
    Int disk_gb = 10
    Int preemptible_attempts = 3
  }

  output {
    File out = "sv_qc.~{name}.tsv"
  }
  command <<<

    set -eu
    svqc ~{metrics} ~{qc_definitions} raw_qc.tsv
    grep -vw "NA" raw_qc.tsv > sv_qc.~{name}.tsv

  >>>
  runtime {
    cpu: 1
    memory: "~{mem_gib} GiB"
    disks: "local-disk ~{disk_gb} HDD"
    bootDiskSizeGb: 10
    docker: sv_pipeline_docker
    preemptible: preemptible_attempts
    maxRetries: 1
  }
}

task GetFilteredSubsampledIndices {
  input {
    File all_strings
    File? exclude_strings
    Int? subset_size
    Int seed
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  String include_indices_filename = "~{prefix}.filtered_subsampled_indices.list"
  String include_strings_filename = "~{prefix}.filtered_subsampled_strings.list"

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 1,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<

    set -euo pipefail
    python3 <<CODE
    import random
    string_array = [line.rstrip() for line in open("~{all_strings}", 'r')]
    exclude_set = set()
    if "~{exclude_strings}" != "":
      exclude_set = set([line.rstrip() for line in open("~{exclude_strings}", 'r')])
    filtered_array = [s for s in string_array if s not in exclude_set]
    subset_size = len(filtered_array) if ~{subset_size} is None else min(~{subset_size}, len(filtered_array))
    result_array = filtered_array
    if subset_size < len(filtered_array):
        random.seed(~{seed})
        selected_indices = random.sample(range(len(filtered_array)), k=subset_size)
        selected_indices.sort()
        result_array = [filtered_array[i] for i in selected_indices]
    with open("~{include_indices_filename}", 'w') as indices, open("~{include_strings_filename}", 'w') as strings:
        for _, sample in enumerate(result_array):
            indices.write(f"{string_array.index(sample)}\n")
            strings.write(f"{sample}\n")
    CODE

  >>>

  output {
    File include_indices_file = include_indices_filename
    Array[Int] include_indices_array = read_lines(include_indices_filename)
    File include_strings_file = include_strings_filename
    Array[String] include_strings_array = read_lines(include_strings_filename)
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

task RandomSubsampleStringArray {
  input {
    File strings
    Int seed
    Int subset_size
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  String subsample_indices_filename = "~{prefix}.subsample_indices.list"
  String subsampled_strings_filename = "~{prefix}.subsampled_strings.list"

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 1,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<

    set -euo pipefail
    python3 <<CODE
    import random
    string_array = [line.rstrip() for line in open("~{strings}", 'r')]
    array_len = len(string_array)
    if ~{subset_size} > array_len:
      raise ValueError("Subsample quantity ~{subset_size} cannot > array length %d" % array_len)
    random.seed(~{seed})
    numbers = random.sample(range(0, array_len), k=~{subset_size})
    numbers.sort()
    with open("~{subsample_indices_filename}", 'w') as indices, open("~{subsampled_strings_filename}", 'w') as strings:
      for num in numbers:
        indices.write(f"{num}\n")
        strings.write(string_array[num] + "\n")
    CODE

  >>>

  output {
    File subsample_indices_file = subsample_indices_filename
    Array[Int] subsample_indices_array = read_lines(subsample_indices_filename)
    File subsampled_strings_file = subsampled_strings_filename
    Array[String] subsampled_strings_array = read_lines(subsampled_strings_filename)
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

task GetSubsampledIndices {
  input {
    File all_strings
    File subset_strings
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  String subsample_indices_filename = "~{prefix}.subsample_indices.list"

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 1,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<

    set -euo pipefail
    python3 <<CODE
    all_strings = [line.rstrip() for line in open("~{all_strings}", 'r')]
    subset_strings = {line.rstrip() for line in open("~{subset_strings}", 'r')}
    if not subset_strings.issubset(set(all_strings)):
      raise ValueError("Subset list must be a subset of full list")
    with open("~{subsample_indices_filename}", 'w') as indices:
      for i, string in enumerate(all_strings):
        if string in subset_strings:
          indices.write(f"{i}\n")
    CODE

  >>>

  output {
    Array[Int] subsample_indices_array = read_lines(subsample_indices_filename)
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


task SubsetPedFile {
  input {
    File ped_file
    File sample_list
    String subset_name = "subset"
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  String ped_subset_filename = basename(ped_file, ".ped") + ".~{subset_name}.ped"

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 3.75,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<

    set -euo pipefail
    awk 'FNR==NR {a[$1]; next}; $2 in a' ~{sample_list} ~{ped_file} > ~{ped_subset_filename}

  >>>

  output {
    File ped_subset_file = ped_subset_filename
  }

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


task ValidatePedFile {
  input {
    File ped_file
    File sample_list
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

  command <<<

    set -euo pipefail
    python /opt/sv-pipeline/scripts/validate_ped.py -p ~{ped_file} -s ~{sample_list}

  >>>

  output {
    File output_ped = ped_file
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


task LocalizeCloudFileWithCredentials {
  input {
    String cloud_file_path
    String service_account_json
    Int disk_size
    String cloud_sdk_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 0.9,
    disk_gb: disk_size,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command {
    set -euo pipefail

    gsutil cp '~{service_account_json}' local.service_account.json
    gcloud auth activate-service-account --key-file='local.service_account.json'

    gsutil cp '~{cloud_file_path}' .
  }

  output {
    File output_file = basename(cloud_file_path)
  }

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: cloud_sdk_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

task GetVcfSize {
    input {
        File vcf
        File vcf_index
        String samtools_cloud_docker
    }

    Int disk_gb = round(10 + size(vcf, "GiB"))
    String num_records_file = "num_records.txt"
    String num_samples_file = "num_samples.txt"
    # if vcf_index is not supplied, try this path automatically:
    String automatic_vcf_index = vcf + ".tbi"

    runtime {
        docker: samtools_cloud_docker
        cpu: 1
        preemptible: 3
        max_retries: 1
        memory: "2 GiB"
        disks: "local-disk " + disk_gb + " HDD"
    }

    command <<<
        set -euo pipefail

        # symlink vcf_index to current working dir
        ln -s ~{vcf_index} .

        export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
        bcftools query -l ~{vcf} | wc -w > ~{num_samples_file}
        # get num records from index.
        {
            bcftools index --nrecords ~{vcf} || {
                # indices built by GATK are broken and won't work. If that happens:
                #   rm symlink to vcf_index
                #   build index on the fly
                #   get records from good index
                #   delete the good index (in case this is run locally)
                rm ~{vcf_index}
                bcftools index -f -t -o "$(basename ~{vcf}).tbi" ~{vcf}
                bcftools index --nrecords ~{vcf}
                rm "$(basename ~{vcf}).tbi"
            }
        } > ~{num_records_file}
    >>>

    output {
        Int num_records = read_int(num_records_file)
        Int num_samples = read_int(num_samples_file)
        Float num_entries = read_float(num_records_file) * num_samples
    }
}


task MaxInts {
    input {
        Array[Int] ints
    }

    command <<<
        awk 'BEGIN {z=0;} {z=(z>=$0?z:$0);} END {print z;}' "~{write_lines(ints)}"
    >>>

    output {
        Int max_int = read_int(stdout())
    }

    runtime {
        docker: "ubuntu:latest"
        cpu: 1
        preemptible: 3
        max_retries: 1
        memory: "1 GiB"
        disks: "local-disk 10 HDD"
    }
}


task WriteLines {
  input {
    Array[String] lines
    String output_filename
    String linux_docker
  }

  command <<<
    cat ~{write_lines(lines)} > ~{output_filename}
  >>>

  output {
    File out = "~{output_filename}"
  }

  runtime {
    cpu: 1
    memory: "0.9 GiB"
    disks: "local-disk 10 HDD"
    docker: linux_docker
    preemptible: 3
    maxRetries: 1
  }
}

task TarFiles {
  input {
    String prefix
    Array[File] files
    String linux_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 1.0,
                               disk_gb: ceil(10 + 2 * size(files, "GB")),
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail
    tar czf ~{prefix}.tar.gz -T ~{write_lines(files)}
  >>>

  output {
    File out = "~{prefix}.tar.gz"
  }

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: linux_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

task UntarFiles {
  input {
    File tar
    String? glob_suffix
    String linux_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 1.0,
                               disk_gb: ceil(10 + 2 * size(tar, "GB")),
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  String glob_arg = "out/*" + select_first([glob_suffix, ""])

  command <<<
    set -euo pipefail
    mkdir out
    tar xzf ~{tar} -C out/
  >>>

  output {
    Array[File] out = glob("~{glob_arg}")
  }

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: linux_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

task CombineTars {
  input {
    File tar1
    File tar2
    String linux_docker
    RuntimeAttr? runtime_attr_override
  }

  String output_filename = basename(tar2)

  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 1.0,
                               disk_gb: ceil(10 + 3 * size([tar1, tar2], "GB")),
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -e
    # Not the most efficient space-wise, but concatenating compressed tars requires that -i be used
    # when decompressing, so this is safer.
    mkdir tmp
    tar xzf ~{tar1} -C tmp/
    tar xzf ~{tar2} -C tmp/
    tar czf ~{output_filename} -C tmp/ .
  >>>

  output {
    File out = "~{output_filename}"
  }

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: linux_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

# Subset a VCF to a specific subset of samples
task SubsetVcfBySamplesList {
  input {
    File vcf
    File? vcf_idx
    File list_of_samples  # List of samples to keep (default, remove_samples = false) or remove (remove_samples = true)
    String? outfile_name
    Boolean remove_samples = false  # If false (default), keep samples in provided list. If true, remove them.
    Boolean remove_private_sites = true  # If true (default), remove sites that are private to excluded samples. If false, keep sites even if no remaining samples are non-ref.
    Boolean keep_af = true  # If true (default), do not recalculate allele frequencies (AC/AF/AN)
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  String vcf_subset_filename = select_first([outfile_name, basename(vcf, ".vcf.gz") + ".subset.vcf.gz"])
  String vcf_subset_idx_filename = vcf_subset_filename + ".tbi"

  String remove_private_sites_flag = if remove_private_sites then " | bcftools view -e 'SVTYPE!=\"CNV\" && COUNT(GT=\"alt\")==0' " else ""
  String keep_af_flag = if keep_af then "--no-update" else ""
  String complement_flag = if remove_samples then "^" else ""

  # Disk must be scaled proportionally to the size of the VCF
  Float input_size = size(vcf, "GiB")
  RuntimeAttr default_attr = object {
    mem_gb: 3.75,
    disk_gb: ceil(10.0 + (input_size * 2)),
    cpu_cores: 1,
    preemptible_tries: 3,
    max_retries: 1,
    boot_disk_gb: 10
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<

    set -euo pipefail

    bcftools view \
      -S ~{complement_flag}~{list_of_samples} \
      --force-samples \
      ~{keep_af_flag} \
      ~{vcf} \
      ~{remove_private_sites_flag} \
      -O z \
      -o ~{vcf_subset_filename}

    tabix -f -p vcf ~{vcf_subset_filename}
    
  >>>

  output {
    File vcf_subset = vcf_subset_filename
    File vcf_subset_index = vcf_subset_idx_filename
  }

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

# Subset a VCF to a specific sample
task SubsetVcfToSample {
  input {
    File vcf
    File? vcf_idx
    String sample
    String? outfile_name
    Boolean remove_sample = false  # If false (default), keep the sample If true, remove it.
    Boolean remove_private_sites = true  # If true (default), remove sites that are private to excluded samples. If false, keep sites even if no remaining samples are non-ref.
    Boolean keep_af = true  # If true (default), do not recalculate allele frequencies (AC/AF/AN)
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  String vcf_subset_filename = select_first([outfile_name, basename(vcf, ".vcf.gz") + ".subset.vcf.gz"])
  String vcf_subset_idx_filename = vcf_subset_filename + ".tbi"

  String remove_private_sites_flag = if remove_private_sites then " | bcftools view -e 'SVTYPE!=\"CNV\" && COUNT(GT=\"alt\")==0' " else ""
  String keep_af_flag = if keep_af then "--no-update" else ""
  String complement_flag = if remove_sample then "^" else ""

  # Disk must be scaled proportionally to the size of the VCF
  Float input_size = size(vcf, "GiB")
  RuntimeAttr default_attr = object {
                               mem_gb: 3.75,
                               disk_gb: ceil(10.0 + (input_size * 2)),
                               cpu_cores: 1,
                               preemptible_tries: 3,
                               max_retries: 1,
                               boot_disk_gb: 10
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<

    set -euo pipefail

    bcftools view \
      -s ~{complement_flag}~{sample} \
      --force-samples \
      ~{keep_af_flag} \
      ~{vcf} \
      ~{remove_private_sites_flag} \
      -O z \
      -o ~{vcf_subset_filename}

    tabix -f -p vcf ~{vcf_subset_filename}

  >>>

  output {
    File vcf_subset = vcf_subset_filename
    File vcf_subset_index = vcf_subset_idx_filename
  }

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

task VcfToBed {

    input {
        File vcf_file
        String variant_interpretation_docker
        String? args
        RuntimeAttr? runtime_attr_override
    }

    Float vcf_size = size(vcf_file, "GB")

    RuntimeAttr default_attr = object {
        mem_gb: 3.75,
        disk_gb: ceil(10 + vcf_size * 1.5),
        cpu_cores: 1,
        preemptible_tries: 2,
        max_retries: 1,
        boot_disk_gb: 8
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output {
        File bed_output = "~{basename}.bed.gz"
    }

    String basename = basename(vcf_file, ".vcf.gz")
    command {
        set -exuo pipefail

        svtk vcf2bed ~{vcf_file} ~{args} ~{basename}.bed
        bgzip ~{basename}.bed
    }

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: variant_interpretation_docker
    }
}

task GetSampleSex {
  input {
    File ped_file
    String sample_id
    String unknown_sex
    String linux_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 0.9,
                               disk_gb: 10,
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<

    set -euo pipefail
    awk -F '\t' '{if ($2=="~{sample_id}") { if ($5 == 1) {print "male"} else if ($5 == 2) {print "female"} else {print "~{unknown_sex}"}}}' ~{ped_file} > ~{sample_id}.sex.txt

    # Fail if the sample id wasn't found
    if ! [ -s ~{sample_id}.sex.txt ]; then
      echo "ERROR: Sample ~{sample_id} not found in ped file ~{ped_file}"
      exit 1
    fi

  >>>

  output {
    File out_file = "${sample_id}.sex.txt"
    String out_string = read_string("${sample_id}.sex.txt")
  }

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: linux_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}