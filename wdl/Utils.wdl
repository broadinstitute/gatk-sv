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
                               disk_gb: 10 + ceil(size(vcf_tar, "GiB")),
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<

    set -euo pipefail
    # Using a named pipe keeps pipefail from triggering
    mkfifo tmppipe
    while read f
    do
      tar -Ozxf ~{vcf_tar} $f &>2 /dev/null > tmppipe &
      bcftools query -l tmppipe
    done < <(tar -tzf ~{vcf_tar} | grep '\.vcf\.gz$') | sort -u > ~{prefix}.txt

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
    String sv_pipeline_base_docker
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
    docker: sv_pipeline_base_docker
    preemptible: preemptible_attempts
    maxRetries: 1
  }

}

task RandomSubsampleStringArray {
  input {
    File strings
    Int seed
    Int subset_size
    String prefix
    String sv_pipeline_base_docker
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
    docker: sv_pipeline_base_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

task GetSubsampledIndices {
  input {
    File all_strings
    File subset_strings
    String prefix
    String sv_pipeline_base_docker
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
    docker: sv_pipeline_base_docker
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

    parameter_meta {
        vcf: {
          localization_optional: true
        }
    }

    Int disk_gb = round(10 + size(vcf_index, "GiB"))
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
        Int num_entries = num_records * num_samples
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

  String glob_arg = "out/*" + glob_suffix

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
    File list_of_samples_to_keep
    String subset_name = "subset"
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  String vcf_subset_filename = basename(vcf, ".vcf.gz") + ".~{subset_name}.vcf.gz"
  String vcf_subset_idx_filename = vcf_subset_filename + ".tbi"

  # Disk must be scaled proportionally to the size of the VCF
  Float input_size = size(vcf, "GiB")
  RuntimeAttr default_attr = object {
    mem_gb: 3.75,
    disk_gb: ceil(10.0 + (input_size * 1.5)),
    cpu_cores: 1,
    preemptible_tries: 3,
    max_retries: 1,
    boot_disk_gb: 10
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<

    set -euo pipefail

    bcftools view \
      -S ~{list_of_samples_to_keep} \
      --force-samples \
      ~{vcf} \
    | bcftools view \
      --min-ac 1 \
      -O z \
      -o ~{vcf_subset_filename}

    tabix -f -p vcf ~{vcf_subset_filename}
    
  >>>

  output {
    File vcf_subset = vcf_subset_filename
    File vcf_subset_idx = vcf_subset_idx_filename
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


task TransferVcfAnnotations {
  input {
    File vcf_to_annotate
    File vcf_to_annotate_index
    File vcf_with_annotations
    File vcf_with_annotations_index
    Array[String] annotations_to_transfer
    String samtools_cloud_docker
    String output_file_name = sub(sub(basename(vcf_to_annotate), ".gz$", ""), ".vcf$", "_annotated.vcf.gz")
  }

  parameter_meta {
    vcf_to_annotate: {
      localization_optional: true
    }
    vcf_with_annotations: {
      localization_optional: true
    }
  }

  Int disk_gb = round(100 + size([vcf_to_annotate, vcf_to_annotate_index,
                                 vcf_with_annotations, vcf_with_annotations_index], "GiB"))

  runtime {
      docker: samtools_cloud_docker
      cpu: 1
      preemptible: 3
      max_retries: 1
      memory: "2 GiB"
      disks: "local-disk 10 HDD"
  }

  command <<<
    # if running in a local mode, this will fail, but it also won't be *needed*
    export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`

    set -euo pipefail

    bcftools annotate \
      -a ~{vcf_with_annotations} \
      -c ~{sep=',' annotations_to_transfer} \
      -Oz -o "~{output_file_name}" \
      --threads 2 \
      ~{vcf_to_annotate}

    bcftools index --tbi "~{output_file_name}" -o "~{output_file_name}.tbi"
  >>>

  output {
    File annotated_vcf = output_file_name
    File annotated_vcf_index = output_file_name + ".tbi"
  }
}
