## Generate BAF file ##

version 1.0

import "Structs.wdl"

workflow BAFFromGVCFs {
  input {
    Array[File?] gvcfs
    Array[String] samples
    Boolean ignore_missing_gvcfs
    File unpadded_intervals_file
    File dbsnp_vcf
    File? dbsnp_vcf_index
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File inclusion_bed
    String batch
    String? gcs_project_for_requester_pays
    String gatk_docker
    String sv_base_mini_docker
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_merge_vcfs
    RuntimeAttr? runtime_attr_baf_gen
    RuntimeAttr? runtime_attr_merge_baf
  }

  Int num_of_original_intervals = length(read_lines(unpadded_intervals_file))

  # Make a 2.5:1 interval number to samples in callset ratio interval list
  Int possible_merge_count = floor(num_of_original_intervals / num_gvcfs / 2.5)
  Int merge_count = if possible_merge_count > 1 then possible_merge_count else 1

  File dbsnp_vcf_index_ = if defined(dbsnp_vcf_index) then select_first([dbsnp_vcf_index]) else dbsnp_vcf + ".idx"

  scatter (i in range(length(samples))) {
    if (defined(gvcfs[i])) {
      String defined_samples_optional_ = samples[i]
      File defined_gvcfs_optional_ = select_first([gvcfs[i]])
      File gvcf_indexes_optional_ = select_first([gvcfs[i]]) + ".tbi"
    }
  }
  Array[String] defined_samples_ = select_all(defined_samples_optional_)
  Array[File] defined_gvcfs_ = select_all(defined_gvcfs_optional_)
  Array[File] gvcf_indexes_ = select_all(gvcf_indexes_optional_)
  Int num_gvcfs = length(defined_gvcfs_)

  call DynamicallyCombineIntervals {
    input:
      intervals = unpadded_intervals_file,
      merge_count = merge_count,
      preemptible = 3
  }

  Array[String] unpadded_intervals = read_lines(DynamicallyCombineIntervals.output_intervals)

  Int disk_size_gb = 10 + ceil((size(defined_gvcfs_, "GB") + size(gvcf_indexes_, "GB")) * 1.5)

  scatter (idx in range(length(unpadded_intervals))) {
    call ImportGVCFs {
      input:
        sample_names = defined_samples_,
        input_gvcfs = defined_gvcfs_,
        input_gvcfs_indices = gvcf_indexes_,
        interval = unpadded_intervals[idx],
        workspace_dir_name = "genomicsdb",
        disk_size = disk_size_gb,
        batch_size = 50,
        docker = gatk_docker,
        gatk_path = "/gatk/gatk",
        gcs_project_for_requester_pays = gcs_project_for_requester_pays,
        preemptible = 3
    }
    call GenotypeGVCFs {
      input:
        workspace_tar = ImportGVCFs.output_genomicsdb,
        interval = unpadded_intervals[idx],
        output_vcf_filename = "~{batch}.~{idx}.vcf.gz",
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict,
        dbsnp_vcf = dbsnp_vcf,
        dbsnp_vcf_index = dbsnp_vcf_index_,
        disk_size = 200,
        docker = gatk_docker,
        gatk_path = "/gatk/gatk",
        preemptible = 3
    }
    call GenerateBAF {
      input:
        vcf = GenotypeGVCFs.output_vcf,
        samples = samples,
        ignore_missing_vcf_samples = ignore_missing_gvcfs,
        batch = batch,
        shard = "~{idx}",
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_baf_gen,
    }
  }

  call MergeEvidenceFiles {
      input:
        files = GenerateBAF.out,
        indexes = GenerateBAF.out_index,
        batch = batch,
        evidence = "BAF",
        inclusion_bed = inclusion_bed,
        sv_base_mini_docker = sv_base_mini_docker,
        runtime_attr_override = runtime_attr_merge_baf
  }

  output {
    File out = MergeEvidenceFiles.merged
    File out_index = MergeEvidenceFiles.merged_idx
  }
}

task DynamicallyCombineIntervals {
  input {
    File intervals
    Int merge_count
    Int preemptible
  }

  command {
    python << CODE
    def parse_interval(interval):
        colon_split = interval.split(":")
        chromosome = colon_split[0]
        dash_split = colon_split[1].split("-")
        start = int(dash_split[0])
        end = int(dash_split[1])
        return chromosome, start, end

    def add_interval(chr, start, end):
        lines_to_write.append(chr + ":" + str(start) + "-" + str(end))
        return chr, start, end

    count = 0
    chain_count = ~{merge_count}
    l_chr, l_start, l_end = "", 0, 0
    lines_to_write = []
    with open("~{intervals}") as f:
        with open("out.intervals", "w") as f1:
            for line in f.readlines():
                # initialization
                if count == 0:
                    w_chr, w_start, w_end = parse_interval(line)
                    count = 1
                    continue
                # reached number to combine, so spit out and start over
                if count == chain_count:
                    l_char, l_start, l_end = add_interval(w_chr, w_start, w_end)
                    w_chr, w_start, w_end = parse_interval(line)
                    count = 1
                    continue

                c_chr, c_start, c_end = parse_interval(line)
                # if adjacent keep the chain going
                if c_chr == w_chr and c_start == w_end + 1:
                    w_end = c_end
                    count += 1
                    continue
                # not adjacent, end here and start a new chain
                else:
                    l_char, l_start, l_end = add_interval(w_chr, w_start, w_end)
                    w_chr, w_start, w_end = parse_interval(line)
                    count = 1
            if l_char != w_chr or l_start != w_start or l_end != w_end:
                add_interval(w_chr, w_start, w_end)
            f1.writelines("\n".join(lines_to_write))
    CODE
  }

  runtime {
    memory: "3 GB"
    preemptible: preemptible
    docker: "python:2.7"
    maxRetries: "1"
  }

  output {
    File output_intervals = "out.intervals"
  }
}

task GenotypeGVCFs {
  input {
    File workspace_tar
    String interval
    String output_vcf_filename
    String gatk_path

    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File dbsnp_vcf
    File? dbsnp_vcf_index

    String docker
    Int disk_size
    Int preemptible
  }

  parameter_meta {
    dbsnp_vcf: {
      localization_optional: true
    }
    dbsnp_vcf_index: {
      localization_optional: true
    }
  }

  command <<<
    set -e

    tar -xf ~{workspace_tar}
    WORKSPACE=$( basename ~{workspace_tar} .tar)

    ~{gatk_path} --java-options "-Xmx5g -Xms5g" \
     GenotypeGVCFs \
     -R ~{ref_fasta} \
     -O ~{output_vcf_filename} \
     -D ~{dbsnp_vcf} \
     -G StandardAnnotation \
     --only-output-calls-starting-in-intervals \
     --allow-old-rms-mapping-quality-annotation-data \
     -V gendb://$WORKSPACE \
     -L ~{interval}
  >>>
  runtime {
    docker: docker
    memory: "7 GB"
    cpu: "2"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: preemptible
    maxRetries: "1"
  }
  output {
    File output_vcf = "~{output_vcf_filename}"
    File output_vcf_index = "~{output_vcf_filename}.tbi"
  }
}

task ImportGVCFs {
  input {
    Array[String] sample_names
    Array[File] input_gvcfs
    Array[File] input_gvcfs_indices
    String interval

    String workspace_dir_name

    String gatk_path
    String? gcs_project_for_requester_pays
    String docker
    Int disk_size
    Int preemptible
    Int batch_size
  }
  parameter_meta {
    input_gvcfs: {
      localization_optional: true
    }
    input_gvcfs_indices: {
      localization_optional: true
    }
  }

  command <<<
    set -e
    set -o pipefail

    python << CODE
    gvcfs = ['~{sep="','" input_gvcfs}']
    sample_names = ['~{sep="','" sample_names}']

    if len(gvcfs)!= len(sample_names):
      exit(1)

    with open("inputs.list", "w") as fi:
      for i in range(len(gvcfs)):
        fi.write(sample_names[i] + "\t" + gvcfs[i] + "\n")

    CODE

    rm -rf ~{workspace_dir_name}

    # The memory setting here is very important and must be several GB lower
    # than the total memory allocated to the VM because this tool uses
    # a significant amount of non-heap memory for native libraries.
    # Also, testing has shown that the multithreaded reader initialization
    # does not scale well beyond 5 threads, so don't increase beyond that.
    ~{gatk_path} --java-options "-Xmx4g -Xms4g" \
    GenomicsDBImport \
    --genomicsdb-workspace-path ~{workspace_dir_name} \
    --batch-size ~{batch_size} \
    -L ~{interval} \
    --sample-name-map inputs.list \
    --reader-threads 5 \
    -ip 500 \
     ~{"--gcs-project-for-requester-pays " + gcs_project_for_requester_pays}

    tar -cf ~{workspace_dir_name}.tar ~{workspace_dir_name}

  >>>
  runtime {
    docker: docker
    memory: "7 GB"
    cpu: "2"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: preemptible
    maxRetries: "1"
  }
  output {
    File output_genomicsdb = "~{workspace_dir_name}.tar"
  }
}

task GenerateBAF {
  input {
    File vcf
    File? vcf_header
    Array[String] samples
    Boolean ignore_missing_vcf_samples
    String batch
    String shard
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  String ignore_missing_vcf_samples_flag = if (ignore_missing_vcf_samples) then "--ignore-missing-vcf-samples" else ""

  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 3.75,
                               disk_gb: 10,
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File out = "~{batch}.~{shard}.txt.gz"
    File out_index = "~{batch}.~{shard}.txt.gz.tbi"
  }
  command <<<

    set -euo pipefail
    python /opt/sv-pipeline/02_evidence_assessment/02d_baftest/scripts/Filegenerate/generate_baf.py \
      --unfiltered \
      --samples-list ~{write_lines(samples)} \
      ~{ignore_missing_vcf_samples_flag} \
      ~{if defined(vcf_header) then "<(cat ~{vcf_header} ~{vcf})" else vcf} \
      | bgzip \
      > ~{batch}.~{shard}.txt.gz

    tabix -s1 -b2 -e2 ~{batch}.~{shard}.txt.gz
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

task MergeEvidenceFiles {
  input {
    Array[File] files
    Array[File] indexes
    String batch
    String evidence
    File inclusion_bed
    Int? disk_gb_override  # Overrides runtime_attr
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  Int disk_size_gb = 10 + ceil(size(files, "GB") * 15)

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 3.75,
    disk_gb: disk_size_gb,
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File merged = "~{batch}.~{evidence}.txt.gz"
    File merged_idx = "~{batch}.~{evidence}.txt.gz.tbi"
  }
  command <<<

    set -euxo pipefail

    mkdir data
    while read file; do
      filename=`basename $file`
      tabix -h -R ~{inclusion_bed} $file > data/$filename.txt
      rm $file
    done < ~{write_lines(files)}

    mkdir tmp
    sort -m -k1,1V -k2,2n -T tmp data/*.txt | bgzip -c > ~{batch}.~{evidence}.txt.gz
    tabix -s1 -b2 -e2 ~{batch}.~{evidence}.txt.gz

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