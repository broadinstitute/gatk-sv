version 1.0

# Gets sample ID from a BAM/CRAM file and also generates a version of the ID that is safe for use in the pipeline

workflow GetSampleID {
  input {
    File bam_or_cram_file

    # Use only for crams in requester pays buckets
    Boolean requester_pays_crams = false

    Int sample_id_hash_length = 6

    # Docker
    String sv_pipeline_docker
    String samtools_cloud_docker
  }

  if (requester_pays_crams) {
    call GetBamIDRequesterPays {
      input:
        bam_or_cram_file = bam_or_cram_file,
        samtools_cloud_docker = samtools_cloud_docker
    }
  }

  if (!requester_pays_crams) {
    call GetBamID {
      input:
        bam_or_cram_file = bam_or_cram_file,
        samtools_cloud_docker = samtools_cloud_docker
    }
  }

  String bam_id_ = select_first([GetBamIDRequesterPays.out, GetBamID.out])
  call InternalSampleID {
    input:
      external_sample_id = bam_id_,
      bam_or_cram_path = bam_or_cram_file,
      hash_length = sample_id_hash_length,
      sv_pipeline_docker = sv_pipeline_docker
  }

  output {
    String bam_id = bam_id_
    String safe_id = InternalSampleID.out
  }
}

task GetBamID {
  input {
    File bam_or_cram_file
    String samtools_cloud_docker
  }

  parameter_meta {
    bam_or_cram_file: {
      localization_optional: true
    }
  }

  output {
    String out = read_lines("sample.txt")[0]
  }
  command <<<
    set -euo pipefail
    export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
    samtools view -H ~{bam_or_cram_file} \
      | grep "^@RG" \
      | awk -F "\t" '{for(i=2;i<=NF;i++){if($i~/^SM:/){a=$i}} print substr(a,4)}' \
      | sort \
      | uniq \
      > sample.txt
    NUM_IDS=$(wc -l < sample.txt)
    if [[ ${NUM_IDS} -eq 0 ]]; then
      echo "No sample IDs were found in the BAM header"
      exit 1
    fi
    if [[ ${NUM_IDS} -gt 1 ]]; then
      echo "Multiple sample IDs were found in the BAM header"
      exit 1
    fi
  >>>
  runtime {
    docker: samtools_cloud_docker
    memory: "1 GB"
    cpu: "1"
    disks: "local-disk 10 HDD"
    preemptible: "3"
    maxRetries: "1"
  }
}

task GetBamIDRequesterPays {
  input {
    File bam_or_cram_file
    String samtools_cloud_docker
  }

  Int disk_size_gb = 10 + ceil(size(bam_or_cram_file, "GB"))

  output {
    String out = read_lines("sample.txt")[0]
  }
  command <<<
    set -euo pipefail
    export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
    samtools view -H ~{bam_or_cram_file} \
      | grep "^@RG" \
      | awk -F "\t" '{for(i=2;i<=NF;i++){if($i~/^SM:/){a=$i}} print substr(a,4)}' \
      | sort \
      | uniq \
      > sample.txt
    NUM_IDS=$(wc -l < sample.txt)
    if [[ ${NUM_IDS} -eq 0 ]]; then
      echo "No sample IDs were found in the BAM header"
      exit 1
    fi
    if [[ ${NUM_IDS} -gt 1 ]]; then
      echo "Multiple sample IDs were found in the BAM header"
      exit 1
    fi
  >>>
  runtime {
    docker: samtools_cloud_docker
    memory: "1 GB"
    cpu: "1"
    disks: "local-disk ~{disk_size_gb} HDD"
    preemptible: "3"
    maxRetries: "1"
  }
}

task InternalSampleID {
  input {
    String external_sample_id
    String bam_or_cram_path
    Int hash_length
    String sv_pipeline_docker
  }

  output {
    String out = read_lines("external_id.txt")[0]
  }
  command <<<
    set -euo pipefail
    HASH=$(echo -n "~{external_sample_id}~{bam_or_cram_path}" | openssl sha1 | awk '{print substr($2,0,~{hash_length})}')
    SAFE_ID=$(echo -n "~{external_sample_id}" | sed 's/[^a-zA-Z0-9]/_/g')
    echo "__${SAFE_ID}__${HASH}" > external_id.txt
  >>>
  runtime {
    docker: sv_pipeline_docker
    memory: "1 GB"
    cpu: "1"
    disks: "local-disk 10 HDD"
    preemptible: "3"
    maxRetries: "1"
  }
}
