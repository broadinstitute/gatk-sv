version 1.0

workflow ExtractFiltersFromVCFs {

  input {
    Array[String] sample_ids
    Array[File] vcf_files

    String sv_pipeline_docker
  }

  scatter (idx in range(length(vcf_files))) {
    call ExtractSampleAndVariants { 
      input:
        vcf_file = vcf_files[idx],
        sample_id = sample_ids[idx],
        sv_pipeline_docker = sv_pipeline_docker
    }
  }

  call MergeJSONs {
    input: 
      json_files = ExtractSampleAndVariants.output_json,
      sv_pipeline_docker = sv_pipeline_docker
  }

  output {
    File filters_json = MergeJSONs.merged_json
  }
}

task ExtractSampleAndVariants {
  input {
    File vcf_file
    String sample_id
    String sv_pipeline_docker
  }

  command <<<
    set -euo pipefail
    python3 <<CODE

    import json
    import gzip

    sample_id = "~{sample_id}"
    vcf_file = "~{vcf_file}"
    output_json = f"{sample_id}.json"

    filter_dict = {}

    with gzip.open(vcf_file, 'rt') as f:
      for line in f:
        if line.startswith("#"):
          continue
        columns = line.strip().split('\t')
        variant_id = columns[2].replace(":", "_").replace("-", "_")
        filters = columns[6].split(';') if columns[6] != "PASS" else []
        
        for filt in filters:
          if filt not in filter_dict:
            filter_dict[filt] = []
          filter_dict[filt].append(variant_id)

    with open(output_json, "w") as f:
      json.dump({sample_id: filter_dict}, f, indent=4)
    CODE
  >>>

  output {
    File output_json = "${sample_id}.json"
  }

  runtime {
    docker: sv_pipeline_docker
  }
}

task MergeJSONs {
  input {
    Array[File] json_files
    String sv_pipeline_docker
  }

  command <<<
    set -euo pipefail
    python3 <<CODE
    import json
    import glob

    output_file = "merged.json"
    merged_dict = {}

    files_str = "~{sep=" " json_files}"
    files = files_str.split()

    for file in files:
        with open(file, "r") as f:
            data = json.load(f)
            merged_dict.update(data)
        
    with open(output_file, "w") as f:
        json.dump(merged_dict, f, indent=4)
    CODE
  >>>

  output {
    File merged_json = "merged.json"
  }

  runtime {
    docker: sv_pipeline_docker
  }
}
