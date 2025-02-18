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
    set -e

    output_json="${sample_id}.json"
    echo "{ \"$sample_id\": {" > "$output_json"

    bcftools query -f '%ID\t%FILTER\n' ~{vcf_file} | \
    awk -v sample="$sample_id" '
    BEGIN { first=1 }
    {
      if ($2 != "PASS") {
        gsub(":", "_", $1);
        if (!seen[$2]) {
          if (!first) print ",";
          first=0;
          printf "\"%s\": [\"%s\"]", $2, $1;
          seen[$2]=1;
        } else {
          printf ", \"%s\"", $1;
        }
      }
    }
    END { print "}}" }' >> "$output_json"
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
    set -e
    echo "{" > merged.json

    first=1
    for f in ~{sep=' ' json_files}; do
      if [[ "$first" -eq 0 ]]; then echo "," >> merged.json; fi
      cat "$f" | jq -c . >> merged.json
      first=0
    done

    echo "}" >> merged.json
  >>>

  output {
    File merged_json = "merged.json"
  }

  runtime {
    docker: sv_pipeline_docker
  }
}
