version 1.0

import "Structs.wdl"
import "TasksMakeCohortVcf.wdl" as tasks_cohort

workflow MakeSlCutoffOptimizationTable {

  input {
    Array[File] vcfs
    File truth_json
    String output_prefix

    File? sample_list
    Boolean labeled_only = true

    Int? subsample_rows
    Int subsample_seed = 42

    String? make_vcf_table_args

    String sv_pipeline_docker
    String linux_docker

    RuntimeAttr? runtime_attr_make_vcf_table
    RuntimeAttr? runtime_attr_concat_vcf_tables
    RuntimeAttr? runtime_attr_concat_class_label_sl_tables
    RuntimeAttr? runtime_attr_subsample_vcf_table
  }

  scatter (i in range(length(vcfs))) {
    call MakeVcfTable {
      input:
        vcf = vcfs[i],
        truth_json = truth_json,
        sample_list = sample_list,
        labeled_only = labeled_only,
        args = make_vcf_table_args,
        output_prefix = "~{output_prefix}.vcf_table_shard_~{i}",
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_make_vcf_table
    }
  }

  call tasks_cohort.ConcatHeaderedTextFiles as ConcatVcfTables {
    input:
      text_files = MakeVcfTable.out,
      gzipped = true,
      output_filename = "~{output_prefix}.vcf_table.tsv.gz",
      linux_docker = linux_docker,
      runtime_attr_override = runtime_attr_concat_vcf_tables
  }

  call tasks_cohort.ConcatHeaderedTextFiles as ConcatClassLabelSlTables {
    input:
      text_files = MakeVcfTable.class_label_sl,
      gzipped = true,
      output_filename = "~{output_prefix}.vcf_table.CLASS_LABEL_SL.tsv.gz",
      linux_docker = linux_docker,
      runtime_attr_override = runtime_attr_concat_class_label_sl_tables
  }

  if (defined(subsample_rows)) {
    call SubsampleVcfTable {
      input:
        table = ConcatVcfTables.out,
        rows = select_first([subsample_rows]),
        seed = subsample_seed,
        output_prefix = output_prefix,
        linux_docker = linux_docker,
        runtime_attr_override = runtime_attr_subsample_vcf_table
    }
  }

  output {
    File vcf_table = ConcatVcfTables.out
    File class_label_sl_table = ConcatClassLabelSlTables.out
    File? vcf_table_subsample = SubsampleVcfTable.out
    Array[File] vcf_table_shards = MakeVcfTable.out
  }
}

task MakeVcfTable {
  input {
    File vcf
    File truth_json
    File? sample_list
    Boolean labeled_only
    String? args
    String output_prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 7.5,
                               disk_gb: ceil(10 + size(vcf, "GB") * 2 + size(truth_json, "GB") * 2),
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File out = "~{output_prefix}.tsv.gz"
    File class_label_sl = "~{output_prefix}.CLASS_LABEL_SL.tsv.gz"
  }
  command <<<
    set -euo pipefail

    TRUTH_JSON="~{truth_json}"
    if ~{defined(sample_list)}; then
      python3 <<CODE
    import json
    with open("~{sample_list}") as f:
        keep = set(line.strip() for line in f if line.strip())
    with open("~{truth_json}") as f:
        data = json.load(f)
    with open("truth.subset.json", "w") as f:
        json.dump({k: v for k, v in data.items() if k in keep}, f)
    CODE
      TRUTH_JSON="truth.subset.json"
    fi

    python /opt/sv-pipeline/scripts/make_sl_table.py \
      --vcf ~{vcf} \
      --truth-json $TRUTH_JSON \
      --out ~{output_prefix}.tsv.gz \
      ~{if labeled_only then "--labeled-only" else ""} \
      ~{args}

    zcat ~{output_prefix}.tsv.gz \
      | awk -F'\t' -v OFS='\t' '
          NR==1 {
            for (i=1; i<=NF; i++) {
              if ($i=="FILTER_CLASS" && !c) c=i
              if ($i=="LABEL" && !l) l=i
              if ($i=="SL" && !s) s=i
            }
            if (!c || !l || !s) { print "Missing required column" > "/dev/stderr"; exit 1 }
            print $c, $l, $s
            next
          }
          { v=$s; if (v=="None" || v=="") v="NA"; print $c, $l, v }' \
      | gzip > ~{output_prefix}.CLASS_LABEL_SL.tsv.gz
  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " SSD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

task SubsampleVcfTable {
  input {
    File table
    Int rows
    Int seed
    String output_prefix
    String linux_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 3.75,
                               disk_gb: ceil(10 + size(table, "GB") * 12),
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  String output_filename = "~{output_prefix}.vcf_table.random_~{rows}.headered.tsv.gz"

  output {
    File out = output_filename
  }
  command <<<
    set -euo pipefail
    zcat ~{table} > table.tsv
    {
      head -n 1 table.tsv
      tail -n +2 table.tsv | shuf -n ~{rows} --random-source=<(yes ~{seed})
    } | gzip > ~{output_filename}
  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " SSD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: linux_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}
