version 1.0

import "Structs.wdl"
import "TasksMakeCohortVcf.wdl" as tasks_cohort
import "GenerateBatchMetrics.wdl" as gbm

workflow MakeSlCutoffOptimizationTable {

  input {
    Array[File] vcfs
    Array[File]? vcf_indexes  # Defaults to "<vcf>.tbi"; only used if region_names is non-empty
    File truth_json
    String output_prefix

    File? sample_list
    Boolean labeled_only = true

    Int? downsample_factor

    String? make_vcf_table_args

    # Genomic context annotation with GATK SVRegionOverlap. Skipped if region_names is empty.
    # For each track NAME, columns NUM_END_OVERLAPS_<NAME> and OVERLAP_FRAC_<NAME> are appended to the tables
    # (NAME is converted to upper-case by the tool).
    Array[File] region_files = []
    Array[File] region_file_indexes = []
    Array[String] region_names = []
    File? reference_dict  # Required if region_names is non-empty
    String? gatk_docker   # Required if region_names is non-empty

    String sv_pipeline_docker
    String linux_docker

    RuntimeAttr? runtime_attr_region_overlap
    RuntimeAttr? runtime_attr_make_vcf_table
    RuntimeAttr? runtime_attr_concat_vcf_tables
    RuntimeAttr? runtime_attr_concat_class_label_sl_tables
    RuntimeAttr? runtime_attr_concat_downsampled_vcf_tables
  }

  Boolean annotate_regions = length(region_names) > 0

  scatter (i in range(length(vcfs))) {
    if (annotate_regions) {
      File default_vcf_index = vcfs[i] + ".tbi"
      File vcf_index_ = if defined(vcf_indexes) then select_first([vcf_indexes])[i] else default_vcf_index
      call gbm.SVRegionOverlap {
        input:
          vcf = vcfs[i],
          vcf_index = vcf_index_,
          reference_dict = select_first([reference_dict]),
          output_prefix = "~{output_prefix}.region_overlap_shard_~{i}",
          region_files = region_files,
          region_file_indexes = region_file_indexes,
          region_names = region_names,
          gatk_docker = select_first([gatk_docker]),
          runtime_attr_override = runtime_attr_region_overlap
      }
    }

    call MakeVcfTable {
      input:
        vcf = select_first([SVRegionOverlap.out, vcfs[i]]),
        truth_json = truth_json,
        sample_list = sample_list,
        labeled_only = labeled_only,
        downsample_factor = downsample_factor,
        args = make_vcf_table_args,
        region_names = region_names,
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

  if (defined(downsample_factor)) {
    call tasks_cohort.ConcatHeaderedTextFiles as ConcatDownsampledVcfTables {
      input:
        text_files = select_all(MakeVcfTable.downsampled),
        gzipped = true,
        output_filename = "~{output_prefix}.vcf_table.downsample_~{select_first([downsample_factor])}.tsv.gz",
        linux_docker = linux_docker,
        runtime_attr_override = runtime_attr_concat_downsampled_vcf_tables
    }
  }

  output {
    File vcf_table = ConcatVcfTables.out
    File class_label_sl_table = ConcatClassLabelSlTables.out
    File? vcf_table_downsampled = ConcatDownsampledVcfTables.out
    Array[File] vcf_table_shards = MakeVcfTable.out
  }
}

task MakeVcfTable {
  input {
    File vcf
    File truth_json
    File? sample_list
    Boolean labeled_only
    Int? downsample_factor
    String? args
    Array[String] region_names = []
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
    File? downsampled = "~{output_prefix}.downsampled.tsv.gz"
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

    # Appends SVRegionOverlap INFO fields to a table, matched on the VID column
    cat > add_region_overlap_columns.py <<'CODE'
    import gzip
    import math
    import sys
    import pysam

    vcf_path, table_path, out_path = sys.argv[1:4]
    names = [n.upper() for n in sys.argv[4:]]
    fields = []
    for n in names:
        fields += ["NUM_END_OVERLAPS_" + n, "OVERLAP_FRAC_" + n]

    def fmt(val):
        if val is None:
            return "None"
        if isinstance(val, float) and math.isnan(val):
            return "NaN"
        if isinstance(val, tuple):
            return ",".join(str(x) for x in val)
        return str(val)

    values = dict()
    with pysam.VariantFile(vcf_path) as vcf:
        missing = [f for f in fields if f not in vcf.header.info]
        if missing:
            sys.stderr.write("INFO fields not found in VCF header: " + ", ".join(missing) + "\n")
            sys.exit(1)
        for r in vcf:
            values[r.id] = [fmt(r.info.get(f, None)) for f in fields]

    default = ["None"] * len(fields)
    with gzip.open(table_path, "rt") as fin, gzip.open(out_path, "wt") as fout:
        header = fin.readline().rstrip("\n").split("\t")
        vid_col = header.index("VID")
        fout.write("\t".join(header + fields) + "\n")
        for line in fin:
            tokens = line.rstrip("\n").split("\t")
            fout.write("\t".join(tokens + values.get(tokens[vid_col], default)) + "\n")
    CODE

    function add_region_overlap_columns() {
      if ~{length(region_names) > 0}; then
        python3 add_region_overlap_columns.py ~{vcf} "$1" "$1.tmp.gz" ~{sep=" " region_names}
        mv "$1.tmp.gz" "$1"
      fi
    }

    python /opt/sv-pipeline/scripts/make_sl_table.py \
      --vcf ~{vcf} \
      --truth-json $TRUTH_JSON \
      --out ~{output_prefix}.tsv.gz \
      ~{if labeled_only then "--labeled-only" else ""} \
      ~{args}
    add_region_overlap_columns ~{output_prefix}.tsv.gz

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

    if ~{defined(downsample_factor)}; then
      python /opt/sv-pipeline/scripts/make_sl_table.py \
        --vcf ~{vcf} \
        --truth-json $TRUTH_JSON \
        --out ~{output_prefix}.downsampled.tsv.gz \
        --downsample-factor ~{default="1" downsample_factor} \
        ~{if labeled_only then "--labeled-only" else ""} \
        ~{args}
      add_region_overlap_columns ~{output_prefix}.downsampled.tsv.gz
    fi
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
