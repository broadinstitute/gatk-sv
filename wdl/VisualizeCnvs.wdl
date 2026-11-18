version 1.0

import "Structs.wdl"

# Plots CNV depth profiles across batches
workflow VisualizeCnvs {
  input {
    Array[File] vcfs
    # tab-delimited files of CNVs
    # chr,start,end,id,svtype,samples
    # samples must be comma-delimited
    Array[File]? cnvs
    # variant IDs to plot, one per line
    # the variants will still be restricted to DELs and DUPs
    File? variant_ids
    Boolean one_sample_per_plot = false
    String plot_prefix
    Int min_size
    Int variants_per_shard = 40
    # padding to use around CNV when plotting, as fraction of CNV length
    Float padding = 0.5

    # Batch (sample set) IDs and their sample lists, one sample ID per line.
    # All arrays of batch data must match the order of batch_name_list.
    Array[String] batch_name_list
    Array[File] batch_sample_lists
    Array[File] medians_files
    Array[File] bincov_files
    Array[File] bincov_index_files

    String sv_base_mini_docker
    String sv_pipeline_docker

    RuntimeAttr? runtime_attr_extract_cnvs
    RuntimeAttr? runtime_attr_make_bincov_query_manifests
    RuntimeAttr? runtime_attr_make_bincov_map
    RuntimeAttr? runtime_attr_make_bincov_index_map
    RuntimeAttr? runtime_attr_make_medians_map
    RuntimeAttr? runtime_attr_subset_bincov_matrix
    RuntimeAttr? runtime_attr_make_bincov_subset_map
    RuntimeAttr? runtime_attr_shard_variants
    RuntimeAttr? runtime_attr_make_plots
    RuntimeAttr? runtime_attr_merge_plot_tars
  }

  if (!defined(cnvs)) {
    scatter (vcf in vcfs) {
      call ExtractCnvs {
        input:
          vcf = vcf,
          variant_ids = variant_ids,
          min_size = min_size,
          sv_base_mini_docker = sv_base_mini_docker,
          runtime_attr_override = runtime_attr_extract_cnvs
      }
    }
  }

  call MakeBincovQueryManifests {
    input:
      cnvs = select_first([cnvs, ExtractCnvs.cnvs]),
      padding = padding,
      batch_name_list = batch_name_list,
      batch_sample_lists = batch_sample_lists,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_make_bincov_query_manifests
  }

  call MapFromArrays as MakeBincovMap {
    input:
      keys = batch_name_list,
      values = bincov_files,
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_make_bincov_map
  }

  call MapFromArrays as MakeBincovIndexMap {
    input:
      keys = batch_name_list,
      values = bincov_index_files,
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_make_bincov_index_map
  }

  call MapFromArrays as MakeMediansMap {
    input:
      keys = batch_name_list,
      values = medians_files,
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_make_medians_map
  }

  scatter (f in MakeBincovQueryManifests.manifests) {
    String batch = basename(f, ".rdx")
    call SubsetBincovMatrix {
      input:
        batch = batch,
        manifest = f,
        bincov = MakeBincovMap.m[batch],
        bincov_index = MakeBincovIndexMap.m[batch],
        medians = MakeMediansMap.m[batch],
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_subset_bincov_matrix
    }
  }

  call MapFromArrays as MakeBincovSubsetMap {
    input:
      keys = SubsetBincovMatrix.batch_out,
      values = SubsetBincovMatrix.subsets,
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_make_bincov_subset_map
  }

  call ShardVariants {
    input:
      variants = MakeBincovQueryManifests.merged_cnvs,
      sample_table = MakeBincovQueryManifests.sample_table,
      variants_per_shard = variants_per_shard,
      bincov_subset_map = write_map(MakeBincovSubsetMap.m),
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_shard_variants
  }

  scatter (i in range(length(ShardVariants.shards))) {
    call MakePlots {
      input:
        variants = ShardVariants.shards[i],
        batches = ShardVariants.batches[i],
        one_sample_per_plot = one_sample_per_plot,
        bincov_tars = read_lines(ShardVariants.bincov_paths[i]),
        sample_table = MakeBincovQueryManifests.sample_table,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_make_plots
    }
  }

  call MergePlotTars {
    input:
      plot_tars = MakePlots.plots,
      plot_prefix = plot_prefix,
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_merge_plot_tars
  }

  output {
    File cnv_plots = MergePlotTars.plots
  }
}

task ExtractCnvs {
  input {
    File vcf
    File? variant_ids
    Int min_size
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  Float disk_size = size(vcf, "GB") * 3 + 16

  RuntimeAttr default_attr = object {
    cpu_cores: 2,
    mem_gb: 4.0,
    disk_gb: ceil(disk_size),
    boot_disk_gb: 8,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_base_mini_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    vcf='~{vcf}'
    min_size='~{min_size}'
    variant_ids='~{if defined(variant_ids) then variant_ids else ""}'

    format='%CHROM\t%POS\t%INFO/END\t%ID\t%INFO/SVTYPE\t[%SAMPLE,]\n'
    svtype_filter='(INFO/SVTYPE == "DEL" || INFO/SVTYPE == "DUP")'
    gt_filter='GT = "alt"'
    if [[ -n "${variant_ids:-}" ]]; then
      mv "${variant_ids}" variant_ids.list
      filter="${svtype_filter} & ID=@variant_ids.list & ${gt_filter}"
    else
      filter="${svtype_filter} & FILTER ~ \"PASS\" & INFO/SVLEN >= ${min_size} & ${gt_filter}"
    fi

    bcftools query --include "${filter}" --format "${format}" "${vcf}" \
      | awk -F'\t' '{sub(/,$/, "", $6); print}' OFS='\t' > cnvs.tsv
  >>>

  output {
    File cnvs = "cnvs.tsv"
  }
}

task MakeBincovQueryManifests {
  input {
    Array[File] cnvs
    Float padding
    Array[String] batch_name_list
    Array[File] batch_sample_lists
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float disk_size = (size(cnvs, "GB") + size(batch_sample_lists, "GB")) * 3 + 30

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 4.0,
    disk_gb: ceil(disk_size),
    boot_disk_gb: 8,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    cnvs='~{write_lines(cnvs)}'
    batch_names='~{write_lines(batch_name_list)}'
    batch_sample_lists='~{write_lines(batch_sample_lists)}'
    padding=~{padding}

    if [[ $(wc -l < "${batch_names}") -ne $(wc -l < "${batch_sample_lists}") ]]; then
      echo "batch_name_list and batch_sample_lists must have the same length" >&2
      exit 1
    fi

    # Build a headerless two-column TSV of sample ID and batch ID
    sample_table='sample_table.tsv'
    : > "${sample_table}"
    while IFS=$'\t' read -r batch samples_file; do
      awk -F'\t' -v batch="${batch}" 'NF > 0 { print $1 "\t" batch }' "${samples_file}" >> "${sample_table}"
    done < <(paste "${batch_names}" "${batch_sample_lists}")

    cat "${cnvs}" | xargs cat > merged_cnvs.tsv

    mkdir manifests
    Rscript /opt/sv-pipeline/scripts/batch_variants.R merged_cnvs.tsv \
      "${sample_table}" manifests "${padding}"
  >>>

  output {
    File merged_cnvs = "merged_cnvs.tsv"
    File sample_table = "sample_table.tsv"
    Array[File] manifests = glob("manifests/*.rdx")
  }
}

task SubsetBincovMatrix {
  input {
    String batch
    File manifest
    File bincov
    File bincov_index
    File medians
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float disk_size = size(bincov, "GB") * 2
    + size(manifest, "GB")
    + size(medians, "GB")
    + 16

  RuntimeAttr default_attr = object {
    cpu_cores: 2,
    mem_gb: 4.0,
    disk_gb: ceil(disk_size),
    boot_disk_gb: 8,
    preemptible_tries: 2,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    batch='~{batch}'
    manifest='~{manifest}'
    bincov='~{bincov}'
    medians='~{medians}'

    mkdir "${batch}"
    Rscript /opt/sv-pipeline/scripts/subset_bincov.R "${manifest}" \
      "${bincov}" "${medians}" "${batch}"

    tar -cvf "${batch}.tar" "${batch}"
  >>>

  output {
    String batch_out = batch
    File subsets = batch + ".tar"
  }
}

task ShardVariants {
  input {
    File variants
    File sample_table
    Int variants_per_shard
    File bincov_subset_map
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  Float disk_size = size(variants, "GB") * 3
    + size(sample_table, "GB")
    + size(bincov_subset_map, "GB")
    + 16

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 2.0,
    disk_gb: ceil(disk_size),
    boot_disk_gb: 8,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_base_mini_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    variants='~{variants}'
    sample_table='~{sample_table}'
    variants_per_shard=~{variants_per_shard}
    bincov_subset_map='~{bincov_subset_map}'

    mkdir shards batches bincov
    split -l "${variants_per_shard}" "${variants}" shards/cnvs_

    # The order of the batch ids in the list for each shard must match the
    # order in the lists of the other data files or the batches could
    # map to the wrong files.
    gawk -F'\t' '
    ARGIND == 1 {
      Sample_arr[$1] = $2
    }
    ARGIND == 2 {
      Bincov_arr[$1] = $2
    }
    ARGIND > 2 {
      split($6, b, /,/)
      for (i in b) {
        Batches[Sample_arr[b[i]]]
      }
    }
    BEGINFILE {
      if (ARGIND > 2) {
        n = split(FILENAME, p, /\//)
        Shard = p[n]
        delete Batches
      }
    }
    ENDFILE {
      if (ARGIND > 2) {
        batches_out = "batches/" Shard
        bincovs_out = "bincov/" Shard
        for (id in Batches) {
          print id > batches_out
          print Bincov_arr[id] > bincovs_out
        }
        close(batches_out)
        close(bincovs_out)
      }
    }' "${sample_table}" "${bincov_subset_map}" shards/cnvs_*
  >>>

  output {
    # Each shard will have the CNVs to plot and each other output with the
    # same basename will have the list of files needed so it is critical
    # that the relative ordering of the files from the glob is the same.
    Array[File] shards = glob("shards/cnvs_*")
    Array[File] batches = glob("batches/cnvs_*")
    Array[File] bincov_paths = glob("bincov/cnvs_*")
  }
}

task MakePlots {
  input {
    File variants
    File batches
    Boolean one_sample_per_plot
    Array[File] bincov_tars
    File sample_table
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Int variant_count = length(read_lines(variants))
  Float input_size = size(bincov_tars, "GB") * 2
    + size(sample_table, "GB")
    + size(variants, "GB")
    + size(batches, "GB")
  Float disk_size = input_size + variant_count * 0.01 + 16

  RuntimeAttr default_attr = object {
    cpu_cores: 2,
    mem_gb: 4.0,
    disk_gb: ceil(disk_size),
    boot_disk_gb: 8,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    variants='~{variants}'
    batches='~{batches}'
    bincov_tars='~{write_lines(bincov_tars)}'
    sample_table='~{sample_table}'
    one_sample_per_plot=~{if one_sample_per_plot then 1 else 0}

    while IFS=$'\t' read -r batch bincov; do
      bn="$(basename "${bincov}" .tar)"
      tar -xf "${bincov}"
      printf '%s\t%s\n' "${batch}" "${bn}" >> bincov_map.tsv
    done < <(paste "${batches}" "${bincov_tars}")

    Rscript /opt/sv-pipeline/scripts/visualize_cnvs.R \
      "${variants}" \
      "${sample_table}" \
      bincov_map.tsv \
      plots \
      "${one_sample_per_plot}"

    tar -cvzf plots.tar.gz plots
  >>>

  output {
    File plots = "plots.tar.gz"
  }
}

task MergePlotTars {
  input {
    Array[File] plot_tars
    String plot_prefix
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  Float disk_size = size(plot_tars, "GB") * 3 + 16

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 1.0,
    disk_gb: ceil(disk_size),
    boot_disk_gb: 8,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_base_mini_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    plot_tars='~{write_lines(plot_tars)}'
    plot_prefix='~{plot_prefix}'

    mkdir temp
    while read -r f; do
      tar -xzf "${f}"
      mv -t temp plots/*.jpg
      rm -r plots
    done < "${plot_tars}"

    mv temp "${plot_prefix}"
    tar -czf "${plot_prefix}.tar.gz" "${plot_prefix}"
  >>>

  output {
    File plots = "${plot_prefix}.tar.gz"
  }
}

task MapFromArrays {
  input {
    Array[String] keys
    Array[String] values
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 1.0,
    disk_gb: 16,
    boot_disk_gb: 8,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_base_mini_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    keys='~{write_lines(keys)}'
    values='~{write_lines(values)}'

    n_keys=$(wc -l < "${keys}")
    n_values=$(wc -l < "${values}")
    if [[ "${n_keys}" -ne "${n_values}" ]]; then
      echo "ERROR: got ${n_keys} keys but ${n_values} values; the arrays must have the same length and order" >&2
      exit 1
    fi

    paste "${keys}" "${values}" > map.tsv
  >>>

  output {
    Map[String, String] m = read_map("map.tsv")
  }
}