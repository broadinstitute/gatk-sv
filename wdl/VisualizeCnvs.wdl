version 1.0

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

    File sample_table # TSV with sample_id, sample_set_id
    Array[String] sample_set_ids
    Array[File] medians_files
    Array[File] bincov_files
    Array[File] bincov_index_files

    String base_docker
    String r_docker
  }

  if (!defined(cnvs)) {
    scatter (vcf in vcfs) {
      call ExtractCnvs {
        input:
          vcf = vcf,
          variant_ids = variant_ids,
          min_size = min_size,
          base_docker = base_docker
      }
    }
  }

  call MakeBincovQueryManifests {
    input:
      cnvs = select_first([cnvs, ExtractCnvs.cnvs]),
      padding = padding,
      sample_table = sample_table,
      r_docker = r_docker
  }

  call MapFromArrays as MakeBincovMap {
    input:
      keys = sample_set_ids,
      values = bincov_files,
      base_docker = base_docker
  }

  call MapFromArrays as MakeBincovIndexMap {
    input:
      keys = sample_set_ids,
      values = bincov_index_files,
      base_docker = base_docker
  }

  call MapFromArrays as MakeMediansMap {
    input:
      keys = sample_set_ids,
      values = medians_files,
      base_docker = base_docker
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
        r_docker = r_docker
    }
  }

  call MapFromArrays as MakeBincovSubsetMap {
    input:
      keys = SubsetBincovMatrix.batch_out,
      values = SubsetBincovMatrix.subsets,
      base_docker = base_docker
  }

  call ShardVariants {
    input:
      variants = MakeBincovQueryManifests.merged_cnvs,
      sample_table = sample_table,
      variants_per_shard = variants_per_shard,
      bincov_subset_map = write_map(MakeBincovSubsetMap.m),
      base_docker = base_docker
  }

  scatter (i in range(length(ShardVariants.shards))) {
    call MakePlots {
      input:
        variants = ShardVariants.shards[i],
        batches = ShardVariants.batches[i],
        one_sample_per_plot = one_sample_per_plot,
        bincov_tars = read_lines(ShardVariants.bincov_paths[i]),
        sample_table = sample_table,
        r_docker = r_docker
    }
  }

  call MergePlotTars {
    input:
      plot_tars = MakePlots.plots,
      plot_prefix = plot_prefix,
      base_docker = base_docker
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
    String base_docker
  }

  Float disk_size = size(vcf, "GB") * 3 + 16

  runtime {
    bootDiskSizeGb: 8
    cpu: 2
    disks: "local-disk ${ceil(disk_size)} HDD"
    docker: base_docker
    maxRetries: 1
    memory: "4 GiB"
    preemptible: 3
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
    File sample_table
    String r_docker
  }

  Float disk_size = (size(cnvs, "GB") + size(sample_table, "GB")) * 3 + 30

  runtime {
    bootDiskSizeGb: 8
    cpu: 1
    disks: "local-disk ${ceil(disk_size)} HDD"
    docker: r_docker
    maxRetries: 1
    memory: "4 GiB"
    preemptible: 3
  }

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    cnvs='~{write_lines(cnvs)}'
    sample_table='~{sample_table}'
    padding=~{padding}

    cat "${cnvs}" | xargs cat > merged_cnvs.tsv

    mkdir manifests
    Rscript /opt/gatk-sv-utils/scripts/batch_variants.R merged_cnvs.tsv \
      "${sample_table}" manifests "${padding}"
  >>>

  output {
    File merged_cnvs = "merged_cnvs.tsv"
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
    String r_docker
  }

  Float disk_size = size(bincov, "GB") * 2
    + size(manifest, "GB")
    + size(medians, "GB")
    + 16

  runtime {
    bootDiskSizeGb: 8
    cpu: 2
    disks: "local-disk ${ceil(disk_size)} HDD"
    docker: r_docker
    maxRetries: 1
    memory: "4 GiB"
    preemptible: 2
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
    Rscript /opt/gatk-sv-utils/scripts/subset_bincov.R "${manifest}" \
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
    String base_docker
  }

  Float disk_size = size(variants, "GB") * 3
    + size(sample_table, "GB")
    + size(bincov_subset_map, "GB")
    + 16

  runtime {
    bootDiskSizeGb: 8
    cpu: 1
    disks: "local-disk ${ceil(disk_size)} HDD"
    docker: base_docker
    maxRetries: 1
    memory: "2 GiB"
    preemptible: 3
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
    String r_docker
  }

  Int variant_count = length(read_lines(variants))
  Float input_size = size(bincov_tars, "GB") * 2
    + size(sample_table, "GB")
    + size(variants, "GB")
    + size(batches, "GB")
  Float disk_size = input_size + variant_count * 0.01 + 16

  runtime {
    bootDiskSizeGb: 8
    cpu: 2
    disks: "local-disk ${ceil(disk_size)} HDD"
    docker: r_docker
    maxRetries: 1
    memory: "4 GiB"
    preemptible: 3
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

    Rscript /opt/gatk-sv-utils/scripts/visualize_cnvs.R \
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
    String base_docker
  }

  Float disk_size = size(plot_tars, "GB") * 3 + 16

  runtime {
    bootDiskSizeGb: 8
    cpu: 1
    disks: "local-disk ${ceil(disk_size)} HDD"
    docker: base_docker
    maxRetries: 1
    memory: "1 GiB"
    preemptible: 3
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
    String base_docker
  }

  runtime {
    bootDiskSizeGb: 8
    cpu: 1
    disks: "local-disk 16 HDD"
    docker: base_docker
    maxRetries: 1
    memory: "1 GiB"
    preemptible: 3
  }

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    paste '~{write_lines(keys)}' '~{write_lines(values)}' > map.tsv
  >>>

  output {
    Map[String, String] m = read_map("map.tsv")
  }
}